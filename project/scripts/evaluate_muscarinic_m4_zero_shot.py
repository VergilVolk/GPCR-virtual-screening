#!/usr/bin/env python3
"""Evaluate frozen DrugCLIP transfer to an unseen M4 pocket and unseen molecules."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from torch import nn
from torch.nn import functional as F


class Projection(nn.Module):
    def __init__(self, state: dict):
        super().__init__()
        self.linear1 = nn.Linear(state["linear1.weight"].shape[1], state["linear1.weight"].shape[0])
        self.linear2 = nn.Linear(state["linear2.weight"].shape[1], state["linear2.weight"].shape[0])
        self.load_state_dict(state)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return F.normalize(self.linear2(F.relu(self.linear1(x))), dim=-1)


def matrix(molecule: np.ndarray, pocket: np.ndarray, state: dict) -> np.ndarray:
    m = Projection(copy.deepcopy(state["mol_project"])).eval()
    p = Projection(copy.deepcopy(state["pocket_project"])).eval()
    with torch.inference_mode():
        return (m(torch.as_tensor(molecule)) @ p(torch.as_tensor(pocket)).T).numpy()


def deltas(scores: np.ndarray, test: pd.DataFrame, mi: dict, pi: dict,
           molecule_override: dict[str, str] | None = None) -> np.ndarray:
    values = []
    for row in test.itertuples():
        molecule = molecule_override.get(row.canonical_smiles, row.canonical_smiles) if molecule_override else row.canonical_smiles
        values.append(scores[mi[molecule], pi[row.positive_subtype]] -
                      scores[mi[molecule], pi[row.negative_subtype]])
    return np.asarray(values)


def metrics(test: pd.DataFrame, delta: np.ndarray) -> dict:
    by_direction = {}
    for direction, group in test.groupby("m4_direction"):
        by_direction[direction] = float((delta[group.index] > 0).mean())
    return {
        "pair_accuracy": float((delta > 0).mean()),
        "direction_accuracy": by_direction,
        "direction_balanced_accuracy": float(np.mean(list(by_direction.values()))),
        "delta_vs_pchembl_gap_spearman": float(spearmanr(delta, test.delta_pchembl).statistic),
    }


def cluster_delta_ci(test: pd.DataFrame, official: np.ndarray, tuned: np.ndarray,
                     n: int, seed: int) -> list[float]:
    molecules = test.canonical_smiles.unique(); rng = np.random.default_rng(seed); values = []
    for _ in range(n):
        sampled = rng.choice(molecules, len(molecules), True)
        index = np.concatenate([np.flatnonzero(test.canonical_smiles.to_numpy() == value) for value in sampled])
        values.append((tuned[index] > 0).mean() - (official[index] > 0).mean())
    return list(map(float, np.quantile(values, [0.025, 0.5, 0.975])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--official-projection", type=Path, required=True)
    parser.add_argument("--tuned-projections", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--permutations", type=int, default=5000)
    args = parser.parse_args()
    arc = np.load(args.representations, allow_pickle=False)
    molecule_ids = list(map(str, arc["molecule_ids"])); pocket_ids = list(map(str, arc["pocket_ids"]))
    test = pd.read_csv(args.test)
    # Retrieval dataset exposes canonical SMILES as molecule IDs.
    mi = {value: i for i, value in enumerate(molecule_ids)}
    pi = {value.split("_")[0]: i for i, value in enumerate(pocket_ids)}
    official_state = torch.load(args.official_projection, map_location="cpu")
    tuned_states = [torch.load(path, map_location="cpu") for path in args.tuned_projections]
    official_matrix = matrix(arc["molecule_representations"].astype(np.float32),
                             arc["pocket_representations"].astype(np.float32), official_state)
    tuned_matrix = np.mean([
        matrix(arc["molecule_representations"].astype(np.float32),
               arc["pocket_representations"].astype(np.float32), state)
        for state in tuned_states
    ], axis=0)
    official = deltas(official_matrix, test, mi, pi); tuned = deltas(tuned_matrix, test, mi, pi)
    rng = np.random.default_rng(20260924); molecules = test.canonical_smiles.unique(); null = []
    for _ in range(args.permutations):
        shuffled = rng.permutation(molecules); mapping = dict(zip(molecules, shuffled))
        null.append(float((deltas(tuned_matrix, test, mi, pi, mapping) > 0).mean()))
    actual = float((tuned > 0).mean())
    report = {
        "protocol": "M4 pocket and molecules absent from triplet training; frozen posthoc zero-shot test",
        "pairs": int(len(test)), "molecules": int(test.canonical_smiles.nunique()),
        "tuned_checkpoints": list(map(str, args.tuned_projections)),
        "official_drugclip": metrics(test, official),
        "gpcr_triplet": metrics(test, tuned),
        "gpcr_triplet_minus_official_molecule_cluster_95ci": cluster_delta_ci(
            test, official, tuned, args.bootstrap, 20260924),
        "molecule_permutation_control": {
            "actual_accuracy": actual,
            "null_median": float(np.median(null)),
            "null_95_interval": list(map(float, np.quantile(null, [0.025, 0.975]))),
            "one_sided_p": float((1 + np.sum(np.asarray(null) >= actual)) / (1 + len(null))),
        },
        "claim_boundary": "M4-vs-subtype activity transfer, not allostery, PAM efficacy, or prospective validation.",
    }
    predictions = test.copy(); predictions["official_delta"] = official; predictions["gpcr_triplet_delta"] = tuned
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    predictions.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
