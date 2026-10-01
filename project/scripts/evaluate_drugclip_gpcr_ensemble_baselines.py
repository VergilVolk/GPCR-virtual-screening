#!/usr/bin/env python3
"""Evaluate frozen DrugCLIP over published GaMD cluster ensembles."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from finetune_drugclip_gpcr_ensemble_screening import aggregate, all_cluster_scores
from finetune_drugclip_gpcr_screening import TARGETS, screening_metrics
from finetune_drugclip_muscarinic_triplet import Proj


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--pocket-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.1)
    args = parser.parse_args()
    archive = np.load(args.representations, allow_pickle=False)
    initial = torch.load(args.projection, map_location="cpu")
    index = {value: i for i, value in enumerate(map(str, archive["molecule_ids"]))}
    pairs_raw = pd.read_csv(args.pairs)
    pairs = pairs_raw[pairs_raw.canonical_smiles.astype(str).isin(index)].copy().reset_index(drop=True)
    pair_mol = torch.as_tensor([index[value] for value in pairs.canonical_smiles.astype(str)], dtype=torch.long)
    pair_target = torch.as_tensor([TARGETS.index(value) for value in pairs.target], dtype=torch.long)
    expected = [f"{target}_cluster{cluster}" for target in TARGETS for cluster in range(10)]
    raw_ids = list(map(str, archive["pocket_ids"]))
    if set(raw_ids) != set(expected) or len(raw_ids) != len(expected):
        raise ValueError("Pocket ID set mismatch")
    pocket_order = [raw_ids.index(value) for value in expected]
    metadata = pd.read_csv(args.pocket_metadata).set_index("pocket_id").loc[expected]
    populations = torch.as_tensor(metadata.population.to_numpy(np.float32).reshape(4, 10))
    populations /= populations.sum(1, keepdim=True)
    mol = Proj(copy.deepcopy(initial["mol_project"])).eval()
    pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode():
        zm = mol(torch.as_tensor(archive["molecule_representations"].astype(np.float32)))
        pocket_rep = archive["pocket_representations"][pocket_order].astype(np.float32)
        zp = pocket(torch.as_tensor(pocket_rep))
        clusters = all_cluster_scores(zm, zp, pair_mol, pair_target)
        tau = torch.full((4,), args.tau)
        scores = {
            "cluster0": aggregate(clusters, populations, pair_target, "cluster0").numpy(),
            "population_mean": aggregate(clusters, populations, pair_target, "population_mean").numpy(),
            "max": aggregate(clusters, populations, pair_target, "max").numpy(),
            "population_lse": aggregate(clusters, populations, pair_target, "population_lse", tau).numpy(),
        }
    report = {"pairs_evaluated": int(len(pairs)), "input_pairs_excluded": int(len(pairs_raw) - len(pairs)),
              "metrics": {name: screening_metrics(pairs, value) for name, value in scores.items()},
              "tau": args.tau,
              "claim_boundary": "Frozen retrospective DrugCLIP ensemble baselines; no efficacy claim."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame({"pair_id": pairs.pair_id, **scores}).to_csv(args.output.with_suffix(".csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
