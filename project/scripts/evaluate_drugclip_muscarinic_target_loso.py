#!/usr/bin/env python3
"""Strict leave-one-muscarinic-target-out DrugCLIP triplet evaluation."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import torch

from finetune_drugclip_muscarinic_triplet import Model, fit, clustered_bootstrap


def pair_scores(model, molecules, pockets, frame, mi, pi):
    scores = model.matrix(molecules, pockets).detach().numpy()
    return np.asarray([scores[mi[r.canonical_smiles], pi[r.positive_subtype]] -
                       scores[mi[r.canonical_smiles], pi[r.negative_subtype]]
                       for r in frame.itertuples()])


def direction_metrics(frame, score, held):
    preferred = frame.positive_subtype.to_numpy() == held
    disfavored = frame.negative_subtype.to_numpy() == held
    a = float((score[preferred] > 0).mean()); b = float((score[disfavored] > 0).mean())
    return {"accuracy": float((score > 0).mean()), "held_preferred_accuracy": a,
            "held_disfavored_accuracy": b, "direction_balanced_accuracy": (a + b) / 2}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", default="20260924,20260925,20260926")
    parser.add_argument("--epochs", type=int, default=500)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--margin", type=float, default=0.1)
    parser.add_argument("--preserve-weight", type=float, default=0.2)
    args = parser.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    arc = np.load(args.representations, allow_pickle=False); initial = torch.load(args.projection, map_location="cpu")
    mids = list(map(str, arc["molecule_ids"])); pids = list(map(str, arc["pocket_ids"]))
    mi = {value: i for i, value in enumerate(mids)}; pi = {value.split("_")[0]: i for i, value in enumerate(pids)}
    molecules = torch.as_tensor(arc["molecule_representations"], dtype=torch.float32)
    pockets = torch.as_tensor(arc["pocket_representations"], dtype=torch.float32)
    frozen = Model(copy.deepcopy(initial)).eval()
    cfg = SimpleNamespace(lr=args.lr, weight_decay=1e-3, epochs=args.epochs,
                          margin=args.margin, preserve_weight=args.preserve_weight,
                          train_scope="pocket")
    rows, reports = [], {}
    for held in ["M1", "M2", "M3", "M4", "M5"]:
        train_raw = pd.read_csv(args.benchmark_dir / f"{held}_train.csv")
        test_raw = pd.read_csv(args.benchmark_dir / f"{held}_test.csv")
        train = train_raw[train_raw.canonical_smiles.isin(mi)].reset_index(drop=True)
        test = test_raw[test_raw.canonical_smiles.isin(mi)].reset_index(drop=True)
        base = pair_scores(frozen, molecules, pockets, test, mi, pi)
        targeted, controls = [], []
        for seed in map(int, args.seeds.split(",")):
            target_model = fit(molecules, pockets, train, mi, pi, initial, seed, cfg, False)
            random_model = fit(molecules, pockets, train, mi, pi, initial, seed, cfg, True)
            target_score = pair_scores(target_model, molecules, pockets, test, mi, pi)
            random_score = pair_scores(random_model, molecules, pockets, test, mi, pi)
            targeted.append(target_score); controls.append(random_score)
            for i, row in test.iterrows():
                rows.append({"held_subtype": held, "canonical_smiles": row.canonical_smiles,
                             "positive_subtype": row.positive_subtype, "negative_subtype": row.negative_subtype,
                             "seed": seed, "frozen_delta": float(base[i]),
                             "random_control_delta": float(random_score[i]),
                             "targeted_delta": float(target_score[i])})
        target_mean = np.mean(targeted, axis=0); control_mean = np.mean(controls, axis=0)
        reports[held] = {
            "train_pairs": int(len(train)), "train_molecules": int(train.canonical_smiles.nunique()),
            "test_pairs": int(len(test)), "test_molecules": int(test.canonical_smiles.nunique()),
            "conformer_failure_train_pairs_removed": int(len(train_raw) - len(train)),
            "conformer_failure_test_pairs_removed": int(len(test_raw) - len(test)),
            "frozen": direction_metrics(test, base, held),
            "random_control": direction_metrics(test, control_mean, held),
            "targeted": direction_metrics(test, target_mean, held),
            "targeted_minus_frozen_molecule_cluster_95ci": clustered_bootstrap(test, base, target_mean),
            "targeted_minus_random_molecule_cluster_95ci": clustered_bootstrap(test, control_mean, target_mean),
            "per_seed_targeted_accuracy": [float((value > 0).mean()) for value in targeted],
        }
    macro = {}
    for method in ["frozen", "random_control", "targeted"]:
        macro[method] = {
            "accuracy": float(np.mean([reports[h][method]["accuracy"] for h in reports])),
            "direction_balanced_accuracy": float(np.mean([
                reports[h][method]["direction_balanced_accuracy"] for h in reports])),
        }
    report = {
        "method": "DrugCLIP pocket projection triplet; strict leave-one-target-out",
        "folds": reports, "macro_average": macro,
        "hyperparameters": {"epochs": args.epochs, "lr": args.lr, "margin": args.margin,
                            "preserve_weight": args.preserve_weight, "seeds": args.seeds},
        "claim_boundary": "Posthoc heterogeneous ChEMBL target-transfer benchmark; not M4 PAM efficacy.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(rows).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
