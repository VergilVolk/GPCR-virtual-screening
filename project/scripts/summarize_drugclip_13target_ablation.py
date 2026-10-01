#!/usr/bin/env python3
"""Paired scaffold-cluster bootstrap for 13-target adapter ablations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from finetune_drugclip_gpcr_screening import binary_metrics


METRICS = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]


def macro(table, scores, targets):
    values = {metric: [] for metric in METRICS}
    for target in targets:
        keep = table.target.eq(target).to_numpy()
        result = binary_metrics(table.loc[keep, "label"].to_numpy(int), scores[keep])
        for metric in METRICS:
            values[metric].append(result[metric])
    return {metric: float(np.mean(items)) for metric, items in values.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, action="append", required=True)
    parser.add_argument("--challenger-predictions", type=Path, action="append")
    parser.add_argument("--reference", default="dual_standard")
    parser.add_argument("--challenger", default="dual_no_preserve")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    frames = [pd.read_csv(path) for path in args.predictions]
    challenger_frames = ([pd.read_csv(path) for path in args.challenger_predictions]
                         if args.challenger_predictions else frames)
    if len(frames) != len(challenger_frames):
        raise ValueError("Reference and challenger seed counts differ")
    table = frames[0]
    identity = ["target", "label", "canonical_smiles", "murcko_scaffold"]
    for frame in frames[1:] + challenger_frames:
        if not frame[identity].equals(table[identity]):
            raise ValueError("Prediction row mismatch across seeds")
    if any(args.reference not in frame for frame in frames):
        raise ValueError(f"Missing reference score column: {args.reference}")
    if any(args.challenger not in frame for frame in challenger_frames):
        raise ValueError(f"Missing challenger score column: {args.challenger}")
    targets = list(table.target.drop_duplicates())
    groups = {(target, scaffold): np.asarray(list(rows), dtype=int)
              for (target, scaffold), rows in table.groupby(["target", "murcko_scaffold"]).groups.items()}
    by_target = {target: [key for key in groups if key[0] == target] for target in targets}
    ref_by_seed = [frame[args.reference].to_numpy(float) for frame in frames]
    challenger_by_seed = [frame[args.challenger].to_numpy(float) for frame in challenger_frames]
    ref = np.mean(ref_by_seed, axis=0)
    challenger = np.mean(challenger_by_seed, axis=0)
    observed_ref = macro(table, ref, targets)
    observed_challenger = macro(table, challenger, targets)
    rng = np.random.default_rng(args.seed)
    differences = {metric: [] for metric in METRICS}
    for _ in range(args.bootstrap):
        sampled = []
        for target in targets:
            keys = by_target[target]
            sampled.extend(np.concatenate([groups[keys[i]] for i in rng.integers(0, len(keys), len(keys))]))
        sampled = np.asarray(sampled, dtype=int)
        sample = table.iloc[sampled].reset_index(drop=True)
        ref_metric = macro(sample, ref[sampled], targets)
        challenger_metric = macro(sample, challenger[sampled], targets)
        for metric in METRICS:
            differences[metric].append(challenger_metric[metric] - ref_metric[metric])
    report = {
        "protocol": "paired target-stratified Murcko-scaffold cluster bootstrap",
        "n_seeds": len(frames),
        "reference": args.reference,
        "challenger": args.challenger,
        "observed": {"reference": observed_ref, "challenger": observed_challenger},
        "per_seed_challenger_minus_reference": [
            {metric: macro(table, challenger_seed, targets)[metric] - macro(table, ref_seed, targets)[metric]
             for metric in METRICS}
            for ref_seed, challenger_seed in zip(ref_by_seed, challenger_by_seed)
        ],
        "challenger_minus_reference_95ci": {
            metric: list(map(float, np.quantile(values, [0.025, 0.5, 0.975])))
            for metric, values in differences.items()
        },
        "claim_boundary": "Retrospective multi-seed paired uncertainty; external confirmation remains required.",
    }
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
