#!/usr/bin/env python3
"""Paired target-bootstrap for the nested CGDA shrinkage result."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


NAMES = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nested", type=Path, required=True)
    parser.add_argument("--ensemble-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()
    nested = json.loads(args.nested.read_text(encoding="utf-8"))
    ensemble = json.loads(args.ensemble_summary.read_text(encoding="utf-8"))
    targets = sorted(nested["folds"])
    delta = np.asarray([
        [nested["folds"][target]["held_metrics"][name] -
         ensemble["per_target"][target]["reference"][name] for name in NAMES]
        for target in targets
    ])
    rng = np.random.default_rng(args.seed)
    samples = np.empty((args.bootstrap, len(NAMES)))
    for index in range(args.bootstrap):
        samples[index] = delta[rng.integers(0, len(targets), len(targets))].mean(axis=0)
    paired = {name: {"mean_delta": float(delta[:, i].mean()),
                     "ci95": [float(v) for v in np.quantile(samples[:, i], [0.025, 0.975])],
                     "targets_improved": int((delta[:, i] > 0).sum())}
              for i, name in enumerate(NAMES)}
    report = {"protocol": "paired outer-target bootstrap for nested CGDA shrinkage",
              "n_targets": len(targets), "paired": paired,
              "claim_boundary": "Retrospective target-level uncertainty; no PAM efficacy claim."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
