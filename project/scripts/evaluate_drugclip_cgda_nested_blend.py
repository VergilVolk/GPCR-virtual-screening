#!/usr/bin/env python3
"""Nested target-LOO shrinkage between two frozen CGDA score channels."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from evaluate_drugclip_official_litpcba_embeddings import metrics


ALPHAS = [0.0, 0.25, 0.5, 0.75, 1.0]
NAMES = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--selection-metric", default="bedroc_alpha80_5")
    parser.add_argument("--left-score", choices=["pocket", "reference", "cgda"], default="reference")
    parser.add_argument("--right-score", choices=["pocket", "reference", "cgda"], default="cgda")
    args = parser.parse_args()
    directories = [path.with_suffix("") for path in args.runs]
    targets = sorted(path.stem for path in directories[0].glob("*.npz"))
    fixed = {}
    for target in targets:
        records = [np.load(directory / f"{target}.npz") for directory in directories]
        labels = records[0]["labels"]
        left = np.mean([record[args.left_score] for record in records], axis=0)
        right = np.mean([record[args.right_score] for record in records], axis=0)
        fixed[target] = {
            str(alpha): metrics(labels, (1 - alpha) * left + alpha * right)
            for alpha in ALPHAS
        }
    folds, counts = {}, {str(alpha): 0 for alpha in ALPHAS}
    for held in targets:
        seen = [target for target in targets if target != held]
        validation = {str(alpha): float(np.mean([fixed[target][str(alpha)][args.selection_metric]
                                                 for target in seen])) for alpha in ALPHAS}
        selected = max(ALPHAS, key=lambda alpha: validation[str(alpha)])
        counts[str(selected)] += 1
        folds[held] = {"selected_alpha": selected, "seen_validation": validation,
                       "held_metrics": fixed[held][str(selected)]}
    nested = {name: float(np.mean([folds[t]["held_metrics"][name] for t in targets])) for name in NAMES}
    fixed_macro = {str(alpha): {name: float(np.mean([fixed[t][str(alpha)][name] for t in targets]))
                                for name in NAMES} for alpha in ALPHAS}
    report = {"protocol": "Outer target-LOO blend selection", "selection_metric": args.selection_metric,
              "left_score": args.left_score, "right_score": args.right_score,
              "nested_macro": nested, "fixed_macro": fixed_macro,
              "selected_counts": counts, "folds": folds}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"nested_macro": nested, "fixed_macro": fixed_macro,
                      "selected_counts": counts}, indent=2))


if __name__ == "__main__":
    main()
