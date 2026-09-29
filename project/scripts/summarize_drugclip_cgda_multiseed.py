#!/usr/bin/env python3
"""Three-seed CGDA ensemble and paired target-bootstrap summary."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from evaluate_drugclip_official_litpcba_embeddings import metrics


METRICS = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    run_dirs = [path.with_suffix("") for path in args.runs]
    targets = sorted(path.stem for path in run_dirs[0].glob("*.npz"))
    per_target = {}
    for target in targets:
        records = [np.load(directory / f"{target}.npz") for directory in run_dirs]
        labels = records[0]["labels"]
        reference = records[0]["reference"]
        ensemble = np.mean([record["cgda"] for record in records], axis=0)
        per_target[target] = {
            "reference": metrics(labels, reference),
            "cgda_ensemble": metrics(labels, ensemble),
        }

    macro = {
        method: {name: float(np.mean([per_target[target][method][name] for target in targets])) for name in METRICS}
        for method in ("reference", "cgda_ensemble")
    }
    differences = np.asarray([
        [per_target[target]["cgda_ensemble"][name] - per_target[target]["reference"][name] for name in METRICS]
        for target in targets
    ])
    rng = np.random.default_rng(args.seed)
    bootstrap = np.empty((args.bootstrap, len(METRICS)), dtype=np.float64)
    for index in range(args.bootstrap):
        sampled = rng.integers(0, len(targets), size=len(targets))
        bootstrap[index] = differences[sampled].mean(axis=0)
    paired = {
        name: {
            "mean_delta": float(differences[:, metric_index].mean()),
            "ci95": [float(value) for value in np.quantile(bootstrap[:, metric_index], [0.025, 0.975])],
            "targets_improved": int((differences[:, metric_index] > 0).sum()),
            "targets_total": len(targets),
        }
        for metric_index, name in enumerate(METRICS)
    }
    report = {
        "protocol": "Mean-score ensemble across independent CGDA seeds; paired outer-target bootstrap",
        "runs": [str(path) for path in args.runs],
        "macro": macro,
        "paired": paired,
        "per_target": per_target,
        "claim_boundary": "Target-level uncertainty only; reference-ligand-assisted retrospective binding screen.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"macro": macro, "paired": paired}, indent=2))


if __name__ == "__main__":
    main()
