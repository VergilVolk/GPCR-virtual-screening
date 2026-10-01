#!/usr/bin/env python3
"""Paired target-bootstrap comparison of PACER projection LoRA and ECFP."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


METRICS = ("roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pacer", type=Path, required=True)
    parser.add_argument("--ecfp", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20260928)
    args = parser.parse_args()

    pacer = json.loads(args.pacer.read_text(encoding="utf-8"))["metrics"]["pacer_loto"]["per_target"]
    ecfp = json.loads(args.ecfp.read_text(encoding="utf-8"))["per_target"]
    targets = sorted(set(pacer) & set(ecfp))
    if not targets:
        raise RuntimeError("No common targets")

    rng = np.random.default_rng(args.seed)
    draws = rng.integers(0, len(targets), size=(args.bootstrap, len(targets)))
    comparisons = {}
    for metric in METRICS:
        delta = np.asarray([pacer[target][metric] - ecfp[target][metric] for target in targets])
        sampled = delta[draws].mean(axis=1)
        comparisons[metric] = {
            "pacer_macro": float(np.mean([pacer[target][metric] for target in targets])),
            "ecfp_macro": float(np.mean([ecfp[target][metric] for target in targets])),
            "delta": float(delta.mean()),
            "target_bootstrap_95ci": list(map(float, np.quantile(sampled, [0.025, 0.5, 0.975]))),
            "improved_targets": int((delta > 0).sum()),
            "tied_targets": int((delta == 0).sum()),
            "degraded_targets": int((delta < 0).sum()),
        }
    report = {
        "comparison": "PACER Science-2026 dual-projection triplet LoRA minus ECFP4 logistic",
        "targets": targets,
        "n_targets": len(targets),
        "bootstrap": args.bootstrap,
        "metrics": comparisons,
        "claim_boundary": (
            "Nine-target development comparison. A positive point estimate is not a stable advantage "
            "unless the paired target-bootstrap interval is above zero."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
