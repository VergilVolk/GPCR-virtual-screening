#!/usr/bin/env python
"""Frozen synthetic benchmark for PACER-FKG estimator behavior.

This benchmark validates mathematics, not biological performance.  The
nonlinear arm changes only the shape of the C+A distribution while preserving
its expected mean, so a mean-difference baseline should not solve it.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pacer_factorial_kernel_graph import linear_norm, signed_kernel_stat


SIGNS = {"candidate_probe": 1, "probe_only": -1, "candidate_no_probe": -1, "apo": 1}


def make(seed, n, d, mode, effect):
    rng = np.random.default_rng(seed)
    groups = {name: rng.normal(size=(n, d)) for name in SIGNS}
    if mode == "mean_interaction":
        groups["candidate_probe"][:, :3] += effect
    elif mode == "shape_interaction":
        direction = rng.choice((-1.0, 1.0), size=(n, 1))
        groups["candidate_probe"][:, :3] += direction * effect
    elif mode != "null":
        raise ValueError(mode)
    return groups


def auc(negative, positive):
    wins = sum(p > n for p in positive for n in negative)
    ties = sum(p == n for p in positive for n in negative)
    return (wins + 0.5 * ties) / (len(positive) * len(negative))


def summarize(values):
    return {"mean": float(np.mean(values)), "sd": float(np.std(values)),
            "q05": float(np.quantile(values, 0.05)), "q95": float(np.quantile(values, 0.95))}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seeds", type=int, default=100)
    p.add_argument("--frames", type=int, default=80)
    p.add_argument("--dimensions", type=int, default=12)
    p.add_argument("--effect", type=float, default=1.5)
    args = p.parse_args()

    scores = {mode: {"kernel_u2": [], "linear": []}
              for mode in ("null", "mean_interaction", "shape_interaction")}
    for mode in scores:
        for seed in range(args.seeds):
            groups = make(seed, args.frames, args.dimensions, mode, args.effect)
            scores[mode]["kernel_u2"].append(signed_kernel_stat(groups, SIGNS)["unbiased_squared"])
            scores[mode]["linear"].append(linear_norm(groups, SIGNS))

    result = {
        "method": "PACER-FKG synthetic estimator benchmark v01",
        "evidence_level": "mathematical_synthetic_control",
        "settings": vars(args) | {"output": str(args.output)},
        "summary": {mode: {metric: summarize(v) for metric, v in metrics.items()}
                    for mode, metrics in scores.items()},
        "auc_signal_vs_null": {
            mode: {metric: auc(scores["null"][metric], scores[mode][metric])
                   for metric in scores[mode]}
            for mode in ("mean_interaction", "shape_interaction")
        },
        "claim_boundary": "Tests estimator behavior on generated distributions only; not evidence of M4 PAM prediction.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
