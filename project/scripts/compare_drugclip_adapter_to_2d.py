#!/usr/bin/env python3
"""Paired uncertainty analysis: DrugCLIP adapter versus frozen 2D OOF baseline."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


def summarize(values: list[float]) -> dict:
    arr = np.asarray(values, dtype=float)
    return {
        "mean": float(arr.mean()),
        "ci95": list(map(float, np.quantile(arr, [0.025, 0.975]))),
        "probability_delta_gt_zero": float((arr > 0).mean()),
        "n_bootstrap_valid": int(len(arr)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--adapter-predictions", type=Path, required=True)
    ap.add_argument("--baseline-predictions", type=Path, required=True)
    ap.add_argument("--benchmark", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--bootstrap", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=20260924)
    args = ap.parse_args()

    adapter = pd.read_csv(args.adapter_predictions)
    # Repeated training seeds are an ensemble at inference, not independent molecules.
    adapter = adapter.groupby(["canonical_molecule_id", "target"], as_index=False)[
        ["bce", "triplet", "triplet_domain"]
    ].mean()
    baseline = pd.read_csv(args.baseline_predictions)[
        ["canonical_molecule_id", "target", "prediction"]
    ].rename(columns={"target": "baseline_target", "prediction": "baseline"})
    benchmark = pd.read_csv(args.benchmark)[["canonical_molecule_id", "source_component"]]
    frame = adapter.merge(baseline, on="canonical_molecule_id", validate="one_to_one")
    frame = frame.merge(benchmark, on="canonical_molecule_id", validate="one_to_one")
    if not np.array_equal(frame.target.to_numpy(), frame.baseline_target.to_numpy()):
        raise ValueError("Target mismatch between adapter and baseline")
    y = frame.target.to_numpy(int)
    rng = np.random.default_rng(args.seed)
    models = ["bce", "triplet", "triplet_domain"]
    observed_baseline = float(roc_auc_score(y, frame.baseline))
    report = {
        "evidence_level": "retrospective_paired_oof_comparison",
        "n": int(len(frame)),
        "baseline_auc": observed_baseline,
        "models": {},
        "claim_boundary": (
            "Intervals quantify retrospective OOF uncertainty. Source-cluster bootstrap has only "
            f"{frame.source_component.nunique()} clusters and does not establish prospective transfer."
        ),
    }
    sources = frame.source_component.unique()
    for model in models:
        observed = float(roc_auc_score(y, frame[model]))
        molecule_delta, cluster_delta = [], []
        for _ in range(args.bootstrap):
            idx = rng.integers(0, len(frame), len(frame))
            if len(np.unique(y[idx])) == 2:
                molecule_delta.append(float(roc_auc_score(y[idx], frame[model].to_numpy()[idx]) -
                                            roc_auc_score(y[idx], frame.baseline.to_numpy()[idx])))
            selected = rng.choice(sources, len(sources), replace=True)
            pieces = [np.flatnonzero(frame.source_component.to_numpy() == source) for source in selected]
            cidx = np.concatenate(pieces)
            if len(np.unique(y[cidx])) == 2:
                cluster_delta.append(float(roc_auc_score(y[cidx], frame[model].to_numpy()[cidx]) -
                                           roc_auc_score(y[cidx], frame.baseline.to_numpy()[cidx])))
        report["models"][model] = {
            "adapter_auc": observed,
            "observed_delta_vs_2d": observed - observed_baseline,
            "paired_molecule_bootstrap_delta": summarize(molecule_delta),
            "source_cluster_bootstrap_delta": summarize(cluster_delta),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
