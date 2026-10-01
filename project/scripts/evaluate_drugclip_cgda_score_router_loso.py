#!/usr/bin/env python3
"""Strict target-LOO, low-capacity router between frozen pocket DrugCLIP and CGDA.

The held target contributes neither labels nor router fitting.  Scores are converted
to within-target percentile ranks so that a global router cannot exploit target
score scale.  This is an exploratory, preregistered-capacity test, not a tuned
stacking sweep.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression

from evaluate_drugclip_cgda_external import METRIC_NAMES, bootstrap
from evaluate_drugclip_official_litpcba_embeddings import metrics


def percentile(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=np.float32)
    ranks[order] = np.arange(len(values), dtype=np.float32)
    return ranks / max(1, len(values) - 1)


def features(pocket: np.ndarray, cgda: np.ndarray) -> np.ndarray:
    p, c = percentile(pocket), percentile(cgda)
    return np.column_stack([p, c, p - c, p * c, np.abs(p - c)]).astype(np.float32)


def balanced_indices(labels: np.ndarray, maximum_negatives: int, rng: np.random.Generator) -> np.ndarray:
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    count = min(len(negative), maximum_negatives, max(len(positive) * 20, 1000))
    return np.concatenate([positive, rng.choice(negative, size=count, replace=False)])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-negatives-per-target", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--bootstrap", type=int, default=20000)
    args = parser.parse_args()

    directories = [path.with_suffix("") for path in args.runs]
    targets = sorted(path.stem for path in directories[0].glob("*.npz"))
    data = {}
    for target in targets:
        records = [np.load(directory / f"{target}.npz") for directory in directories]
        labels = records[0]["labels"].astype(np.int64)
        data[target] = {"labels": labels, "pocket": records[0]["pocket"].astype(np.float32),
                        "cgda": np.mean([record["cgda"] for record in records], axis=0).astype(np.float32)}
        data[target]["features"] = features(data[target]["pocket"], data[target]["cgda"])

    folds = {}
    for fold_index, held in enumerate(targets):
        rng = np.random.default_rng(args.seed + fold_index)
        train_x, train_y = [], []
        for target in targets:
            if target == held:
                continue
            idx = balanced_indices(data[target]["labels"], args.max_negatives_per_target, rng)
            train_x.append(data[target]["features"][idx])
            train_y.append(data[target]["labels"][idx])
        # Fixed L2 logistic router: five score-derived features, no target ID.
        model = LogisticRegression(C=1.0, penalty="l2", solver="lbfgs", max_iter=300,
                                   class_weight="balanced", random_state=args.seed + fold_index)
        model.fit(np.concatenate(train_x), np.concatenate(train_y))
        routed = model.predict_proba(data[held]["features"])[:, 1]
        labels = data[held]["labels"]
        folds[held] = {"n": int(len(labels)), "positives": int(labels.sum()),
                       "pocket": metrics(labels, data[held]["pocket"]),
                       "cgda": metrics(labels, data[held]["cgda"]),
                       "router": metrics(labels, routed),
                       "coefficients": model.coef_[0].tolist(), "intercept": float(model.intercept_[0])}
        print(f"completed {fold_index + 1}/{len(targets)} {held}", flush=True)

    macro = {method: {name: float(np.mean([folds[t][method][name] for t in targets]))
                      for name in METRIC_NAMES} for method in ("pocket", "cgda", "router")}
    report = {"protocol": "Outer target-LOO frozen-score router; held-target labels absent from fitting",
              "features": ["pocket_percentile", "cgda_percentile", "difference", "product", "absolute_difference"],
              "model": "fixed L2 logistic regression C=1.0; no target identifiers; within-target score percentiles",
              "macro": macro,
              "paired_router_vs_pocket": bootstrap({t: {"reference": folds[t]["pocket"], "cgda": folds[t]["router"]}
                                                       for t in targets}, args.bootstrap, args.seed),
              "paired_router_vs_cgda": bootstrap({t: {"reference": folds[t]["cgda"], "cgda": folds[t]["router"]}
                                                     for t in targets}, args.bootstrap, args.seed + 1),
              "folds": folds,
              "parameters": vars(args),
              "claim_boundary": "Exploratory score-routing benchmark; requires independent external confirmation before a method claim."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(json.dumps({"macro": macro, "vs_pocket": report["paired_router_vs_pocket"],
                      "vs_cgda": report["paired_router_vs_cgda"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
