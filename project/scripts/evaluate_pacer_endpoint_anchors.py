# -*- coding: utf-8 -*-
"""Evaluate endpoint-spanning anchors before any later-generation external test.

The existing predicted-span policy samples the 25/50/75% predicted quantiles.
This diagnostic compares a frozen min/median/max policy that explicitly covers
the predicted potency range.  All comparisons remain outer-series held out.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_pacer_fs_baselines import (
    DATA, MIN_GROUP, adapted, base_model, choose_predicted_span,
    features, safe_rho, train_delta,
)


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_endpoint_anchor_v01"


def endpoints(target, base):
    ordered = target[np.argsort(base[target])]
    return np.asarray([ordered[0], ordered[len(ordered) // 2], ordered[-1]], int)


def bootstrap(delta, seed=20260830, n=100_000):
    values = np.asarray(delta, float)
    rng = np.random.default_rng(seed)
    sampled = values[rng.integers(0, len(values), size=(n, len(values)))].mean(1)
    return {"estimate": float(values.mean()),
            "ci95": [float(x) for x in np.quantile(sampled, [.025, .975])],
            "probability_gt_zero": float(np.mean(sampled > 0))}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bvs = features(data.canonical_smiles)
    y = data.pEC50.to_numpy(float); groups = data.source_component.astype(str).to_numpy()
    evaluated = [g for g, n in data.source_component.value_counts().items() if n >= MIN_GROUP]
    rows = []
    for oi, held in enumerate(evaluated):
        target = np.where(groups == held)[0]; train = np.where(groups != held)[0]
        base = np.full(len(data), np.nan)
        base[target] = base_model(X, D, y, train, target)
        delta = train_delta(train, X, D, y, groups, bvs, seed=9100 + oi)
        policies = {
            "quartile_span_25_50_75": choose_predicted_span(target, 3, base),
            "endpoint_span_min_med_max": endpoints(target, base),
        }
        for policy, support in policies.items():
            support_set = set(map(int, support))
            query = np.asarray([i for i in target if int(i) not in support_set], int)
            pred = adapted("DeltaSARHybrid", base, support, query, y, X, D, bvs, delta)
            rows.append({
                "group": held, "policy": policy, "n_group": len(target), "n_query": len(query),
                "Spearman": safe_rho(y[query], pred),
                "support_ids": "|".join(data.loc[support, "canonical_molecule_id"].astype(str)),
                "support_true_range": float(np.ptp(y[support])),
                "support_predicted_range": float(np.ptp(base[support])),
            })
        print(held, flush=True)
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "series_results.csv", index=False)
    summary = result.groupby("policy").agg(
        macro_Spearman=("Spearman", "mean"), median_Spearman=("Spearman", "median"),
        worst_Spearman=("Spearman", "min"), positive_series=("Spearman", lambda x: int((x > 0).sum())),
        mean_anchor_true_range=("support_true_range", "mean"), n_series=("group", "size"),
    ).reset_index()
    summary.to_csv(OUT / "summary.csv", index=False)
    pivot = result.pivot(index="group", columns="policy", values="Spearman")
    comparison = bootstrap(
        pivot["endpoint_span_min_med_max"] - pivot["quartile_span_25_50_75"]
    )
    audit = {
        "protocol": "11 outer-series-held-out functional SAR episodes; fixed 3 anchors",
        "comparison": "endpoint min/median/max minus prior 25/50/75 predicted-span",
        "paired_series_bootstrap": comparison,
        "promote_for_future_holdout": bool(comparison["ci95"][0] > 0),
        "external_2026_table1_labels_used": False,
        "claim_boundary": "Internal anchor-policy evaluation only.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(summary.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
