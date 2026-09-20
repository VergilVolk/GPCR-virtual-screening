# -*- coding: utf-8 -*-
"""Evaluate whether predicted-span anchors improve PACER-FS calibration.

Ranking is intentionally not the endpoint because a constant series offset
cannot change ranks.  The primary endpoint is error in the estimated series
offset; query MAE is secondary.  Policies see no outer labels.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from sklearn.metrics import mean_absolute_error

from run_pacer_fs_baselines import (
    DATA, MIN_GROUP, choose_diverse, choose_facility, choose_predicted_span,
    features,
)


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_anchor_policy_v01"
OUT.mkdir(parents=True, exist_ok=True)
N_RANDOM = 200


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bv = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float)
    g = d.source_component.astype(str).to_numpy()
    groups = [x for x, n in d.source_component.value_counts().items() if n >= MIN_GROUP]
    rows = []
    rng = np.random.default_rng(20260829)

    for oi, group in enumerate(groups):
        target = np.where(g == group)[0]
        train = np.where(g != group)[0]
        centered = np.zeros(len(d), float)
        for source in set(g[train]):
            idx = train[g[train] == source]
            centered[idx] = y[idx] - y[idx].mean()
        sel = np.argsort(X[train].var(0))[-1024:]
        xx = np.c_[X[:, sel], D]
        score = np.full(len(d), np.nan)
        score[target] = LGBMRegressor(
            n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
            min_child_samples=12, reg_alpha=.1, reg_lambda=1,
            random_state=42, verbosity=-1, n_jobs=6,
        ).fit(xx[train], centered[train]).predict(xx[target])
        oracle_offset = float(np.mean(y[target] - score[target]))

        policies = {
            "diverse": [choose_diverse(target, 3, bv)],
            "facility": [choose_facility(target, 3, bv)],
            "joint_facility": [choose_facility(target, 3, bv, score)],
            "predicted_span": [choose_predicted_span(target, 3, score)],
            "random": [rng.choice(target, 3, replace=False) for _ in range(N_RANDOM)],
        }
        for policy, supports in policies.items():
            for rep, support in enumerate(supports):
                query = np.asarray([i for i in target if i not in set(support)], int)
                offset = float(np.mean(y[support] - score[support]))
                pred = score[query] + offset
                rows.append({
                    "group": group, "policy": policy, "repeat": rep,
                    "offset_abs_error": abs(offset - oracle_offset),
                    "query_MAE": float(mean_absolute_error(y[query], pred)),
                    "anchor_true_range": float(np.ptp(y[support])),
                    "support_ids": "|".join(d.loc[support, "canonical_molecule_id"].astype(str)),
                })
        print(group, flush=True)

    raw = pd.DataFrame(rows)
    raw.to_csv(OUT / "episodes.csv", index=False)
    per = raw.groupby(["policy", "group"])[
        ["offset_abs_error", "query_MAE", "anchor_true_range"]
    ].mean().reset_index()
    summary = per.groupby("policy").agg(
        macro_offset_abs_error=("offset_abs_error", "mean"),
        macro_query_MAE=("query_MAE", "mean"),
        macro_anchor_true_range=("anchor_true_range", "mean"),
        n_series=("group", "size"),
    ).sort_values("macro_offset_abs_error")
    summary.to_csv(OUT / "summary.csv")

    pivot = per.pivot(index="group", columns="policy", values="offset_abs_error")
    delta = pivot.predicted_span - pivot.random
    boot = []
    brng = np.random.default_rng(829)
    vals = delta.to_numpy(float)
    for _ in range(50000):
        boot.append(float(np.mean(brng.choice(vals, len(vals), replace=True))))
    payload = {
        "primary_endpoint": "series_offset_absolute_error",
        "predicted_span_minus_random": float(delta.mean()),
        "ci95": [float(np.quantile(boot, .025)), float(np.quantile(boot, .975))],
        "probability_predicted_span_better": float(np.mean(np.asarray(boot) < 0)),
        "interpretation": (
            "Negative delta favors predicted-span. If CI crosses zero, anchor selection remains a heuristic."
        ),
    }
    (OUT / "audit.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(summary.to_string())
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
