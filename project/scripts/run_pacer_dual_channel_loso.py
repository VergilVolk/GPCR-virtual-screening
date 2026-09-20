# -*- coding: utf-8 -*-
"""PACER-FS dual-channel LOSO: relative ranking plus mean-preserving calibration."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from sklearn.metrics import mean_absolute_error

from run_pacer_fs_baselines import (
    DATA, MIN_GROUP, adapted, choose_predicted_span, features, safe_rho, train_delta,
)


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_dual_channel_v01"
OUT.mkdir(parents=True, exist_ok=True)


def fit(xx, target, train, test):
    return LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1, random_state=42,
        verbosity=-1, n_jobs=6,
    ).fit(xx[train], target[train]).predict(xx[test])


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bv = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float); g = d.source_component.astype(str).to_numpy()
    groups = [x for x, n in d.source_component.value_counts().items() if n >= MIN_GROUP]
    rows = []
    for oi, group in enumerate(groups):
        target = np.where(g == group)[0]; train = np.where(g != group)[0]
        sel = np.argsort(X[train].var(0))[-1024:]; xx = np.c_[X[:, sel], D]
        absolute = np.full(len(d), np.nan); absolute[target] = fit(xx, y, train, target)
        z = np.zeros(len(d))
        for source in set(g[train]):
            idx = train[g[train] == source]; z[idx] = y[idx] - y[idx].mean()
        centered = np.full(len(d), np.nan); centered[target] = fit(xx, z, train, target)
        support = choose_predicted_span(target, 3, absolute)
        query = np.asarray([i for i in target if i not in set(support)], int)
        absolute[target] += np.mean(y[support] - absolute[support])
        centered[target] += np.mean(y[support] - centered[support])
        delta = train_delta(train, X, D, y, g, bv, 8000 + oi, cliff_weight=0)
        p_abs = adapted("DeltaSARHybrid", absolute, support, query, y, X, D, bv, delta)
        p_ctr = adapted("DeltaSARHybrid", centered, support, query, y, X, D, bv, delta)
        predictions = {
            "AbsoluteCalibrated": absolute[query],
            "CenteredCalibrated": centered[query],
            "DeltaSAR-Absolute": p_abs,
            "DeltaSAR-Centered": p_ctr,
            "DualChannel-Absolute": p_abs - np.mean(p_abs - absolute[query]),
            "DualChannel-Centered": p_ctr - np.mean(p_ctr - centered[query]),
        }
        for method, pred in predictions.items():
            rows.append({"group": group, "method": method, "n_query": len(query),
                         "Spearman": safe_rho(y[query], pred),
                         "MAE": float(mean_absolute_error(y[query], pred))})
        print(group, flush=True)
    out = pd.DataFrame(rows); out.to_csv(OUT / "outer_episodes.csv", index=False)
    summary = out.groupby("method").agg(
        macro_Spearman=("Spearman", "mean"), median_Spearman=("Spearman", "median"),
        worst_Spearman=("Spearman", "min"), positive_series=("Spearman", lambda x: int((x > 0).sum())),
        macro_MAE=("MAE", "mean"), n_series=("group", "size"),
    ).sort_values("macro_Spearman", ascending=False)
    summary.to_csv(OUT / "summary.csv")
    (OUT / "metrics.json").write_text(json.dumps(summary.reset_index().to_dict("records"), indent=2), encoding="utf-8")
    print(summary.to_string())


if __name__ == "__main__":
    main()
