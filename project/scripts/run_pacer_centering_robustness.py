# -*- coding: utf-8 -*-
"""Hyperparameter-region robustness test for series-centered LightGBM."""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor

from run_pacer_fs_baselines import DATA, MIN_GROUP, choose_predicted_span, features, safe_rho


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_centering_robustness_v01"
OUT.mkdir(parents=True, exist_ok=True)
CONFIGS = list(itertools.product([7, 15, 31], [4, 6], [8, 20]))


def centered(y, g, train):
    z = np.zeros(len(y))
    for group in set(g[train]):
        idx = train[g[train] == group]
        z[idx] = y[idx] - y[idx].mean()
    return z


def model(leaves, depth, child):
    return LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=leaves, max_depth=depth,
        min_child_samples=child, reg_alpha=.1, reg_lambda=1,
        random_state=42, verbosity=-1, n_jobs=6,
    )


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, _ = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float); g = d.source_component.astype(str).to_numpy()
    groups = [x for x, n in d.source_component.value_counts().items() if n >= MIN_GROUP]
    rows = []
    for ci, (leaves, depth, child) in enumerate(CONFIGS):
        for group in groups:
            outer = np.where(g == group)[0]; train = np.where(g != group)[0]
            sel = np.argsort(X[train].var(0))[-1024:]; xx = np.c_[X[:, sel], D]
            pa = model(leaves, depth, child).fit(xx[train], y[train]).predict(xx[outer])
            z = centered(y, g, train)
            pc = model(leaves, depth, child).fit(xx[train], z[train]).predict(xx[outer])
            full = np.full(len(d), np.nan); full[outer] = pa
            support = choose_predicted_span(outer, 3, full)
            query = np.asarray([i for i in outer if i not in set(support)], int)
            loc = {idx: j for j, idx in enumerate(outer)}
            ql = np.asarray([loc[i] for i in query])
            rows.append({"config": ci, "leaves": leaves, "depth": depth, "child": child,
                         "group": group, "absolute_rho": safe_rho(y[query], pa[ql]),
                         "centered_rho": safe_rho(y[query], pc[ql])})
        print(f"config {ci + 1}/{len(CONFIGS)}", flush=True)
    out = pd.DataFrame(rows); out["delta_rho"] = out.centered_rho - out.absolute_rho
    out.to_csv(OUT / "per_series.csv", index=False)
    summary = out.groupby(["config", "leaves", "depth", "child"]).agg(
        absolute_macro=("absolute_rho", "mean"), centered_macro=("centered_rho", "mean"),
        delta_macro=("delta_rho", "mean"), improved_series=("delta_rho", lambda x: int((x > 0).sum())),
        centered_worst=("centered_rho", "min"),
    ).reset_index()
    summary.to_csv(OUT / "summary.csv", index=False)
    metrics = {
        "n_configs": len(summary),
        "configs_with_positive_macro_delta": int((summary.delta_macro > 0).sum()),
        "median_delta_macro": float(summary.delta_macro.median()),
        "min_delta_macro": float(summary.delta_macro.min()),
        "max_delta_macro": float(summary.delta_macro.max()),
        "median_centered_macro": float(summary.centered_macro.median()),
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(summary.to_string(index=False)); print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
