# -*- coding: utf-8 -*-
"""Strong model-family controls for the PACER series-centering hypothesis."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import rankdata
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from run_pacer_fs_baselines import DATA, MIN_GROUP, choose_predicted_span, features, safe_rho


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_centered_suite_v01"
OUT.mkdir(parents=True, exist_ok=True)


def models(seed):
    return {
        "LightGBM": LGBMRegressor(
            n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
            min_child_samples=12, reg_alpha=.1, reg_lambda=1,
            random_state=seed, verbosity=-1, n_jobs=6,
        ),
        "RandomForest": RandomForestRegressor(
            n_estimators=500, min_samples_leaf=2, max_features=.35,
            random_state=seed, n_jobs=6,
        ),
        "ExtraTrees": ExtraTreesRegressor(
            n_estimators=500, min_samples_leaf=2, max_features=.35,
            random_state=seed, n_jobs=6,
        ),
        "Ridge": make_pipeline(StandardScaler(), Ridge(alpha=20.0)),
    }


def centered_labels(y, g, train):
    z = np.zeros(len(y), float)
    for group in set(g[train]):
        idx = train[g[train] == group]
        z[idx] = y[idx] - y[idx].mean()
    return z


def mean_offset(score, support_local, query_local, y_outer):
    return score[query_local] + np.mean(y_outer[support_local] - score[support_local])


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, _ = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float)
    g = d.source_component.astype(str).to_numpy()
    groups = [x for x, n in d.source_component.value_counts().items() if n >= MIN_GROUP]
    rows = []

    for oi, group in enumerate(groups):
        outer = np.where(g == group)[0]
        train = np.where(g != group)[0]
        sel = np.argsort(X[train].var(0))[-1024:]
        xx = np.c_[X[:, sel], D]
        z = centered_labels(y, g, train)

        # Common label-free support split from the fixed absolute LightGBM.
        absolute = models(42)["LightGBM"].fit(xx[train], y[train]).predict(xx[outer])
        full = np.full(len(d), np.nan); full[outer] = absolute
        support = choose_predicted_span(outer, 3, full)
        query = np.asarray([i for i in outer if i not in set(support)], int)
        loc = {idx: j for j, idx in enumerate(outer)}
        sl = np.asarray([loc[i] for i in support]); ql = np.asarray([loc[i] for i in query])

        centered_scores = {}
        for name, model in models(100 + oi).items():
            score = model.fit(xx[train], z[train]).predict(xx[outer])
            centered_scores[name] = score
            pred = mean_offset(score, sl, ql, y[outer])
            rows.append({"group": group, "method": f"Centered-{name}",
                         "Spearman": safe_rho(y[query], pred),
                         "MAE": float(mean_absolute_error(y[query], pred))})

        for ensemble_name, members in {
            "Centered-EnsembleTrees": ["LightGBM", "RandomForest", "ExtraTrees"],
            "Centered-EnsembleAll": list(centered_scores),
        }.items():
            # Rank averaging is scale invariant across heterogeneous learners.
            score = np.mean([
                rankdata(centered_scores[m], method="average") / len(outer) for m in members
            ], axis=0)
            pred = mean_offset(score, sl, ql, y[outer])
            rows.append({"group": group, "method": ensemble_name,
                         "Spearman": safe_rho(y[query], pred),
                         "MAE": float(mean_absolute_error(y[query], pred))})
        print(group, flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "outer_episodes.csv", index=False)
    summary = out.groupby("method").agg(
        macro_Spearman=("Spearman", "mean"), median_Spearman=("Spearman", "median"),
        worst_Spearman=("Spearman", "min"), positive_series=("Spearman", lambda x: int((x > 0).sum())),
        macro_MAE=("MAE", "mean"), n_series=("Spearman", "size"),
    ).sort_values("macro_Spearman", ascending=False)
    summary.to_csv(OUT / "summary.csv")
    (OUT / "metrics.json").write_text(json.dumps(summary.reset_index().to_dict("records"), indent=2), encoding="utf-8")
    print(summary.to_string())


if __name__ == "__main__":
    main()
