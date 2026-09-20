# -*- coding: utf-8 -*-
"""Test whether assay/series normalization improves unseen-series SAR ranking.

This is a diagnostic ablation, not a final benchmark.  Every outer series is
excluded from fitting.  Three predicted-span anchors are used only to restore
the absolute pEC50 offset; they cannot alter the reported ranking.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRanker, LGBMRegressor
from scipy.stats import rankdata
from sklearn.metrics import mean_absolute_error

from run_pacer_fs_baselines import DATA, DESC, MIN_GROUP, choose_predicted_span, features, safe_rho


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_series_normalized_v01"
OUT.mkdir(parents=True, exist_ok=True)


def matrix(X, D, train):
    sel = np.argsort(X[train].var(0))[-1024:]
    return np.c_[X[:, sel], D]


def regressor(seed=42):
    return LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1, random_state=seed,
        verbosity=-1, n_jobs=6,
    )


def group_targets(y, g, train, mode):
    target = np.zeros(len(y), float)
    for group in sorted(set(g[train])):
        idx = train[g[train] == group]
        yy = y[idx]
        if mode == "center":
            target[idx] = yy - yy.mean()
        elif mode == "zscore":
            target[idx] = (yy - yy.mean()) / max(yy.std(), .15)
        elif mode == "percentile":
            target[idx] = (rankdata(yy, method="average") - .5) / len(yy)
        else:
            raise ValueError(mode)
    return target


def offset_calibrate(score, support, query, y):
    offset = float(np.mean(y[support] - score[support]))
    return score[query] + offset


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, _ = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float)
    g = d.source_component.astype(str).to_numpy()
    groups = [x for x, n in d.source_component.value_counts().items() if n >= MIN_GROUP]
    rows = []

    for oi, outer_group in enumerate(groups):
        outer = np.where(g == outer_group)[0]
        train = np.where(g != outer_group)[0]
        xx = matrix(X, D, train)
        scores = {}

        scores["Absolute"] = regressor(10 + oi).fit(xx[train], y[train]).predict(xx[outer])
        for mode in ["center", "zscore", "percentile"]:
            t = group_targets(y, g, train, mode)
            scores[mode.capitalize()] = regressor(100 + oi).fit(xx[train], t[train]).predict(xx[outer])

        # Ranking objective: each source series is one query group.  Five
        # within-series potency relevance levels remove assay-scale offsets.
        order = np.argsort(g[train], kind="stable")
        sorted_train = train[order]
        rel = np.zeros(len(sorted_train), dtype=int)
        group_sizes = []
        pos = 0
        for source in sorted(set(g[train])):
            idx = sorted_train[g[sorted_train] == source]
            vals = y[idx]
            pct = (rankdata(vals, method="average") - .5) / len(vals)
            rel[pos:pos + len(idx)] = np.minimum((pct * 5).astype(int), 4)
            group_sizes.append(len(idx))
            pos += len(idx)
        ranker = LGBMRanker(
            objective="lambdarank", metric="ndcg", label_gain=[0, 1, 3, 7, 15],
            n_estimators=300, learning_rate=.03, num_leaves=15, max_depth=5,
            min_child_samples=12, reg_alpha=.1, reg_lambda=1,
            random_state=500 + oi, verbosity=-1, n_jobs=6,
        )
        ranker.fit(xx[sorted_train], rel, group=group_sizes)
        scores["LambdaRank"] = ranker.predict(xx[outer])

        # One common support/query split selected from the absolute zero-shot
        # model; no labels are used in anchor selection.
        base_full = np.full(len(d), np.nan)
        base_full[outer] = scores["Absolute"]
        support = choose_predicted_span(outer, 3, base_full)
        query = np.asarray([i for i in outer if i not in set(support)], dtype=int)
        local = {idx: j for j, idx in enumerate(outer)}
        s_local = np.asarray([local[i] for i in support])
        q_local = np.asarray([local[i] for i in query])

        for method, score in scores.items():
            pred = score[q_local] if method == "Absolute" else offset_calibrate(
                score, s_local, q_local, y[outer]
            )
            # Absolute already has a global scale, but apply the same permitted
            # offset for a fair absolute-error comparison.
            if method == "Absolute":
                pred = offset_calibrate(score, s_local, q_local, y[outer])
            rows.append({
                "group": outer_group, "method": method, "n_query": len(query),
                "Spearman": safe_rho(y[query], pred),
                "MAE": float(mean_absolute_error(y[query], pred)),
            })
        print(outer_group, flush=True)

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
