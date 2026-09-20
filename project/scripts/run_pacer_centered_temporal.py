# -*- coding: utf-8 -*-
"""Retrospective 2022+ test of the frozen series-centered PACER-SAR model."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from sklearn.metrics import mean_absolute_error

from run_pacer_fs_baselines import choose_predicted_span, features, safe_rho


P = Path(__file__).resolve().parents[1]
DATA = P / "data" / "benchmarks" / "m4_pam_v1" / "temporal" / "potency_molecules_temporal.csv"
OUT = P / "results" / "pacer_centered_temporal_v01"
OUT.mkdir(parents=True, exist_ok=True)


def model():
    return LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1, random_state=42,
        verbosity=-1, n_jobs=6,
    )


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, _ = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float); g = d.source_component.astype(str).to_numpy()
    test_mask = d.temporal_split_enriched.eq("test_ge_2022").to_numpy()
    train = np.where(~test_mask)[0]
    # The only evaluable 2022+ medicinal-chemistry series.
    test_groups = [x for x in sorted(set(g[test_mask])) if np.sum(test_mask & (g == x)) >= 8]
    rows = []
    for group in test_groups:
        test = np.where(test_mask & (g == group))[0]
        sel = np.argsort(X[train].var(0))[-1024:]; xx = np.c_[X[:, sel], D]
        z = np.zeros(len(d))
        for source in set(g[train]):
            idx = train[g[train] == source]; z[idx] = y[idx] - y[idx].mean()
        pa = model().fit(xx[train], y[train]).predict(xx[test])
        pc = model().fit(xx[train], z[train]).predict(xx[test])
        full = np.full(len(d), np.nan); full[test] = pa
        support = choose_predicted_span(test, 3, full)
        query = np.asarray([i for i in test if i not in set(support)], int)
        loc = {idx: j for j, idx in enumerate(test)}; ql = np.asarray([loc[i] for i in query])
        for name, pred in [("Absolute", pa[ql]), ("Centered", pc[ql])]:
            # Offset from three anchors for absolute pEC50; ranking is unchanged.
            sl = np.asarray([loc[i] for i in support])
            calibrated = pred + np.mean(y[support] - (pa if name == "Absolute" else pc)[sl])
            rows.append({"group": group, "method": name, "n_query": len(query),
                         "Spearman": safe_rho(y[query], calibrated),
                         "MAE": float(mean_absolute_error(y[query], calibrated))})
    out = pd.DataFrame(rows); out.to_csv(OUT / "metrics.csv", index=False)
    payload = {"status": "retrospective_single_series_not_independent_prospective_test",
               "rows": out.to_dict("records")}
    (OUT / "metrics.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
