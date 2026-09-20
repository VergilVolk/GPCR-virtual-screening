# -*- coding: utf-8 -*-
"""Three-anchor adaptation on the single retrospective 2022+ M4 PAM series."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from sklearn.metrics import mean_absolute_error

from run_pacer_fs_baselines import (
    adapted, choose_predicted_span, features, safe_rho, train_delta,
)


P = Path(__file__).resolve().parents[1]
DATA = P / "data" / "benchmarks" / "m4_pam_v1" / "temporal" / "potency_molecules_temporal.csv"
OUT = P / "results" / "pacer_temporal_fewshot_v01"
OUT.mkdir(parents=True, exist_ok=True)


def fit_predict(xx, target, train, test):
    m = LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1, random_state=42,
        verbosity=-1, n_jobs=6,
    )
    return m.fit(xx[train], target[train]).predict(xx[test])


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bv = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float); g = d.source_component.astype(str).to_numpy()
    future = d.temporal_split_enriched.eq("test_ge_2022").to_numpy()
    train = np.where(~future)[0]
    group = "SCOMP0002"; target = np.where(future & (g == group))[0]
    sel = np.argsort(X[train].var(0))[-1024:]; xx = np.c_[X[:, sel], D]

    absolute = np.full(len(d), np.nan)
    absolute[target] = fit_predict(xx, y, train, target)
    z = np.zeros(len(d))
    for source in set(g[train]):
        idx = train[g[train] == source]; z[idx] = y[idx] - y[idx].mean()
    centered = np.full(len(d), np.nan)
    centered[target] = fit_predict(xx, z, train, target)

    support = choose_predicted_span(target, 3, absolute)
    query = np.asarray([i for i in target if i not in set(support)], int)
    # Restore an absolute scale without changing centered ranking.
    centered[target] += np.mean(y[support] - centered[support])
    absolute[target] += np.mean(y[support] - absolute[support])
    delta = train_delta(train, X, D, y, g, bv, 7000, cliff_weight=0)

    delta_abs = adapted("DeltaSARHybrid", absolute, support, query, y, X, D, bv, delta)
    delta_ctr = adapted("DeltaSARHybrid", centered, support, query, y, X, D, bv, delta)
    # ΔSAR is permitted to change ordering but not the anchor-calibrated batch
    # location.  This separates relative SAR from absolute potency calibration.
    delta_abs_mp = delta_abs - np.mean(delta_abs - absolute[query])
    delta_ctr_mp = delta_ctr - np.mean(delta_ctr - centered[query])
    predictions = {
        "Absolute": absolute[query],
        "Centered": centered[query],
        "KernelResidual-Absolute": adapted("KernelResidual", absolute, support, query, y, X, D, bv, delta),
        "TanimotoGP-Absolute": adapted("TanimotoGP", absolute, support, query, y, X, D, bv, delta),
        "DeltaSARHybrid-Absolute": delta_abs,
        "DeltaSARHybrid-Centered": delta_ctr,
        "DualChannel-Absolute": delta_abs_mp,
        "DualChannel-Centered": delta_ctr_mp,
    }
    rows = []
    for method, pred in predictions.items():
        rows.append({"method": method, "group": group, "n_query": len(query),
                     "Spearman": safe_rho(y[query], pred),
                     "MAE": float(mean_absolute_error(y[query], pred))})
    out = pd.DataFrame(rows).sort_values("Spearman", ascending=False)
    out.to_csv(OUT / "metrics.csv", index=False)
    payload = {"status": "retrospective_single_future_series_exploratory",
               "support_ids": d.loc[support, "canonical_molecule_id"].astype(str).tolist(),
               "rows": out.to_dict("records")}
    (OUT / "metrics.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
