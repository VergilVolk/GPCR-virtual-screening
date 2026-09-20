# -*- coding: utf-8 -*-
"""PACER-FS v0.2 preregistered reliability-gated three-shot adaptation."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error

from run_pacer_fs_baselines import (
    DATA, MIN_GROUP, base_model, choose_predicted_span, features, pair_matrix,
    safe_rho, sims, train_delta,
)


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_fs_gated_v02"
OUT.mkdir(parents=True, exist_ok=True)
SHOT = 3
SERIES_TEMP = 0.15
LOCAL_SIGMA = 0.50
COVERAGE_SCALE = 0.50
MAX_BLEND = 0.60


def anchor_absolute(delta_pack, support, query, y, X, D, bv):
    model, sel, mean, std, _ = delta_pack
    aa = np.repeat(support, len(query))
    bb = np.tile(query, len(support))
    z = pair_matrix(aa, bb, X, D, bv, sel, mean, std)
    dp = model.predict(z).reshape(len(support), len(query))
    return y[support, None] + dp


def weighted_summary(values, similarities):
    w = np.maximum(similarities, 1e-4) ** 3
    mean = (values * w).sum(0) / (w.sum(0) + 1e-8)
    var = ((values - mean[None, :]) ** 2 * w).sum(0) / (w.sum(0) + 1e-8)
    return mean, np.sqrt(np.maximum(var, 0))


def support_reliability(delta_pack, support, base, y, X, D, bv):
    reconstructed = []
    for j in range(len(support)):
        source = np.delete(support, j)
        query = support[j:j + 1]
        vals = anchor_absolute(delta_pack, source, query, y, X, D, bv)
        s = sims(query, source, bv).T
        pred, _ = weighted_summary(vals, s)
        reconstructed.append(pred[0])
    delta_mae = float(mean_absolute_error(y[support], reconstructed))
    base_mae = float(mean_absolute_error(y[support], base[support]))
    reliability = float(1 / (1 + np.exp((delta_mae - base_mae) / SERIES_TEMP)))
    return reliability, base_mae, delta_mae


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bv = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float)
    g = d.source_component.astype(str).to_numpy()
    groups = [x for x, n in d.source_component.value_counts().items() if n >= MIN_GROUP]
    rows = []

    for oi, group in enumerate(groups):
        target = np.where(g == group)[0]
        train = np.where(g != group)[0]
        base = np.full(len(d), np.nan)
        base[target] = base_model(X, D, y, train, target)
        delta = train_delta(train, X, D, y, g, bv, 5000 + oi, cliff_weight=0)
        support = choose_predicted_span(target, SHOT, base)
        query = np.asarray([i for i in target if i not in set(support)], dtype=int)

        vals = anchor_absolute(delta, support, query, y, X, D, bv)
        similarity = sims(query, support, bv).T
        delta_pred, disagreement = weighted_summary(vals, similarity)
        r_series, support_base_mae, support_delta_mae = support_reliability(
            delta, support, base, y, X, D, bv
        )
        coverage = similarity.max(0)
        r_local = np.exp(-(disagreement / LOCAL_SIGMA) ** 2) * np.clip(
            coverage / COVERAGE_SCALE, 0, 1
        )
        blend = MAX_BLEND * r_series * r_local
        pred = base[query] + blend * (delta_pred - base[query])

        rows.append({
            "group": group, "n_group": len(target), "n_query": len(query),
            "Spearman": safe_rho(y[query], pred),
            "Base_Spearman": safe_rho(y[query], base[query]),
            "MAE": float(mean_absolute_error(y[query], pred)),
            "Base_MAE": float(mean_absolute_error(y[query], base[query])),
            "series_reliability": r_series,
            "support_base_MAE": support_base_mae,
            "support_delta_MAE": support_delta_mae,
            "mean_local_reliability": float(r_local.mean()),
            "mean_blend": float(blend.mean()),
            "support_ids": "|".join(d.loc[support, "canonical_molecule_id"].astype(str)),
        })
        print(group, f"r={r_series:.3f} blend={blend.mean():.3f}", flush=True)

    out = pd.DataFrame(rows)
    out["delta_Spearman"] = out.Spearman - out.Base_Spearman
    out["delta_MAE"] = out.MAE - out.Base_MAE
    out.to_csv(OUT / "outer_episodes.csv", index=False)

    rng = np.random.default_rng(20260829)
    draw = rng.integers(0, len(out), size=(50_000, len(out)))
    dr = out.delta_Spearman.to_numpy()
    dm = out.delta_MAE.to_numpy()
    metrics = {
        "protocol": "preregistered_fixed_3shot_loso",
        "n_series": len(out),
        "macro_Spearman": float(out.Spearman.mean()),
        "base_macro_Spearman_same_queries": float(out.Base_Spearman.mean()),
        "delta_macro_Spearman": float(dr.mean()),
        "delta_Spearman_ci95": [float(x) for x in np.quantile(dr[draw].mean(1), [.025, .975])],
        "probability_delta_Spearman_gt_0": float((dr[draw].mean(1) > 0).mean()),
        "worst_Spearman": float(out.Spearman.min()),
        "base_worst_Spearman_same_queries": float(out.Base_Spearman.min()),
        "positive_series": int((out.Spearman > 0).sum()),
        "macro_MAE": float(out.MAE.mean()),
        "base_macro_MAE_same_queries": float(out.Base_MAE.mean()),
        "delta_macro_MAE": float(dm.mean()),
        "delta_MAE_ci95": [float(x) for x in np.quantile(dm[draw].mean(1), [.025, .975])],
        "success_primary": bool(np.quantile(dr[draw].mean(1), .025) > 0),
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
