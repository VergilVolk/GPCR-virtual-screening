# -*- coding: utf-8 -*-
"""Strict nested evaluation of PACER-FS at a fixed three-shot budget.

For each outer series, the complete method configuration is selected using
leave-one-series-out episodes among the remaining development series.  The
outer labels are never used for model fitting or configuration selection;
only the three selected anchor labels are revealed at adaptation time.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error

from run_pacer_fs_baselines import (
    DATA, MIN_GROUP, adapted, base_model, choose_diverse, choose_facility,
    choose_predicted_span, features, safe_rho, train_delta,
)


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_fs_nested_v01"
OUT.mkdir(parents=True, exist_ok=True)
SHOT = 3
POLICIES = ("diverse", "facility", "joint_facility", "predicted_span")
CONFIGS = (
    ("Base", "predicted_span", 0.0),
    ("KernelResidual", "predicted_span", 0.0),
    ("TanimotoGP", "predicted_span", 0.0),
    ("DeltaSARHybrid", "diverse", 0.0),
    ("DeltaSARHybrid", "facility", 0.0),
    ("DeltaSARHybrid", "joint_facility", 0.0),
    ("DeltaSARHybrid", "predicted_span", 0.0),
    ("DeltaSARHybridCliff", "predicted_span", 3.0),
)


def support_for(policy, target, base, bv):
    if policy == "diverse":
        return choose_diverse(target, SHOT, bv)
    if policy == "facility":
        return choose_facility(target, SHOT, bv)
    if policy == "joint_facility":
        return choose_facility(target, SHOT, bv, base)
    if policy == "predicted_span":
        return choose_predicted_span(target, SHOT, base)
    raise ValueError(policy)


def evaluate_episode(target, base, y, X, D, bv, delta_plain, delta_cliff, method, policy):
    support = support_for(policy, target, base, bv)
    query = np.asarray([i for i in target if i not in set(support)], dtype=int)
    if method == "Base":
        pred = base[query]
    else:
        pack = delta_cliff if method == "DeltaSARHybridCliff" else delta_plain
        adapted_name = "DeltaSARHybrid" if method == "DeltaSARHybridCliff" else method
        pred = adapted(adapted_name, base, support, query, y, X, D, bv, pack)
    return {
        "Spearman": safe_rho(y[query], pred),
        "MAE": float(mean_absolute_error(y[query], pred)),
        "support": support,
        "query": query,
    }


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bv = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float)
    g = d.source_component.astype(str).to_numpy()
    groups = [x for x, n in d.source_component.value_counts().items() if n >= MIN_GROUP]
    outer_rows, selection_rows = [], []

    for oi, outer_group in enumerate(groups):
        outer = np.where(g == outer_group)[0]
        inner_groups = [x for x in groups if x != outer_group]
        inner_scores = {c: [] for c in CONFIGS}

        # True inner validation: both the outer series and current inner series
        # are absent from every fitted model.
        for ii, inner_group in enumerate(inner_groups):
            inner = np.where(g == inner_group)[0]
            train = np.where((g != outer_group) & (g != inner_group))[0]
            base = np.full(len(d), np.nan)
            base[inner] = base_model(X, D, y, train, inner)
            plain = train_delta(train, X, D, y, g, bv, 1000 + oi * 100 + ii, cliff_weight=0)
            cliff = train_delta(train, X, D, y, g, bv, 2000 + oi * 100 + ii, cliff_weight=3)
            for config in CONFIGS:
                method, policy, _ = config
                z = evaluate_episode(inner, base, y, X, D, bv, plain, cliff, method, policy)
                inner_scores[config].append((z["Spearman"], z["MAE"]))

        ranked = []
        for config, values in inner_scores.items():
            rho = float(np.mean([z[0] for z in values]))
            mae = float(np.mean([z[1] for z in values]))
            ranked.append((rho, -mae, config))
            selection_rows.append({
                "outer_group": outer_group, "method": config[0], "policy": config[1],
                "cliff_weight": config[2], "inner_macro_Spearman": rho, "inner_macro_MAE": mae,
            })
        ranked.sort(reverse=True)
        selected = ranked[0][2]

        train = np.where(g != outer_group)[0]
        base = np.full(len(d), np.nan)
        base[outer] = base_model(X, D, y, train, outer)
        plain = train_delta(train, X, D, y, g, bv, 3000 + oi, cliff_weight=0)
        cliff = train_delta(train, X, D, y, g, bv, 4000 + oi, cliff_weight=3)
        method, policy, cliff_weight = selected
        z = evaluate_episode(outer, base, y, X, D, bv, plain, cliff, method, policy)
        base_rho = safe_rho(y[z["query"]], base[z["query"]])
        base_mae = float(mean_absolute_error(y[z["query"]], base[z["query"]]))
        outer_rows.append({
            "group": outer_group, "n_group": len(outer), "n_query": len(z["query"]),
            "selected_method": method, "selected_policy": policy, "cliff_weight": cliff_weight,
            "inner_macro_Spearman": ranked[0][0], "Spearman": z["Spearman"],
            "Base_Spearman": base_rho, "delta_Spearman": z["Spearman"] - base_rho,
            "MAE": z["MAE"], "Base_MAE": base_mae, "delta_MAE": z["MAE"] - base_mae,
            "support_ids": "|".join(d.loc[z["support"], "canonical_molecule_id"].astype(str)),
        })
        print(outer_group, "selected", method, policy,
              f"rho={z['Spearman']:.3f} base={base_rho:.3f}", flush=True)

    out = pd.DataFrame(outer_rows)
    sel = pd.DataFrame(selection_rows)
    out.to_csv(OUT / "outer_episodes.csv", index=False)
    sel.to_csv(OUT / "inner_selection.csv", index=False)
    metrics = {
        "protocol": "nested_leave_one_series_out_3shot",
        "n_outer_series": len(out),
        "macro_Spearman": float(out.Spearman.mean()),
        "base_macro_Spearman_same_queries": float(out.Base_Spearman.mean()),
        "delta_macro_Spearman": float(out.delta_Spearman.mean()),
        "median_Spearman": float(out.Spearman.median()),
        "worst_Spearman": float(out.Spearman.min()),
        "positive_series": int((out.Spearman > 0).sum()),
        "macro_MAE": float(out.MAE.mean()),
        "base_macro_MAE_same_queries": float(out.Base_MAE.mean()),
        "delta_macro_MAE": float(out.delta_MAE.mean()),
        "selection_counts": out.groupby(["selected_method", "selected_policy"]).size().astype(int).to_dict(),
    }
    # JSON cannot encode tuple keys.
    metrics["selection_counts"] = {"|".join(k): v for k, v in metrics["selection_counts"].items()}
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
