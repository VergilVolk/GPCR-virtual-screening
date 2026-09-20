# -*- coding: utf-8 -*-
"""Leakage and chemical-domain audit for centered PACER-SAR."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from rdkit import DataStructs

from run_pacer_fs_baselines import DATA, MIN_GROUP, choose_predicted_span, features, safe_rho


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_centered_domain_audit_v01"
OUT.mkdir(parents=True, exist_ok=True)


def fit(xx, target, train):
    return LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1, random_state=42,
        verbosity=-1, n_jobs=6,
    ).fit(xx[train], target[train])


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bv = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float); g = d.source_component.astype(str).to_numpy()
    groups = [x for x, n in d.source_component.value_counts().items() if n >= MIN_GROUP]
    molecule_rows = []
    for group in groups:
        outer = np.where(g == group)[0]; train = np.where(g != group)[0]
        sel = np.argsort(X[train].var(0))[-1024:]; xx = np.c_[X[:, sel], D]
        z = np.zeros(len(d))
        for source in set(g[train]):
            idx = train[g[train] == source]; z[idx] = y[idx] - y[idx].mean()
        pa = fit(xx, y, train).predict(xx[outer]); pc = fit(xx, z, train).predict(xx[outer])
        full = np.full(len(d), np.nan); full[outer] = pa
        support = choose_predicted_span(outer, 3, full); support_set = set(support)
        train_scaffolds = set(d.loc[train, "murcko_scaffold"])
        for j, idx in enumerate(outer):
            if idx in support_set: continue
            mx = max(DataStructs.BulkTanimotoSimilarity(bv[idx], [bv[t] for t in train]))
            molecule_rows.append({
                "group": group, "canonical_molecule_id": d.at[idx, "canonical_molecule_id"],
                "pEC50": y[idx], "absolute_prediction": pa[j], "centered_prediction": pc[j],
                "max_train_tanimoto": mx,
                "novel_scaffold": d.at[idx, "murcko_scaffold"] not in train_scaffolds,
            })
    mol = pd.DataFrame(molecule_rows); mol.to_csv(OUT / "molecule_predictions.csv", index=False)
    rows = []
    domains = {
        "all_query": np.ones(len(mol), dtype=bool),
        "tanimoto_lt_0.70": mol.max_train_tanimoto.to_numpy() < .70,
        "novel_murcko_scaffold": mol.novel_scaffold.to_numpy(bool),
        "novel_scaffold_and_tanimoto_lt_0.70": mol.novel_scaffold.to_numpy(bool) & (mol.max_train_tanimoto.to_numpy() < .70),
    }
    for domain, mask in domains.items():
        sub = mol[mask]
        for group, x in sub.groupby("group"):
            if len(x) < 5: continue
            rows.append({"domain": domain, "group": group, "n": len(x),
                         "absolute_rho": safe_rho(x.pEC50.to_numpy(), x.absolute_prediction.to_numpy()),
                         "centered_rho": safe_rho(x.pEC50.to_numpy(), x.centered_prediction.to_numpy())})
    per = pd.DataFrame(rows); per["delta_rho"] = per.centered_rho - per.absolute_rho
    per.to_csv(OUT / "per_series.csv", index=False)
    summary = per.groupby("domain").agg(
        n_series=("group", "size"), n_molecules=("n", "sum"),
        absolute_macro=("absolute_rho", "mean"), centered_macro=("centered_rho", "mean"),
        delta_macro=("delta_rho", "mean"), centered_worst=("centered_rho", "min"),
    ).reset_index()
    summary.to_csv(OUT / "summary.csv", index=False)
    (OUT / "metrics.json").write_text(json.dumps(summary.to_dict("records"), indent=2), encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
