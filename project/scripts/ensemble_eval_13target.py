#!/usr/bin/env python3
"""3-seed ensemble evaluation for the 13-target LOSO migration:
mean of per-seed tuned scores per row -> macro metrics + paired
target-stratified scaffold-cluster bootstrap (same protocol as the
per-seed summarizer). Also evaluates mean of official columns as sanity."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd


def macro(table, col, targets):
    from sklearn.metrics import roc_auc_score, average_precision_score
    rows = []
    for t in targets:
        sub = table[table.target == t]
        y, s = sub.label.values, sub[col].values
        order = np.argsort(-s, kind="mergesort")
        rows.append({
            "roc_auc": roc_auc_score(y, s),
            "pr_auc": average_precision_score(y, s),
            "bedroc_alpha20": bedroc(y[order], 20.0),
            "ef1pct": (y[order[:max(1, len(y) // 100)]].mean() / y.mean()) if y.mean() else np.nan,
        })
    return {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}, rows


def bedroc(y_sorted, alpha):
    n = len(y_sorted)
    r = np.arange(1, n + 1)
    return float((y_sorted * np.exp(-alpha * r / n)).sum() / (y_sorted.sum() * (1 - np.exp(-alpha) / (1 - np.exp(-alpha / n))) ) )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--predictions", type=Path, nargs="+", required=True)
    p.add_argument("--ecfp-predictions", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=20000)
    p.add_argument("--seed", type=int, default=20260926)
    a = p.parse_args()
    tables = [pd.read_csv(f) for f in a.predictions]
    keys = ["target", "label", "canonical_smiles", "murcko_scaffold"]
    base = tables[0][keys].copy()
    base["official"] = tables[0]["official"]
    base["tuned"] = np.mean([t["tuned"].values for t in tables], axis=0)
    targets = sorted(base.target.unique())
    ens, _ = macro(base, "tuned", targets)
    off, _ = macro(base, "official", targets)
    result = {"protocol": "ensemble = mean of 3 per-seed tuned scores (identical rows)",
              "n_targets": len(targets), "ensemble_macro": ens, "official_macro": off}
    if a.ecfp_predictions:
        e = pd.read_csv(a.ecfp_predictions)
        m = base.merge(e[keys + ["ecfp4_logistic"]], on=keys, how="inner")
        ecfp, _ = macro(m, "ecfp4_logistic", targets)
        result["ecfp4_macro"] = ecfp
    # paired bootstrap: target-stratified scaffold-cluster resampling (vectorized)
    rng = np.random.default_rng(a.seed)
    base["cluster"] = base.target + "|" + base.murcko_scaffold.astype(str)
    clusters = sorted(base.cluster.unique())
    cl_idx = {c: i for i, c in enumerate(clusters)}
    row_cluster = base.cluster.map(cl_idx).values
    row_y = base.label.values.astype(int)
    row_ens, row_off = base.tuned.values, base.official.values
    tgt_clusters = {t: np.unique(row_cluster[base.target.values == t]) for t in targets}
    cl_rows = {i: np.flatnonzero(row_cluster == i) for i in range(len(clusters))}

    def fast_auc(y, s):
        n1, n0 = int(y.sum()), int((1 - y).sum())
        if not n1 or not n0:
            return None
        r = np.argsort(np.argsort(s))
        return (r[y == 1].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)

    d = []
    for _ in range(a.bootstrap):
        be, bo = [], []
        for t, cls in tgt_clusters.items():
            pick = rng.choice(cls, size=len(cls), replace=True)
            sel = np.concatenate([cl_rows[c] for c in pick])
            ya = row_y[sel]
            if ya.sum() == 0 or ya.min() == ya.max():
                continue
            ue, uo = fast_auc(ya, row_ens[sel]), fast_auc(ya, row_off[sel])
            if ue is not None:
                be.append(ue); bo.append(uo)
        if be:
            d.append(np.mean(be) - np.mean(bo))
    d = np.array(d)
    result["bootstrap_95ci"] = {"ens_minus_off_roc": [float(np.quantile(d, .025)), float(np.quantile(d, .975))]}
    a.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "protocol"}, indent=2))


if __name__ == "__main__":
    main()
