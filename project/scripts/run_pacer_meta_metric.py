# -*- coding: utf-8 -*-
"""PACER-MetaMetric: nested episodic metric learning for 3-shot M4 PAM SAR."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

from run_pacer_delta_residual import (
    DATA, MIN_GROUP, N_ANCHOR, features, fit_base, predict_base,
    inner_oos_base, predicted_span, similarities, bootstrap_delta,
)
from run_pacer_fs_baselines import train_delta, adapted

P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_meta_metric_v01"
OUT.mkdir(parents=True, exist_ok=True)
EPOCHS = 220
SEED = 830


class MetricNet(torch.nn.Module):
    def __init__(self, n_features):
        super().__init__()
        self.encoder = torch.nn.Sequential(
            torch.nn.Linear(n_features, 128), torch.nn.ReLU(),
            torch.nn.Linear(128, 32),
        )
        self.log_temperature = torch.nn.Parameter(torch.tensor(-1.0))

    def embed(self, x):
        return torch.nn.functional.normalize(self.encoder(x), dim=-1)

    def residual(self, xq, xs, support_residual):
        zq, zs = self.embed(xq), self.embed(xs)
        distance = ((zq[:, None, :] - zs[None, :, :]) ** 2).sum(-1)
        temperature = torch.exp(self.log_temperature).clamp(.05, 2.0)
        weight = torch.softmax(-distance / temperature, dim=1)
        return weight @ support_residual


def prepare_matrix(X, D, base_oos, train, target):
    sel = np.argsort(X[train].var(0))[-256:]
    mean = D[train].mean(0); std = D[train].std(0); std[std < 1e-8] = 1
    # Base prediction is a label-free covariate and is standardized only on
    # outer-training OOS values.
    bmean = float(np.nanmean(base_oos[train])); bstd = float(np.nanstd(base_oos[train])) or 1.0
    matrix = np.c_[X[:, sel], (D - mean) / std, (base_oos - bmean) / bstd]
    return matrix.astype(np.float32), sel, mean, std, bmean, bstd


def safe_rho(y, p):
    value = spearmanr(y, p).statistic
    return float(value) if np.isfinite(value) else 0.0


def train_metric(matrix, y, groups, base_oos, train, seed):
    torch.manual_seed(seed); np.random.seed(seed)
    model = MetricNet(matrix.shape[1])
    optimizer = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
    tx = torch.from_numpy(matrix)
    ty = torch.from_numpy(y.astype(np.float32))
    tb = torch.from_numpy(base_oos.astype(np.float32))
    train_groups = sorted(set(groups[train]))
    rng = np.random.default_rng(seed)
    losses = []
    for epoch in range(EPOCHS):
        optimizer.zero_grad(); total = torch.tensor(0.0); used = 0
        for group in train_groups:
            idx = train[groups[train] == group]
            if len(idx) <= N_ANCHOR + 2:
                continue
            if epoch % 5 == 0:
                support = predicted_span(idx, N_ANCHOR, base_oos)
            else:
                support = rng.choice(idx, N_ANCHOR, replace=False)
            support_set = set(map(int, support))
            query = np.asarray([i for i in idx if int(i) not in support_set], int)
            sr = ty[support] - tb[support]
            pred = tb[query] + model.residual(tx[query], tx[support], sr)
            mse = torch.mean((pred - ty[query]) ** 2)
            # Deterministic adjacent-by-potency pairs provide a low-variance
            # rank objective without using outer-test labels.
            order = query[np.argsort(y[query])]
            if len(order) > 1:
                low, high = order[:-1], order[1:]
                pair_pred_low = tb[low] + model.residual(tx[low], tx[support], sr)
                pair_pred_high = tb[high] + model.residual(tx[high], tx[support], sr)
                rank_loss = torch.nn.functional.softplus(-(pair_pred_high - pair_pred_low)).mean()
            else:
                rank_loss = torch.tensor(0.0)
            total = total + mse + .2 * rank_loss; used += 1
        loss = total / max(used, 1)
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step(); losses.append(float(loss.detach()))
    return model, losses


def main():
    torch.set_num_threads(6)
    data = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bvs = features(data.canonical_smiles)
    y = data.pEC50.to_numpy(float)
    groups = data.source_component.astype(str).to_numpy()
    eval_groups = [g for g, n in data.source_component.value_counts().items() if n >= MIN_GROUP]
    rows, learning, molecule_rows = [], [], []
    for oi, held in enumerate(eval_groups):
        target = np.flatnonzero(groups == held); train = np.flatnonzero(groups != held)
        model, sel = fit_base(X, D, y, train, 5000 + oi)
        base = np.full(len(y), np.nan); base[target] = predict_base(model, sel, X, D, target)
        base_oos = inner_oos_base(X, D, y, groups, train, 50 + oi)
        # Fill target base before building the feature matrix; labels remain hidden.
        base_oos[target] = base[target]
        matrix, *_ = prepare_matrix(X, D, base_oos, train, target)
        metric, losses = train_metric(matrix, y, groups, base_oos, train, SEED + oi)
        support = predicted_span(target, N_ANCHOR, base)
        query = np.asarray([i for i in target if i not in set(support)], int)
        with torch.no_grad():
            tx = torch.from_numpy(matrix)
            sr = torch.from_numpy((y[support] - base[support]).astype(np.float32))
            meta = base[query] + metric.residual(tx[query], tx[support], sr).numpy()
        delta_pack = train_delta(train, X, D, y, groups, bvs, seed=42 + oi)
        fs = adapted("DeltaSARHybrid", base, support, query, y, X, D, bvs, delta_pack)
        hybrid = .5 * meta + .5 * fs
        sim = similarities(query, support, bvs); w = np.maximum(sim, 1e-4) ** 3
        fixed = base[query] + (w * (y[support] - base[support])[None, :]).sum(1) / (w.sum(1) + 1.5)
        predictions = {
            "Absolute": base[query], "AnchorResidual": fixed, "PACER-FS": fs,
            "PACER-MetaMetric": meta, "PACER-Hybrid": hybrid,
        }
        for name, pred in predictions.items():
            rows.append({
                "group": held, "method": name, "n_group": len(target), "n_query": len(query),
                "Spearman": safe_rho(y[query], pred), "MAE": float(mean_absolute_error(y[query], pred)),
                "support_ids": "|".join(data.loc[support, "canonical_molecule_id"]),
            })
            for idx, value in zip(query, pred):
                molecule_rows.append({
                    "group": held, "method": name,
                    "canonical_molecule_id": data.loc[idx, "canonical_molecule_id"],
                    "true_pEC50": y[idx], "prediction": float(value),
                })
        learning.append({"group": held, "initial_loss": losses[0], "final_loss": losses[-1],
                         "learned_temperature": float(torch.exp(metric.log_temperature).detach())})
        print(held, "rho", safe_rho(y[query], meta), "loss", losses[-1], flush=True)

    result = pd.DataFrame(rows); result.to_csv(OUT / "series_results.csv", index=False)
    pd.DataFrame(molecule_rows).to_csv(OUT / "molecule_predictions.csv", index=False)
    pd.DataFrame(learning).to_csv(OUT / "training_audit.csv", index=False)
    summary = result.groupby("method").agg(
        macro_Spearman=("Spearman", "mean"), median_Spearman=("Spearman", "median"),
        worst_Spearman=("Spearman", "min"), positive_series=("Spearman", lambda x: int((x > 0).sum())),
        macro_MAE=("MAE", "mean"), n_series=("Spearman", "size"),
    ).reset_index().sort_values("macro_Spearman", ascending=False)
    summary.to_csv(OUT / "summary.csv", index=False)
    delta = bootstrap_delta(result, "PACER-MetaMetric")
    hybrid_delta = bootstrap_delta(result, "PACER-Hybrid")
    score = float(summary.set_index("method").loc["PACER-MetaMetric", "macro_Spearman"])
    audit = {
        "protocol": "nested 11-series LOSO; 3 anchors; 220 fixed episodic epochs",
        "outer_query_labels_used_for_training_or_selection": False,
        "PACER_MetaMetric_vs_Absolute": delta,
        "PACER_Hybrid_vs_Absolute": hybrid_delta,
        "macro_Spearman": score,
        "PACER_Hybrid_macro_Spearman": float(
            summary.set_index("method").loc["PACER-Hybrid", "macro_Spearman"]
        ),
        "preregistered_pass": bool(
            hybrid_delta["ci95"][0] > 0
            and summary.set_index("method").loc["PACER-Hybrid", "macro_Spearman"] >= .263
        ),
        "claim_boundary": "internal few-shot functional potency benchmark only",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(summary.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
