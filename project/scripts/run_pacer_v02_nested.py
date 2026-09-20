# -*- coding: utf-8 -*-
"""PACER-M4 v0.2: nested selection, calibration and cross-model fusion.

The final three negative-bearing medicinal-chemistry series are outer tests.
For each outer test series:
  1. inner scaffold-group CV selects multitask loss weights;
  2. inner OOF predictions fit Platt calibration and a decision threshold;
  3. inner OOF predictions select fusion weight with Extra Trees;
  4. models are retrained on all outer-training molecules;
  5. the untouched outer series is predicted once.

No choice is made using outer-test labels.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score, balanced_accuracy_score, brier_score_loss,
    matthews_corrcoef, roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold

import run_pacer_multitask as core

PROJECT = Path(__file__).resolve().parents[1]
BENCH = PROJECT / "data" / "benchmarks" / "m4_pam_v1"
OUT = PROJECT / "results" / "pacer_v02_nested"
SEEDS = [42, 43, 44]
INNER_SEED = 1701
INNER_EPOCHS = 70
FINAL_EPOCHS = 120
CONFIGS = [
    {"name": "reg025_pair010", "reg": 0.25, "pair": 0.10},
    {"name": "reg050_pair010", "reg": 0.50, "pair": 0.10},
    {"name": "reg050_pair025", "reg": 0.50, "pair": 0.25},
    {"name": "reg050_pair050", "reg": 0.50, "pair": 0.50},
    {"name": "reg100_pair025", "reg": 1.00, "pair": 0.25},
]


def mcc_threshold(y: np.ndarray, p: np.ndarray) -> tuple[float, float]:
    best = (-2.0, 0.5)
    for threshold in np.linspace(0.05, 0.95, 181):
        score = matthews_corrcoef(y, p >= threshold)
        if score > best[0]:
            best = (float(score), float(threshold))
    return best[1], best[0]


def enrichment(y: np.ndarray, score: np.ndarray, fraction: float) -> float:
    k = max(1, int(np.ceil(len(y) * fraction)))
    idx = np.argsort(score)[::-1][:k]
    return float(y[idx].mean() / y.mean())


def metrics(y: np.ndarray, p: np.ndarray, threshold: float) -> dict:
    pred = p >= threshold
    return {
        "ROC_AUC": float(roc_auc_score(y, p)),
        "PR_AUC": float(average_precision_score(y, p)),
        "MCC": float(matthews_corrcoef(y, pred)),
        "BalancedAcc": float(balanced_accuracy_score(y, pred)),
        "Brier": float(brier_score_loss(y, p)),
        "threshold": float(threshold),
        "positive_prediction_rate": float(pred.mean()),
        "EF_5pct": enrichment(y, p, 0.05),
        "EF_10pct": enrichment(y, p, 0.10),
    }


def fit_platt(y: np.ndarray, p: np.ndarray) -> LogisticRegression:
    eps = 1e-6
    logits = np.log(np.clip(p, eps, 1-eps) / np.clip(1-p, eps, 1-eps)).reshape(-1, 1)
    model = LogisticRegression(C=1.0, class_weight="balanced", random_state=42)
    model.fit(logits, y)
    return model


def apply_platt(model: LogisticRegression, p: np.ndarray) -> np.ndarray:
    eps = 1e-6
    logits = np.log(np.clip(p, eps, 1-eps) / np.clip(1-p, eps, 1-eps)).reshape(-1, 1)
    return model.predict_proba(logits)[:, 1]


def train_predict_neural(train_idx, test_idx, config, seed, df, fps, desc, pairs, epochs):
    train_x, test_x = core.preprocess(fps[train_idx], fps[test_idx], desc[train_idx], desc[test_idx])
    model, mean, std = core.train_model(
        train_x, df.loc[train_idx, "target"].to_numpy(int),
        df.loc[train_idx, "pEC50"].to_numpy(float),
        df.loc[train_idx, "canonical_molecule_id"].tolist(), pairs,
        "MultiTaskPair", seed, epochs, config["reg"], config["pair"],
    )
    model.eval()
    with torch.no_grad():
        logits, reg = model(torch.from_numpy(test_x))
    return torch.sigmoid(logits).numpy(), reg.numpy() * std + mean


def train_predict_tree(train_idx, test_idx, fps, desc, y, seed):
    train_x, test_x = core.preprocess(fps[train_idx], fps[test_idx], desc[train_idx], desc[test_idx])
    model = ExtraTreesClassifier(
        n_estimators=600, class_weight="balanced", min_samples_leaf=2,
        random_state=seed, n_jobs=-1,
    )
    model.fit(train_x, y[train_idx])
    return model.predict_proba(test_x)[:, 1]


def inner_predictions(outer_train: np.ndarray, config: dict, df, fps, desc, pairs):
    y = df["target"].to_numpy(int)
    groups = df["murcko_scaffold"].to_numpy(str)
    local = np.where(outer_train)[0]
    splitter = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=INNER_SEED)
    neural = np.full(len(df), np.nan); tree = np.full(len(df), np.nan)
    for inner_train_local, inner_val_local in splitter.split(
        np.zeros((len(local), 1)), y[local], groups[local]
    ):
        tr, va = local[inner_train_local], local[inner_val_local]
        neural[va], _ = train_predict_neural(
            tr, va, config, 42, df, fps, desc, pairs, INNER_EPOCHS
        )
        tree[va] = train_predict_tree(tr, va, fps, desc, y, 42)
    return local, neural[local], tree[local]


def choose_config(outer_train, df, fps, desc, pairs):
    y = df["target"].to_numpy(int)
    candidates = []
    cache = {}
    for config in CONFIGS:
        local, neural, tree = inner_predictions(outer_train, config, df, fps, desc, pairs)
        # Ranking-first nested criterion with a small Brier penalty.
        auc = roc_auc_score(y[local], neural)
        brier = brier_score_loss(y[local], neural)
        score = auc - 0.15 * brier
        candidates.append({**config, "inner_AUC": float(auc), "inner_Brier": float(brier),
                           "selection_score": float(score)})
        cache[config["name"]] = (local, neural, tree)
    best = max(candidates, key=lambda x: x["selection_score"])
    return best, candidates, cache[best["name"]]


def choose_fusion(y, neural, tree):
    best = None
    for alpha in np.linspace(0, 1, 21):
        fused = alpha * neural + (1-alpha) * tree
        auc = roc_auc_score(y, fused); brier = brier_score_loss(y, fused)
        score = auc - 0.15 * brier
        item = {"alpha_neural": float(alpha), "AUC": float(auc),
                "Brier": float(brier), "selection_score": float(score)}
        if best is None or item["selection_score"] > best["selection_score"]:
            best = item
    return best


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    strict = pd.read_csv(BENCH / "pam_vs_inactive.csv")
    potency = pd.read_csv(BENCH / "potency_molecules.csv")[["canonical_molecule_id", "pEC50"]]
    df = strict.merge(potency, on="canonical_molecule_id", how="left").reset_index(drop=True)
    pairs = pd.read_csv(BENCH / "activity_cliffs" / "near_neighbor_pairs.csv")
    fps, desc = core.raw_features(df["canonical_smiles"].tolist())
    y = df["target"].to_numpy(int)
    eval_mask = df["series_holdout_fold"].to_numpy() >= 0
    predictions = {name: np.full(len(df), np.nan) for name in ["neural_raw", "neural_cal", "tree_raw", "fused"]}
    potency_pred = np.full(len(df), np.nan)
    thresholds = {}; fold_reports = []

    for fold in sorted(x for x in df["series_holdout_fold"].unique() if x >= 0):
        test = df["series_holdout_fold"].to_numpy() == fold
        train = ~test
        best, candidates, (inner_idx, inner_neural, inner_tree) = choose_config(
            train, df, fps, desc, pairs
        )
        calibrator = fit_platt(y[inner_idx], inner_neural)
        inner_cal = apply_platt(calibrator, inner_neural)
        fusion = choose_fusion(y[inner_idx], inner_cal, inner_tree)
        inner_fused = fusion["alpha_neural"] * inner_cal + (1-fusion["alpha_neural"]) * inner_tree
        threshold, inner_mcc = mcc_threshold(y[inner_idx], inner_fused)

        neural_seeds, reg_seeds, tree_seeds = [], [], []
        tr_idx, te_idx = np.where(train)[0], np.where(test)[0]
        for seed in SEEDS:
            p, reg = train_predict_neural(
                tr_idx, te_idx, best, seed, df, fps, desc, pairs, FINAL_EPOCHS
            )
            neural_seeds.append(p); reg_seeds.append(reg)
            tree_seeds.append(train_predict_tree(tr_idx, te_idx, fps, desc, y, seed))
        neural_raw = np.vstack(neural_seeds).mean(0)
        neural_cal = apply_platt(calibrator, neural_raw)
        tree_raw = np.vstack(tree_seeds).mean(0)
        fused = fusion["alpha_neural"] * neural_cal + (1-fusion["alpha_neural"]) * tree_raw
        predictions["neural_raw"][test] = neural_raw
        predictions["neural_cal"][test] = neural_cal
        predictions["tree_raw"][test] = tree_raw
        predictions["fused"][test] = fused
        potency_pred[test] = np.vstack(reg_seeds).mean(0)
        thresholds[int(fold)] = threshold
        fold_reports.append({
            "fold": int(fold), "n_test": int(test.sum()), "selected_config": best,
            "all_inner_configs": candidates, "fusion": fusion,
            "inner_threshold": threshold, "inner_MCC": inner_mcc,
            "outer_fused": metrics(y[test], fused, threshold),
        })
        print("fold", fold, "config", best["name"], "fusion", fusion,
              "outer", fold_reports[-1]["outer_fused"], flush=True)

    eval_df = df[eval_mask].copy()
    eval_y = y[eval_mask]
    threshold_vector = np.asarray([thresholds[int(f)] for f in eval_df["series_holdout_fold"]])
    report = {"folds": fold_reports, "aggregate": {}}
    for name, vector in predictions.items():
        p = vector[eval_mask]
        # Fold-specific nested thresholds; metrics computed explicitly here.
        pred = p >= threshold_vector
        item = {
            "ROC_AUC": float(roc_auc_score(eval_y, p)),
            "PR_AUC": float(average_precision_score(eval_y, p)),
            "MCC": float(matthews_corrcoef(eval_y, pred)),
            "BalancedAcc": float(balanced_accuracy_score(eval_y, pred)),
            "Brier": float(brier_score_loss(eval_y, p)),
            "positive_prediction_rate": float(pred.mean()),
            "EF_5pct": enrichment(eval_y, p, 0.05),
            "EF_10pct": enrichment(eval_y, p, 0.10),
        }
        report["aggregate"][name] = item
    report["macro_outer_fused"] = {
        metric: float(np.nanmean([fold["outer_fused"].get(metric, np.nan) for fold in fold_reports]))
        for metric in ["ROC_AUC", "PR_AUC", "MCC", "BalancedAcc", "Brier", "EF_5pct", "EF_10pct"]
    }
    report["worst_outer_fused"] = {
        "ROC_AUC": float(min(fold["outer_fused"]["ROC_AUC"] for fold in fold_reports)),
        "BalancedAcc": float(min(fold["outer_fused"]["BalancedAcc"] for fold in fold_reports)),
        "EF_5pct": float(min(fold["outer_fused"]["EF_5pct"] for fold in fold_reports)),
    }
    report["validity_warning"] = (
        "Pooled metrics combine series with different class prevalence and score scales. "
        "Use macro/worst outer-series metrics as the primary domain-generalization result."
    )
    eval_df["pEC50_prediction"] = potency_pred[eval_mask]
    for name, vector in predictions.items(): eval_df[name] = vector[eval_mask]
    eval_df["nested_threshold"] = threshold_vector
    eval_df.to_csv(OUT / "series_nested_oof.csv", index=False)
    with open(OUT / "nested_metrics.json", "w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
    print(json.dumps({"pooled": report["aggregate"], "macro": report["macro_outer_fused"],
                      "worst": report["worst_outer_fused"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
