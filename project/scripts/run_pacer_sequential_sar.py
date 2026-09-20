# -*- coding: utf-8 -*-
"""Frozen within-paper sequential SAR challenge: Table 1 -> Generation 3."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, roc_auc_score

from run_pacer_fs_baselines import adapted, features, sims, train_delta


P = Path(__file__).resolve().parents[1]
HIST = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
EXT = P / "data" / "benchmarks" / "m4_pam_v1" / "external_2026_vu6025733"
OUT = P / "results" / "pacer_sequential_sar_v01"
SEED = 20260830


def regressor(seed):
    return LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1,
        random_state=seed, verbosity=-1, n_jobs=6,
    )


def fit_predict(X, D, y, train, test, seed):
    selected = np.argsort(X[train].var(0))[-1024:]
    xx = np.c_[X[:, selected], D]
    return regressor(seed).fit(xx[train], y[train]).predict(xx[test])


def rho(y, p):
    value = spearmanr(y, p).statistic
    return float(value) if np.isfinite(value) else 0.0


def concordance(y, p):
    correct = comparable = 0
    for i in range(len(y)):
        for j in range(i + 1, len(y)):
            if y[i] == y[j]:
                continue
            comparable += 1
            correct += int(np.sign(y[i] - y[j]) == np.sign(p[i] - p[j]))
    return float(correct / comparable), comparable


def top3_recall(y, p):
    true = set(np.argsort(y)[-3:]); predicted = set(np.argsort(p)[-3:])
    return float(len(true & predicted) / 3)


def bootstrap_rho(y, p, n=20_000):
    rng = np.random.default_rng(SEED)
    values = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        values.append(rho(y[idx], p[idx]))
    return [float(x) for x in np.quantile(values, [.025, .975])]


def bootstrap_rho_delta(y, candidate, reference, n=50_000):
    rng = np.random.default_rng(SEED + 1)
    values = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        values.append(rho(y[idx], candidate[idx]) - rho(y[idx], reference[idx]))
    values = np.asarray(values, float)
    return {
        "estimate": float(rho(y, candidate) - rho(y, reference)),
        "ci95": [float(x) for x in np.quantile(values, [.025, .975])],
        "probability_gt_zero": float(np.mean(values > 0)),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    hist = pd.read_csv(HIST).reset_index(drop=True)
    table1 = pd.read_csv(EXT / "external_potency.csv").reset_index(drop=True)
    gen3 = pd.read_csv(EXT / "generation3_potency.csv").reset_index(drop=True)
    smiles = pd.concat([hist.canonical_smiles, table1.canonical_smiles,
                        gen3.canonical_smiles], ignore_index=True)
    X, D, bvs = features(smiles)
    nh, n1, n3 = len(hist), len(table1), len(gen3)
    ih = np.arange(nh); i1 = np.arange(nh, nh + n1); i3 = np.arange(nh + n1, nh + n1 + n3)
    eligible_local = np.where(gen3.eligible_sequential_query.astype(bool))[0]
    controls_local = np.where(~gen3.eligible_sequential_query.astype(bool))[0]
    query = i3[eligible_local]; controls = i3[controls_local]
    y = np.r_[hist.pEC50.to_numpy(float), table1.pEC50.to_numpy(float),
              gen3.pEC50_eval_ceiling.to_numpy(float)]
    groups = np.r_[hist.source_component.astype(str).to_numpy(),
                   np.repeat("VU2026_TABLE1", n1), np.repeat("VU2026_GEN3", n3)]

    historical_absolute = fit_predict(X, D, y, ih, i3, 101)
    augmented = np.r_[ih, i1]
    augmented_absolute = fit_predict(X, D, y, augmented, i3, 102)

    similarity = sims(i3, i1, bvs)
    weight = np.maximum(similarity, 1e-6) ** 3
    table1_knn = (weight @ y[i1]) / weight.sum(axis=1)

    delta_historical = train_delta(ih, X, D, y, groups, bvs, seed=103)
    delta_augmented = train_delta(augmented, X, D, y, groups, bvs, seed=104)
    hist_delta = adapted("DeltaSAR", np.full(len(y), np.nan), i1, i3,
                         y, X, D, bvs, delta_historical)
    aug_delta = adapted("DeltaSAR", np.full(len(y), np.nan), i1, i3,
                        y, X, D, bvs, delta_augmented)

    all_predictions = {
        "Historical-Absolute": historical_absolute,
        "Augmented-Absolute": augmented_absolute,
        "Table1-kNN": table1_knn,
        "Historical-DeltaSAR": hist_delta,
        "Augmented-DeltaSAR": aug_delta,
    }
    true_query = y[query]
    exact_query = gen3.iloc[eligible_local].censoring.eq("exact").to_numpy()
    potent = (true_query >= 7.0).astype(int)
    metrics = []
    for name, full_pred in all_predictions.items():
        pred = full_pred[eligible_local]
        c, pairs = concordance(true_query, pred)
        rho_ci = bootstrap_rho(true_query, pred)
        metrics.append({
            "method": name, "n_query": len(query),
            "ceiling_Spearman": rho(true_query, pred),
            "rho_bootstrap_low": rho_ci[0],
            "rho_bootstrap_high": rho_ci[1],
            "pairwise_concordance": c, "n_comparable_pairs": pairs,
            "potent_100nM_ROC_AUC": float(roc_auc_score(potent, pred)),
            "top3_recall": top3_recall(true_query, pred),
            "MAE_exact_only": float(mean_absolute_error(true_query[exact_query], pred[exact_query])),
        })
    metric_frame = pd.DataFrame(metrics).sort_values("ceiling_Spearman", ascending=False)
    metric_frame.to_csv(OUT / "metrics.csv", index=False)

    predictions = gen3[["compound_id", "canonical_smiles", "pEC50", "censoring",
                        "eligible_sequential_query", "structural_note"]].copy()
    for name, pred in all_predictions.items():
        predictions[name] = pred
    predictions.to_csv(OUT / "predictions.csv", index=False)

    audit = {
        "status": "retrospective_within_paper_sequential_sar_challenge",
        "n_round1_labels": n1, "n_main_query": len(query), "n_overlap_controls": len(controls),
        "n_censored_inactive_query": int((~exact_query).sum()),
        "generation3_labels_used_for_training_or_selection": False,
        "rules_frozen_before_first_prediction": True,
        "analyst_blinded_to_generation3_labels": False,
        "best_method_by_ceiling_spearman": str(metric_frame.iloc[0].method),
        "augmented_absolute_vs_historical_absolute": bootstrap_rho_delta(
            true_query, augmented_absolute[eligible_local], historical_absolute[eligible_local]
        ),
        "claim_boundary": "Small within-paper sequential challenge; not independent or prospective SOTA evidence.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(metric_frame.to_string(index=False))
    print("\nPredictions")
    print(predictions[["compound_id", "pEC50", "censoring", "eligible_sequential_query"]
                      + list(all_predictions)].to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
