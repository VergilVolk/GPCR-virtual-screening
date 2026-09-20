# -*- coding: utf-8 -*-
"""Frozen zero-/three-shot evaluation on the independent 2026 VU6025733 SAR."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

from run_pacer_fs_baselines import (
    adapted, choose_predicted_span, features, sims, train_delta,
)


P = Path(__file__).resolve().parents[1]
TRAIN = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
EXTERNAL = P / "data" / "benchmarks" / "m4_pam_v1" / "external_2026_vu6025733" / "external_potency.csv"
OUT = P / "results" / "pacer_external_vu6025733_v01"
SEED = 20260830


def model(seed: int) -> LGBMRegressor:
    return LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1,
        random_state=seed, verbosity=-1, n_jobs=6,
    )


def rho(y, p) -> float:
    value = spearmanr(y, p).statistic
    return float(value) if np.isfinite(value) else 0.0


def bootstrap_rho_delta(y, candidate, reference, n=50_000):
    rng = np.random.default_rng(SEED)
    delta = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        a, b = rho(y[idx], candidate[idx]), rho(y[idx], reference[idx])
        delta.append(a - b)
    values = np.asarray(delta, float)
    return {
        "estimate": float(rho(y, candidate) - rho(y, reference)),
        "ci95": [float(x) for x in np.quantile(values, [.025, .975])],
        "probability_gt_zero": float(np.mean(values > 0)),
        "bootstrap_replicates": n,
    }


def pair_direction_accuracy(frame, prediction):
    mapping = dict(zip(frame.compound_id, prediction))
    rows = []
    for suffix in sorted(set(x[1:] for x in frame.compound_id)):
        a, b = "5" + suffix, "6" + suffix
        true_delta = float(frame.loc[frame.compound_id == b, "pEC50"].iloc[0]
                           - frame.loc[frame.compound_id == a, "pEC50"].iloc[0])
        pred_delta = float(mapping[b] - mapping[a])
        rows.append({"pair": f"{a}->{b}", "true_delta": true_delta,
                     "predicted_delta": pred_delta,
                     "direction_correct": bool(np.sign(true_delta) == np.sign(pred_delta))})
    return rows, float(np.mean([r["direction_correct"] for r in rows]))


def top_quartile_recall(y, p):
    k = max(1, int(np.ceil(len(y) * .25)))
    truth = set(np.argsort(y)[-k:]); predicted = set(np.argsort(p)[-k:])
    return float(len(truth & predicted) / k)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    train = pd.read_csv(TRAIN).reset_index(drop=True)
    ext = pd.read_csv(EXTERNAL)
    ext = ext[ext.eligible_zero_shot.astype(bool)].reset_index(drop=True)
    combined_smiles = pd.concat([train.canonical_smiles, ext.canonical_smiles], ignore_index=True)
    X, D, bvs = features(combined_smiles)
    n_train = len(train)
    tr = np.arange(n_train); te = np.arange(n_train, n_train + len(ext))
    y_train = train.pEC50.to_numpy(float); y_ext = ext.pEC50.to_numpy(float)
    groups_train = train.source_component.astype(str).to_numpy()

    selected = np.argsort(X[tr].var(0))[-1024:]
    matrix = np.c_[X[:, selected], D]
    absolute = model(11).fit(matrix[tr], y_train).predict(matrix[te])

    centered_target = np.zeros(n_train, float)
    for group in sorted(set(groups_train)):
        idx = np.where(groups_train == group)[0]
        centered_target[idx] = y_train[idx] - y_train[idx].mean()
    centered_relative = model(12).fit(matrix[tr], centered_target).predict(matrix[te])
    centered = centered_relative + float(np.mean(y_train))

    similarity = sims(te, tr, bvs)
    nearest_idx = np.argmax(similarity, axis=1)
    nearest = y_train[nearest_idx]
    weights = np.maximum(similarity, 1e-6) ** 3
    knn = (weights @ y_train) / weights.sum(axis=1)

    zero = {
        "Absolute-QSAR": absolute,
        "PACER-Centered": centered,
        "1NN": nearest,
        "Similarity-kNN": knn,
    }
    zero_rows = []
    pair_rows = []
    for method_name, pred in zero.items():
        pairs, pair_acc = pair_direction_accuracy(ext, pred)
        pair_rows.extend({"protocol": "zero_shot", "method": method_name, **r} for r in pairs)
        zero_rows.append({
            "protocol": "zero_shot", "method": method_name, "n_query": len(ext),
            "Spearman": rho(y_ext, pred), "MAE": float(mean_absolute_error(y_ext, pred)),
            "top_quartile_recall": top_quartile_recall(y_ext, pred),
            "matched_pair_direction_accuracy": pair_acc,
        })

    # Predicted-span anchors are selected from label-blind absolute predictions.
    absolute_full = np.full(n_train + len(ext), np.nan); absolute_full[te] = absolute
    centered_full = np.full(n_train + len(ext), np.nan); centered_full[te] = centered
    support = choose_predicted_span(te, 3, absolute_full)
    support_set = set(map(int, support))
    query = np.asarray([i for i in te if int(i) not in support_set], int)
    y_all = np.r_[y_train, y_ext]
    g_all = np.r_[groups_train, np.repeat("EXTERNAL_2026", len(ext))]
    absolute_full[te] += float(np.mean(y_all[support] - absolute_full[support]))
    centered_full[te] += float(np.mean(y_all[support] - centered_full[support]))
    delta = train_delta(tr, X, D, y_all, g_all, bvs, seed=SEED, cliff_weight=0)
    fs_centered = adapted("DeltaSARHybrid", centered_full, support, query,
                          y_all, X, D, bvs, delta)
    fs_absolute = adapted("DeltaSARHybrid", absolute_full, support, query,
                          y_all, X, D, bvs, delta)
    local_q = query - n_train
    few = {
        "Absolute-3shot-offset": absolute_full[query],
        "Centered-3shot-offset": centered_full[query],
        "PACER-FS-Centered-3shot": fs_centered,
        "DeltaSAR-Absolute-3shot": fs_absolute,
    }
    few_rows = []
    for method_name, pred in few.items():
        few_rows.append({
            "protocol": "three_shot", "method": method_name, "n_query": len(query),
            "Spearman": rho(y_all[query], pred),
            "MAE": float(mean_absolute_error(y_all[query], pred)),
            "top_quartile_recall": top_quartile_recall(y_all[query], pred),
            "matched_pair_direction_accuracy": np.nan,
        })

    metrics = pd.DataFrame(zero_rows + few_rows)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(OUT / "matched_pair_directions.csv", index=False)

    predictions = ext[["compound_id", "canonical_smiles", "EC50_nM", "pEC50",
                       "max_train_tanimoto_ecfp4", "nearest_train_id"]].copy()
    for method_name, pred in zero.items():
        predictions[method_name] = pred
    predictions["three_shot_role"] = ["support" if i in support_set else "query" for i in te]
    for method_name, pred in few.items():
        predictions[method_name] = np.nan
        predictions.loc[local_q, method_name] = pred
    predictions.to_csv(OUT / "predictions.csv", index=False)

    z_abs = zero["Absolute-QSAR"]; z_ctr = zero["PACER-Centered"]
    q_abs = few["Absolute-3shot-offset"]; q_fs = few["PACER-FS-Centered-3shot"]
    audit = {
        "status": "independent_temporal_external_evaluation",
        "dataset": "M4_VU6025733_2026_TABLE1_HIGHCONF_V01",
        "n_external": len(ext),
        "n_exact_training_overlap": 0,
        "anchor_compound_ids": ext.iloc[support - n_train].compound_id.tolist(),
        "anchor_selection": "three predicted-span quantiles from zero-shot Absolute-QSAR; labels hidden",
        "zero_shot_centered_vs_absolute": bootstrap_rho_delta(y_ext, z_ctr, z_abs),
        "three_shot_pacer_fs_vs_absolute": bootstrap_rho_delta(y_all[query], q_fs, q_abs),
        "external_labels_used_for_training_or_model_selection": False,
        "model_configuration_frozen_before_prediction": True,
        "claim_boundary": "External functional SAR ranking only; not prospective PAM identity or wet-lab validation.",
    }
    audit["preregistered_external_support"] = bool(
        audit["three_shot_pacer_fs_vs_absolute"]["ci95"][0] > 0
    )
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(metrics.sort_values(["protocol", "Spearman"], ascending=[True, False]).to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
