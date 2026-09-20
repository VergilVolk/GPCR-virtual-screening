"""Retrospective cross-campaign replication of assay-conditional adaptation."""

from __future__ import annotations

import json
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import mean_absolute_error

from evaluate_multidomain_sar import features, fit_predict, inner_select, local_fit_predict


P = Path(__file__).resolve().parents[1]
HIST = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
EXT = P / "data" / "benchmarks" / "m4_pam_v1" / "external_2026_vu6025733"
OLD = P / "results" / "pacer_sequential_sar_v01" / "predictions.csv"
OUT = P / "results" / "pacer_multidomain_vu_sequential_v01"


def rho(y, prediction):
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def exact_permutation_p(y, prediction):
    left = rankdata(y).astype(float); left -= left.mean()
    right = rankdata(prediction).astype(float); right -= right.mean()
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator == 0:
        return 1.0
    observed = abs(float(left @ right / denominator))
    extreme, total = 0, 0
    for order in itertools.permutations(range(len(right))):
        total += 1
        value = abs(float(left @ right[list(order)] / denominator))
        extreme += int(value >= observed - 1e-12)
    return float(extreme / total)


def bootstrap_delta(y, candidate, reference, n=50_000):
    rng = np.random.default_rng(20260830)
    values = []
    for _ in range(n):
        index = rng.integers(0, len(y), len(y))
        values.append(rho(y[index], candidate[index]) - rho(y[index], reference[index]))
    values = np.asarray(values, float)
    return {
        "estimate": rho(y, candidate) - rho(y, reference),
        "ci95": [float(x) for x in np.quantile(values, [.025, .975])],
        "probability_gt_zero": float(np.mean(values > 0)),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    history = pd.read_csv(HIST)
    support = pd.read_csv(EXT / "external_potency.csv")
    query = pd.read_csv(EXT / "generation3_potency.csv")
    query = query[query.eligible_sequential_query.astype(bool)].reset_index(drop=True)
    combined = pd.concat([
        support[["canonical_smiles", "pEC50"]],
        query[["canonical_smiles", "pEC50"]],
    ], ignore_index=True)
    history_X = features(history.canonical_smiles)
    local_X = features(combined.canonical_smiles)
    history_y = history.pEC50.to_numpy(float)
    local_y = combined.pEC50.to_numpy(float)
    train = np.arange(len(support))
    test = np.arange(len(support), len(combined))

    choices = {
        mode: inner_select(history_X, history_y, local_X, local_y, train, mode)
        for mode in ("local", "pooled", "conditional")
    }
    predictions = {
        "Local-Molecular-Ridge": local_fit_predict(
            local_X, local_y, train, test, choices["local"][0]
        ),
        "Pooled-Ridge": fit_predict(
            history_X, history_y, local_X, local_y, train, test, False,
            choices["pooled"][0], choices["pooled"][1],
        ),
        "PACER-AssayConditional": fit_predict(
            history_X, history_y, local_X, local_y, train, test, True,
            choices["conditional"][0], choices["conditional"][1],
        ),
    }
    old = pd.read_csv(OLD)
    old = old[old.eligible_sequential_query.astype(bool)].set_index("compound_id")
    for name in ["Historical-Absolute", "Augmented-Absolute", "Table1-kNN"]:
        predictions[name] = old.loc[query.compound_id, name].to_numpy(float)

    y = query.pEC50_eval_ceiling.to_numpy(float)
    metrics = pd.DataFrame([{
        "method": name, "n": len(query), "ceiling_Spearman": rho(y, prediction),
        "MAE_vs_eval_ceiling": float(mean_absolute_error(y, prediction)),
        "exact_permutation_p_two_sided": exact_permutation_p(y, prediction),
    } for name, prediction in predictions.items()]).sort_values("ceiling_Spearman", ascending=False)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    result = query[["compound_id", "canonical_smiles", "pEC50", "pEC50_eval_ceiling", "censoring"]].copy()
    for name, prediction in predictions.items():
        result[name] = prediction
    result.to_csv(OUT / "predictions.csv", index=False)
    audit = {
        "status": "retrospective_second_campaign_replication",
        "n_support_round1": len(support),
        "n_query_generation3": len(query),
        "inner_selected_hyperparameters": choices,
        "generation3_labels_used_for_training_or_selection": False,
        "pooled_vs_augmented_absolute_spearman_delta": bootstrap_delta(
            y, predictions["Pooled-Ridge"], predictions["Augmented-Absolute"]
        ),
        "independent_claim": False,
        "warning": "Generation-3 labels had been examined in earlier project work; this is replication, not untouched validation.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(metrics.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
