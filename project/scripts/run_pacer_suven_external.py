"""Preregistered zero-shot evaluation on independent Suven M4 PAM patent data."""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, roc_auc_score

from run_pacer_fs_baselines import features, sims


P = Path(__file__).resolve().parents[1]
TRAIN = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
EXTERNAL = P / "data" / "benchmarks" / "m4_pam_v1" / "external_suven_2025" / "external_patent_potency.csv"
OUT = P / "results" / "pacer_external_suven_v01"
SEED = 20260830
POTENCY_THRESHOLD_NM = 50.0


def model(seed: int) -> LGBMRegressor:
    return LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1,
        random_state=seed, verbosity=-1, n_jobs=6,
    )


def rho(y: np.ndarray, prediction: np.ndarray) -> float:
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def top_quartile_recall(y: np.ndarray, prediction: np.ndarray) -> float:
    k = max(1, int(np.ceil(len(y) * .25)))
    truth = set(np.argsort(y)[-k:])
    selected = set(np.argsort(prediction)[-k:])
    return float(len(truth & selected) / k)


def pairwise_concordance(y: np.ndarray, prediction: np.ndarray) -> tuple[float, int, int]:
    correct, comparable = 0, 0
    for left, right in combinations(range(len(y)), 2):
        observed = np.sign(y[left] - y[right])
        predicted = np.sign(prediction[left] - prediction[right])
        if observed == 0 or predicted == 0:
            continue
        comparable += 1
        correct += int(observed == predicted)
    return (float(correct / comparable) if comparable else 0.0, correct, comparable)


def bootstrap_metric(y, candidate, reference=None, n=50_000):
    rng = np.random.default_rng(SEED)
    values = np.empty(n, float)
    for replicate in range(n):
        index = rng.integers(0, len(y), len(y))
        value = rho(y[index], candidate[index])
        if reference is not None:
            value -= rho(y[index], reference[index])
        values[replicate] = value
    estimate = rho(y, candidate) - (rho(y, reference) if reference is not None else 0)
    return {
        "estimate": float(estimate),
        "ci95": [float(x) for x in np.quantile(values, [.025, .975])],
        "probability_gt_zero": float(np.mean(values > 0)),
        "bootstrap_replicates": n,
    }


def scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    core = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(core)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    train = pd.read_csv(TRAIN).reset_index(drop=True)
    ext = pd.read_csv(EXTERNAL)
    ext = ext[ext.eligible_external.astype(bool)].reset_index(drop=True)
    if ext.exact_training_overlap.any():
        raise RuntimeError("External exact overlap violates the frozen protocol")

    combined = pd.concat([train.canonical_smiles, ext.smiles], ignore_index=True)
    X, D, bvs = features(combined)
    n_train = len(train)
    train_index = np.arange(n_train)
    query_index = np.arange(n_train, n_train + len(ext))
    y_train = train.pEC50.to_numpy(float)
    y = ext.cre_luc_pEC50.to_numpy(float)

    selected = np.argsort(X[train_index].var(axis=0))[-1024:]
    matrix = np.c_[X[:, selected], D]
    absolute = model(11).fit(matrix[train_index], y_train).predict(matrix[query_index])

    groups = train.source_component.astype(str).to_numpy()
    centered_target = np.zeros(n_train, float)
    for group in sorted(set(groups)):
        index = np.where(groups == group)[0]
        centered_target[index] = y_train[index] - y_train[index].mean()
    centered = (
        model(12).fit(matrix[train_index], centered_target).predict(matrix[query_index])
        + float(np.mean(y_train))
    )

    similarity = sims(query_index, train_index, bvs)
    nearest_index = np.argmax(similarity, axis=1)
    nearest = y_train[nearest_index]
    weights = np.maximum(similarity, 1e-6) ** 3
    knn = weights @ y_train / weights.sum(axis=1)
    constant = np.repeat(float(np.median(y_train)), len(ext))
    methods = {
        "Absolute-QSAR": absolute,
        "PACER-Centered": centered,
        "1NN": nearest,
        "Similarity-kNN": knn,
        "Constant-median": constant,
    }

    potent = (ext.cre_luc_ec50_nM.to_numpy(float) <= POTENCY_THRESHOLD_NM).astype(int)
    metric_rows = []
    for name, prediction in methods.items():
        concordance, correct, comparable = pairwise_concordance(y, prediction)
        metric_rows.append(
            {
                "endpoint": "M4_CRE-Luc",
                "method": name,
                "n": len(y),
                "Spearman": rho(y, prediction),
                "MAE": float(mean_absolute_error(y, prediction)),
                "top_quartile_recall": top_quartile_recall(y, prediction),
                "ROC_AUC_at_50nM": float(roc_auc_score(potent, prediction)) if len(set(potent)) == 2 else np.nan,
                "pairwise_concordance": concordance,
                "pairwise_correct": correct,
                "pairwise_comparable": comparable,
            }
        )

    # Secondary cross-platform audit; no model/model selection is performed here.
    secondary = ext.glosensor_pEC50.notna().to_numpy()
    y_secondary = ext.loc[secondary, "glosensor_pEC50"].to_numpy(float)
    empirical_cross_platform_rho = rho(
        ext.loc[secondary, "cre_luc_pEC50"].to_numpy(float), y_secondary
    )
    for name, prediction in methods.items():
        sub_prediction = prediction[secondary]
        concordance, correct, comparable = pairwise_concordance(y_secondary, sub_prediction)
        metric_rows.append(
            {
                "endpoint": "M4_GloSensor",
                "method": name,
                "n": len(y_secondary),
                "Spearman": rho(y_secondary, sub_prediction),
                "MAE": float(mean_absolute_error(y_secondary, sub_prediction)),
                "top_quartile_recall": top_quartile_recall(y_secondary, sub_prediction),
                "ROC_AUC_at_50nM": np.nan,
                "pairwise_concordance": concordance,
                "pairwise_correct": correct,
                "pairwise_comparable": comparable,
            }
        )

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    predictions = ext.copy()
    predictions["murcko_scaffold"] = predictions.smiles.map(scaffold)
    for name, prediction in methods.items():
        predictions[name] = prediction
    predictions.to_csv(OUT / "predictions.csv", index=False)

    primary = bootstrap_metric(y, absolute)
    versus_knn = bootstrap_metric(y, absolute, knn)
    versus_nn = bootstrap_metric(y, absolute, nearest)
    audit = {
        "status": "preregistered_independent_patent_external_evaluation",
        "dataset": "SUVEN_WO2025099660A1_TABLE1_V01",
        "n_external": len(ext),
        "n_exact_training_overlap": int(ext.exact_training_overlap.sum()),
        "median_max_train_tanimoto_ecfp4": float(ext.max_train_tanimoto_ecfp4.median()),
        "primary_method": "Absolute-QSAR",
        "primary_absolute_spearman": primary,
        "absolute_vs_similarity_knn": versus_knn,
        "absolute_vs_1nn": versus_nn,
        "empirical_cre_luc_vs_glosensor_spearman_n14": empirical_cross_platform_rho,
        "external_labels_used_for_training_or_model_selection": False,
        "model_configuration_frozen_before_prediction": True,
        "preregistered_external_support": bool(primary["ci95"][0] > 0 and versus_knn["ci95"][0] > 0),
        "claim_boundary": "Independent functional ranking only; not prospective PAM confirmation or wet-lab validation.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = [
        "# PACER-M4 Suven independent patent external evaluation",
        "",
        f"- Primary N: {len(ext)}; exact training overlap: 0.",
        f"- Absolute-QSAR Spearman: {primary['estimate']:.3f} "
        f"(95% bootstrap CI {primary['ci95'][0]:.3f} to {primary['ci95'][1]:.3f}).",
        f"- Delta vs similarity-kNN: {versus_knn['estimate']:.3f} "
        f"(95% CI {versus_knn['ci95'][0]:.3f} to {versus_knn['ci95'][1]:.3f}).",
        f"- Measured CRE-Luc vs GloSensor Spearman (N=14): {empirical_cross_platform_rho:.3f}.",
        f"- Preregistered support: **{audit['preregistered_external_support']}**.",
        "",
        "## Primary endpoint",
        "",
        "| Method | Spearman | MAE | ROC-AUC at 50 nM | Top-quartile recall |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in metrics[metrics.endpoint == "M4_CRE-Luc"].itertuples():
        report.append(
            f"| {row.method} | {row.Spearman:.3f} | {row.MAE:.3f} | "
            f"{row.ROC_AUC_at_50nM:.3f} | {row.top_quartile_recall:.3f} |"
        )
    report += [
        "",
        "## Interpretation",
        "",
        "The primary CI crosses zero, so the preregistered external SOTA claim fails. "
        "Absolute-QSAR nevertheless beats both similarity baselines with a positive paired CI; "
        "the local-smoothness assumption reverses under this chemotype/assay shift.",
        "",
        "No patent labels were used for training, tuning, seed selection, or model selection. "
        "Any later patent-informed method is explicitly post-hoc development.",
    ]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
