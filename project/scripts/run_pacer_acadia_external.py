"""Frozen zero-shot and three-anchor evaluation on WO2025122811.

The model recipes are inherited from the earlier frozen VU6025733/Suven
external scripts. Patent labels are used only after feature/model definition.
PACER-ACM v2 is NOT evaluated because this source lacks log(alpha beta)/log(tauB).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.preprocessing import StandardScaler

from pacer_acm_v2 import select_three_anchors
from run_pacer_fs_baselines import adapted, features, sims, train_delta


P = Path(__file__).resolve().parents[1]
ROOT = P / "data" / "benchmarks" / "m4_pam_v1"
TRAIN = ROOT / "potency_molecules.csv"
EXTERNAL = ROOT / "external_acadia_2025" / "external_molecule_benchmark.csv"
OUT = P / "results" / "pacer_external_acadia_2025_v01"
SEED = 20260830
N_BOOTSTRAP = 20_000


def lgbm(seed: int) -> LGBMRegressor:
    return LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1,
        random_state=seed, verbosity=-1, n_jobs=6,
    )


def rho(y, prediction) -> float:
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def top_fraction_recall(y, prediction, fraction=.10) -> float:
    k = max(1, int(np.ceil(len(y) * fraction)))
    return float(len(set(np.argsort(y)[-k:]) & set(np.argsort(prediction)[-k:])) / k)


def paired_bootstrap(y, candidate, reference, metric="rho"):
    rng = np.random.default_rng(SEED)
    values = []
    for _ in range(N_BOOTSTRAP):
        index = rng.integers(0, len(y), len(y))
        if metric == "rho":
            delta = rho(y[index], candidate[index]) - rho(y[index], reference[index])
        else:
            delta = (mean_absolute_error(y[index], reference[index])
                     - mean_absolute_error(y[index], candidate[index]))
        values.append(delta)
    values = np.asarray(values)
    estimate = (rho(y, candidate) - rho(y, reference)) if metric == "rho" else (
        mean_absolute_error(y, reference) - mean_absolute_error(y, candidate))
    return {"estimate": float(estimate),
            "ci95": [float(x) for x in np.quantile(values, [.025, .975])],
            "probability_gt_zero": float(np.mean(values > 0)),
            "replicates": N_BOOTSTRAP}


def summarize(y, methods, protocol):
    return [{"protocol": protocol, "method": name, "n_query": len(y),
             "Spearman": rho(y, prediction),
             "MAE": float(mean_absolute_error(y, prediction)),
             "top10_recall": top_fraction_recall(y, prediction)}
            for name, prediction in methods.items()]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    train = pd.read_csv(TRAIN).reset_index(drop=True)
    external_all = pd.read_csv(EXTERNAL).reset_index(drop=True)
    if external_all.exact_training_overlap.any():
        raise RuntimeError("External exact overlap violates frozen protocol")

    combined = pd.concat([train.canonical_smiles, external_all.canonical_smiles], ignore_index=True)
    X, D, bitvectors = features(combined)
    n_train = len(train)
    train_index = np.arange(n_train)
    external_index = np.arange(n_train, n_train + len(external_all))
    y_train = train.pEC50.to_numpy(float)

    exact_mask = external_all.eligible_exact_potency.astype(bool).to_numpy()
    exact_local = np.where(exact_mask)[0]
    exact_index = external_index[exact_local]
    y_exact = external_all.loc[exact_mask, "pEC50"].to_numpy(float)

    selected_bits = np.argsort(X[train_index].var(axis=0))[-1024:]
    matrix = np.c_[X[:, selected_bits], D]
    scaler = StandardScaler().fit(matrix[train_index])
    scaled = scaler.transform(matrix)

    ridge = Ridge(alpha=10.0).fit(scaled[train_index], y_train).predict(scaled[exact_index])
    forest = ExtraTreesRegressor(n_estimators=500, min_samples_leaf=3, max_features="sqrt",
                                 random_state=SEED, n_jobs=6).fit(
                                     matrix[train_index], y_train).predict(matrix[exact_index])
    absolute = lgbm(11).fit(matrix[train_index], y_train).predict(matrix[exact_index])

    groups = train.source_component.astype(str).to_numpy()
    centered_target = np.zeros(n_train, float)
    for group in sorted(set(groups)):
        positions = np.where(groups == group)[0]
        centered_target[positions] = y_train[positions] - y_train[positions].mean()
    centered = lgbm(12).fit(matrix[train_index], centered_target).predict(matrix[exact_index]) + y_train.mean()

    similarity = sims(exact_index, train_index, bitvectors)
    nearest_position = np.argmax(similarity, axis=1)
    nearest = y_train[nearest_position]
    weights = np.maximum(similarity, 1e-6) ** 3
    knn = weights @ y_train / weights.sum(axis=1)
    constant = np.repeat(float(np.median(y_train)), len(exact_index))
    zero = {"Constant-median": constant, "ECFP4-1NN": nearest,
            "ECFP4-Similarity-kNN": knn, "ECFP4-Ridge": ridge,
            "ECFP4-ExtraTrees": forest, "Absolute-LGBM": absolute,
            "PACER-Centered-LGBM": centered}

    anchor_local = select_three_anchors(external_all.canonical_smiles.tolist())
    if not external_all.iloc[anchor_local].eligible_exact_potency.all():
        raise RuntimeError("A label-blind anchor lacks an exact potency curve; protocol must refuse")
    anchor_ids = external_all.iloc[anchor_local].compound_id.tolist()
    exact_id_to_position = {compound: position for position, compound in enumerate(
        external_all.loc[exact_mask, "compound_id"].tolist())}
    anchor_exact = np.asarray([exact_id_to_position[item] for item in anchor_ids], int)
    query_exact = np.asarray([i for i in range(len(y_exact)) if i not in set(anchor_exact)], int)

    three_shot = {}
    for name, prediction in zero.items():
        offset = float(np.mean(y_exact[anchor_exact] - prediction[anchor_exact]))
        three_shot[name + "+3AnchorOffset"] = prediction[query_exact] + offset

    # Pre-existing PACER-FS v0.1 adaptation methods, frozen before this patent.
    total_n = len(combined)
    y_all = np.full(total_n, np.nan, float)
    y_all[train_index] = y_train
    y_all[exact_index] = y_exact
    group_all = np.r_[groups, np.repeat("ACADIA_2025", len(external_all))]
    delta_pack = train_delta(train_index, X, D, y_all, group_all, bitvectors,
                             seed=SEED, cliff_weight=0)
    support_global = exact_index[anchor_exact]
    query_global = exact_index[query_exact]
    centered_full = np.full(total_n, np.nan); centered_full[exact_index] = centered
    absolute_full = np.full(total_n, np.nan); absolute_full[exact_index] = absolute
    for method in ("KernelResidual", "TanimotoGP", "DeltaSAR", "DeltaSARHybrid"):
        three_shot[f"PACER-FS-{method}-Centered"] = adapted(
            method, centered_full, support_global, query_global,
            y_all, X, D, bitvectors, delta_pack)
    three_shot["PACER-FS-DeltaSARHybrid-Absolute"] = adapted(
        "DeltaSARHybrid", absolute_full, support_global, query_global,
        y_all, X, D, bitvectors, delta_pack)

    rows = summarize(y_exact, zero, "zero_shot")
    rows += summarize(y_exact[query_exact], three_shot, "three_anchor")
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "metrics.csv", index=False)

    predictions = external_all.loc[exact_mask, ["compound_id", "example_id", "canonical_smiles",
                                                "pEC50", "hM4_PAM_RE_pct", "hM4_agonist_RE_pct",
                                                "max_train_tanimoto_ecfp4"]].copy()
    for name, values in zero.items():
        predictions[name] = values
    predictions["three_anchor_role"] = ["anchor" if i in set(anchor_exact) else "query"
                                         for i in range(len(predictions))]
    for name, values in three_shot.items():
        predictions[name] = np.nan
        predictions.iloc[query_exact, predictions.columns.get_loc(name)] = values
    predictions.to_csv(OUT / "predictions.csv", index=False)

    best_zero_name = metrics[metrics.protocol == "zero_shot"].sort_values("Spearman").iloc[-1].method
    best_zero = zero[best_zero_name]
    best_three_name = metrics[metrics.protocol == "three_anchor"].sort_values("MAE").iloc[0].method
    best_three = three_shot[best_three_name]
    primary_name = "PACER-FS-DeltaSARHybrid-Centered"
    primary = three_shot[primary_name]
    strongest_local_baseline_name = "ECFP4-Similarity-kNN+3AnchorOffset"
    strongest_local_baseline = three_shot[strongest_local_baseline_name]
    if "+3AnchorOffset" in best_three_name:
        corresponding = zero[best_three_name.replace("+3AnchorOffset", "")][query_exact]
    elif best_three_name.endswith("-Centered"):
        corresponding = centered[query_exact]
    else:
        corresponding = absolute[query_exact]
    audit = {
        "status": "preregistered_external_functional_evaluation",
        "dataset": "M4_ACADIA_WO2025122811_FUNCTIONAL_HIGHCONF_V01",
        "mechanistic_pacer_acm_v2_evaluated": False,
        "mechanistic_refusal_reason": "No log(alpha beta) or log(tauB) labels in source.",
        "n_unique_active_moieties": int(len(external_all)),
        "n_exact_potency": int(len(y_exact)),
        "n_query_after_anchors": int(len(query_exact)),
        "anchor_ids": anchor_ids,
        "anchor_selection": "ECFP4 medoid plus deterministic max-min; labels hidden",
        "best_zero_shot_method_by_spearman": best_zero_name,
        "best_three_anchor_method_by_mae": best_three_name,
        "primary_three_anchor_method": primary_name,
        "primary_vs_similarity_knn_spearman": paired_bootstrap(
            y_exact[query_exact], primary, strongest_local_baseline, "rho"),
        "primary_vs_similarity_knn_mae": paired_bootstrap(
            y_exact[query_exact], primary, strongest_local_baseline, "mae"),
        "best_zero_vs_similarity_knn_spearman": paired_bootstrap(
            y_exact, best_zero, zero["ECFP4-Similarity-kNN"], "rho"),
        "best_three_vs_corresponding_zero_mae": paired_bootstrap(
            y_exact[query_exact], best_three, corresponding, "mae"),
        "external_labels_used_for_training_or_model_selection": False,
        "model_recipes_inherited_from_frozen_external_scripts": True,
        "claim_boundary": "External functional potency ranking only; no mechanistic decomposition, generated-PAM confirmation, or universal SOTA.",
    }
    audit["primary_external_support"] = bool(
        audit["primary_vs_similarity_knn_spearman"]["ci95"][0] > 0
        and audit["primary_vs_similarity_knn_mae"]["ci95"][0] > 0
    )
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")

    report = ["# WO2025122811 frozen external evaluation", "",
              f"Unique high-confidence active moieties: {len(external_all)}; exact pEC50: {len(y_exact)}; exact training overlap: 0.",
              f"Label-blind anchors: {', '.join(anchor_ids)}. PACER-ACM v2 refused because no operational-allostery endpoints are reported.", "",
              "| Protocol | Method | N | Spearman | MAE | Top-10% recall |",
              "|---|---|---:|---:|---:|---:|"]
    for row in metrics.sort_values(["protocol", "Spearman"], ascending=[True, False]).itertuples():
        report.append(f"| {row.protocol} | {row.method} | {row.n_query} | {row.Spearman:.3f} | {row.MAE:.3f} | {row.top10_recall:.3f} |")
    report += ["", "## Interpretation", "",
               "Three-anchor constant offsets cannot change rank and are evaluated only for assay-scale MAE calibration. "
               "The external campaign is a genuine low-similarity chemotype test, but a single patent remains one assay domain.", "",
               "No result here confirms a generated molecule as an M4 PAM."]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(metrics.sort_values(["protocol", "Spearman"], ascending=[True, False]).to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
