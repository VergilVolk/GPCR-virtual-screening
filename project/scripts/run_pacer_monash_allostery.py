"""Run the preregistered LY2033298 mechanistic allostery stress test."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import brier_score_loss, mean_absolute_error, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler

from evaluate_multidomain_sar import features as ridge_features
from run_pacer_fs_baselines import features as frozen_features


P = Path(__file__).resolve().parents[1]
HISTORY = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
DATA = P / "data" / "benchmarks" / "m4_pam_v1" / "external_monash_ly2033298" / "external_allostery.csv"
FUNCTIONAL = P / "results" / "m4_gamd_ensemble" / "validation_set_primary.csv"
OUT = P / "results" / "pacer_monash_allostery_v01"
ENDPOINTS = ("functional_pKB", "log_alpha_beta", "log_tauB")
TEST_FAMILIES = ("O-alkyl", "5-halogen", "N-alkyl", "core-variant")
RIDGE_ALPHAS = (0.1, 1.0, 10.0, 100.0)
HISTORY_WEIGHTS = (0.03, 0.1, 0.3, 1.0)
SEED = 20260830


def rho(y, prediction) -> float:
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def frozen_absolute_predictions(history: pd.DataFrame, frame: pd.DataFrame):
    combined = pd.concat([history.canonical_smiles, frame.canonical_smiles], ignore_index=True)
    X, D, bitvectors = frozen_features(combined)
    n_history = len(history)
    train = np.arange(n_history)
    query = np.arange(n_history, len(combined))
    selected = np.argsort(X[train].var(axis=0))[-1024:]
    design = np.c_[X[:, selected], D]
    estimator = LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1,
        random_state=11, verbosity=-1, n_jobs=6,
    ).fit(design[train], history.pEC50.to_numpy(float))
    absolute = estimator.predict(design[query])
    similarity = np.asarray([
        DataStructs.BulkTanimotoSimilarity(bitvectors[i], bitvectors[:n_history])
        for i in query
    ])
    nearest = history.pEC50.to_numpy(float)[np.argmax(similarity, axis=1)]
    weights = np.maximum(similarity, 1e-6) ** 3
    knn = weights @ history.pEC50.to_numpy(float) / weights.sum(axis=1)
    return absolute, nearest, knn


def conditional_design(history_X, local_X):
    history = np.c_[history_X, np.zeros_like(history_X), np.zeros(len(history_X))]
    local = np.c_[local_X, local_X, np.ones(len(local_X))]
    return history, local


def fit_single(history_X, history_y, local_X, local_y, train, query, mode, alpha, weight=0.0):
    if mode == "local":
        scaler = StandardScaler().fit(local_X[train])
        model = Ridge(alpha=alpha, solver="lsqr").fit(scaler.transform(local_X[train]), local_y[train])
        return model.predict(scaler.transform(local_X[query]))
    scaler = StandardScaler().fit(np.vstack([history_X, local_X[train]]))
    history_scaled = scaler.transform(history_X)
    train_scaled = scaler.transform(local_X[train])
    query_scaled = scaler.transform(local_X[query])
    if mode == "conditional":
        dh, dt = conditional_design(history_scaled, train_scaled)
        _, dq = conditional_design(np.empty((0, history_X.shape[1])), query_scaled)
    else:
        dh, dt, dq = history_scaled, train_scaled, query_scaled
    design = np.vstack([dh, dt])
    target = np.r_[history_y, local_y[train]]
    sample_weight = np.r_[np.repeat(weight, len(history_y)), np.ones(len(train))]
    model = Ridge(alpha=alpha, solver="lsqr").fit(design, target, sample_weight=sample_weight)
    return model.predict(dq)


def inner_splits(train, groups):
    unique = np.unique(groups[train])
    for held in unique:
        validation = train[groups[train] == held]
        fitting = train[groups[train] != held]
        if len(fitting) >= 3 and len(validation):
            yield fitting, validation


def select_single(history_X, history_y, local_X, local_y, train, groups, mode):
    candidates = [(a, 0.0) for a in RIDGE_ALPHAS] if mode == "local" else [
        (a, w) for a in RIDGE_ALPHAS for w in HISTORY_WEIGHTS
    ]
    best, best_loss = candidates[0], np.inf
    for alpha, weight in candidates:
        observed, predicted = [], []
        for fitting, validation in inner_splits(train, groups):
            pred = fit_single(
                history_X, history_y, local_X, local_y, fitting, validation,
                mode, alpha, weight,
            )
            observed.extend(local_y[validation]); predicted.extend(pred)
        loss = mean_absolute_error(observed, predicted) if observed else np.inf
        if loss < best_loss:
            best, best_loss = (alpha, weight), loss
    return best


def fit_pls(local_X, historical_prior, targets, train, query, components):
    design = np.c_[local_X, historical_prior]
    x_scaler = StandardScaler().fit(design[train])
    x_train = x_scaler.transform(design[train])
    x_query = x_scaler.transform(design[query])
    components = min(components, len(train) - 1, targets.shape[1])
    model = PLSRegression(n_components=max(1, components), scale=True, max_iter=1000)
    model.fit(x_train, targets[train])
    return model.predict(x_query)


def select_pls(local_X, historical_prior, targets, train, groups):
    scale = np.std(targets[train], axis=0)
    scale[scale < 1e-8] = 1.0
    best, best_loss = 1, np.inf
    for components in (1, 2, 3):
        losses = []
        for fitting, validation in inner_splits(train, groups):
            prediction = fit_pls(local_X, historical_prior, targets, fitting, validation, components)
            losses.extend(np.mean(np.abs(prediction - targets[validation]) / scale, axis=1))
        loss = float(np.mean(losses)) if losses else np.inf
        if loss < best_loss:
            best, best_loss = components, loss
    return best


def local_knn(local_smiles, target, train, query):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [generator.GetFingerprint(Chem.MolFromSmiles(item)) for item in local_smiles]
    output = []
    for index in query:
        similarity = np.asarray(DataStructs.BulkTanimotoSimilarity(fps[index], [fps[i] for i in train]))
        top = np.argsort(similarity)[-min(3, len(train)):]
        weights = np.maximum(similarity[top], 1e-6) ** 3
        output.append(float(np.average(target[train[top]], weights=weights)))
    return np.asarray(output)


def family_bootstrap_delta(predictions, candidate, reference, n=20000):
    rng = np.random.default_rng(SEED)
    families = predictions.modification_family.unique()
    values = []
    for _ in range(n):
        selected = rng.choice(families, len(families), replace=True)
        blocks = [predictions.index[predictions.modification_family == item].to_numpy() for item in selected]
        index = np.concatenate(blocks)
        values.append(rho(predictions.loc[index, "observed"], predictions.loc[index, candidate]) -
                      rho(predictions.loc[index, "observed"], predictions.loc[index, reference]))
    values = np.asarray(values)
    estimate = rho(predictions.observed, predictions[candidate]) - rho(predictions.observed, predictions[reference])
    return {"estimate": float(estimate), "ci95": [float(x) for x in np.quantile(values, [.025, .975])],
            "probability_gt_zero": float(np.mean(values > 0)), "n_boot": n}


def zero_shot(frame, absolute, nearest, similarity_knn):
    rows = []
    eligible = (frame.pam_label == 1) & (~frame.exact_history_overlap)
    methods = {"Historical-Absolute": absolute, "Historical-1NN": nearest,
               "Historical-Similarity-kNN": similarity_knn}
    predictions = frame[["compound_id", "modification_family", *ENDPOINTS,
                         "exact_history_overlap"]].copy()
    for name, values in methods.items():
        predictions[name] = values
    for endpoint in ENDPOINTS:
        mask = eligible & frame[endpoint].notna()
        for name, values in methods.items():
            rows.append({"task": "zero_shot", "endpoint": endpoint, "method": name,
                         "n": int(mask.sum()), "Spearman": rho(frame.loc[mask, endpoint], values[mask])})
    return pd.DataFrame(rows), predictions


def leave_family(frame, history, absolute):
    active = frame[(frame.pam_label == 1)].reset_index(drop=True)
    local_X = ridge_features(active.canonical_smiles)
    history_X = ridge_features(history.canonical_smiles)
    history_y = history.pEC50.to_numpy(float)
    targets = active[list(ENDPOINTS)].to_numpy(float)
    groups = active.modification_family.to_numpy(str)
    prior = absolute[frame.pam_label.to_numpy(bool)]
    prediction = {endpoint: {name: np.full(len(active), np.nan) for name in
                  ("Local-Ridge", "Pooled-Ridge", "PACER-AssayConditional", "Local-3NN", "PACER-Mechanism-PLS")}
                  for endpoint in ENDPOINTS}
    fold_rows = []
    for family in TEST_FAMILIES:
        query = np.where(groups == family)[0]
        train = np.where(groups != family)[0]
        if not len(query):
            continue
        pls_components = select_pls(local_X, prior, targets, train, groups)
        pls_prediction = fit_pls(local_X, prior, targets, train, query, pls_components)
        for endpoint_index, endpoint in enumerate(ENDPOINTS):
            y = targets[:, endpoint_index]
            for mode, name in (("local", "Local-Ridge"), ("pooled", "Pooled-Ridge"),
                               ("conditional", "PACER-AssayConditional")):
                alpha, weight = select_single(history_X, history_y, local_X, y, train, groups, mode)
                prediction[endpoint][name][query] = fit_single(
                    history_X, history_y, local_X, y, train, query, mode, alpha, weight
                )
                fold_rows.append({"family": family, "endpoint": endpoint, "method": name,
                                  "alpha": alpha, "history_weight": weight})
            prediction[endpoint]["Local-3NN"][query] = local_knn(active.canonical_smiles, y, train, query)
            prediction[endpoint]["PACER-Mechanism-PLS"][query] = pls_prediction[:, endpoint_index]
            fold_rows.append({"family": family, "endpoint": endpoint, "method": "PACER-Mechanism-PLS",
                              "pls_components": pls_components})

    rows, metric_rows = [], []
    eligible = active.modification_family.isin(TEST_FAMILIES) & (~active.exact_history_overlap)
    for endpoint in ENDPOINTS:
        for index in np.where(eligible)[0]:
            row = {"compound_id": active.loc[index, "compound_id"],
                   "modification_family": groups[index], "endpoint": endpoint,
                   "observed": targets[index, ENDPOINTS.index(endpoint)]}
            for name, values in prediction[endpoint].items():
                row[name] = values[index]
            rows.append(row)
        endpoint_rows = pd.DataFrame([item for item in rows if item["endpoint"] == endpoint])
        for name in prediction[endpoint]:
            family_rhos = []
            for _, block in endpoint_rows.groupby("modification_family"):
                if len(block) >= 3:
                    family_rhos.append(rho(block.observed, block[name]))
            metric_rows.append({"task": "leave_family_out", "endpoint": endpoint, "method": name,
                                "n": len(endpoint_rows), "Spearman": rho(endpoint_rows.observed, endpoint_rows[name]),
                                "MAE": float(mean_absolute_error(endpoint_rows.observed, endpoint_rows[name])),
                                "worst_family_Spearman_n_ge_3": min(family_rhos) if family_rhos else np.nan,
                                "prediction_range": float(np.ptp(endpoint_rows[name]))})
    return pd.DataFrame(metric_rows), pd.DataFrame(rows), pd.DataFrame(fold_rows)


def activity_cliff_classification(frame):
    train = pd.read_csv(FUNCTIONAL)
    combined = pd.concat([train.canonical_smiles, frame.canonical_smiles], ignore_index=True)
    X = ridge_features(combined)
    n = len(train)
    scaler = StandardScaler().fit(X[:n])
    train_X, query_X = scaler.transform(X[:n]), scaler.transform(X[n:])
    y = train.target.to_numpy(int)
    methods = {
        "A-tier-Logistic": LogisticRegression(C=0.1, class_weight="balanced", max_iter=5000,
                                               random_state=SEED).fit(train_X, y).predict_proba(query_X)[:, 1],
        "A-tier-ExtraTrees": ExtraTreesClassifier(n_estimators=1000, min_samples_leaf=2,
                                                   class_weight="balanced", random_state=SEED,
                                                   n_jobs=-1).fit(X[:n], y).predict_proba(X[n:])[:, 1],
    }
    output = frame[["compound_id", "modification_family", "pam_label", "exact_history_overlap"]].copy()
    for name, values in methods.items():
        output[name] = values
    metrics = []
    for name, values in methods.items():
        metrics.append({"task": "zero_shot_active_inactive", "method": name, "n": len(frame),
                        "n_inactive": int((frame.pam_label == 0).sum()),
                        "ROC_AUC_descriptive_only": float(roc_auc_score(frame.pam_label, values)),
                        "Brier": float(brier_score_loss(frame.pam_label, values)),
                        "mean_probability_active": float(np.mean(values[frame.pam_label == 1])),
                        "mean_probability_inactive": float(np.mean(values[frame.pam_label == 0]))})
    return pd.DataFrame(metrics), output


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(DATA)
    history = pd.read_csv(HISTORY)
    absolute, nearest, similarity_knn = frozen_absolute_predictions(history, frame)
    zero_metrics, zero_predictions = zero_shot(frame, absolute, nearest, similarity_knn)
    family_metrics, family_predictions, folds = leave_family(frame, history, absolute)
    cliff_metrics, cliff_predictions = activity_cliff_classification(frame)
    metrics = pd.concat([zero_metrics, family_metrics, cliff_metrics], ignore_index=True, sort=False)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    zero_predictions.to_csv(OUT / "zero_shot_predictions.csv", index=False)
    family_predictions.to_csv(OUT / "leave_family_predictions.csv", index=False)
    folds.to_csv(OUT / "folds.csv", index=False)
    cliff_predictions.to_csv(OUT / "activity_cliff_predictions.csv", index=False)

    primary = family_predictions[family_predictions.endpoint == "log_alpha_beta"].reset_index(drop=True)
    primary_metrics = family_metrics[family_metrics.endpoint == "log_alpha_beta"].sort_values("Spearman", ascending=False)
    best_baseline = primary_metrics[primary_metrics.method != "PACER-Mechanism-PLS"].iloc[0].method
    delta = family_bootstrap_delta(primary, "PACER-Mechanism-PLS", best_baseline)
    mechanism_row = primary_metrics[primary_metrics.method == "PACER-Mechanism-PLS"].iloc[0]
    baseline_row = primary_metrics[primary_metrics.method == best_baseline].iloc[0]
    supported = bool(delta["ci95"][0] > 0 and
                     mechanism_row.worst_family_Spearman_n_ge_3 >= baseline_row.worst_family_Spearman_n_ge_3)
    audit = {
        "status": "preregistered_external_endpoint_transfer_and_mechanistic_stress_test",
        "primary_endpoint": "log_alpha_beta",
        "primary_method": "PACER-Mechanism-PLS",
        "best_preregistered_baseline": best_baseline,
        "paired_family_bootstrap_delta": delta,
        "mechanism_support": supported,
        "external_labels_used_for_zero_shot_models": False,
        "external_labels_used_within_leave_family_training_only": True,
        "n_strict_continuous_queries": int(len(primary)),
        "n_exact_overlap_excluded_from_metrics": int(frame.exact_history_overlap.sum()),
        "claim_boundary": "Small same-scaffold endpoint-transfer stress test. It can test mechanism decomposition, not establish unseen-chemotype SOTA or confirm a PAM.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = [
        "# PACER-M4 Monash LY2033298 mechanistic stress test", "",
        f"- Strict continuous queries: {len(primary)}; one exact historical overlap excluded.",
        f"- Primary endpoint: log(alpha*beta); primary method: PACER-Mechanism-PLS.",
        f"- Best baseline: {best_baseline}.",
        f"- Paired family-bootstrap delta: {delta['estimate']:+.3f} "
        f"(95% CI {delta['ci95'][0]:+.3f} to {delta['ci95'][1]:+.3f}).",
        f"- Preregistered mechanism support: **{supported}**.", "",
        "## Leave-family-out continuous endpoints", "",
        "| Endpoint | Method | N | Spearman | MAE | Worst family (N>=3) |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in family_metrics.sort_values(["endpoint", "Spearman"], ascending=[True, False]).itertuples():
        report.append(f"| {row.endpoint} | {row.method} | {int(row.n)} | {row.Spearman:.3f} | "
                      f"{row.MAE:.3f} | {row.worst_family_Spearman_n_ge_3:.3f} |")
    report += ["", "## Activity-cliff classification", "",
               "Only three inactive N-acyl analogues exist; AUC is descriptive and is not a SOTA claim.", "",
               "| Method | AUC (descriptive) | Brier | P(active), active | P(active), inactive |",
               "|---|---:|---:|---:|---:|"]
    for row in cliff_metrics.itertuples():
        report.append(f"| {row.method} | {row.ROC_AUC_descriptive_only:.3f} | {row.Brier:.3f} | "
                      f"{row.mean_probability_active:.3f} | {row.mean_probability_inactive:.3f} |")
    report += ["", "## Claim boundary", "",
               "The benchmark is an independent functional-endpoint transfer test but remains in the LY2033298 scaffold neighborhood. "
               "No score in this report confirms a compound as a PAM; confirmation requires an ACh concentration-response matrix and operational-model fitting."]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
