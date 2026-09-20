"""Nested assay-conditional domain adaptation for M4 PAM potency.

Historical calcium/ACh data define a shared molecular prior.  Suven CRE-Luc
data define a domain-specific residual branch.  Hyperparameters and historical
sample weight are selected only inside each outer patent training partition.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold, LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler

from analyze_suven_factor_sar import core_state, headgroup


P = Path(__file__).resolve().parents[1]
HIST = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
PATENT = P / "results" / "pacer_external_suven_v01" / "predictions.csv"
OUT = P / "results" / "pacer_multidomain_sar_v01"
ALPHAS = (0.1, 1.0, 10.0)
HISTORY_WEIGHTS = (0.03, 0.1, 0.3, 1.0)
SEED = 20260830


def features(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=512)
    rows = []
    for item in smiles:
        mol = Chem.MolFromSmiles(item)
        fingerprint = generator.GetFingerprintAsNumPy(mol).astype(float)
        descriptor = np.asarray([
            Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol),
            Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
            Descriptors.NumRotatableBonds(mol), Descriptors.RingCount(mol),
            Descriptors.FractionCSP3(mol), Descriptors.HeavyAtomCount(mol),
        ], float)
        rows.append(np.r_[fingerprint, descriptor])
    return np.asarray(rows, float)


def conditional_design(history_X, local_X):
    history = np.c_[history_X, np.zeros_like(history_X), np.zeros(len(history_X))]
    local = np.c_[local_X, local_X, np.ones(len(local_X))]
    return history, local


def fit_predict(history_X, history_y, local_X, local_y, train, query, conditional, alpha, weight):
    scaler = StandardScaler().fit(np.vstack([history_X, local_X[train]]))
    history_scaled = scaler.transform(history_X)
    train_scaled = scaler.transform(local_X[train])
    query_scaled = scaler.transform(local_X[query])
    if conditional:
        design_history, design_train = conditional_design(history_scaled, train_scaled)
        _, design_query = conditional_design(np.empty((0, history_X.shape[1])), query_scaled)
    else:
        design_history, design_train, design_query = history_scaled, train_scaled, query_scaled
    design = np.vstack([design_history, design_train])
    target = np.r_[history_y, local_y[train]]
    sample_weight = np.r_[np.repeat(weight, len(history_y)), np.ones(len(train))]
    estimator = Ridge(alpha=alpha, solver="lsqr").fit(design, target, sample_weight=sample_weight)
    return estimator.predict(design_query)


def local_fit_predict(local_X, local_y, train, query, alpha):
    scaler = StandardScaler().fit(local_X[train])
    estimator = Ridge(alpha=alpha, solver="lsqr").fit(scaler.transform(local_X[train]), local_y[train])
    return estimator.predict(scaler.transform(local_X[query]))


def inner_select(history_X, history_y, local_X, local_y, outer_train, mode):
    splitter = KFold(n_splits=min(4, len(outer_train)), shuffle=True, random_state=SEED)
    candidates = [(alpha, None) for alpha in ALPHAS] if mode == "local" else [
        (alpha, weight) for alpha in ALPHAS for weight in HISTORY_WEIGHTS
    ]
    best, best_mae = None, np.inf
    for alpha, weight in candidates:
        observed, predicted = [], []
        for inner_train_position, inner_query_position in splitter.split(outer_train):
            inner_train = outer_train[inner_train_position]
            inner_query = outer_train[inner_query_position]
            if mode == "local":
                value = local_fit_predict(local_X, local_y, inner_train, inner_query, alpha)
            else:
                value = fit_predict(
                    history_X, history_y, local_X, local_y, inner_train, inner_query,
                    conditional=(mode == "conditional"), alpha=alpha, weight=float(weight),
                )
            observed.extend(local_y[inner_query])
            predicted.extend(value)
        loss = mean_absolute_error(observed, predicted)
        if loss < best_mae:
            best, best_mae = (alpha, weight), loss
    return best


def rho(y, prediction):
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def evaluate(protocol, groups, history_X, history_y, local_X, frame):
    y = frame.cre_luc_pEC50.to_numpy(float)
    methods = {
        "Historical-Absolute": frame["Absolute-QSAR"].to_numpy(float).copy(),
        "Local-Molecular-Ridge": np.full(len(frame), np.nan),
        "Pooled-Ridge": np.full(len(frame), np.nan),
        "PACER-AssayConditional": np.full(len(frame), np.nan),
    }
    folds = []
    for fold, (train, query) in enumerate(LeaveOneGroupOut().split(local_X, y, groups)):
        local_choice = inner_select(history_X, history_y, local_X, y, train, "local")
        pooled_choice = inner_select(history_X, history_y, local_X, y, train, "pooled")
        conditional_choice = inner_select(history_X, history_y, local_X, y, train, "conditional")
        methods["Local-Molecular-Ridge"][query] = local_fit_predict(
            local_X, y, train, query, local_choice[0]
        )
        methods["Pooled-Ridge"][query] = fit_predict(
            history_X, history_y, local_X, y, train, query, False,
            pooled_choice[0], pooled_choice[1],
        )
        methods["PACER-AssayConditional"][query] = fit_predict(
            history_X, history_y, local_X, y, train, query, True,
            conditional_choice[0], conditional_choice[1],
        )
        folds.append({
            "protocol": protocol, "fold": fold,
            "query_ids": ";".join(frame.iloc[query].compound_id.astype(str)),
            "local_alpha": local_choice[0],
            "pooled_alpha": pooled_choice[0], "pooled_history_weight": pooled_choice[1],
            "conditional_alpha": conditional_choice[0],
            "conditional_history_weight": conditional_choice[1],
        })
    metrics = [{
        "protocol": protocol, "method": name, "n": len(y),
        "Spearman": rho(y, prediction),
        "MAE": float(mean_absolute_error(y, prediction)),
        "prediction_range": float(np.ptp(prediction)),
    } for name, prediction in methods.items()]
    predictions = frame[["compound_id", "cre_luc_ec50_nM", "cre_luc_pEC50"]].copy()
    predictions["protocol"] = protocol
    for name, prediction in methods.items():
        predictions[name] = prediction
    return metrics, folds, predictions


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    history = pd.read_csv(HIST)
    frame = pd.read_csv(PATENT)
    frame["headgroup"] = frame.name.map(headgroup)
    states = frame.name.map(core_state)
    frame["core_family"] = [item[0] for item in states]
    frame["core_state"] = [item[1] for item in states]
    history_X = features(history.canonical_smiles)
    local_X = features(frame.smiles)
    history_y = history.pEC50.to_numpy(float)

    metric_rows, fold_rows, prediction_frames = [], [], []
    for protocol, groups in [
        ("leave_core_state_out", frame.core_state.to_numpy()),
        ("leave_headgroup_out", frame.headgroup.to_numpy()),
        ("leave_core_family_out", frame.core_family.to_numpy()),
    ]:
        metrics, folds, predictions = evaluate(
            protocol, groups, history_X, history_y, local_X, frame
        )
        metric_rows.extend(metrics)
        fold_rows.extend(folds)
        prediction_frames.append(predictions)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    pd.DataFrame(fold_rows).to_csv(OUT / "folds.csv", index=False)
    pd.concat(prediction_frames, ignore_index=True).to_csv(OUT / "predictions.csv", index=False)
    primary = metrics[(metrics.protocol == "leave_core_state_out") &
                      (metrics.method == "PACER-AssayConditional")].iloc[0]
    local = metrics[(metrics.protocol == "leave_core_state_out") &
                    (metrics.method == "Local-Molecular-Ridge")].iloc[0]
    audit = {
        "status": "post_hoc_nested_multidomain_development",
        "primary_method": "PACER-AssayConditional",
        "primary_protocol": "leave_core_state_out",
        "primary_spearman": float(primary.Spearman),
        "local_ridge_spearman": float(local.Spearman),
        "delta": float(primary.Spearman - local.Spearman),
        "all_hyperparameters_selected_inside_outer_training": True,
        "independent_external_claim": False,
        "claim_boundary": "Patent-informed nested development; requires a new untouched external series.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(metrics.sort_values(["protocol", "Spearman"], ascending=[True, False]).to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
