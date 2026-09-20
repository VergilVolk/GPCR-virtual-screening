"""Nested developmental evaluation of factorized SAR on the Suven matrix.

This is explicitly post-hoc algorithm development.  The independent zero-shot
result in pacer_external_suven_v01 remains immutable and primary.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold, LeaveOneGroupOut, LeaveOneOut
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from analyze_suven_factor_sar import core_state, headgroup


P = Path(__file__).resolve().parents[1]
DATA = P / "results" / "pacer_external_suven_v01" / "predictions.csv"
OUT = P / "results" / "pacer_factor_sar_nested_v01"
ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0)
SEED = 20260830


def rho(y, prediction):
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def molecular_features(smiles: pd.Series):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mols = [Chem.MolFromSmiles(item) for item in smiles]
    fingerprints = np.asarray([generator.GetFingerprintAsNumPy(mol) for mol in mols], float)
    descriptors = np.asarray([
        [
            Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol),
            Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
            Descriptors.NumRotatableBonds(mol), Descriptors.RingCount(mol),
            Descriptors.FractionCSP3(mol), Descriptors.HeavyAtomCount(mol),
        ] for mol in mols
    ], float)
    return np.c_[fingerprints, descriptors], [generator.GetFingerprint(mol) for mol in mols]


def choose_alpha(X, y):
    if len(y) < 6:
        return 1.0
    splitter = KFold(n_splits=min(5, len(y)), shuffle=True, random_state=SEED)
    best_alpha, best_loss = None, np.inf
    for alpha in ALPHAS:
        predictions = np.full(len(y), np.nan)
        for train, query in splitter.split(X):
            scaler = StandardScaler(with_mean=False).fit(X[train])
            predictions[query] = Ridge(alpha=alpha).fit(
                scaler.transform(X[train]), y[train]
            ).predict(scaler.transform(X[query]))
        loss = mean_absolute_error(y, predictions)
        if loss < best_loss:
            best_alpha, best_loss = alpha, loss
    return float(best_alpha)


def ridge_predict(X, y, train, query, residual_base=None):
    target = y[train] if residual_base is None else y[train] - residual_base[train]
    alpha = choose_alpha(X[train], target)
    scaler = StandardScaler(with_mean=False).fit(X[train])
    prediction = Ridge(alpha=alpha).fit(
        scaler.transform(X[train]), target
    ).predict(scaler.transform(X[query]))
    if residual_base is not None:
        prediction += residual_base[query]
    return prediction, alpha


def nearest_predict(fp, y, train, query):
    result = []
    for index in query:
        similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(fp[index], [fp[i] for i in train]))
        result.append(y[train[int(np.argmax(similarities))]])
    return np.asarray(result, float)


def evaluate_protocol(frame, factor_X, molecular_X, fp, protocol, groups=None):
    y = frame.cre_luc_pEC50.to_numpy(float)
    historical = frame["Absolute-QSAR"].to_numpy(float)
    if groups is None:
        splits = LeaveOneOut().split(frame)
    else:
        splits = LeaveOneGroupOut().split(frame, y, groups)
    methods = {
        "Historical-Absolute": historical.copy(),
        "Local-1NN": np.full(len(frame), np.nan),
        "Molecular-Ridge": np.full(len(frame), np.nan),
        "PACER-MolecularResidual": np.full(len(frame), np.nan),
        "FactorSAR": np.full(len(frame), np.nan),
        "PACER-FactorResidual": np.full(len(frame), np.nan),
    }
    fold_rows = []
    for fold, (train, query) in enumerate(splits):
        methods["Local-1NN"][query] = nearest_predict(fp, y, train, query)
        molecular, molecular_alpha = ridge_predict(molecular_X, y, train, query)
        molecular_residual, molecular_residual_alpha = ridge_predict(
            molecular_X, y, train, query, historical
        )
        factor, factor_alpha = ridge_predict(factor_X, y, train, query)
        residual, residual_alpha = ridge_predict(factor_X, y, train, query, historical)
        methods["Molecular-Ridge"][query] = molecular
        methods["PACER-MolecularResidual"][query] = molecular_residual
        methods["FactorSAR"][query] = factor
        methods["PACER-FactorResidual"][query] = residual
        fold_rows.append({
            "protocol": protocol, "fold": fold, "n_train": len(train), "n_query": len(query),
            "query_ids": ";".join(frame.iloc[query].compound_id.astype(str)),
            "molecular_alpha": molecular_alpha, "factor_alpha": factor_alpha,
            "molecular_residual_alpha": molecular_residual_alpha,
            "residual_alpha": residual_alpha,
        })
    rows = []
    for name, prediction in methods.items():
        rows.append({
            "protocol": protocol, "method": name, "n": len(y),
            "Spearman": rho(y, prediction),
            "MAE": float(mean_absolute_error(y, prediction)),
            "prediction_range": float(np.ptp(prediction)),
        })
    prediction_frame = frame[["compound_id", "cre_luc_ec50_nM", "cre_luc_pEC50"]].copy()
    prediction_frame["protocol"] = protocol
    for name, prediction in methods.items():
        prediction_frame[name] = prediction
    return rows, fold_rows, prediction_frame


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(DATA)
    frame["headgroup"] = frame.name.map(headgroup)
    states = frame.name.map(core_state)
    frame["core_family"] = [item[0] for item in states]
    frame["core_state"] = [item[1] for item in states]
    encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    factor_X = encoder.fit_transform(frame[["headgroup", "core_state"]])
    molecular_X, fp = molecular_features(frame.smiles)

    all_metrics, all_folds, all_predictions = [], [], []
    protocols = [
        ("edge_LOO", None),
        ("leave_core_state_out", frame.core_state.to_numpy()),
        ("leave_headgroup_out", frame.headgroup.to_numpy()),
        ("leave_core_family_out", frame.core_family.to_numpy()),
    ]
    for protocol, groups in protocols:
        metrics, folds, predictions = evaluate_protocol(
            frame, factor_X, molecular_X, fp, protocol, groups
        )
        all_metrics.extend(metrics)
        all_folds.extend(folds)
        all_predictions.append(predictions)
    metrics = pd.DataFrame(all_metrics)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    pd.DataFrame(all_folds).to_csv(OUT / "folds.csv", index=False)
    pd.concat(all_predictions, ignore_index=True).to_csv(OUT / "predictions.csv", index=False)

    primary = metrics[(metrics.protocol == "leave_core_state_out") &
                      (metrics.method == "Molecular-Ridge")].iloc[0]
    reference = metrics[(metrics.protocol == "leave_core_state_out") &
                        (metrics.method == "Local-1NN")].iloc[0]
    audit = {
        "status": "post_hoc_nested_developmental_benchmark",
        "primary_development_protocol": "leave_core_state_out",
        "primary_method": "Molecular-Ridge",
        "primary_spearman": float(primary.Spearman),
        "local_1nn_spearman": float(reference.Spearman),
        "delta_spearman": float(primary.Spearman - reference.Spearman),
        "patent_labels_used_for_algorithm_development": True,
        "independent_external_claim": False,
        "claim_boundary": "Developmental cross-validation only; requires a new untouched external series.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(metrics.sort_values(["protocol", "Spearman"], ascending=[True, False]).to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
