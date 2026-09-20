# -*- coding: utf-8 -*-
"""Differential SAR baselines for PACER-M4 activity-cliff prediction.

Compares:
  - AbsoluteQSAR: delta inferred from molecule-level RandomForest OOF predictions
  - DeltaRidge: antisymmetric linear model on molecular feature differences
  - DeltaRF: nonlinear context-aware pair model

For every fold, pair models are trained only on pairs whose two molecules are
both in the molecule-level training partition. Test metrics use only pairs whose
two molecules are both in the held-out partition.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(__file__).resolve().parents[1]
BENCH = PROJECT / "data" / "benchmarks" / "m4_pam_v1"
PAIRS = BENCH / "activity_cliffs" / "near_neighbor_pairs.csv"
BASELINES = PROJECT / "results" / "pacer_baselines_v1"
OUT = PROJECT / "results" / "pacer_delta_sar_v1"
SEED = 42

DESC = [
    Descriptors.MolWt, Descriptors.MolLogP, Descriptors.NumHDonors,
    Descriptors.NumHAcceptors, Descriptors.TPSA,
    Descriptors.NumRotatableBonds, Descriptors.NumAromaticRings,
    Descriptors.FractionCSP3, Descriptors.HeavyAtomCount,
]


def features(smiles: list[str]) -> np.ndarray:
    output = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros(2048, dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)
        desc = np.asarray([fn(mol) for fn in DESC], dtype=np.float32)
        output.append(np.concatenate([arr, desc]))
    return np.vstack(output)


def pair_features(pair_frame: pd.DataFrame, matrix: np.ndarray,
                  id_to_idx: dict[str, int], context: bool) -> np.ndarray:
    ia = np.asarray([id_to_idx[x] for x in pair_frame["mol_a"]], dtype=int)
    ib = np.asarray([id_to_idx[x] for x in pair_frame["mol_b"]], dtype=int)
    delta = matrix[ib] - matrix[ia]
    if not context:
        return delta
    midpoint = (matrix[ib] + matrix[ia]) * 0.5
    return np.hstack([delta, midpoint])


def metrics(frame: pd.DataFrame, prediction: np.ndarray) -> dict:
    true = frame["delta_pEC50_b_minus_a"].to_numpy(float)
    abs_true = np.abs(true)
    evaluable = abs_true >= 0.30
    cliffs = frame["is_cliff"].to_numpy(int) == 1
    result = {
        "n_pairs": int(len(frame)),
        "n_cliffs": int(cliffs.sum()),
        "delta_MAE": float(mean_absolute_error(true, prediction)),
        "delta_Spearman": float(pd.Series(true).corr(pd.Series(prediction), method="spearman")),
        "direction_accuracy_delta_ge_0.3": float(
            np.mean(np.sign(true[evaluable]) == np.sign(prediction[evaluable]))
        ) if evaluable.any() else float("nan"),
        "cliff_direction_accuracy": float(
            np.mean(np.sign(true[cliffs]) == np.sign(prediction[cliffs]))
        ) if cliffs.any() else float("nan"),
        "cliff_delta_MAE": float(mean_absolute_error(true[cliffs], prediction[cliffs]))
        if cliffs.any() else float("nan"),
    }
    binary = frame["pair_class"].isin(["cliff", "smooth"]).to_numpy()
    labels = frame.loc[binary, "is_cliff"].to_numpy(int)
    if len(np.unique(labels)) == 2:
        result["cliff_vs_smooth_ROC_AUC"] = float(
            roc_auc_score(labels, np.abs(prediction[binary]))
        )
    else:
        result["cliff_vs_smooth_ROC_AUC"] = float("nan")
    return result


def absolute_qsar_delta(split: str, pairs: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    pred = pd.read_csv(BASELINES / f"potency_{split}_RandomForest_oof.csv")
    pred_map = pred.set_index("canonical_molecule_id")["prediction"].to_dict()
    fold_col = f"{split}_pair_fold"
    eligible = pairs[pairs[fold_col] >= 0].copy()
    delta = np.asarray([
        pred_map[b] - pred_map[a] for a, b in zip(eligible["mol_a"], eligible["mol_b"])
    ])
    return eligible, delta


def run_pair_model(split: str, model_name: str, mols: pd.DataFrame,
                   pairs: pd.DataFrame, X: np.ndarray, id_to_idx: dict[str, int]) -> tuple[pd.DataFrame, np.ndarray, list]:
    mol_fold_col = f"{split}_fold"
    pair_fold_col = f"{split}_pair_fold"
    outputs = []
    fold_reports = []
    for fold in sorted(mols[mol_fold_col].unique()):
        train_ids = set(mols.loc[mols[mol_fold_col] != fold, "canonical_molecule_id"])
        train_pairs = pairs[pairs["mol_a"].isin(train_ids) & pairs["mol_b"].isin(train_ids)].copy()
        test_pairs = pairs[pairs[pair_fold_col] == fold].copy()
        if test_pairs.empty:
            continue

        # Add reversed orientation to force the training set to expose antisymmetry.
        reverse = train_pairs.copy()
        reverse[["mol_a", "mol_b"]] = reverse[["mol_b", "mol_a"]].to_numpy()
        reverse["delta_pEC50_b_minus_a"] *= -1
        train_aug = pd.concat([train_pairs, reverse], ignore_index=True)

        if model_name == "DeltaRidge":
            train_X = pair_features(train_aug, X, id_to_idx, context=False)
            test_X = pair_features(test_pairs, X, id_to_idx, context=False)
            model = make_pipeline(StandardScaler(), Ridge(alpha=10.0))
        elif model_name == "DeltaRF":
            train_X = pair_features(train_aug, X, id_to_idx, context=True)
            test_X = pair_features(test_pairs, X, id_to_idx, context=True)
            model = RandomForestRegressor(
                n_estimators=500, min_samples_leaf=2, max_features=0.35,
                random_state=SEED, n_jobs=-1,
            )
        else:
            raise ValueError(model_name)
        model.fit(train_X, train_aug["delta_pEC50_b_minus_a"].to_numpy(float))
        prediction = model.predict(test_X)
        test_pairs = test_pairs.copy()
        test_pairs["prediction"] = prediction
        outputs.append(test_pairs)
        fold_reports.append({"fold": int(fold), "n_train_pairs": int(len(train_pairs)),
                             **metrics(test_pairs, prediction)})
    combined = pd.concat(outputs, ignore_index=True)
    return combined, combined["prediction"].to_numpy(float), fold_reports


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mols = pd.read_csv(BENCH / "potency_molecules.csv")
    pairs = pd.read_csv(PAIRS)
    X = features(mols["canonical_smiles"].tolist())
    id_to_idx = {mid: i for i, mid in enumerate(mols["canonical_molecule_id"])}
    report = {}

    for split in ["scaffold", "source"]:
        eligible, absolute_pred = absolute_qsar_delta(split, pairs)
        report[f"{split}|AbsoluteQSAR"] = {
            "aggregate": metrics(eligible, absolute_pred), "folds": []
        }
        absolute_out = eligible.copy()
        absolute_out["prediction"] = absolute_pred
        absolute_out.to_csv(OUT / f"{split}_AbsoluteQSAR_pairs.csv", index=False)
        print(split, "AbsoluteQSAR", report[f"{split}|AbsoluteQSAR"]["aggregate"])

        for model_name in ["DeltaRidge", "DeltaRF"]:
            combined, prediction, folds = run_pair_model(
                split, model_name, mols, pairs, X, id_to_idx
            )
            report[f"{split}|{model_name}"] = {
                "aggregate": metrics(combined, prediction), "folds": folds
            }
            combined.to_csv(OUT / f"{split}_{model_name}_pairs.csv", index=False)
            print(split, model_name, report[f"{split}|{model_name}"]["aggregate"])

    with open(OUT / "delta_sar_metrics.json", "w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
    print(f"Saved differential SAR baselines to {OUT}")


if __name__ == "__main__":
    main()

