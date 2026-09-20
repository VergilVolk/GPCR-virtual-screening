# -*- coding: utf-8 -*-
"""Leakage-resistant PACER-M4 baseline suite on frozen benchmark folds.

Baselines:
  1. Tanimoto kNN: local chemical-neighborhood baseline
  2. Random Forest: ECFP4 + physicochemical descriptors
  3. Extra Trees: strong nonlinear ligand-only baseline

The script consumes preassigned benchmark folds and never optimizes on test
folds. It writes molecule-level OOF predictions, fold metrics and aggregate
bootstrap confidence intervals. These are the minimum baselines every PACER-M4
module must beat.
"""
from __future__ import annotations

import json
import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    matthews_corrcoef,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)

RDLogger.DisableLog("rdApp.*")
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PROJECT = Path(__file__).resolve().parents[1]
BENCH = PROJECT / "data" / "benchmarks" / "m4_pam_v1"
OUT = PROJECT / "results" / "pacer_baselines_v1"
SEED = 42
N_BOOT = 1000

DESC = [
    Descriptors.MolWt, Descriptors.MolLogP, Descriptors.NumHDonors,
    Descriptors.NumHAcceptors, Descriptors.TPSA,
    Descriptors.NumRotatableBonds, Descriptors.NumAromaticRings,
    Descriptors.FractionCSP3, Descriptors.HeavyAtomCount,
]


def featurize(smiles: list[str]) -> tuple[np.ndarray, list]:
    rows, fps = [], []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid benchmark SMILES: {smi}")
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros(2048, dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)
        desc = np.asarray([fn(mol) for fn in DESC], dtype=np.float32)
        rows.append(np.concatenate([arr, desc]))
        fps.append(fp)
    return np.vstack(rows), fps


def knn_predict(train_fps, test_fps, y_train: np.ndarray, k: int = 5) -> np.ndarray:
    pred = []
    for fp in test_fps:
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(fp, train_fps), dtype=float)
        top = np.argsort(sims)[-min(k, len(sims)):]
        weights = np.maximum(sims[top], 1e-6) ** 2
        pred.append(float(np.average(y_train[top], weights=weights)))
    return np.asarray(pred)


def reg_metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "MAE": float(mean_absolute_error(y, pred)),
        "RMSE": float(np.sqrt(mean_squared_error(y, pred))),
        "R2": float(r2_score(y, pred)),
        "Spearman": float(pd.Series(y).corr(pd.Series(pred), method="spearman")),
    }


def clf_metrics(y: np.ndarray, prob: np.ndarray) -> dict[str, float]:
    pred = (prob >= 0.5).astype(int)
    return {
        "ROC_AUC": float(roc_auc_score(y, prob)),
        "PR_AUC": float(average_precision_score(y, prob)),
        "MCC": float(matthews_corrcoef(y, pred)),
        "BalancedAcc": float(balanced_accuracy_score(y, pred)),
        "Brier": float(brier_score_loss(y, prob)),
    }


def bootstrap_ci(y: np.ndarray, pred: np.ndarray, metric_fn, keys: list[str]) -> dict:
    rng = np.random.RandomState(SEED)
    values = {key: [] for key in keys}
    n = len(y)
    accepted = 0
    for _ in range(N_BOOT * 2):
        idx = rng.randint(0, n, n)
        if len(np.unique(y[idx])) < 2 and "ROC_AUC" in keys:
            continue
        metrics = metric_fn(y[idx], pred[idx])
        for key in keys:
            value = metrics[key]
            if np.isfinite(value):
                values[key].append(value)
        accepted += 1
        if accepted >= N_BOOT:
            break
    return {
        key: {
            "estimate": metric_fn(y, pred)[key],
            "ci95": [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))],
        }
        for key, vals in values.items() if vals
    }


def models(classification: bool):
    if classification:
        return {
            "RandomForest": RandomForestClassifier(
                n_estimators=600, class_weight="balanced_subsample",
                min_samples_leaf=2, random_state=SEED, n_jobs=-1,
            ),
            "ExtraTrees": ExtraTreesClassifier(
                n_estimators=600, class_weight="balanced",
                min_samples_leaf=2, random_state=SEED, n_jobs=-1,
            ),
        }
    return {
        "RandomForest": RandomForestRegressor(
            n_estimators=600, min_samples_leaf=2, random_state=SEED, n_jobs=-1,
        ),
        "ExtraTrees": ExtraTreesRegressor(
            n_estimators=600, min_samples_leaf=2, random_state=SEED, n_jobs=-1,
        ),
    }


def run_view(name: str, filename: str, target: str, classification: bool,
             requested_splits: list[str] | None = None) -> dict:
    df = pd.read_csv(BENCH / filename)
    X, fps = featurize(df["canonical_smiles"].tolist())
    y = df[target].to_numpy(dtype=int if classification else float)
    report = {}
    split_specs = [("scaffold", "scaffold_fold"), ("source", "source_fold")]
    if classification and "series_holdout_fold" in df.columns:
        split_specs.append(("series", "series_holdout_fold"))
    if requested_splits:
        split_specs = [item for item in split_specs if item[0] in requested_splits]
    for split_name, fold_col in split_specs:
        evaluation_mask = df[fold_col].to_numpy() >= 0
        for model_name in ["TanimotoKNN", "RandomForest", "ExtraTrees"]:
            oof = np.full(len(df), np.nan, dtype=float)
            fold_reports = []
            for fold in sorted(x for x in df[fold_col].unique() if x >= 0):
                test = df[fold_col].to_numpy() == fold
                train = ~test
                if model_name == "TanimotoKNN":
                    train_fps = [fps[i] for i in np.where(train)[0]]
                    test_fps = [fps[i] for i in np.where(test)[0]]
                    prediction = knn_predict(train_fps, test_fps, y[train])
                else:
                    model = models(classification)[model_name]
                    model.fit(X[train], y[train])
                    prediction = (model.predict_proba(X[test])[:, 1]
                                  if classification else model.predict(X[test]))
                oof[test] = prediction
                if classification and len(np.unique(y[test])) < 2:
                    fold_metrics = {"warning": "single-class test fold; rank metrics undefined"}
                else:
                    fold_metrics = clf_metrics(y[test], prediction) if classification \
                        else reg_metrics(y[test], prediction)
                fold_reports.append({"fold": int(fold), "n": int(test.sum()), **fold_metrics})

            if np.isnan(oof[evaluation_mask]).any():
                raise RuntimeError(f"Missing OOF predictions for {name}/{split_name}/{model_name}")
            metric_fn = clf_metrics if classification else reg_metrics
            keys = (["ROC_AUC", "PR_AUC", "MCC", "BalancedAcc", "Brier"]
                    if classification else ["MAE", "RMSE", "R2", "Spearman"])
            eval_y = y[evaluation_mask]
            eval_oof = oof[evaluation_mask]
            aggregate = metric_fn(eval_y, eval_oof)
            ci = bootstrap_ci(eval_y, eval_oof, metric_fn, keys)
            key = f"{split_name}|{model_name}"
            report[key] = {"aggregate_oof": aggregate, "bootstrap": ci, "folds": fold_reports}

            pred_out = df.loc[evaluation_mask, ["canonical_molecule_id", "canonical_smiles", fold_col]].copy()
            pred_out["target"] = eval_y
            pred_out["prediction"] = eval_oof
            pred_out.to_csv(OUT / f"{name}_{split_name}_{model_name}_oof.csv", index=False)
            print(name, key, aggregate)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--views", default="potency,pam_vs_inactive,threshold_10uM")
    parser.add_argument("--splits", default="", help="Optional comma-separated split filter")
    args = parser.parse_args()
    views = {x.strip() for x in args.views.split(",") if x.strip()}
    splits = [x.strip() for x in args.splits.split(",") if x.strip()] or None
    OUT.mkdir(parents=True, exist_ok=True)
    metrics_path = OUT / "baseline_metrics.json"
    all_results = json.load(open(metrics_path, encoding="utf-8")) if metrics_path.exists() else {}
    if "potency" in views:
        all_results.setdefault("potency", {}).update(
            run_view("potency", "potency_molecules.csv", "pEC50", False, splits))
    if "pam_vs_inactive" in views:
        all_results.setdefault("pam_vs_inactive", {}).update(
            run_view("pam_vs_inactive", "pam_vs_inactive.csv", "target", True, splits))
    if "threshold_10uM" in views:
        all_results.setdefault("threshold_10uM", {}).update(
            run_view("threshold_10uM", "threshold_10uM.csv", "target", True, splits))
    with open(metrics_path, "w", encoding="utf-8") as handle:
        json.dump(all_results, handle, ensure_ascii=False, indent=2)
    print(f"Saved baseline benchmark to {OUT}")


if __name__ == "__main__":
    main()
