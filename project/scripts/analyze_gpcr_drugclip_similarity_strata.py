#!/usr/bin/env python3
"""Audit DrugCLIP GPCR retrieval by distance from the training chemistry.

The script reconstructs the ECFP OOF baseline under the frozen scaffold folds,
computes each test molecule's maximum Tanimoto similarity to its training fold,
and reports retrieval metrics in pre-defined similarity strata.  This prevents
an apparently strong average from hiding chemotype memorisation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.multiclass import OneVsRestClassifier


TARGETS = ["B2AR", "CCR2", "M2R", "M4R"]
STRATA = [
    ("remote_lt_0.30", -np.inf, 0.30),
    ("low_0.30_0.50", 0.30, 0.50),
    ("medium_0.50_0.70", 0.50, 0.70),
    ("close_ge_0.70", 0.70, np.inf),
]


def fingerprints(smiles: list[str]):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return [generator.GetFingerprint(Chem.MolFromSmiles(value)) for value in smiles]


def matrix(fps) -> np.ndarray:
    out = np.zeros((len(fps), 2048), dtype=np.float32)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, out[i])
    return out


def metrics(labels: np.ndarray, scores: np.ndarray) -> dict:
    order = np.argsort(-scores, axis=1)
    ranks = np.asarray([int(np.flatnonzero(order[i] == y)[0]) + 1 for i, y in enumerate(labels)])
    per_target = {}
    for j, target in enumerate(TARGETS):
        keep = labels == j
        per_target[target] = {
            "n": int(keep.sum()),
            "recall_at_1": None if not keep.any() else float((ranks[keep] == 1).mean()),
            "mrr": None if not keep.any() else float((1.0 / ranks[keep]).mean()),
        }
    present = [v for v in per_target.values() if v["n"]]
    one_hot = np.eye(len(TARGETS), dtype=int)[labels]
    return {
        "n": int(len(labels)),
        "recall_at_1": float((ranks == 1).mean()),
        "recall_at_2": float((ranks <= 2).mean()),
        "mrr": float((1.0 / ranks).mean()),
        "macro_recall_at_1_present_targets": float(np.mean([v["recall_at_1"] for v in present])),
        "pair_roc_auc": float(roc_auc_score(one_hot.ravel(), scores.ravel())),
        "pair_pr_auc": float(average_precision_score(one_hot.ravel(), scores.ravel())),
        "per_target": per_target,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.benchmark)
    predictions = pd.read_csv(args.predictions)
    available = set(predictions.canonical_molecule_id.astype(str))
    table = table[table.canonical_molecule_id.astype(str).isin(available)].copy().reset_index(drop=True)
    labels = np.asarray([TARGETS.index(value) for value in table.target], dtype=int)
    folds = table.scaffold_fold.to_numpy(int)
    fps = fingerprints(table.canonical_smiles.astype(str).tolist())
    x = matrix(fps)

    ecfp_scores = np.full((len(table), len(TARGETS)), np.nan)
    max_train_similarity = np.full(len(table), np.nan)
    for fold in sorted(np.unique(folds)):
        train_idx = np.flatnonzero(folds != fold)
        test_idx = np.flatnonzero(folds == fold)
        clf = OneVsRestClassifier(LogisticRegression(
            C=1.0, class_weight="balanced", max_iter=3000,
            solver="liblinear", random_state=20260925))
        clf.fit(x[train_idx], labels[train_idx])
        ecfp_scores[test_idx] = clf.predict_proba(x[test_idx])[:, np.argsort(clf.classes_)]
        train_fps = [fps[i] for i in train_idx]
        for i in test_idx:
            max_train_similarity[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], train_fps))

    grouped = predictions.groupby("canonical_molecule_id", sort=False)
    score_sets = {"ecfp_logistic": ecfp_scores}
    for mode in ["official", "random", "ce", "hard"]:
        values = grouped[[f"{mode}_{target}" for target in TARGETS]].mean()
        values = values.reindex(table.canonical_molecule_id.astype(str))
        score_sets[mode] = values.to_numpy(float)

    report = {
        "protocol": "five-fold exact-Murcko-scaffold OOF; max ECFP4 Tanimoto computed only to each test fold's training molecules",
        "predefined_similarity_strata": [name for name, _, _ in STRATA],
        "overall": {name: metrics(labels, values) for name, values in score_sets.items()},
        "strata": {},
        "claim_boundary": "Retrospective target retrieval. Similarity-stratum analysis is diagnostic, not a prospective binding or efficacy validation.",
    }
    for name, lower, upper in STRATA:
        keep = (max_train_similarity >= lower) & (max_train_similarity < upper)
        report["strata"][name] = {
            "n": int(keep.sum()),
            "similarity_min": None if not keep.any() else float(max_train_similarity[keep].min()),
            "similarity_median": None if not keep.any() else float(np.median(max_train_similarity[keep])),
            "similarity_max": None if not keep.any() else float(max_train_similarity[keep].max()),
            "metrics": {method: metrics(labels[keep], values[keep]) for method, values in score_sets.items()} if keep.any() else {},
        }

    audit = table[["canonical_molecule_id", "target", "scaffold_fold"]].copy()
    audit["max_train_ecfp4_tanimoto"] = max_train_similarity
    for method, values in score_sets.items():
        audit[f"{method}_predicted_target"] = [TARGETS[i] for i in np.argmax(values, axis=1)]
        audit[f"{method}_correct"] = np.argmax(values, axis=1) == labels
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    audit.to_csv(args.output.with_suffix(".csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
