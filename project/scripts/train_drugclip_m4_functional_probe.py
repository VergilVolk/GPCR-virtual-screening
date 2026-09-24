#!/usr/bin/env python3
"""CPU-only functional probes on frozen DrugCLIP M4 representations.

This is deliberately not called full DrugCLIP fine-tuning: the pretrained
encoders remain frozen.  The experiment asks whether a low-capacity supervised
adapter can recover PAM-function signal under predefined scaffold/source
holdouts, beyond DrugCLIP's zero-shot binding cosine.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, matthews_corrcoef, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)


def metric_block(y: np.ndarray, p: np.ndarray) -> dict:
    pred = p >= 0.5
    return {
        "n": int(len(y)),
        "n_positive": int(y.sum()),
        "n_negative": int((y == 0).sum()),
        "prevalence": float(y.mean()),
        "roc_auc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
        "ap_lift_over_prevalence": float(average_precision_score(y, p) - y.mean()),
        "balanced_accuracy_0p5": float(balanced_accuracy_score(y, pred)),
        "mcc_0p5": float(matthews_corrcoef(y, pred)),
    }


def bootstrap_auc(y: np.ndarray, p: np.ndarray, repeats: int, seed: int) -> list[float]:
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    values = []
    for _ in range(repeats):
        idx = np.concatenate((rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)))
        values.append(float(roc_auc_score(y[idx], p[idx])))
    return list(map(float, np.quantile(values, [0.025, 0.975])))


def oof_probe(x: np.ndarray, y: np.ndarray, folds: np.ndarray, c: float) -> tuple[np.ndarray, list[dict]]:
    valid_folds = sorted(int(v) for v in np.unique(folds) if int(v) >= 0)
    prediction = np.full(len(y), np.nan, dtype=float)
    audits = []
    for fold in valid_folds:
        test = folds == fold
        train = (folds >= 0) & ~test
        if len(np.unique(y[train])) != 2 or len(np.unique(y[test])) != 2:
            audits.append({"fold": fold, "status": "not_evaluable_single_class",
                           "n_train": int(train.sum()), "n_test": int(test.sum())})
            continue
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(C=c, class_weight="balanced", max_iter=5000,
                               solver="liblinear", random_state=20260924),
        )
        model.fit(x[train], y[train])
        prediction[test] = model.predict_proba(x[test])[:, 1]
        audits.append({"fold": fold, "status": "evaluated",
                       "n_train": int(train.sum()), "n_test": int(test.sum()),
                       "test_positive": int(y[test].sum()), "test_negative": int((y[test] == 0).sum())})
    return prediction, audits


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings", type=Path, required=True)
    ap.add_argument("--benchmark", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--c", type=float, default=0.1,
                    help="Preregistered inverse L2 strength; no test-fold tuning")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260924)
    args = ap.parse_args()

    archive = np.load(args.embeddings, allow_pickle=False)
    ids = list(map(str, archive["molecule_ids"]))
    molecule = archive["molecule_embeddings"].astype(np.float64)
    scores = archive["scores"].astype(np.float64).T
    pocket_ids = list(map(str, archive["pocket_ids"]))
    if molecule.shape != (len(ids), 128) or scores.shape != (len(ids), len(pocket_ids)):
        raise ValueError("Unexpected DrugCLIP archive dimensions")
    if not np.isfinite(molecule).all() or not np.isfinite(scores).all():
        raise ValueError("Non-finite DrugCLIP values")

    table = pd.read_csv(args.benchmark)
    required = {"canonical_molecule_id", "canonical_smiles", "target", "scaffold_fold", "source_fold"}
    if missing := required - set(table.columns):
        raise ValueError(f"Missing columns: {sorted(missing)}")
    lookup = {s: i for i, s in enumerate(ids)}
    smiles = table.canonical_smiles.map(canonical)
    if missing := [s for s in smiles if s not in lookup]:
        raise ValueError(f"Missing embedded molecules: {missing[:5]}")
    order = np.asarray([lookup[s] for s in smiles], dtype=int)
    molecule, scores = molecule[order], scores[order]
    y = table.target.to_numpy(np.int64)

    # Low-capacity arms. State cosines test binding compatibility; the frozen
    # molecule embedding tests whether pretrained chemical geometry transfers.
    arms = {
        "drugclip_state_cosines": scores,
        "drugclip_molecule_embedding": molecule,
        "drugclip_embedding_plus_state_cosines": np.concatenate((molecule, scores), axis=1),
    }
    protocols = {"scaffold_holdout": table.scaffold_fold.to_numpy(int),
                 "source_holdout": table.source_fold.to_numpy(int)}
    if "series_holdout_fold" in table:
        protocols["series_holdout_assigned_subset"] = table.series_holdout_fold.to_numpy(int)

    report = {
        "evidence_level": "retrospective_frozen_drugclip_functional_probe",
        "encoder_update": "none_frozen",
        "classifier": {"type": "L2_logistic_regression", "C": args.c,
                       "class_weight": "balanced", "test_fold_tuning": False},
        "benchmark": str(args.benchmark),
        "n_molecules": int(len(y)),
        "n_positive": int(y.sum()),
        "n_negative": int((y == 0).sum()),
        "pocket_ids": pocket_ids,
        "results": {},
        "claim_boundary": (
            "This evaluates a frozen DrugCLIP functional adapter retrospectively. "
            "It is not experimental PAM validation and is not full backbone fine-tuning."
        ),
    }
    predictions = table[["canonical_molecule_id", "target", "scaffold_fold", "source_fold"]].copy()
    for pidx, (protocol, folds) in enumerate(protocols.items()):
        report["results"][protocol] = {}
        for aidx, (arm, x) in enumerate(arms.items()):
            pred, audits = oof_probe(x, y, folds, args.c)
            keep = np.isfinite(pred)
            if len(np.unique(y[keep])) != 2:
                result = {"not_evaluable": "OOF predictions contain fewer than two classes"}
            else:
                result = metric_block(y[keep], pred[keep])
                result["roc_auc_stratified_bootstrap_95ci"] = bootstrap_auc(
                    y[keep], pred[keep], args.bootstrap, args.seed + 100 * pidx + aidx)
            result["fold_audit"] = audits
            result["n_oof"] = int(keep.sum())
            report["results"][protocol][arm] = result
            predictions[f"{protocol}__{arm}"] = pred

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    predictions.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
