#!/usr/bin/env python3
"""Leakage-audited 2D baselines for PubChem CHRM4 PAM AID 624126."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def matrix(smiles: list[str], fpgen) -> tuple[np.ndarray, list]:
    array = np.zeros((len(smiles), 2048), dtype=np.uint8)
    fps = []
    for i, value in enumerate(smiles):
        fp = fpgen.GetFingerprint(Chem.MolFromSmiles(value)); fps.append(fp)
        DataStructs.ConvertToNumpyArray(fp, array[i])
    return array, fps


def metrics(y: np.ndarray, score: np.ndarray, inactive_weight: float = 1.0) -> dict:
    weights = np.where(y == 0, inactive_weight, 1.0)
    order = np.argsort(-score)
    out = {
        "roc_auc": float(roc_auc_score(y, score)),
        "average_precision_sampled": float(average_precision_score(y, score)),
        "average_precision_prevalence_corrected": float(average_precision_score(y, score, sample_weight=weights)),
    }
    prevalence = float(y.mean())
    for fraction in (0.005, 0.01, 0.05):
        k = max(1, int(np.ceil(len(y) * fraction)))
        out[f"EF@{100*fraction:g}pct"] = float(y[order[:k]].mean() / prevalence)
        out[f"hits@{100*fraction:g}pct"] = int(y[order[:k]].sum())
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--screen", type=Path, required=True)
    p.add_argument("--training", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    screen = pd.read_csv(a.screen)
    train = pd.read_csv(a.training)
    train_set = set(train.canonical_smiles.astype(str))
    screen["exact_training_overlap"] = screen.canonical_smiles.isin(train_set)
    test = screen[~screen.exact_training_overlap].copy().reset_index(drop=True)
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    x_train, train_fps = matrix(train.canonical_smiles.astype(str).tolist(), fpgen)
    x_test, test_fps = matrix(test.canonical_smiles.astype(str).tolist(), fpgen)
    y_train = train.target.to_numpy(int); y = test.pam_primary_screen_label.to_numpy(int)

    ridge = LogisticRegression(max_iter=3000, class_weight="balanced", C=1.0,
                               solver="liblinear", random_state=42).fit(x_train, y_train)
    rf = RandomForestClassifier(n_estimators=500, class_weight="balanced_subsample",
                                min_samples_leaf=2, max_features="sqrt", n_jobs=-1,
                                random_state=42).fit(x_train, y_train)
    test["ecfp_ridge"] = ridge.predict_proba(x_test)[:, 1]
    test["ecfp_rf"] = rf.predict_proba(x_test)[:, 1]
    test["max_train_similarity"] = [max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)) for fp in test_fps]

    descriptors = []
    for s in test.canonical_smiles:
        mol = Chem.MolFromSmiles(s)
        descriptors.append([Descriptors.MolWt(mol), Crippen.MolLogP(mol),
                            Descriptors.TPSA(mol), Descriptors.NumHDonors(mol),
                            Descriptors.NumHAcceptors(mol)])
    desc_test = np.asarray(descriptors, dtype=float)
    desc_train = []
    for s in train.canonical_smiles:
        mol = Chem.MolFromSmiles(s)
        desc_train.append([Descriptors.MolWt(mol), Crippen.MolLogP(mol),
                           Descriptors.TPSA(mol), Descriptors.NumHDonors(mol),
                           Descriptors.NumHAcceptors(mol)])
    desc_model = make_pipeline(StandardScaler(), LogisticRegression(
        max_iter=3000, class_weight="balanced", random_state=42)).fit(np.asarray(desc_train), y_train)
    test["physchem_logistic"] = desc_model.predict_proba(desc_test)[:, 1]
    test["murcko_scaffold"] = test.canonical_smiles.map(
        lambda s: MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s)))

    # Secondary property-matched panel: up to five unique nearby inactives per active.
    scaler = StandardScaler().fit(desc_test)
    z = scaler.transform(desc_test)
    pos_idx = np.flatnonzero(y == 1); neg_idx = np.flatnonzero(y == 0)
    nn = NearestNeighbors(n_neighbors=min(20, len(neg_idx))).fit(z[neg_idx])
    neighbors = nn.kneighbors(z[pos_idx], return_distance=False)
    picked = set()
    for row in neighbors:
        for local in row:
            picked.add(int(neg_idx[local]))
            if len(picked) >= 5 * len(pos_idx):
                break
    matched_idx = np.asarray(sorted(set(pos_idx) | picked), dtype=int)
    inactive_weight = 362354 / max(1, int((y == 0).sum()))
    methods = ["ecfp_ridge", "ecfp_rf", "max_train_similarity", "physchem_logistic"]
    result = {
        "dataset": "PubChem AID 624126 CHRM4 PAM primary screen",
        "n_after_exact_overlap_exclusion": int(len(test)),
        "active": int(y.sum()), "inactive": int((y == 0).sum()),
        "exact_training_overlap_removed": int(screen.exact_training_overlap.sum()),
        "unique_scaffolds": int(test.murcko_scaffold.nunique()),
        "full_deterministic_subset": {m: metrics(y, test[m].to_numpy(float), inactive_weight) for m in methods},
        "property_matched_secondary": {
            "n": int(len(matched_idx)), "active": int(y[matched_idx].sum()),
            "inactive": int((y[matched_idx] == 0).sum()),
            "methods": {m: metrics(y[matched_idx], test[m].to_numpy(float)[matched_idx]) for m in methods},
        },
        "claim_boundary": (
            "AID 624126 is a singlicate primary functional PAM screen. It measures large-scale "
            "screen transfer and enrichment, not confirmatory PAM pharmacology."
        ),
    }
    a.output_dir.mkdir(parents=True, exist_ok=True)
    (a.output_dir / "baseline_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    test.to_csv(a.output_dir / "baseline_predictions.csv", index=False)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
