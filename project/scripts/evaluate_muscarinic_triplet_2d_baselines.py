#!/usr/bin/env python3
"""Frozen 2D baselines for the molecule-disjoint muscarinic triplet test.

Each subtype receives an independent pChEMBL model trained only on the expanded
training molecules.  The benchmark asks whether predicted activity is higher
for the experimentally preferred subtype of the same molecule.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge


def fp_matrix(smiles: list[str]) -> np.ndarray:
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    matrix = np.zeros((len(smiles), 2048), dtype=np.float32)
    for i, value in enumerate(smiles):
        DataStructs.ConvertToNumpyArray(gen.GetFingerprint(Chem.MolFromSmiles(value)), matrix[i])
    return matrix


def clustered_delta_ci(test: pd.DataFrame, baseline: np.ndarray, target: np.ndarray,
                       n: int = 5000, seed: int = 20260924) -> list[float]:
    molecules = test.canonical_smiles.unique(); rng = np.random.default_rng(seed); values = []
    for _ in range(n):
        sampled = rng.choice(molecules, len(molecules), replace=True)
        index = np.concatenate([np.flatnonzero(test.canonical_smiles.to_numpy() == value)
                                for value in sampled])
        values.append(float((target[index] > 0).mean() - (baseline[index] > 0).mean()))
    return list(map(float, np.quantile(values, [0.025, 0.5, 0.975])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--drugclip-predictions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    train = pd.read_csv(args.train); test = pd.read_csv(args.test)
    if set(train.canonical_smiles) & set(test.canonical_smiles):
        raise ValueError("Training and strict test molecules overlap")
    long = pd.concat([
        train[["canonical_smiles", "positive_subtype", "positive_pchembl"]].rename(
            columns={"positive_subtype": "subtype", "positive_pchembl": "pchembl"}),
        train[["canonical_smiles", "negative_subtype", "negative_pchembl"]].rename(
            columns={"negative_subtype": "subtype", "negative_pchembl": "pchembl"}),
    ]).groupby(["canonical_smiles", "subtype"], as_index=False).pchembl.mean()
    test_x = fp_matrix(test.canonical_smiles.astype(str).tolist())
    methods = {"subtype_prior": np.zeros(len(test)), "ecfp_ridge": np.zeros(len(test)),
               "ecfp_random_forest": np.zeros(len(test))}
    models = {}
    for subtype, frame in long.groupby("subtype"):
        x = fp_matrix(frame.canonical_smiles.astype(str).tolist()); y = frame.pchembl.to_numpy(float)
        models[subtype] = {
            "subtype_prior": float(y.mean()),
            "ecfp_ridge": Ridge(alpha=10.0).fit(x, y),
            "ecfp_random_forest": RandomForestRegressor(
                n_estimators=600, min_samples_leaf=2, max_features=0.35,
                random_state=42, n_jobs=-1).fit(x, y),
        }
    for i, row in enumerate(test.itertuples()):
        for method in methods:
            positive = models[row.positive_subtype][method]
            negative = models[row.negative_subtype][method]
            if method != "subtype_prior":
                positive = positive.predict(test_x[i:i + 1])[0]
                negative = negative.predict(test_x[i:i + 1])[0]
            methods[method][i] = positive - negative
    drugclip = pd.read_csv(args.drugclip_predictions).groupby(
        ["canonical_smiles", "positive_subtype", "negative_subtype"], as_index=False
    ).targeted_delta.mean()
    merged = test.merge(drugclip, on=["canonical_smiles", "positive_subtype", "negative_subtype"],
                        validate="one_to_one")
    targeted = merged.targeted_delta.to_numpy(float)
    report = {
        "protocol": "Molecule-disjoint strict test; no test-label tuning",
        "train_molecules": int(train.canonical_smiles.nunique()),
        "test_molecules": int(test.canonical_smiles.nunique()),
        "test_pairs": int(len(test)),
        "accuracy": {name: float((score > 0).mean()) for name, score in methods.items()},
        "drugclip_gpcr_triplet_accuracy": float((targeted > 0).mean()),
        "drugclip_minus_2d_molecule_cluster_bootstrap_95ci": {
            name: clustered_delta_ci(test, score, targeted) for name, score in methods.items()
        },
        "claim_boundary": "Subtype activity ranking, not M4 PAM efficacy or prospective screening.",
    }
    predictions = test[["canonical_smiles", "positive_subtype", "negative_subtype"]].copy()
    for name, score in methods.items(): predictions[name + "_delta"] = score
    predictions["drugclip_gpcr_triplet_delta"] = targeted
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    predictions.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
