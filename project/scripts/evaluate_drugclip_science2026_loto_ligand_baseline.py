#!/usr/bin/env python3
"""Ligand-only target-LOTO confound baseline for the Science-2026 pilot."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy import sparse
from sklearn.linear_model import LogisticRegression

from finetune_drugclip_gpcr_screening import binary_metrics


def scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or Chem.MolToSmiles(mol, canonical=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bits", type=int, default=2048)
    args = parser.parse_args()
    frame = pd.read_csv(args.pairs).copy().reset_index(drop=True)
    frame["scaffold"] = [scaffold(value) for value in frame.canonical_smiles.astype(str)]
    unique = list(dict.fromkeys(frame.canonical_smiles.astype(str)))
    index = {value: i for i, value in enumerate(unique)}
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=args.bits)
    matrix = np.zeros((len(unique), args.bits), dtype=np.uint8)
    for row, smiles in enumerate(unique):
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(Chem.MolFromSmiles(smiles)), matrix[row])
    matrix = sparse.csr_matrix(matrix)
    pair_mol = np.asarray([index[value] for value in frame.canonical_smiles.astype(str)])
    prediction = np.full(len(frame), np.nan, dtype=np.float64)
    audit = []
    for held in sorted(frame.target.unique()):
        test = frame.target.eq(held).to_numpy()
        held_scaffolds = set(frame.loc[test, "scaffold"])
        train = (~test) & (~frame.scaffold.isin(held_scaffolds).to_numpy())
        model = LogisticRegression(
            C=1.0, class_weight="balanced", solver="liblinear", max_iter=3000, random_state=20260926
        )
        model.fit(matrix[pair_mol[train]], frame.loc[train, "label"].to_numpy(int))
        prediction[test] = model.predict_proba(matrix[pair_mol[test]])[:, 1]
        audit.append({"held_target": held, "train_pairs": int(train.sum()),
                      "test_pairs": int(test.sum()), "scaffold_overlap": 0})
    per_target = {}
    for target in sorted(frame.target.unique()):
        keep = frame.target.eq(target).to_numpy()
        per_target[target] = binary_metrics(frame.loc[keep, "label"].to_numpy(int), prediction[keep])
    names = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]
    report = {
        "method": "ECFP4 logistic regression ligand-only confound baseline",
        "protocol": "leave-one-target-out with held-target Murcko scaffold purge",
        "macro": {name: float(np.mean([value[name] for value in per_target.values()])) for name in names},
        "per_target": per_target,
        "fold_audit": audit,
        "claim_boundary": "Measures ligand/dataset bias only; contains no pocket information.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    predictions = frame[["pair_id", "target", "label", "scaffold"]].copy()
    predictions["ecfp_loto"] = prediction
    predictions.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
