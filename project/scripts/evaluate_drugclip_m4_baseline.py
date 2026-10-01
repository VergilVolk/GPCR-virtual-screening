#!/usr/bin/env python3
"""Evaluate frozen DrugCLIP scores on the frozen M4 PAM/inactive benchmark.

This is an unsupervised binding-compatibility baseline, not a PAM classifier
training result. Molecule labels are used only after all embeddings are frozen.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from sklearn.metrics import average_precision_score, roc_auc_score


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"Invalid benchmark SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)


def bootstrap_auc(y: np.ndarray, score: np.ndarray, repeats: int, seed: int) -> list[float]:
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    values = []
    for _ in range(repeats):
        idx = np.concatenate([rng.choice(pos, len(pos), replace=True),
                              rng.choice(neg, len(neg), replace=True)])
        values.append(float(roc_auc_score(y[idx], score[idx])))
    return values


def metrics(y: np.ndarray, score: np.ndarray, repeats: int, seed: int) -> dict:
    aucs = bootstrap_auc(y, score, repeats, seed)
    return {
        "roc_auc": float(roc_auc_score(y, score)),
        "roc_auc_stratified_bootstrap_95ci": list(map(float, np.quantile(aucs, [0.025, 0.975]))),
        "average_precision": float(average_precision_score(y, score)),
        "prevalence": float(y.mean()),
        "n": int(len(y)), "n_positive": int(y.sum()), "n_negative": int((y == 0).sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260922)
    args = parser.parse_args()

    archive = np.load(args.embeddings, allow_pickle=False)
    molecule_ids = list(map(str, archive["molecule_ids"]))
    pocket_ids = list(map(str, archive["pocket_ids"]))
    scores = archive["scores"].astype(float)
    if scores.shape != (len(pocket_ids), len(molecule_ids)):
        raise ValueError("Score matrix shape does not match identifiers")
    if len(set(molecule_ids)) != len(molecule_ids):
        raise ValueError("DrugCLIP molecule identifiers are not unique")

    benchmark = pd.read_csv(args.benchmark)
    required = {"canonical_molecule_id", "canonical_smiles", "target", "source_component"}
    if missing := required - set(benchmark.columns):
        raise ValueError(f"Missing benchmark columns: {sorted(missing)}")
    benchmark["drugclip_smiles"] = benchmark.canonical_smiles.map(canonical)
    index = {smiles: i for i, smiles in enumerate(molecule_ids)}
    missing = benchmark.loc[~benchmark.drugclip_smiles.isin(index), "canonical_molecule_id"].tolist()
    if missing:
        raise ValueError(f"Benchmark molecules missing from embeddings: {missing[:10]}")
    order = np.asarray([index[s] for s in benchmark.drugclip_smiles], dtype=int)
    ordered = scores[:, order]
    y = benchmark.target.to_numpy(int)

    methods = {pocket: ordered[i] for i, pocket in enumerate(pocket_ids)}
    if len(pocket_ids) > 1:
        methods["state_mean"] = ordered.mean(axis=0)
        methods["state_max"] = ordered.max(axis=0)
    report = {
        "evidence_level": "retrospective_unsupervised_binding_baseline",
        "benchmark": str(args.benchmark),
        "methods": {},
        "claim_boundary": (
            "DrugCLIP scores binding compatibility. Labels were not used to fit or select a score, "
            "but this benchmark does not experimentally confirm PAM function."
        ),
    }
    score_table = benchmark[["canonical_molecule_id", "target", "source_component"]].copy()
    for offset, (name, score) in enumerate(methods.items()):
        score_table[name] = score
        result = metrics(y, score, args.bootstrap, args.seed + offset)
        per_source = {}
        for source, group in benchmark.groupby("source_component"):
            idx = group.index.to_numpy()
            if len(np.unique(y[idx])) == 2:
                per_source[str(source)] = metrics(y[idx], score[idx], min(args.bootstrap, 500),
                                                  args.seed + offset + len(per_source) + 1)
        result["per_source_with_both_classes"] = per_source
        report["methods"][name] = result
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    score_table.to_csv(args.output.with_suffix(".scores.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
