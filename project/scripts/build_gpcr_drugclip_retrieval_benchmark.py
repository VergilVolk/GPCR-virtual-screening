#!/usr/bin/env python3
"""Build a four-target GPCR active-ligand retrieval benchmark.

Exact multi-target molecules are removed from the single-label primary task.
Murcko scaffolds are assigned wholly to one of five folds.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def canonical(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        return None
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def scaffold(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    value = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
    return value or Chem.MolToSmiles(molecule, canonical=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--libraries", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()
    rows, failures = [], []
    for target in ["B2AR", "CCR2", "M2R", "M4R"]:
        path = args.libraries / target / f"{target}_actives.smi"
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            parts = line.split()
            if not parts:
                continue
            smiles = canonical(parts[0])
            if smiles is None:
                failures.append({"target": target, "line": line_number, "raw": line})
                continue
            rows.append({"target": target, "canonical_smiles": smiles,
                         "source_id": parts[1] if len(parts) > 1 else f"{target}_{line_number}"})
    raw = pd.DataFrame(rows).drop_duplicates(["target", "canonical_smiles"]).reset_index(drop=True)
    target_count = raw.groupby("canonical_smiles").target.nunique()
    ambiguous = set(target_count[target_count > 1].index)
    primary = raw[~raw.canonical_smiles.isin(ambiguous)].copy().reset_index(drop=True)
    primary["murcko_scaffold"] = primary.canonical_smiles.map(scaffold)
    # Greedy vector balancing is more reliable than StratifiedGroupKFold for
    # the highly imbalanced 16/64/170/2254 target counts while still keeping a
    # scaffold indivisible.
    targets = ["B2AR", "CCR2", "M2R", "M4R"]
    vectors = primary.groupby(["murcko_scaffold", "target"]).size().unstack(fill_value=0)
    vectors = vectors.reindex(columns=targets, fill_value=0)
    totals = vectors.sum(0).to_numpy(float); ideal = totals / args.folds
    assigned = [dict() for _ in range(args.folds)]; counts = [pd.Series(0.0, index=targets) for _ in range(args.folds)]
    ordered = sorted(vectors.index, key=lambda value: (
        -float((vectors.loc[value].to_numpy(float) / ideal).max()),
        -int(vectors.loc[value].sum()), str(value)))
    for value in ordered:
        vector = vectors.loc[value].astype(float)
        def cost(fold: int):
            after = counts[fold] + vector
            target_cost = float(((after.to_numpy() / ideal) ** 2).sum())
            size_cost = float((after.sum() / (len(primary) / args.folds)) ** 2)
            return target_cost + 0.05 * size_cost, float(counts[fold].sum()), fold
        fold = min(range(args.folds), key=cost)
        assigned[fold][value] = True; counts[fold] += vector
    fold_map = {value: fold for fold, values in enumerate(assigned) for value in values}
    primary["scaffold_fold"] = primary.murcko_scaffold.map(fold_map).astype(int)
    primary.insert(0, "canonical_molecule_id", [f"GPCRRET{i:05d}" for i in range(len(primary))])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    primary.to_csv(args.output_dir / "retrieval_benchmark.csv", index=False)
    primary[["canonical_molecule_id", "canonical_smiles"]].to_csv(
        args.output_dir / "molecules.csv", index=False)
    raw[raw.canonical_smiles.isin(ambiguous)].to_csv(
        args.output_dir / "excluded_multitarget_molecules.csv", index=False)
    fold_counts = primary.groupby(["scaffold_fold", "target"]).size().unstack(fill_value=0)
    audit = {
        "raw_unique_target_molecule_pairs": int(len(raw)),
        "primary_single_label_molecules": int(len(primary)),
        "exact_multitarget_molecules_excluded": int(len(ambiguous)),
        "invalid_smiles": failures,
        "unique_scaffolds": int(primary.murcko_scaffold.nunique()),
        "target_counts": primary.target.value_counts().sort_index().astype(int).to_dict(),
        "fold_target_counts": {str(i): {k: int(v) for k, v in row.items()}
                               for i, row in fold_counts.to_dict("index").items()},
        "split_rule": "Five-fold deterministic greedy target-vector balancing; Murcko scaffold is indivisible.",
        "primary_metrics": ["macro Recall@1", "Recall@2", "MRR", "pair ROC-AUC", "pair PR-AUC"],
        "claim_boundary": "Retrieval of literature-curated GPCR allosteric-modulator target; not efficacy or prospective binding.",
    }
    (args.output_dir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
