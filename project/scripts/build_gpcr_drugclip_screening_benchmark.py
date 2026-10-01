#!/usr/bin/env python3
"""Build a leakage-controlled four-GPCR active/decoy screening benchmark.

The same canonical molecule and Murcko scaffold are assigned globally to one
fold even when they occur for multiple targets.  Decoys are selected by a
stable hash, never by a model score.  The resulting pair table supports
target-specific virtual-screening evaluation and direct joins to the published
docking-score files through source_id.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


TARGETS = ["B2AR", "CCR2", "M2R", "M4R"]


def canonical(value: str) -> str | None:
    molecule = Chem.MolFromSmiles(value)
    return None if molecule is None else Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def murcko(value: str) -> str:
    molecule = Chem.MolFromSmiles(value)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
    return scaffold or Chem.MolToSmiles(molecule, canonical=True)


def stable_key(*values: str) -> str:
    return hashlib.sha256("|".join(values).encode("utf-8")).hexdigest()


def read_library(path: Path, target: str, label: int) -> tuple[pd.DataFrame, list[dict]]:
    rows, failures = [], []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        parts = line.split()
        if not parts:
            continue
        value = canonical(parts[0])
        if value is None:
            failures.append({"target": target, "label": label, "line": line_number, "raw": line})
            continue
        rows.append({"target": target, "label": label, "canonical_smiles": value,
                     "source_id": parts[1] if len(parts) > 1 else f"{target}_{label}_{line_number}"})
    return pd.DataFrame(rows).drop_duplicates(["target", "label", "canonical_smiles"]), failures


def assign_folds(pairs: pd.DataFrame, folds: int) -> dict[str, int]:
    targets = [f"{target}_{label}" for target in TARGETS for label in [0, 1]]
    vectors = pairs.groupby(["murcko_scaffold", "target", "label"]).size().unstack([1, 2], fill_value=0)
    vectors.columns = [f"{target}_{label}" for target, label in vectors.columns]
    vectors = vectors.reindex(columns=targets, fill_value=0)
    totals = vectors.sum(0).to_numpy(float)
    ideal = np.maximum(totals / folds, 1.0)
    counts = [pd.Series(0.0, index=targets) for _ in range(folds)]
    assignments = [{} for _ in range(folds)]
    ordered = sorted(vectors.index, key=lambda value: (
        -float((vectors.loc[value].to_numpy(float) / ideal).max()),
        -int(vectors.loc[value].sum()), str(value)))
    for value in ordered:
        vector = vectors.loc[value].astype(float)
        def cost(fold):
            # Minimise the imbalance of the complete fold-by-class table, not
            # merely the load of the candidate fold.  The latter can starve a
            # rare target/label stratum when a scaffold spans other strata.
            trial = np.stack([item.to_numpy(float) for item in counts])
            trial[fold] += vector.to_numpy(float)
            imbalance = float((((trial - ideal[None, :]) / ideal[None, :]) ** 2).sum())
            return imbalance, float(counts[fold].sum()), fold
        fold = min(range(folds), key=cost)
        assignments[fold][value] = True
        counts[fold] += vector
    return {value: fold for fold, values in enumerate(assignments) for value in values}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--libraries", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--decoys-per-active", type=int, default=10)
    parser.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()

    selected, invalid, conflicts, source_counts = [], [], [], {}
    for target in TARGETS:
        active, failed = read_library(args.libraries / target / f"{target}_actives.smi", target, 1)
        decoy, failed_decoy = read_library(args.libraries / target / f"{target}_decoys.smi", target, 0)
        invalid.extend(failed + failed_decoy)
        source_counts[target] = {"actives": int(len(active)), "decoys": int(len(decoy))}
        overlap = set(active.canonical_smiles) & set(decoy.canonical_smiles)
        if overlap:
            conflicts.extend({"target": target, "canonical_smiles": value} for value in sorted(overlap))
            active = active[~active.canonical_smiles.isin(overlap)]
            decoy = decoy[~decoy.canonical_smiles.isin(overlap)]
        keep_n = min(len(decoy), args.decoys_per_active * len(active))
        decoy = decoy.assign(_key=[stable_key(target, value) for value in decoy.canonical_smiles])
        decoy = decoy.sort_values("_key").head(keep_n).drop(columns="_key")
        selected.extend([active, decoy])

    pairs = pd.concat(selected, ignore_index=True)
    pairs["murcko_scaffold"] = pairs.canonical_smiles.map(murcko)
    fold_map = assign_folds(pairs, args.folds)
    pairs["scaffold_fold"] = pairs.murcko_scaffold.map(fold_map).astype(int)
    unique = pairs[["canonical_smiles"]].drop_duplicates().sort_values("canonical_smiles").reset_index(drop=True)
    unique.insert(0, "canonical_molecule_id", [f"GPCRVS{i:06d}" for i in range(len(unique))])
    pairs = pairs.merge(unique, on="canonical_smiles", how="left", validate="many_to_one")
    pairs.insert(0, "pair_id", [f"GPCRPAIR{i:06d}" for i in range(len(pairs))])
    pairs = pairs[["pair_id", "canonical_molecule_id", "target", "label", "source_id",
                   "canonical_smiles", "murcko_scaffold", "scaffold_fold"]]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.output_dir / "screening_pairs.csv", index=False)
    unique.to_csv(args.output_dir / "molecules.csv", index=False)
    pd.DataFrame(conflicts).to_csv(args.output_dir / "excluded_label_conflicts.csv", index=False)
    fold_counts = pairs.groupby(["scaffold_fold", "target", "label"]).size().unstack([1, 2], fill_value=0)
    audit = {
        "source_counts": source_counts,
        "decoys_per_active": args.decoys_per_active,
        "selected_pairs": int(len(pairs)),
        "unique_molecules": int(len(unique)),
        "unique_scaffolds": int(pairs.murcko_scaffold.nunique()),
        "selected_counts": {target: {
            "actives": int(((pairs.target == target) & (pairs.label == 1)).sum()),
            "decoys": int(((pairs.target == target) & (pairs.label == 0)).sum())}
            for target in TARGETS},
        "invalid_rows": invalid,
        "within_target_label_conflicts_excluded": len(conflicts),
        "fold_counts": {str(fold): {f"{target}_{label}": int(row.get((target, label), 0))
                                     for target in TARGETS for label in [0, 1]}
                        for fold, row in fold_counts.iterrows()},
        "selection_rule": "All non-conflicting actives; deterministic SHA256 decoy subset; no model-dependent selection.",
        "split_rule": "Global Murcko scaffold is indivisible across every target and label.",
        "primary_metrics": ["macro ROC-AUC", "macro PR-AUC", "macro BEDROC20", "macro EF1%", "macro EF5%"],
        "claim_boundary": "Retrospective active/decoy screening; not affinity, PAM efficacy, or prospective validation.",
    }
    (args.output_dir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
