#!/usr/bin/env python3
"""Build strict leave-one-muscarinic-subtype-out transfer folds.

For each held subtype, every test molecule is removed from that fold's training
triplets. The held pocket is never used in supervision.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--activity", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-positive", type=float, default=6.0)
    parser.add_argument("--min-delta", type=float, default=1.0)
    args = parser.parse_args()
    data = pd.read_csv(args.activity)
    data = data[~data.overlaps_pacer_m4.astype(bool)].copy()
    rows = []
    for smiles, group in data.groupby("canonical_smiles", sort=False):
        for positive, negative in itertools.permutations(group.itertuples(), 2):
            delta = float(positive.pchembl_median - negative.pchembl_median)
            if positive.pchembl_median >= args.min_positive and delta >= args.min_delta:
                rows.append({
                    "canonical_smiles": smiles,
                    "positive_subtype": positive.subtype,
                    "negative_subtype": negative.subtype,
                    "positive_pchembl": float(positive.pchembl_median),
                    "negative_pchembl": float(negative.pchembl_median),
                    "delta_pchembl": delta,
                })
    all_pairs = pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fold_audit = {}
    used_molecules = set()
    for held in ["M1", "M2", "M3", "M4", "M5"]:
        test = all_pairs[(all_pairs.positive_subtype == held) | (all_pairs.negative_subtype == held)].copy()
        test_molecules = set(test.canonical_smiles)
        train = all_pairs[
            (all_pairs.positive_subtype != held) & (all_pairs.negative_subtype != held) &
            ~all_pairs.canonical_smiles.isin(test_molecules)
        ].copy()
        train.to_csv(args.output_dir / f"{held}_train.csv", index=False)
        test.to_csv(args.output_dir / f"{held}_test.csv", index=False)
        used_molecules.update(train.canonical_smiles); used_molecules.update(test.canonical_smiles)
        fold_audit[held] = {
            "train_pairs": int(len(train)), "train_molecules": int(train.canonical_smiles.nunique()),
            "test_pairs": int(len(test)), "test_molecules": int(test.canonical_smiles.nunique()),
            "molecule_overlap": int(len(set(train.canonical_smiles) & test_molecules)),
            "held_pocket_in_training": False,
        }
    molecules = pd.DataFrame({"canonical_smiles": sorted(used_molecules)})
    molecules.insert(0, "canonical_molecule_id", [f"MUSCLOSO{i:05d}" for i in range(len(molecules))])
    molecules.to_csv(args.output_dir / "molecules.csv", index=False)
    audit = {
        "all_qualifying_pairs": int(len(all_pairs)),
        "all_qualifying_molecules": int(all_pairs.canonical_smiles.nunique()),
        "encoded_molecules": int(len(molecules)),
        "folds": fold_audit,
        "rules": f"PACER-overlap excluded; positive pChEMBL >= {args.min_positive}; delta >= {args.min_delta}; held-target test molecules removed from training.",
        "claim_boundary": "Heterogeneous ChEMBL muscarinic activity, not PAM-specific.",
    }
    (args.output_dir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
