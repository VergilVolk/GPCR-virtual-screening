#!/usr/bin/env python3
"""Freeze an M4-held-out, molecule-disjoint muscarinic activity test.

The DrugCLIP pocket projection was trained only with M1/M2/M3/M5 triplets.
This test asks whether it transfers to an unseen M4 pocket. All molecules seen
in that triplet training set and all PACER-M4 overlaps are excluded.
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
    parser.add_argument("--training-triplets", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-positive", type=float, default=6.0)
    parser.add_argument("--min-delta", type=float, default=1.0)
    args = parser.parse_args()
    activity = pd.read_csv(args.activity)
    training = pd.read_csv(args.training_triplets)
    training_molecules = set(training.canonical_smiles)
    activity = activity[~activity.overlaps_pacer_m4.astype(bool)].copy()
    rows = []
    for smiles, group in activity.groupby("canonical_smiles", sort=False):
        if smiles in training_molecules or "M4" not in set(group.subtype):
            continue
        for positive, negative in itertools.permutations(group.itertuples(), 2):
            if "M4" not in {positive.subtype, negative.subtype}:
                continue
            delta = float(positive.pchembl_median - negative.pchembl_median)
            if positive.pchembl_median >= args.min_positive and delta >= args.min_delta:
                rows.append({
                    "canonical_smiles": smiles,
                    "positive_subtype": positive.subtype,
                    "negative_subtype": negative.subtype,
                    "positive_pchembl": float(positive.pchembl_median),
                    "negative_pchembl": float(negative.pchembl_median),
                    "delta_pchembl": delta,
                    "m4_direction": "M4_preferred" if positive.subtype == "M4" else "M4_disfavored",
                })
    test = pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)
    molecules = test[["canonical_smiles"]].drop_duplicates().reset_index(drop=True)
    molecules.insert(0, "canonical_molecule_id", [f"M4ZERO{i:04d}" for i in range(len(molecules))])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    test.to_csv(args.output_dir / "m4_zero_shot_pairs.csv", index=False)
    molecules.to_csv(args.output_dir / "molecules.csv", index=False)
    audit = {
        "pairs": int(len(test)),
        "molecules": int(test.canonical_smiles.nunique()),
        "m4_preferred_pairs": int((test.m4_direction == "M4_preferred").sum()),
        "m4_disfavored_pairs": int((test.m4_direction == "M4_disfavored").sum()),
        "training_molecule_overlap": int(len(set(test.canonical_smiles) & training_molecules)),
        "pacer_m4_overlap": 0,
        "rules": f"M4-containing pair; positive pChEMBL >= {args.min_positive}; delta >= {args.min_delta}; triplet-training molecules excluded.",
        "evidence_level": "posthoc_frozen_zero_shot_test",
        "claim_boundary": "Heterogeneous ChEMBL subtype activity; not M4 PAM-specific.",
    }
    (args.output_dir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
