# -*- coding: utf-8 -*-
"""Build the PACER-M4 activity-cliff and local-SAR pair benchmark.

Definitions (frozen defaults):
  near-neighbor pair : ECFP4 Tanimoto >= 0.70
  activity cliff     : abs(delta pEC50) >= 1.00 (>=10-fold potency change)
  strong cliff       : abs(delta pEC50) >= 1.50
  smooth pair        : abs(delta pEC50) <= 0.30

Pairs retain both molecule-level folds. A pair is eligible for a held-out pair
test only when both molecules belong to the same held-out fold. This prevents a
training molecule from appearing on one side of a test pair.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(__file__).resolve().parents[1]
INPUT = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
OUT = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "activity_cliffs"
SIM_THRESHOLD = 0.70
CLIFF_DELTA = 1.00
STRONG_CLIFF_DELTA = 1.50
SMOOTH_DELTA = 0.30


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT).sort_values("canonical_molecule_id").reset_index(drop=True)
    mols = [Chem.MolFromSmiles(s) for s in df["canonical_smiles"]]
    if any(m is None for m in mols):
        raise ValueError("Invalid SMILES in frozen potency benchmark")
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]

    rows = []
    for i in range(len(df)):
        similarities = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
        for offset, similarity in enumerate(similarities):
            if similarity < SIM_THRESHOLD:
                continue
            j = i + 1 + offset
            left, right = df.iloc[i], df.iloc[j]
            delta = float(right.pEC50 - left.pEC50)
            abs_delta = abs(delta)
            # Conservative lower bound using observed replicate dispersion.
            pooled_sd = float(np.sqrt(float(left.pEC50_sd) ** 2 + float(right.pEC50_sd) ** 2))
            conservative_delta = max(0.0, abs_delta - 1.96 * pooled_sd)
            if abs_delta >= CLIFF_DELTA:
                pair_class = "cliff"
            elif abs_delta <= SMOOTH_DELTA:
                pair_class = "smooth"
            else:
                pair_class = "intermediate"
            scaffold_pair_fold = (int(left.scaffold_fold)
                                  if left.scaffold_fold == right.scaffold_fold else -1)
            source_pair_fold = (int(left.source_fold)
                                if left.source_fold == right.source_fold else -1)
            stronger = (right.canonical_molecule_id if delta > 0 else left.canonical_molecule_id
                        if delta < 0 else "tie")
            rows.append({
                "pair_id": f"{left.canonical_molecule_id}__{right.canonical_molecule_id}",
                "mol_a": left.canonical_molecule_id,
                "mol_b": right.canonical_molecule_id,
                "smiles_a": left.canonical_smiles,
                "smiles_b": right.canonical_smiles,
                "pEC50_a": float(left.pEC50),
                "pEC50_b": float(right.pEC50),
                "delta_pEC50_b_minus_a": delta,
                "abs_delta_pEC50": abs_delta,
                "pooled_observed_sd": pooled_sd,
                "conservative_abs_delta": conservative_delta,
                "tanimoto": float(similarity),
                "pair_class": pair_class,
                "is_cliff": int(abs_delta >= CLIFF_DELTA),
                "is_strong_cliff": int(abs_delta >= STRONG_CLIFF_DELTA),
                "is_conservative_cliff": int(conservative_delta >= CLIFF_DELTA),
                "stronger_molecule": stronger,
                "same_murcko_scaffold": int(left.murcko_scaffold == right.murcko_scaffold),
                "same_primary_source": int(left.primary_source_id == right.primary_source_id),
                "same_source_component": int(left.source_component == right.source_component),
                "scaffold_fold_a": int(left.scaffold_fold),
                "scaffold_fold_b": int(right.scaffold_fold),
                "scaffold_pair_fold": scaffold_pair_fold,
                "source_fold_a": int(left.source_fold),
                "source_fold_b": int(right.source_fold),
                "source_pair_fold": source_pair_fold,
            })

    pairs = pd.DataFrame(rows).sort_values(
        ["is_cliff", "abs_delta_pEC50", "tanimoto"], ascending=[False, False, False]
    )
    pairs.to_csv(OUT / "near_neighbor_pairs.csv", index=False)
    pairs[pairs["pair_class"].isin(["cliff", "smooth"])].to_csv(
        OUT / "cliff_vs_smooth.csv", index=False
    )
    pairs[pairs["is_cliff"] == 1].to_csv(OUT / "activity_cliffs.csv", index=False)

    audit = {
        "definitions": {
            "similarity_threshold": SIM_THRESHOLD,
            "cliff_delta_pEC50": CLIFF_DELTA,
            "strong_cliff_delta_pEC50": STRONG_CLIFF_DELTA,
            "smooth_delta_pEC50": SMOOTH_DELTA,
        },
        "n_molecules": int(len(df)),
        "n_near_pairs": int(len(pairs)),
        "pair_class_counts": {str(k): int(v) for k, v in pairs["pair_class"].value_counts().items()},
        "n_cliffs": int(pairs["is_cliff"].sum()),
        "n_strong_cliffs": int(pairs["is_strong_cliff"].sum()),
        "n_conservative_cliffs": int(pairs["is_conservative_cliff"].sum()),
        "n_same_scaffold_cliffs": int(((pairs["is_cliff"] == 1) &
                                       (pairs["same_murcko_scaffold"] == 1)).sum()),
        "n_same_source_cliffs": int(((pairs["is_cliff"] == 1) &
                                     (pairs["same_primary_source"] == 1)).sum()),
        "pair_test_coverage": {
            "scaffold": {
                "eligible_pairs": int((pairs["scaffold_pair_fold"] >= 0).sum()),
                "eligible_cliffs": int(((pairs["scaffold_pair_fold"] >= 0) &
                                         (pairs["is_cliff"] == 1)).sum()),
            },
            "source": {
                "eligible_pairs": int((pairs["source_pair_fold"] >= 0).sum()),
                "eligible_cliffs": int(((pairs["source_pair_fold"] >= 0) &
                                         (pairs["is_cliff"] == 1)).sum()),
            },
        },
        "leakage_rule": (
            "A pair is test-eligible only when both molecules share the held-out fold; "
            "pairs crossing folds are excluded from pair-level test metrics."
        ),
    }
    with open(OUT / "audit_report.json", "w", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2)
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

