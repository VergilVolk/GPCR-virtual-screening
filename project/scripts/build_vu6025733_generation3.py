# -*- coding: utf-8 -*-
"""Reconstruct the third-generation (29a-j) VU6025733 SAR table."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from build_vu6025733_external import HEAD


P = Path(__file__).resolve().parents[1]
OUT = P / "data" / "benchmarks" / "m4_pam_v1" / "external_2026_vu6025733"
TABLE1 = OUT / "external_potency.csv"

# The paper states that 29b/c are resolved single enantiomers but the absolute
# configuration is unknown.  @/@@ preserve their image-defined distinction;
# they must not be interpreted as assigned R/S configurations.
ROWS = {
    "29a": ("6", "c4ccc5c(c4)OC(F)(F)C(F)(F)O5", 3400, "exact", "tetrafluoro dioxine"),
    "29b": ("5", "c4ccc5c(c4)O[C@H](C)CO5", 75, "exact", "single enantiomer; absolute configuration unknown"),
    "29c": ("6", "c4ccc5c(c4)O[C@@H](C)CO5", 42, "exact", "single enantiomer; absolute configuration unknown"),
    "29d": ("5", "c4ccc5c(c4)OCCO5", 27, "exact", "d4 isotopologue; isotope normalized for 2D model"),
    "29e": ("6", "c4ccc5c(c4)OCCO5", 38, "exact", "d4 isotopologue; isotope normalized for 2D model"),
    "29f": ("6", "c4ccc5c(c4C)OCCO5", 1240, "exact", "5-methyl"),
    "29g": ("6", "c4ccc5c(c4C7CC7)OCCO5", None, "inactive_lt_4.5", "5-cyclopropyl"),
    "29h": ("6", "c4ccc5c(c4C#N)OCCO5", None, "inactive_lt_4.5", "5-cyano"),
    "29i": ("6", "c4ccc5c(c4F)OCCO5", 280, "exact", "5-fluoro"),
    "29j": ("6", "c4cc(F)c5c(c4)OCCO5", 72, "exact", "8-fluoro"),
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    table1 = pd.read_csv(TABLE1)
    known_2d = set(table1.canonical_smiles)
    rows = []
    for compound_id, (head, aryl, ec50, censoring, note) in ROWS.items():
        mol = Chem.MolFromSmiles(HEAD[head].format(AR=aryl))
        if mol is None:
            raise ValueError(compound_id)
        smi = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        exact = ec50 is not None
        p = 9.0 - np.log10(float(ec50)) if exact else 4.5
        rows.append({
            "compound_id": compound_id, "canonical_smiles": smi,
            "EC50_nM": ec50, "pEC50": p,
            "censoring": censoring,
            "pEC50_eval_ceiling": p,
            "head_group": "7-methyl" if head == "5" else "7,8-dimethyl",
            "structural_note": note,
            "molecular_formula_isotope_normalized": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight_isotope_normalized": Descriptors.MolWt(mol),
            "exact_table1_2d_overlap": smi in known_2d,
            "eligible_sequential_query": smi not in known_2d,
            "assay": "hM4/Gqi5-CHO calcium mobilization; ACh EC20",
            "source_doi": "10.1021/acschemneuro.5c00963",
            "structure_source": "Scheme 2/3 + Table 3 image/prose",
            "reconstruction_confidence": "high",
        })
    frame = pd.DataFrame(rows)
    path = OUT / "generation3_potency.csv"
    frame.to_csv(path, index=False)
    audit = {
        "dataset_id": "M4_VU6025733_2026_GENERATION3_V01",
        "n_total": len(frame),
        "n_sequential_query": int(frame.eligible_sequential_query.sum()),
        "n_isotope_normalized_2d_overlap_controls": int(frame.exact_table1_2d_overlap.sum()),
        "n_censored_inactive_queries": int((frame.eligible_sequential_query & frame.censoring.ne("exact")).sum()),
        "evaluation_policy": "Inactives use pEC50<4.5 as an ordering ceiling, not an exact regression label.",
        "analyst_blinding": "not blinded to paper labels; algorithmic rules must be frozen before prediction",
        "claim_boundary": "Within-paper sequential medicinal-chemistry challenge, not independent-publication validation.",
        "csv_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    (OUT / "generation3_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(frame[["compound_id", "EC50_nM", "pEC50", "censoring", "exact_table1_2d_overlap",
                 "eligible_sequential_query", "molecular_formula_isotope_normalized"]].to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
