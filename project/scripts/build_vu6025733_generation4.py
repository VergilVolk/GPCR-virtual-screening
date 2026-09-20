# -*- coding: utf-8 -*-
"""Reconstruct the high-confidence stereochemistry/isotope subset 33h-r."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


P = Path(__file__).resolve().parents[1]
OUT = P / "data" / "benchmarks" / "m4_pam_v1" / "external_2026_vu6025733"
CORE = {
    "mono": "C1=C(C)C({N})=NN3C1=NN=C3",
    "dimethyl": "CC1=C(C)C({N})=NN3C1=NN=C3",
}
AR_H = "Oc4ccc5c(c4)OCCO5"
AR_D4 = "Oc4ccc5c(c4)O[C]([2H])([2H])[C]([2H])([2H])O5"

# @ choices were enumerated and verified with RDKit CIP assignment:
# (@,@)=(3S,4R), (@@,@@)=(3R,4S), (@,@@)=(3S,4S), (@@,@)=(3R,4R).
CONFIG = {
    "3S4R": ("@", "@"), "3R4S": ("@@", "@@"),
    "3S4S": ("@", "@@"), "3R4R": ("@@", "@"),
}
POTENCY = {
    "33h": 23, "33i": 47, "33j": 42, "33k": 81,
    "33l": 147, "33m": 184, "33n": 91, "33o": 52,
    "33p": 25, "33q": 33, "33r": 23,
}
SPEC = {
    "33h": ("dimethyl", "3S4R", False, False),
    "33i": ("dimethyl", "3S4R", True, False),
    "33j": ("dimethyl", "3R4S", False, False),
    "33k": ("dimethyl", "3R4S", True, False),
    "33l": ("dimethyl", "3S4S", False, False),
    "33m": ("dimethyl", "3S4S", True, False),
    "33n": ("dimethyl", "3R4R", False, False),
    "33o": ("dimethyl", "3R4R", True, False),
    "33p": ("dimethyl", None, False, True),
    "33q": ("dimethyl", None, True, True),
    "33r": ("mono", None, True, True),
}


def fluoropiperidine(config, ar):
    a, b = CONFIG[config]
    return f"N2C[C{a}H](F)[C{b}H]({ar})CC2"


def piperidine_d(ar):
    return f"N2CCC({ar})([2H])CC2"


def isotope_stripped_smiles(mol):
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetIsotope(0)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for compound_id, (head, config, d4, piperidine_d4) in SPEC.items():
        ar = AR_D4 if d4 else AR_H
        ring = piperidine_d(ar) if piperidine_d4 else fluoropiperidine(config, ar)
        mol = Chem.MolFromSmiles(CORE[head].format(N=ring))
        if mol is None:
            raise ValueError(compound_id)
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        ec50 = POTENCY[compound_id]
        rows.append({
            "compound_id": compound_id,
            "canonical_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            "isotope_stripped_smiles": isotope_stripped_smiles(mol),
            "EC50_nM": ec50, "pEC50": 9.0 - np.log10(ec50),
            "censoring": "exact", "stereochemical_class": config or "nonfluoro",
            "dioxin_d4": d4, "piperidine_4_d": piperidine_d4,
            "rdkit_chiral_centers": "|".join(f"{i}:{c}" for i, c in centers),
            "primary_stereochemical_query": compound_id in {"33h", "33j", "33l", "33n"},
            "isotope_or_headgroup_control": compound_id not in {"33h", "33j", "33l", "33n"},
            "molecular_formula_rdkit": rdMolDescriptors.CalcMolFormula(mol),
            "assay": "hM4/Gqi5-CHO calcium mobilization; ACh EC20",
            "source_doi": "10.1021/acschemneuro.5c00963",
            "structure_source": "Scheme 4/5 + Table 4 image/prose",
            "reconstruction_confidence": "high",
        })
    frame = pd.DataFrame(rows)
    path = OUT / "generation4_stereo_potency.csv"
    frame.to_csv(path, index=False)
    audit = {
        "dataset_id": "M4_VU6025733_2026_GENERATION4_STEREO_V01",
        "n_total": len(frame), "n_primary_stereoisomers": int(frame.primary_stereochemical_query.sum()),
        "n_isotope_or_headgroup_controls": int(frame.isotope_or_headgroup_control.sum()),
        "primary_question": "Can a model rank the four 3-fluoro-piperidine stereoisomers?",
        "claim_boundary": "Tiny stereochemical stress test; isotope labels are not independent chemotypes.",
        "csv_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    (OUT / "generation4_stereo_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(frame[["compound_id", "EC50_nM", "pEC50", "stereochemical_class", "dioxin_d4",
                 "piperidine_4_d", "rdkit_chiral_centers", "molecular_formula_rdkit"]].to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
