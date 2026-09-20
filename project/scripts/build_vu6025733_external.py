# -*- coding: utf-8 -*-
"""Build the unambiguous 2026 VU6025733 external SAR subset.

The paper reports structures as images.  This first frozen subset contains only
analogues whose structures are uniquely specified by Scheme 1 and Table 1:
the shared piperidin-4-yloxy/triazolopyridazine scaffold, two head-group states
(5: 7-methyl; 6: 7,8-dimethyl), and nine unambiguous aryl ether groups.

Ambiguous image-only entries are deliberately excluded until independently
reconstructed and checked.  Labels are never used to build the structures.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


P = Path(__file__).resolve().parents[1]
OUT = P / "data" / "benchmarks" / "m4_pam_v1" / "external_2026_vu6025733"
TRAIN = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
DOI = "10.1021/acschemneuro.5c00963"

# 5 has R1=H/R2=Me; 6 has R1=R2=Me in the paper's Scheme 1.
# {AR} is the aromatic atom directly bonded to the piperidine ether oxygen.
HEAD = {
    "5": "C1=C(C)C(N2CCC(O{AR})CC2)=NN3C1=NN=C3",
    "6": "CC1=C(C)C(N2CCC(O{AR})CC2)=NN3C1=NN=C3",
}

# Table-1 rows whose substitution pattern is explicit in the prose and drawing.
ARYL = {
    "a": ("4-fluorophenyl", "c4ccc(F)cc4"),
    "b": ("3-fluorophenyl", "c4cc(F)ccc4"),
    "c": ("2-fluorophenyl", "c4c(F)cccc4"),
    "d": ("3,4-difluorophenyl", "c4cc(F)c(F)cc4"),
    "e": ("2,3-difluorophenyl", "c4c(F)c(F)ccc4"),
    "f": ("4-trifluoromethylphenyl", "c4ccc(C(F)(F)F)cc4"),
    "g": ("4-methylphenyl", "c4ccc(C)cc4"),
    "h": ("3-methylphenyl", "c4cc(C)ccc4"),
    "q": ("2,3-dihydrobenzo[b][1,4]dioxin-6-yl", "c4ccc5c(c4)OCCO5"),
}

POTENCY_NM = {
    "5a": 152, "6a": 131, "5b": 253, "6b": 237,
    "5c": 4030, "6c": 1990, "5d": 404, "6d": 327,
    "5e": 1640, "6e": 1400, "5f": 6020, "6f": 4050,
    "5g": 324, "6g": 280, "5h": 160, "6h": 159,
    "5q": 134, "6q": 38,
}


def canonical(smiles: str) -> tuple[str, Chem.Mol]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid reconstructed SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True), mol


def fp(mol: Chem.Mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for suffix, (aryl_name, aryl) in ARYL.items():
        for family in ("5", "6"):
            compound_id = family + suffix
            smi, mol = canonical(HEAD[family].format(AR=aryl))
            ec50 = float(POTENCY_NM[compound_id])
            rows.append({
                "compound_id": compound_id,
                "canonical_smiles": smi,
                "EC50_nM": ec50,
                "pEC50": 9.0 - np.log10(ec50),
                "censoring": "exact",
                "head_group": "7-methyl" if family == "5" else "7,8-dimethyl",
                "aryl_ether": aryl_name,
                "assay": "hM4/Gqi5-CHO calcium mobilization; ACh EC20",
                "functional_mode": "PAM",
                "source_doi": DOI,
                "structure_source": "Scheme 1 + Table 1 image/prose",
                "reconstruction_confidence": "high",
                "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
                "molecular_weight": Descriptors.MolWt(mol),
            })
    external = pd.DataFrame(rows).sort_values("compound_id").reset_index(drop=True)

    train = pd.read_csv(TRAIN).reset_index(drop=True)
    train_mols = [Chem.MolFromSmiles(s) for s in train.canonical_smiles]
    if any(m is None for m in train_mols):
        raise ValueError("invalid SMILES in frozen training benchmark")
    train_canonical = [Chem.MolToSmiles(m, canonical=True, isomericSmiles=True) for m in train_mols]
    train_fps = [fp(m) for m in train_mols]
    exact_set = set(train_canonical)

    exact, max_sim, nearest = [], [], []
    for smi in external.canonical_smiles:
        mol = Chem.MolFromSmiles(smi); query_fp = fp(mol)
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(query_fp, train_fps), float)
        j = int(np.argmax(sims))
        exact.append(smi in exact_set)
        max_sim.append(float(sims[j]))
        nearest.append(str(train.loc[j, "canonical_molecule_id"]))
    external["exact_training_overlap"] = exact
    external["max_train_tanimoto_ecfp4"] = max_sim
    external["nearest_train_id"] = nearest
    external["eligible_zero_shot"] = ~external.exact_training_overlap
    external.to_csv(OUT / "external_potency.csv", index=False)

    eligible = external[external.eligible_zero_shot]
    audit = {
        "dataset_id": "M4_VU6025733_2026_TABLE1_HIGHCONF_V01",
        "source_doi": DOI,
        "publication_date": "2026-01-25",
        "n_reconstructed": int(len(external)),
        "n_exact_training_overlap": int(external.exact_training_overlap.sum()),
        "n_zero_shot_eligible": int(len(eligible)),
        "label_range_pEC50": [float(external.pEC50.min()), float(external.pEC50.max())],
        "median_max_train_tanimoto": float(external.max_train_tanimoto_ecfp4.median()),
        "max_train_tanimoto": float(external.max_train_tanimoto_ecfp4.max()),
        "structure_policy": "Only uniquely reconstructable Scheme-1/Table-1 structures; ambiguous entries excluded.",
        "endpoint_match": "Exact functional endpoint family: hM4/Gqi5-CHO calcium mobilization with ACh EC20.",
        "claim_boundary": "Independent temporal SAR validation, not prospective wet-lab confirmation.",
        "csv_sha256": hashlib.sha256((OUT / "external_potency.csv").read_bytes()).hexdigest(),
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(external[["compound_id", "EC50_nM", "pEC50", "exact_training_overlap",
                    "max_train_tanimoto_ecfp4", "nearest_train_id"]].to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
