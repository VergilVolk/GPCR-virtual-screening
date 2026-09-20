"""Build the Monash LY2033298 mechanistic allostery benchmark.

The source is ChEMBL document CHEMBL3769348 (Szabo et al., 2015).
ChEMBL stores both raw and transformed copies of several measurements; this
builder keeps one explicitly named transformed endpoint per molecule and never
converts missing/not-active entries into invented continuous values.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


P = Path(__file__).resolve().parents[1]
RAW = P / "data" / "papers" / "ly2033298_2015" / "chembl_document_activities.json"
HISTORY = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
OUT = P / "data" / "benchmarks" / "m4_pam_v1" / "external_monash_ly2033298"

COMPOUND_NAMES = {
    "CHEMBL3770346": "8a_LY2033298",
    "CHEMBL3770870": "8b_O-ethyl",
    "CHEMBL3769740": "8c_O-butyl",
    "CHEMBL3770397": "8d_O-hexyl",
    "CHEMBL3770667": "8e_5-bromo",
    "CHEMBL3770686": "8f_5-iodo",
    "CHEMBL3770968": "9a_N-ethyl",
    "CHEMBL3770360": "9b_N-butyl",
    "CHEMBL3770771": "9c_N-hexyl",
    "CHEMBL3769684": "10a_N-acetyl",
    "CHEMBL3770719": "10b_N-butanoyl",
    "CHEMBL3769529": "10c_N-hexanoyl",
    "CHEMBL3770312": "12_core-variant",
    "CHEMBL2313377": "18_core-variant",
    "CHEMBL3769542": "21_core-variant",
    "CHEMBL1212988": "VU10004",
}

FAMILIES = {
    "CHEMBL3770346": "reference",
    "CHEMBL3770870": "O-alkyl",
    "CHEMBL3769740": "O-alkyl",
    "CHEMBL3770397": "O-alkyl",
    "CHEMBL3770667": "5-halogen",
    "CHEMBL3770686": "5-halogen",
    "CHEMBL3770968": "N-alkyl",
    "CHEMBL3770360": "N-alkyl",
    "CHEMBL3770771": "N-alkyl",
    "CHEMBL3769684": "N-acyl-inactive",
    "CHEMBL3770719": "N-acyl-inactive",
    "CHEMBL3769529": "N-acyl-inactive",
    "CHEMBL3770312": "core-variant",
    "CHEMBL2313377": "core-variant",
    "CHEMBL3769542": "core-variant",
    "CHEMBL1212988": "core-variant",
}

ENDPOINTS = {
    "CHEMBL3772172": ("functional_pKB", "pKb"),
    "CHEMBL3772173": ("log_tauB", "log(IA)"),
    "CHEMBL3772174": ("log_alpha_beta", "log(activity)"),
    "CHEMBL3772176": ("log_alpha", "log(activity)"),
}


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def one_value(rows: pd.DataFrame, assay: str, standard_type: str) -> float:
    values = rows.loc[
        (rows.assay_chembl_id == assay) & (rows.standard_type == standard_type),
        "standard_value",
    ].dropna().astype(float).unique()
    if len(values) > 1:
        raise ValueError(f"Conflicting {assay}/{standard_type} values: {values}")
    return float(values[0]) if len(values) == 1 else np.nan


def binding_pki(rows: pd.DataFrame) -> float:
    # ChEMBL stores pKi in the raw `type/value` pair while standardizing the
    # corresponding record to Ki in nM. Select the raw pKi copy explicitly.
    values = rows.loc[
        (rows.assay_chembl_id == "CHEMBL3772175") & (rows.type == "pKi"),
        "value",
    ].dropna().astype(float).unique()
    if len(values) > 1:
        raise ValueError(f"Conflicting binding pKi values: {values}")
    return float(values[0]) if len(values) == 1 else np.nan


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = json.loads(RAW.read_text(encoding="utf-8"))
    activities = pd.DataFrame(payload["activities"])
    if activities.document_chembl_id.nunique() != 1 or activities.document_chembl_id.iloc[0] != "CHEMBL3769348":
        raise ValueError("Unexpected ChEMBL document")

    rows = []
    for molecule_id, group in activities.groupby("molecule_chembl_id", sort=True):
        smiles = group.canonical_smiles.dropna().unique()
        if len(smiles) != 1:
            raise ValueError(f"Expected one structure for {molecule_id}: {smiles}")
        row = {
            "compound_id": COMPOUND_NAMES[molecule_id],
            "molecule_chembl_id": molecule_id,
            "canonical_smiles": canonical(smiles[0]),
            "modification_family": FAMILIES[molecule_id],
        }
        for assay, (name, kind) in ENDPOINTS.items():
            row[name] = one_value(group, assay, kind)
        row["binding_pKi"] = binding_pki(group)
        status = group.loc[group.assay_chembl_id == "CHEMBL3772174", "activity_comment"].dropna()
        row["functional_status"] = "inactive" if (status == "Not Active").any() else "active"
        row["pam_label"] = int(row["functional_status"] == "active")
        rows.append(row)

    frame = pd.DataFrame(rows).sort_values("compound_id").reset_index(drop=True)
    history = pd.read_csv(HISTORY)
    history_smiles = [canonical(item) for item in history.canonical_smiles]
    history_set = set(history_smiles)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    history_fp = [generator.GetFingerprint(Chem.MolFromSmiles(item)) for item in history_smiles]
    exact, similarities, nearest = [], [], []
    for smiles in frame.canonical_smiles:
        fp = generator.GetFingerprint(Chem.MolFromSmiles(smiles))
        sims = DataStructs.BulkTanimotoSimilarity(fp, history_fp)
        index = int(np.argmax(sims))
        exact.append(smiles in history_set)
        similarities.append(float(sims[index]))
        nearest.append(str(history.iloc[index].canonical_molecule_id))
    frame["exact_history_overlap"] = exact
    frame["max_history_tanimoto_ecfp4"] = similarities
    frame["nearest_history_id"] = nearest

    output = OUT / "external_allostery.csv"
    frame.to_csv(output, index=False)
    audit = {
        "status": "frozen_external_endpoint_transfer_benchmark",
        "source_document": "CHEMBL3769348",
        "doi": "10.1039/C5MD00334B",
        "n_molecules": int(len(frame)),
        "n_functional_active": int(frame.pam_label.sum()),
        "n_functional_inactive": int((1 - frame.pam_label).sum()),
        "n_complete_functional_triplet": int(frame[["functional_pKB", "log_alpha_beta", "log_tauB"]].notna().all(axis=1).sum()),
        "n_binding_alpha": int(frame.log_alpha.notna().sum()),
        "n_exact_history_overlap": int(frame.exact_history_overlap.sum()),
        "median_max_history_tanimoto": float(frame.max_history_tanimoto_ecfp4.median()),
        "raw_sha256": hashlib.sha256(RAW.read_bytes()).hexdigest(),
        "csv_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "missing_policy": "Not Active remains a categorical inactive; no continuous endpoint is imputed.",
        "claim_boundary": "Independent publication and endpoint transfer, but not an unseen-scaffold benchmark; one exact historical overlap and strong LY2033298-neighborhood similarity are disclosed.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(frame.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
