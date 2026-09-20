"""Build the 2023/2024 M4 multi-pathway operational-allostery panel."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


P = Path(__file__).resolve().parents[1]
RAW = P / "data" / "benchmarks" / "m4_pam_v1" / "mechanistic_chembl" / "raw_mechanistic_activities.json"
MONASH = P / "data" / "benchmarks" / "m4_pam_v1" / "external_monash_ly2033298" / "external_allostery.csv"
OUT = P / "data" / "benchmarks" / "m4_pam_v1" / "mechanistic_chembl"

ASSAYS = {
    "CHEMBL5623659": "pKB",
    "CHEMBL5623660": "binding_log_alpha",
    "CHEMBL5623661": "cAMP_log_tauB",
    "CHEMBL5623662": "cAMP_log_alpha_beta",
    "CHEMBL5623663": "arrestin_log_tauB",
    "CHEMBL5623664": "arrestin_log_alpha_beta",
    "CHEMBL5526937": "pKB",
    "CHEMBL5526938": "GoB_log_tauB",
    "CHEMBL5526940": "GoB_log_alpha_beta",
    "CHEMBL5526942": "cAMP_log_tauB",
    "CHEMBL5526944": "cAMP_log_alpha_beta",
    "CHEMBL5526946": "arrestin_log_tauB",
    "CHEMBL5526948": "arrestin_log_alpha_beta",
}
DOCUMENTS = {
    "CHEMBL5620424": "Jorg2023",
    "CHEMBL5523758": "Liu2024",
}


def parent_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(smiles)
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    parent = max(fragments, key=lambda item: item.GetNumHeavyAtoms())
    return Chem.MolToSmiles(parent, isomericSmiles=True)


def main() -> None:
    payload = json.loads(RAW.read_text(encoding="utf-8"))
    data = pd.DataFrame(payload["activities"])
    data = data[data.document_chembl_id.isin(DOCUMENTS) & data.assay_chembl_id.isin(ASSAYS)].copy()
    data["endpoint"] = data.assay_chembl_id.map(ASSAYS)
    data["document"] = data.document_chembl_id.map(DOCUMENTS)
    data["value_numeric"] = pd.to_numeric(data.standard_value, errors="coerce")
    if data.value_numeric.isna().any():
        raise ValueError("Unexpected missing operational-allostery value")
    duplicates = data.groupby(["document", "molecule_chembl_id", "endpoint"]).value_numeric.nunique()
    if (duplicates > 1).any():
        raise ValueError(f"Conflicting duplicated endpoints: {duplicates[duplicates > 1]}")

    identity = data.groupby(["document", "document_chembl_id", "molecule_chembl_id"], as_index=False).agg(
        canonical_smiles=("canonical_smiles", "first")
    )
    identity["canonical_smiles"] = identity.canonical_smiles.map(parent_smiles)
    values = data.pivot_table(index=["document", "document_chembl_id", "molecule_chembl_id"],
                              columns="endpoint", values="value_numeric", aggfunc="first").reset_index()
    frame = identity.merge(values, on=["document", "document_chembl_id", "molecule_chembl_id"], validate="one_to_one")
    tau = [item for item in frame.columns if item.endswith("log_tauB")]
    cooperativity = [item for item in frame.columns if item.endswith("log_alpha_beta")]
    frame["consensus_log_tauB"] = frame[tau].mean(axis=1, skipna=True)
    frame["consensus_log_alpha_beta"] = frame[cooperativity].mean(axis=1, skipna=True)

    # Cross-document and Monash similarity audits are label-independent.
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [generator.GetFingerprint(Chem.MolFromSmiles(item)) for item in frame.canonical_smiles]
    exact_other, max_other = [], []
    for index, row in frame.iterrows():
        other = np.where(frame.document.to_numpy() != row.document)[0]
        sims = DataStructs.BulkTanimotoSimilarity(fps[index], [fps[item] for item in other])
        max_other.append(float(max(sims)))
        exact_other.append(any(frame.iloc[item].canonical_smiles == row.canonical_smiles for item in other))
    frame["exact_other_document_overlap"] = exact_other
    frame["max_other_document_tanimoto"] = max_other

    monash = pd.read_csv(MONASH)
    monash_fp = [generator.GetFingerprint(Chem.MolFromSmiles(item)) for item in monash.canonical_smiles]
    frame["max_monash_tanimoto"] = [max(DataStructs.BulkTanimotoSimilarity(fp, monash_fp)) for fp in fps]
    output = OUT / "mechanistic_multidomain.csv"
    frame.sort_values(["document", "molecule_chembl_id"]).to_csv(output, index=False)
    audit = {
        "status": "frozen_multidomain_mechanistic_panel",
        "n_molecules": int(len(frame)),
        "by_document": frame.document.value_counts().to_dict(),
        "n_exact_cross_document_overlap": int(frame.exact_other_document_overlap.sum() // 2),
        "median_cross_document_tanimoto": float(frame.max_other_document_tanimoto.median()),
        "median_monash_tanimoto": float(frame.max_monash_tanimoto.median()),
        "endpoints": sorted(ASSAYS.values()),
        "csv_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }
    (OUT / "mechanistic_multidomain_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(frame.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
