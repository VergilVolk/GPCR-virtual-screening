"""Build the frozen high-confidence WO2025122811 external benchmark."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


P = Path(__file__).resolve().parents[1]
ROOT = P / "data" / "benchmarks" / "m4_pam_v1"
SOURCE = ROOT / "external_acadia_2025"
TRAIN = ROOT / "potency_molecules.csv"


def canonical(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(smiles)
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def collapse_active_moieties(external: pd.DataFrame) -> pd.DataFrame:
    """Collapse salt/free-base or repeated preparations of one active moiety.

    Patent examples 12/13 and 32/33 are the same active structures prepared in
    different solid forms.  Keeping them as independent ML rows would permit
    exact-structure leakage between support and query.  Functional values are
    therefore aggregated at canonical active-moiety level.
    """
    rows = []
    for smiles, part in external.groupby("canonical_smiles", sort=False):
        first = part.iloc[0].copy()
        exact = part[part.hM4_PAM_EC50_qualifier.eq("=") & part.pEC50.notna()]
        first["compound_id"] = "|".join(part.compound_id.astype(str))
        first["example_id"] = "|".join(part.example_id.astype(int).astype(str))
        first["source_example_count"] = len(part)
        first["source_example_ids"] = first["example_id"]
        if len(exact):
            first["pEC50"] = float(exact.pEC50.median())
            first["hM4_PAM_EC50_nM"] = float(10 ** (9.0 - first["pEC50"]))
            first["hM4_PAM_EC50_qualifier"] = "="
        for column in ["hM4_PAM_RE_pct", "hM4_agonist_EC50_nM",
                       "hM4_agonist_RE_pct", "hM2_PAM_EC50_nM", "hM2_PAM_RE_pct"]:
            values = pd.to_numeric(part[column], errors="coerce").dropna()
            first[column] = float(values.median()) if len(values) else np.nan
        first["pam_re_ge_50"] = int(first.hM4_PAM_RE_pct >= 50)
        first["intrinsic_agonist_re_ge_50"] = int(first.hM4_agonist_RE_pct >= 50)
        first["eligible_exact_potency"] = bool(len(exact) and not first.exact_training_overlap)
        first["patent_name"] = " || ".join(part.patent_name.astype(str).unique())
        rows.append(first)
    return pd.DataFrame(rows).reset_index(drop=True)


def main() -> None:
    structures = pd.read_csv(SOURCE / "structures_opsin_audit.csv")
    structures = structures[structures.structure_status == "high_confidence"].copy()
    labels = pd.read_csv(SOURCE / "functional_labels_frozen.csv")
    external = structures.merge(labels, on="example_id", validate="one_to_one")
    external["compound_id"] = external.example_id.map(lambda value: f"ACADIA-{int(value):02d}")
    external["canonical_smiles"] = external.canonical_smiles.map(canonical)
    external["pEC50"] = np.where(
        external.hM4_PAM_EC50_qualifier.eq("="),
        9.0 - np.log10(external.hM4_PAM_EC50_nM.astype(float)), np.nan,
    )
    external["pam_re_ge_50"] = external.hM4_PAM_RE_pct.ge(50).astype(int)
    external["intrinsic_agonist_re_ge_50"] = external.hM4_agonist_RE_pct.ge(50).astype(int)

    train = pd.read_csv(TRAIN)
    train_smiles = train.canonical_smiles.map(canonical).tolist()
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    train_fps = [generator.GetFingerprint(Chem.MolFromSmiles(value)) for value in train_smiles]
    train_set = set(train_smiles)
    max_similarity, nearest_id = [], []
    for smiles in external.canonical_smiles:
        fp = generator.GetFingerprint(Chem.MolFromSmiles(smiles))
        similarities = np.asarray([DataStructs.TanimotoSimilarity(fp, item) for item in train_fps])
        nearest = int(np.argmax(similarities))
        max_similarity.append(float(similarities[nearest]))
        nearest_id.append(str(train.iloc[nearest].get("molecule_id", train.iloc[nearest].get("compound_id", nearest))))
    external["exact_training_overlap"] = external.canonical_smiles.isin(train_set)
    external["max_train_tanimoto_ecfp4"] = max_similarity
    external["nearest_train_id"] = nearest_id
    external["eligible_exact_potency"] = (
        external.hM4_PAM_EC50_qualifier.eq("=") & ~external.exact_training_overlap
    )

    columns = [
        "compound_id", "example_id", "canonical_smiles", "molecular_formula",
        "hM4_PAM_EC50_nM", "hM4_PAM_EC50_qualifier", "pEC50", "hM4_PAM_RE_pct",
        "pam_re_ge_50", "hM4_agonist_EC50_nM", "hM4_agonist_EC50_qualifier",
        "hM4_agonist_RE_pct", "intrinsic_agonist_re_ge_50", "hM2_PAM_EC50_nM",
        "hM2_PAM_EC50_qualifier", "hM2_PAM_RE_pct", "exact_training_overlap",
        "max_train_tanimoto_ecfp4", "nearest_train_id", "eligible_exact_potency",
        "structure_status", "mass_delta_da", "patent_name",
    ]
    output = SOURCE / "external_functional_benchmark.csv"
    external[columns].to_csv(output, index=False)
    molecule = collapse_active_moieties(external)
    molecule_columns = columns + ["source_example_count", "source_example_ids"]
    molecule_output = SOURCE / "external_molecule_benchmark.csv"
    molecule[molecule_columns].to_csv(molecule_output, index=False)
    duplicate_groups = external.groupby("canonical_smiles").filter(lambda x: len(x) > 1)
    duplicate_map = (duplicate_groups.groupby("canonical_smiles").compound_id
                     .apply(lambda x: x.astype(str).tolist()).to_dict())
    audit = {
        "dataset_id": "M4_ACADIA_WO2025122811_FUNCTIONAL_HIGHCONF_V01",
        "n_patent_examples": 49,
        "n_structure_high_confidence": int(len(external)),
        "n_exact_potency": int(external.eligible_exact_potency.sum()),
        "n_unique_active_moieties": int(len(molecule)),
        "n_unique_exact_potency": int(molecule.eligible_exact_potency.sum()),
        "n_duplicate_active_moiety_groups": int(len(duplicate_map)),
        "duplicate_active_moiety_groups": duplicate_map,
        "n_pam_re_below_50": int((external.pam_re_ge_50 == 0).sum()),
        "n_intrinsic_agonist_re_ge_50": int(external.intrinsic_agonist_re_ge_50.sum()),
        "n_exact_training_overlap": int(external.exact_training_overlap.sum()),
        "median_max_train_tanimoto_ecfp4": float(external.max_train_tanimoto_ecfp4.median()),
        "max_train_tanimoto_ecfp4": float(external.max_train_tanimoto_ecfp4.max()),
        "csv_sha256": hashlib.sha256(output.read_bytes()).hexdigest().upper(),
        "molecule_csv_sha256": hashlib.sha256(molecule_output.read_bytes()).hexdigest().upper(),
        "labels_used_for_method_selection": False,
        "claim_boundary": "External functional potency/efficacy benchmark; not mechanistic allostery or prospective PAM confirmation.",
    }
    (SOURCE / "external_benchmark_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
