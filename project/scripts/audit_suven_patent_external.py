"""Freeze chemical-domain and source-QC metadata for the Suven holdout."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Draw, rdFingerprintGenerator


ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "project" / "data" / "benchmarks" / "m4_pam_v1"
EXT_DIR = BENCH / "external_suven_2025"
OCR = ROOT / "project" / "data" / "papers" / "WO2025099660A1" / "compound_tables" / "compound_row_ocr.csv"


def main() -> None:
    train = pd.read_csv(BENCH / "potency_molecules.csv")
    ext = pd.read_csv(EXT_DIR / "external_patent_potency.csv")
    # Idempotence: remove fields produced by a previous QC pass.
    ext = ext.drop(
        columns=[
            "exact_training_overlap", "max_train_tanimoto_ecfp4", "nearest_train_id",
            "eligible_external", "reported_mh_plus", "calculated_mh_plus",
            "mass_delta_Da", "mass_qc",
        ],
        errors="ignore",
    )
    train_mols = [Chem.MolFromSmiles(s) for s in train.canonical_smiles]
    ext_mols = [Chem.MolFromSmiles(s) for s in ext.smiles]
    if any(mol is None for mol in train_mols + ext_mols):
        raise ValueError("Invalid SMILES in frozen train/external set")

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    train_fp = [generator.GetFingerprint(mol) for mol in train_mols]
    ext_fp = [generator.GetFingerprint(mol) for mol in ext_mols]
    train_canonical = {Chem.MolToSmiles(mol, isomericSmiles=True) for mol in train_mols}
    max_similarity, nearest_id, exact = [], [], []
    for mol, fp in zip(ext_mols, ext_fp):
        similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(fp, train_fp), float)
        best = int(np.argmax(similarities))
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        max_similarity.append(float(similarities[best]))
        nearest_id.append(str(train.iloc[best].canonical_molecule_id))
        exact.append(canonical in train_canonical)
    ext["exact_training_overlap"] = exact
    ext["max_train_tanimoto_ecfp4"] = max_similarity
    ext["nearest_train_id"] = nearest_id
    ext["eligible_external"] = ~ext.exact_training_overlap

    ocr = pd.read_csv(OCR)
    mass = (
        ocr.loc[ocr.reported_mh_plus.notna(), ["compound_id", "reported_mh_plus"]]
        .drop_duplicates("compound_id", keep="first")
    )
    # Two values are plainly readable in the source but sit across page/row splits.
    mass = pd.concat(
        [mass, pd.DataFrame({"compound_id": [1, 41], "reported_mh_plus": [368.0, 382.0]})],
        ignore_index=True,
    ).drop_duplicates("compound_id", keep="last")
    ext = ext.merge(mass, on="compound_id", how="left")
    ext["calculated_mh_plus"] = [Descriptors.ExactMolWt(mol) + 1.007276 for mol in ext_mols]
    ext["mass_delta_Da"] = ext.reported_mh_plus - ext.calculated_mh_plus
    ext["mass_qc"] = np.where(
        ext.reported_mh_plus.isna(), "not_machine_recovered",
        np.where(ext.mass_delta_Da.abs() <= 0.7, "pass", "manual_review"),
    )
    # Example 88 is drawn/named as a constitutional isomer of example 87 and
    # therefore has the same formula.  The patent prints "351.(M+H)+" for 88
    # while 87 prints 352.2; retain, do not silently repair, this source typo.
    ext.loc[ext.compound_id == 88, "mass_qc"] = "source_mass_typo_suspected"
    ext.to_csv(EXT_DIR / "external_patent_potency.csv", index=False)

    legends = [f"Ex {row.compound_id} | CRE {row.cre_luc_ec50_nM:g} nM" for row in ext.itertuples()]
    image = Draw.MolsToGridImage(ext_mols, molsPerRow=5, subImgSize=(320, 250), legends=legends)
    image.save(EXT_DIR / "structure_qc_grid.png")

    audit = {
        "status": "source_and_chemical_domain_frozen_before_model_prediction",
        "n_external": len(ext),
        "n_exact_training_overlap": int(ext.exact_training_overlap.sum()),
        "n_eligible_external": int(ext.eligible_external.sum()),
        "median_max_train_tanimoto_ecfp4": float(ext.max_train_tanimoto_ecfp4.median()),
        "max_max_train_tanimoto_ecfp4": float(ext.max_train_tanimoto_ecfp4.max()),
        "n_machine_recovered_mass": int(ext.reported_mh_plus.notna().sum()),
        "n_mass_qc_pass": int((ext.mass_qc == "pass").sum()),
        "n_mass_manual_review": int((ext.mass_qc == "manual_review").sum()),
        "mass_manual_review_ids": ext.loc[ext.mass_qc == "manual_review", "compound_id"].astype(int).tolist(),
        "source_mass_typo_suspected_ids": ext.loc[
            ext.mass_qc == "source_mass_typo_suspected", "compound_id"
        ].astype(int).tolist(),
        "external_labels_used_for_training_or_model_selection": False,
        "model_predictions_run": False,
    }
    (EXT_DIR / "chemical_qc_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
