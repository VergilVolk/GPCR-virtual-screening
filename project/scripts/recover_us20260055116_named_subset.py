"""Recover the explicitly named US20260055116A1 M4 PAM subset.

Only compounds with a complete systematic name in a synthesis heading are
eligible.  Image-only structures are deliberately excluded from the primary
benchmark instead of being hand-transcribed.
"""

from __future__ import annotations

import hashlib
import html
import json
import re
import urllib.parse
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "benchmarks" / "m4_pam_v1"
SOURCE = ROOT / "external_us20260055116"
HTML = SOURCE / "US20260055116A1_google.html"
PDF = SOURCE / "US20260055116A1.pdf"
TRAIN = ROOT / "potency_molecules.csv"
OPSIN = "https://opsin.ch.cam.ac.uk/opsin/{}.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def headings(text: str):
    output = []
    for raw in re.findall(r"<heading[^>]*>(.*?)</heading>", text, re.S):
        value = html.unescape(re.sub(r"<[^>]+>", "", raw)).strip()
        match = re.search(r"Synthesis of Compound (I-\d+):\s*(.+)$", value, re.I)
        if match:
            name = re.sub(r"\s+", " ", match.group(2)).strip()
            name = re.sub(r"\s+as a Racemic Mixture of Stereoisomers$", "", name, flags=re.I)
            output.append((match.group(1).upper(), name))
    return output


def opsin(name: str):
    url = OPSIN.format(urllib.parse.quote(name, safe=""))
    response = requests.get(url, timeout=60)
    if response.status_code != 200:
        return None, f"HTTP_{response.status_code}"
    payload = response.json()
    smiles = payload.get("smiles")
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None, "RDKIT_PARSE_FAILED"
    return Chem.MolToSmiles(mol, isomericSmiles=True), "OK"


def functional_table() -> pd.DataFrame:
    table = pd.read_html(HTML)[49].copy()
    table.columns = ["compound_id", "human_m4_perk", "rat_m4_perk", "human_m4_gtpgs"]
    table = table[table.compound_id.astype(str).str.match(r"I-\d+$")].copy()
    for column in table.columns[1:]:
        table[column] = table[column].astype(str).str.upper().replace({"NAN": "N/A"})
    return table


def main() -> None:
    text = HTML.read_text(encoding="utf-8")
    label = functional_table()
    rows = []
    for compound_id, name in headings(text):
        smiles, status = opsin(name)
        mass = np.nan
        if smiles:
            mass = float(Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles)) + 1.007276)
        rows.append({"compound_id": compound_id, "systematic_name": name,
                     "canonical_smiles": smiles, "opsin_status": status,
                     "calculated_MH_plus": mass})
    structures = pd.DataFrame(rows).drop_duplicates("compound_id").merge(
        label, on="compound_id", how="left", validate="one_to_one")

    train = pd.read_csv(TRAIN)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    train_mols = [Chem.MolFromSmiles(value) for value in train.canonical_smiles]
    train_fps = [generator.GetFingerprint(mol) for mol in train_mols]
    train_set = {Chem.MolToSmiles(mol, isomericSmiles=True) for mol in train_mols}
    max_sim, nearest = [], []
    for smiles in structures.canonical_smiles:
        if not isinstance(smiles, str):
            max_sim.append(np.nan); nearest.append(""); continue
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
        fp = generator.GetFingerprint(Chem.MolFromSmiles(canonical))
        similarities = np.asarray([DataStructs.TanimotoSimilarity(fp, x) for x in train_fps])
        idx = int(np.argmax(similarities))
        max_sim.append(float(similarities[idx]))
        nearest.append(str(train.iloc[idx].canonical_molecule_id))
    structures["max_train_tanimoto_ecfp4"] = max_sim
    structures["nearest_train_id"] = nearest
    structures["exact_train_overlap"] = structures.canonical_smiles.isin(train_set)
    structures["high_confidence_named_structure"] = structures.opsin_status.eq("OK")
    structures.to_csv(SOURCE / "named_structure_functional_subset.csv", index=False)

    valid = structures[structures.high_confidence_named_structure]
    duplicate = valid[valid.duplicated("canonical_smiles", keep=False)]
    counts = {column: valid[column].value_counts(dropna=False).to_dict()
              for column in ["human_m4_perk", "rat_m4_perk", "human_m4_gtpgs"]}
    audit = {
        "dataset_id": "US20260055116_NAMED_FUNCTIONAL_SUBSET_V01",
        "pdf_sha256": sha256(PDF), "html_sha256": sha256(HTML),
        "n_table23_rows": int(len(label)), "n_explicitly_named": int(len(structures)),
        "n_opsin_rdkit_success": int(len(valid)),
        "n_unique_active_moieties": int(valid.canonical_smiles.nunique()),
        "n_duplicate_structure_rows": int(len(duplicate)),
        "n_exact_training_overlap": int(valid.exact_train_overlap.sum()),
        "median_max_train_tanimoto": float(valid.max_train_tanimoto_ecfp4.median()),
        "max_train_tanimoto": float(valid.max_train_tanimoto_ecfp4.max()),
        "endpoint_class_counts": counts,
        "image_only_compounds_excluded": True,
        "labels_visible_during_source_qualification": True,
        "claim_boundary": "Named high-confidence ordinal external stress test; not pristine blind validation or exact-potency evidence.",
    }
    (SOURCE / "named_subset_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(structures.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
