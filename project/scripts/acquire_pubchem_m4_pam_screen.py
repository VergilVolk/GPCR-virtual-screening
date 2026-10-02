#!/usr/bin/env python3
"""Acquire a deterministic, auditable subset of PubChem CHRM4 PAM AID 624126.

All reported actives are retained.  Inactives are selected without structures
or model scores by SHA-256 ordering of CID, so acquisition is label-aware only
at the class-balancing level and cannot cherry-pick easy chemical negatives.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from urllib.request import urlopen

import pandas as pd
from rdkit import Chem


AID = 624126
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def get_json(url: str, attempts: int = 6) -> dict:
    for attempt in range(attempts):
        try:
            with urlopen(url, timeout=120) as response:
                return json.load(response)
        except Exception:
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def assay_cids(outcome: str) -> list[int]:
    data = get_json(f"{PUG}/assay/aid/{AID}/cids/JSON?cids_type={outcome}")
    return list(map(int, data["InformationList"]["Information"][0]["CID"]))


def fetch_properties(cids: list[int], batch_size: int = 200) -> list[dict]:
    rows = []
    for start in range(0, len(cids), batch_size):
        batch = cids[start:start + batch_size]
        joined = ",".join(map(str, batch))
        url = (f"{PUG}/compound/cid/{joined}/property/CanonicalSMILES,IsomericSMILES,"
               "InChIKey,MolecularFormula,MolecularWeight/JSON")
        data = get_json(url)
        rows.extend(data["PropertyTable"]["Properties"])
        if start:
            time.sleep(0.12)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--inactive", type=int, default=20000)
    args = parser.parse_args()
    active = assay_cids("active")
    inactive_all = assay_cids("inactive")
    inactive = sorted(
        inactive_all,
        key=lambda cid: hashlib.sha256(f"AID624126:{cid}".encode()).digest()
    )[:args.inactive]
    selected = [(cid, 1) for cid in active] + [(cid, 0) for cid in inactive]
    properties = {int(r["CID"]): r for r in fetch_properties([cid for cid, _ in selected])}
    rows = []
    for cid, label in selected:
        p = properties.get(cid)
        if not p:
            continue
        smiles = p.get("SMILES") or p.get("ConnectivitySMILES")
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            continue
        rows.append({
            "pubchem_cid": cid,
            "canonical_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            "inchi_key": p.get("InChIKey"), "molecular_formula": p.get("MolecularFormula"),
            "molecular_weight": pd.to_numeric(p.get("MolecularWeight"), errors="coerce"),
            "pam_primary_screen_label": label,
        })
    table = pd.DataFrame(rows).drop_duplicates("canonical_smiles", keep="first")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "aid624126_screen_subset.csv", index=False)
    audit = {
        "dataset_id": "PUBCHEM_AID624126_CHRM4_PAM_PRIMARY_V01",
        "source_url": "https://pubchem.ncbi.nlm.nih.gov/bioassay/624126",
        "assay": "Human CHRM4 PAM primary cell-based calcium screen, singlicate",
        "pubchem_active_cids": len(active),
        "pubchem_inactive_cids": len(inactive_all),
        "requested_inactive_subset": args.inactive,
        "selection": "all active CIDs plus SHA256(AID624126:CID)-ordered inactive CIDs",
        "retrieved_unique_structures": int(len(table)),
        "retrieved_active": int(table.pam_primary_screen_label.sum()),
        "retrieved_inactive": int((table.pam_primary_screen_label == 0).sum()),
        "label_quality": "screening_level",
        "claim_boundary": (
            "Large independent primary-screen benchmark, not confirmatory PAM ground truth. "
            "Singlicate hits may include assay artifacts and false positives."
        ),
    }
    (args.output_dir / "source.audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
