#!/usr/bin/env python3
"""Create DrugCLIP-compatible crystal-ligand LMDBs for an external DUD-E subset.

This script is label-blind: it only converts each supplied crystal ligand into
one reference-ligand record. Candidate labels remain untouched in `mols.lmdb`
and are read only by the final evaluator.
"""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import lmdb
import numpy as np
from rdkit import Chem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    audit = {}
    for folder in sorted(path for path in args.root.iterdir() if path.is_dir()):
        ligand_path = folder / "crystal_ligand.mol2"
        # DUD-E MOL2 atom typing can violate RDKit's valence rules. DrugCLIP's
        # input path consumes only atom symbols and coordinates, so retain the
        # deposited geometry without chemistry re-sanitisation.
        mol = Chem.MolFromMol2File(str(ligand_path), sanitize=False, removeHs=False)
        if mol is None or not mol.GetNumConformers():
            raise RuntimeError(f"Unable to parse crystal ligand: {ligand_path}")
        record = {
            "atoms": [atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() != "H"],
            "coordinates": [np.asarray([
                (mol.GetConformer().GetAtomPosition(atom.GetIdx()).x,
                 mol.GetConformer().GetAtomPosition(atom.GetIdx()).y,
                 mol.GetConformer().GetAtomPosition(atom.GetIdx()).z)
                for atom in mol.GetAtoms() if atom.GetSymbol() != "H"
            ], dtype=np.float32)],
            "smi": f"DUD_E_CRYSTAL_{folder.name}",
            "mol": None,
            "label": 1,
        }
        output = folder / "reference.lmdb"
        env = lmdb.open(str(output), subdir=False, map_size=2**30)
        with env.begin(write=True) as txn:
            txn.put(b"0", pickle.dumps(record))
        env.sync()
        env.close()
        audit[folder.name] = {"reference": str(output), "atoms": len(record["atoms"]),
                              "smiles": record["smi"], "source": str(ligand_path)}
    report = {"targets": sorted(audit), "records": audit,
              "claim_boundary": "Crystal reference conversion only; no DUD-E assay labels read."}
    (args.root / "reference_build.audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
