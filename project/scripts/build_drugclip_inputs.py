#!/usr/bin/env python3
"""Create small DrugCLIP LMDB inputs without downloading the 4.7 GB demo set."""
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import lmdb
import numpy as np
import pandas as pd


def write_lmdb(records: list[dict], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(output), subdir=False, map_size=2**34)
    with env.begin(write=True) as txn:
        for idx, record in enumerate(records):
            txn.put(str(idx).encode("ascii"), pickle.dumps(record))
    env.sync(); env.close()


def molecules(args: argparse.Namespace) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    table = pd.read_csv(args.csv)
    missing = {args.id_column, args.smiles_column} - set(table.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    records = []
    for _, row in table.iterrows():
        candidate_id, smiles = str(row[args.id_column]), str(row[args.smiles_column])
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES for {candidate_id}: {smiles}")
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3(); params.randomSeed = int(args.seed)
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=args.conformers, params=params))
        if not conf_ids:
            raise RuntimeError(f"No conformer generated for {candidate_id}")
        try:
            AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1)
        except Exception:
            pass
        mol = Chem.RemoveHs(mol)
        coordinates = [np.asarray(mol.GetConformer(i).GetPositions(), dtype=np.float32)
                       for i in range(mol.GetNumConformers())]
        records.append({
            "atoms": [atom.GetSymbol() for atom in mol.GetAtoms()],
            "coordinates": coordinates,
            "smi": Chem.MolToSmiles(mol, canonical=True),
            "mol": mol,
            "label": candidate_id,
        })
    write_lmdb(records, args.output)
    print(f"wrote {len(records)} molecules to {args.output}")


def pdb_atoms(path: Path, chain: str | None, resids: set[int]) -> tuple[list[str], list[np.ndarray]]:
    names, coords = [], []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        if chain and line[21].strip() != chain:
            continue
        try:
            resid = int(line[22:26])
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])],
                           dtype=np.float32)
        except ValueError:
            continue
        if resid in resids:
            names.append(line[12:16].strip())
            coords.append(xyz)
    return names, coords


def pocket(args: argparse.Namespace) -> None:
    resids = {int(x) for x in args.resids.split(",") if x.strip()}
    names, coords = pdb_atoms(args.pdb, args.chain or None, resids)
    if not names:
        raise RuntimeError("Pocket selection is empty; check chain and residue numbering")
    if len(names) > args.max_atoms:
        raise RuntimeError(f"Pocket has {len(names)} atoms > max {args.max_atoms}; predeclare a smaller set")
    record = {"pocket": args.name, "pocket_atoms": names, "pocket_coordinates": coords}
    write_lmdb([record], args.output)
    print(f"wrote pocket {args.name}: {len(names)} atoms, {len(resids)} requested residues")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    mol = sub.add_parser("molecules")
    mol.add_argument("--csv", type=Path, required=True)
    mol.add_argument("--output", type=Path, required=True)
    mol.add_argument("--id-column", default="candidate_id")
    mol.add_argument("--smiles-column", default="smiles")
    mol.add_argument("--conformers", type=int, default=10)
    mol.add_argument("--seed", type=int, default=20260922)
    mol.set_defaults(function=molecules)
    poc = sub.add_parser("pocket")
    poc.add_argument("--pdb", type=Path, required=True)
    poc.add_argument("--output", type=Path, required=True)
    poc.add_argument("--name", required=True)
    poc.add_argument("--chain", default="R")
    poc.add_argument("--resids", default="89,92,93,96,184,186,190,423,432,433,435,436,439")
    poc.add_argument("--max-atoms", type=int, default=256)
    poc.set_defaults(function=pocket)
    args = parser.parse_args(); args.function(args)


if __name__ == "__main__":
    main()
