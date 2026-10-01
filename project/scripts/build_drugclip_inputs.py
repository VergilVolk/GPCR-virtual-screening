#!/usr/bin/env python3
"""Create small DrugCLIP LMDB inputs without downloading the 4.7 GB demo set."""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import lmdb
import numpy as np
import pandas as pd


def prepare_molecule_task(task: tuple[str, str, int, int, bool]) -> tuple[dict | None, dict | None]:
    """Build one deterministic RDKit record; safe for multiprocessing."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    candidate_id, smiles, conformers, seed, skip_optimize = task
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, {"id": candidate_id, "smiles": smiles, "reason": "invalid_smiles"}
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3(); params.randomSeed = int(seed)
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=conformers, params=params))
    if not conf_ids:
        return None, {"id": candidate_id, "smiles": smiles, "reason": "no_conformer"}
    if not skip_optimize:
        try:
            AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1)
        except Exception:
            pass
    mol = Chem.RemoveHs(mol)
    coordinates = [np.asarray(mol.GetConformer(i).GetPositions(), dtype=np.float32)
                   for i in range(mol.GetNumConformers())]
    return ({"atoms": [atom.GetSymbol() for atom in mol.GetAtoms()],
             "coordinates": coordinates, "smi": Chem.MolToSmiles(mol, canonical=True),
             "mol": mol, "label": candidate_id}, None)


def write_lmdb(records: list[dict], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(output), subdir=False, map_size=2**34)
    with env.begin(write=True) as txn:
        for idx, record in enumerate(records):
            txn.put(str(idx).encode("ascii"), pickle.dumps(record))
    env.sync(); env.close()


def molecules(args: argparse.Namespace) -> None:
    table = pd.read_csv(args.csv)
    missing = {args.id_column, args.smiles_column} - set(table.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    records, failures = [], []
    tasks = [(str(row[args.id_column]), str(row[args.smiles_column]), args.conformers,
              args.seed, args.skip_optimize) for _, row in table.iterrows()]
    if args.workers == 1:
        iterator = map(prepare_molecule_task, tasks)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn"))
        iterator = executor.map(prepare_molecule_task, tasks, chunksize=args.worker_chunksize)
    for row_number, (record, failure) in enumerate(iterator, 1):
        if failure is not None:
            if not args.skip_failures:
                if failure["reason"] == "invalid_smiles":
                    raise ValueError(f"Invalid SMILES for {failure['id']}: {failure['smiles']}")
                raise RuntimeError(f"No conformer generated for {failure['id']}")
            failures.append(failure)
        else:
            records.append(record)
        if args.progress_every and row_number % args.progress_every == 0:
            print(f"processed={row_number}/{len(table)} written={len(records)} failed={len(failures)}", flush=True)
    if executor is not None:
        executor.shutdown()
    write_lmdb(records, args.output)
    audit = {"input_rows": len(table), "written": len(records), "failed": len(failures),
             "failures": failures}
    args.output.with_suffix(args.output.suffix + ".audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8")
    print(f"wrote {len(records)} molecules to {args.output}; failed={len(failures)}")


def pdb_atoms(path: Path, chain: str | None, resids: set[int], first_model_only: bool = False,
              exclude_hydrogens: bool = False) -> tuple[list[str], list[np.ndarray]]:
    names, coords = [], []
    seen_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            if seen_model and first_model_only:
                break
            seen_model = True
            continue
        if line.startswith("ENDMDL") and first_model_only:
            break
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
        atom_name = line[12:16].strip()
        element = line[76:78].strip() or atom_name[:1]
        if exclude_hydrogens and element.upper() == "H":
            continue
        if resid in resids:
            names.append(atom_name)
            coords.append(xyz)
    return names, coords


def pocket(args: argparse.Namespace) -> None:
    resids = {int(x) for x in args.resids.split(",") if x.strip()}
    names, coords = pdb_atoms(args.pdb, args.chain or None, resids,
                               args.first_model_only, args.exclude_hydrogens)
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
    mol.add_argument("--skip-failures", action="store_true",
                     help="Record invalid/no-conformer rows in an audit instead of aborting the batch")
    mol.add_argument("--skip-optimize", action="store_true",
                     help="Keep ETKDG coordinates without MMFF relaxation")
    mol.add_argument("--progress-every", type=int, default=1000)
    mol.add_argument("--workers", type=int, default=1)
    mol.add_argument("--worker-chunksize", type=int, default=16)
    mol.set_defaults(function=molecules)
    poc = sub.add_parser("pocket")
    poc.add_argument("--pdb", type=Path, required=True)
    poc.add_argument("--output", type=Path, required=True)
    poc.add_argument("--name", required=True)
    poc.add_argument("--chain", default="R")
    poc.add_argument("--resids", default="89,92,93,96,184,186,190,423,432,433,435,436,439")
    poc.add_argument("--max-atoms", type=int, default=256)
    poc.add_argument("--first-model-only", action="store_true")
    poc.add_argument("--exclude-hydrogens", action="store_true")
    poc.set_defaults(function=pocket)
    args = parser.parse_args(); args.function(args)


if __name__ == "__main__":
    main()
