#!/usr/bin/env python3
"""Build a leakage-audited M4R allosteric domain-adaptation benchmark.

The source active/decoy files are the public M4R benchmark used by the local
GaMD ensemble replication.  Functional PACER-M4 molecules are excluded before
sampling so this stage cannot see the downstream functional benchmark.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def read_smi(path: Path, label: int) -> pd.DataFrame:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.strip().split()
        if not fields: continue
        mol = Chem.MolFromSmiles(fields[0])
        if mol is None: continue
        rows.append((Chem.MolToSmiles(mol, canonical=True), fields[1] if len(fields) > 1 else "", label))
    return pd.DataFrame(rows, columns=["canonical_smiles", "source_id", "target"])


def scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    value = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    return value or f"ACYCLIC:{hashlib.sha1(smiles.encode()).hexdigest()[:12]}"


def balanced_scaffold_folds(frame: pd.DataFrame, n_folds: int) -> dict[str, int]:
    groups = frame.groupby("murcko_scaffold").agg(n=("target", "size"), pos=("target", "sum"))
    groups = groups.assign(neg=groups.n-groups.pos).sort_values(["n", "pos"], ascending=False)
    totals = np.zeros((n_folds, 3), dtype=float); assignment = {}
    for name, row in groups.iterrows():
        costs = []
        for fold in range(n_folds):
            trial = totals.copy(); trial[fold] += [row.n, row.pos, row.neg]
            costs.append(np.std(trial[:, 0]) + 2*np.std(trial[:, 1]) + 2*np.std(trial[:, 2]))
        choice = int(np.argmin(costs)); assignment[name] = choice
        totals[choice] += [row.n, row.pos, row.neg]
    return assignment


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--actives", type=Path, required=True)
    ap.add_argument("--decoys", type=Path, required=True)
    ap.add_argument("--functional", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--decoys-per-active", type=int, default=5)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260924)
    args = ap.parse_args()

    active = read_smi(args.actives, 1); decoy = read_smi(args.decoys, 0)
    raw = {"actives": len(active), "decoys": len(decoy)}
    functional = pd.read_csv(args.functional)
    excluded = set(functional.canonical_smiles.map(lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(str(s)), canonical=True)))
    active_overlap = active.canonical_smiles.isin(excluded); decoy_overlap = decoy.canonical_smiles.isin(excluded)
    active = active[~active_overlap].drop_duplicates("canonical_smiles")
    decoy = decoy[~decoy_overlap].drop_duplicates("canonical_smiles")
    conflicts = set(active.canonical_smiles) & set(decoy.canonical_smiles)
    active = active[~active.canonical_smiles.isin(conflicts)]
    decoy = decoy[~decoy.canonical_smiles.isin(conflicts)]
    rng = np.random.default_rng(args.seed)
    n_decoy = min(len(decoy), args.decoys_per_active * len(active))
    decoy = decoy.iloc[np.sort(rng.choice(len(decoy), n_decoy, replace=False))]
    frame = pd.concat([active, decoy], ignore_index=True)
    frame.insert(0, "canonical_molecule_id", [f"M4RDOM{i:06d}" for i in range(len(frame))])
    frame["murcko_scaffold"] = frame.canonical_smiles.map(scaffold)
    assignment = balanced_scaffold_folds(frame, args.folds)
    frame["scaffold_fold"] = frame.murcko_scaffold.map(assignment).astype(int)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)
    audit = {
        "source_counts": raw, "functional_exact_overlap_excluded": {
            "active": int(active_overlap.sum()), "decoy": int(decoy_overlap.sum())},
        "conflicting_labels_excluded": len(conflicts), "output_rows": len(frame),
        "output_actives": int(frame.target.sum()), "output_decoys": int((frame.target == 0).sum()),
        "unique_scaffolds": int(frame.murcko_scaffold.nunique()),
        "folds": frame.groupby("scaffold_fold").target.agg(["size", "sum"]).reset_index().to_dict("records"),
        "rules": "Exact functional molecules excluded; labels deduplicated; deterministic random subset of public property-matched decoys; entire Murcko scaffolds held in one fold.",
        "claim_boundary": "Broad M4R allosteric active-vs-decoy membership, not ACh-dependent functional PAM efficacy."
    }
    args.output.with_suffix(".audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__": main()
