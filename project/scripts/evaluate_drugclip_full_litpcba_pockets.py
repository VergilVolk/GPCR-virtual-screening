#!/usr/bin/env python3
"""Evaluate full LIT-PCBA with independently encoded DrugCLIP pocket vectors."""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import lmdb
import numpy as np

from evaluate_drugclip_official_litpcba_embeddings import metrics


def lmdb_count(path):
    env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
    with env.begin() as txn:
        count = txn.stat()["entries"]
    env.close()
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pocket-root", type=Path, required=True)
    parser.add_argument("--pocket-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    archive = np.load(args.pocket_archive, allow_pickle=False)
    all_pockets = np.asarray(archive["pocket_embeddings"], dtype=np.float32)
    pocket_files = sorted(args.pocket_root.glob("*.pockets.lmdb"))
    counts = [lmdb_count(path) for path in pocket_files]
    if sum(counts) != len(all_pockets):
        raise ValueError(f"Pocket count mismatch: manifest={sum(counts)} archive={len(all_pockets)}")
    offsets = np.cumsum([0] + counts)

    per_target = {}
    total = 0
    audit = {}
    for path, start, end in zip(pocket_files, offsets[:-1], offsets[1:]):
        target = path.name.removesuffix(".pockets.lmdb")
        molecule, _, labels = pickle.load(
            (args.root / target / "drugclip_emb" / "mols.lmdb.pkl").open("rb"))
        molecules = np.asarray(molecule, dtype=np.float32)
        labels = np.asarray(labels, dtype=np.int64)
        pockets = all_pockets[start:end]
        scores = (pockets @ molecules.T).max(axis=0)
        per_target[target] = metrics(labels, scores)
        total += len(labels)
        audit[target] = {"pocket_records": int(end - start), "candidate_rows": int(len(labels))}
        print(target, json.dumps(per_target[target]), flush=True)

    names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    macro = {name: float(np.mean([value[name] for value in per_target.values()])) for name in names}
    report = {
        "protocol": "Official molecule embeddings plus independently encoded official DrugCLIP pocket embeddings; max across pocket structures",
        "checkpoint": "project/tools/DrugCLIP/artifacts/checkpoint_best.pt",
        "n_targets": len(per_target),
        "n_pairs": total,
        "n_pocket_records": int(len(all_pockets)),
        "macro": macro,
        "per_target": per_target,
        "audit": audit,
        "claim_boundary": "Frozen pocket-molecule compatibility benchmark; no adaptation or PAM efficacy claim.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
