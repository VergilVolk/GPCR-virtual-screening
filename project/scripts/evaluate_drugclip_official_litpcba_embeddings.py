#!/usr/bin/env python3
"""Reproduce the official full LIT-PCBA DrugCLIP embedding benchmark."""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
from pathlib import Path

import numpy as np
from rdkit.ML.Scoring.Scoring import CalcBEDROC
from sklearn.metrics import average_precision_score, roc_auc_score


def enrichment(labels, scores, fraction):
    n = max(1, int(np.ceil(len(labels) * fraction)))
    top = np.argsort(-scores, kind="mergesort")[:n]
    return float(labels[top].mean() / labels.mean())


def metrics(labels, scores):
    ranked = [[float(score), int(label)] for score, label in sorted(
        zip(scores, labels), key=lambda item: item[0], reverse=True)]
    result = {
        "n": int(len(labels)),
        "actives": int(labels.sum()),
        "roc_auc": float(roc_auc_score(labels, scores)),
        "pr_auc": float(average_precision_score(labels, scores)),
        "bedroc_alpha80_5": float(CalcBEDROC(ranked, 1, 80.5)),
    }
    for fraction in (0.005, 0.01, 0.02, 0.05):
        result[f"ef{fraction:g}"] = enrichment(labels, scores, fraction)
    return result


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--source-zip", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    per_target = {}
    total_pairs = 0
    for folder in sorted(path for path in args.root.iterdir() if path.is_dir()):
        emb = folder / "drugclip_emb"
        molecule_path = emb / "mols.lmdb.pkl"
        pocket_path = emb / "ligand.lmdb.pkl"
        if not molecule_path.exists() or not pocket_path.exists():
            continue
        molecule_reps, molecule_names, labels = pickle.load(molecule_path.open("rb"))
        pocket_reps, _, _ = pickle.load(pocket_path.open("rb"))
        molecule_reps = np.asarray(molecule_reps, dtype=np.float32)
        pocket_reps = np.asarray(pocket_reps, dtype=np.float32)
        labels = np.asarray(labels, dtype=np.int64)
        if len(molecule_reps) != len(labels) or len(molecule_names) != len(labels):
            raise ValueError(f"Row mismatch for {folder.name}")
        if not np.isfinite(molecule_reps).all() or not np.isfinite(pocket_reps).all():
            raise ValueError(f"Non-finite embedding for {folder.name}")
        scores = (pocket_reps @ molecule_reps.T).max(axis=0)
        per_target[folder.name] = metrics(labels, scores)
        total_pairs += len(labels)
        print(folder.name, json.dumps(per_target[folder.name]), flush=True)

    metric_names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    macro = {name: float(np.mean([value[name] for value in per_target.values()])) for name in metric_names}
    report = {
        "protocol": "Official precomputed DrugCLIP embeddings; max score across target pocket structures",
        "source_zip": str(args.source_zip),
        "source_zip_bytes": args.source_zip.stat().st_size,
        "source_zip_sha256": sha256(args.source_zip),
        "n_targets": len(per_target),
        "n_pairs": total_pairs,
        "macro": macro,
        "per_target": per_target,
        "claim_boundary": "Official embedding reproduction only; no adaptation or PAM efficacy claim.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
