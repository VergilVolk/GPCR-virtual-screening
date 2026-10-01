#!/usr/bin/env python3
"""Nested target-LOO fusion of DrugCLIP pocket and reference-ligand retrieval."""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import lmdb
import numpy as np
from scipy.stats import rankdata

from evaluate_drugclip_official_litpcba_embeddings import metrics


WEIGHTS = [0.0, 0.25, 0.5, 0.75, 1.0]


def pocket_map(root, archive_path):
    archive = np.load(archive_path, allow_pickle=False)
    embeddings = np.asarray(archive["pocket_embeddings"], dtype=np.float32)
    files = sorted(root.glob("*.pockets.lmdb"))
    counts = []
    for path in files:
        env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
        with env.begin() as txn:
            counts.append(txn.stat()["entries"])
        env.close()
    offsets = np.cumsum([0] + counts)
    return {path.name.removesuffix(".pockets.lmdb"): embeddings[start:end]
            for path, start, end in zip(files, offsets[:-1], offsets[1:])}


def normalize(scores, mode):
    if mode == "raw":
        return scores
    if mode == "zscore":
        return (scores - scores.mean()) / max(scores.std(), 1e-8)
    return rankdata(scores, method="average").astype(np.float32) / len(scores)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pocket-root", type=Path, required=True)
    parser.add_argument("--pocket-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--normalization", choices=["raw", "zscore", "rank"], default="rank")
    parser.add_argument("--selection-metric", choices=["bedroc_alpha80_5", "ef0.01", "roc_auc", "pr_auc"], default="bedroc_alpha80_5")
    args = parser.parse_args()

    pockets = pocket_map(args.pocket_root, args.pocket_archive)
    cache, fixed = {}, {}
    for target in sorted(pockets):
        folder = args.root / target / "drugclip_emb"
        molecules, _, labels = pickle.load((folder / "mols.lmdb.pkl").open("rb"))
        references, _, _ = pickle.load((folder / "ligand.lmdb.pkl").open("rb"))
        molecules = np.asarray(molecules, dtype=np.float32)
        labels = np.asarray(labels, dtype=np.int64)
        pocket_score = (pockets[target] @ molecules.T).max(axis=0)
        reference_score = (np.asarray(references, dtype=np.float32) @ molecules.T).max(axis=0)
        pocket_norm = normalize(pocket_score, args.normalization)
        reference_norm = normalize(reference_score, args.normalization)
        cache[target] = {"labels": labels, "pocket": pocket_score, "reference": reference_score, "fusion": {}}
        fixed[target] = {}
        for weight in WEIGHTS:
            score = (1 - weight) * pocket_norm + weight * reference_norm
            cache[target]["fusion"][str(weight)] = score
            fixed[target][str(weight)] = metrics(labels, score)
        print(f"computed {target}", flush=True)

    targets = sorted(cache)
    folds, selected_counts = {}, {str(weight): 0 for weight in WEIGHTS}
    for held in targets:
        seen = [target for target in targets if target != held]
        validation = {
            str(weight): float(np.mean([fixed[target][str(weight)][args.selection_metric] for target in seen]))
            for weight in WEIGHTS
        }
        selected = max(WEIGHTS, key=lambda weight: (validation[str(weight)], -weight))
        selected_counts[str(selected)] += 1
        folds[held] = {
            "selected_reference_weight": selected,
            "seen_target_validation": validation,
            "held_metrics": fixed[held][str(selected)],
        }

    metric_names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    nested_macro = {name: float(np.mean([folds[target]["held_metrics"][name] for target in targets])) for name in metric_names}
    fixed_macro = {
        str(weight): {name: float(np.mean([fixed[target][str(weight)][name] for target in targets])) for name in metric_names}
        for weight in WEIGHTS
    }
    report = {
        "protocol": "Outer leave-one-target-out selection of pocket/reference-ligand fusion weight",
        "normalization": args.normalization,
        "selection_metric": args.selection_metric,
        "n_targets": len(targets),
        "nested_macro": nested_macro,
        "fixed_weight_macro": fixed_macro,
        "selected_weight_counts": selected_counts,
        "folds": folds,
        "claim_boundary": "Reference-ligand-assisted retrospective screening; not ligand-free zero-shot, PAM efficacy, or prospective validation.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
