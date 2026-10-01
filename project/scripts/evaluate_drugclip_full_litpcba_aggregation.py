#!/usr/bin/env python3
"""Nested leave-target-out selection of multi-pocket DrugCLIP aggregation."""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
from rdkit.ML.Scoring.Scoring import CalcBEDROC
from sklearn.metrics import average_precision_score, roc_auc_score


METHODS = ["max", "mean", "median", "top2_mean", "top3_mean", "lse_0.01", "lse_0.05", "lse_0.10"]


def aggregate(matrix, method):
    if method == "max":
        return matrix.max(axis=0)
    if method == "mean":
        return matrix.mean(axis=0)
    if method == "median":
        return np.median(matrix, axis=0)
    if method.startswith("top"):
        k = min(int(method[3]), len(matrix))
        return np.partition(matrix, len(matrix) - k, axis=0)[-k:].mean(axis=0)
    tau = float(method.split("_")[1])
    peak = matrix.max(axis=0)
    return peak + tau * np.log(np.exp((matrix - peak) / tau).mean(axis=0))


def metrics(labels, scores):
    ranked = [[float(score), int(label)] for score, label in sorted(
        zip(scores, labels), key=lambda item: item[0], reverse=True)]
    n = max(1, int(np.ceil(len(labels) * 0.01)))
    top = np.argsort(-scores, kind="mergesort")[:n]
    return {
        "roc_auc": float(roc_auc_score(labels, scores)),
        "pr_auc": float(average_precision_score(labels, scores)),
        "bedroc_alpha80_5": float(CalcBEDROC(ranked, 1, 80.5)),
        "ef0.01": float(labels[top].mean() / labels.mean()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--selection-metric", default="bedroc_alpha80_5",
                        choices=["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.01"])
    args = parser.parse_args()

    all_metrics, score_cache, labels_cache = {}, {}, {}
    folders = sorted(path for path in args.root.iterdir()
                     if path.is_dir() and (path / "drugclip_emb" / "mols.lmdb.pkl").exists())
    for folder in folders:
        molecule, _, labels = pickle.load((folder / "drugclip_emb" / "mols.lmdb.pkl").open("rb"))
        pocket, _, _ = pickle.load((folder / "drugclip_emb" / "ligand.lmdb.pkl").open("rb"))
        matrix = np.asarray(pocket, dtype=np.float32) @ np.asarray(molecule, dtype=np.float32).T
        labels = np.asarray(labels, dtype=np.int64)
        labels_cache[folder.name] = labels
        score_cache[folder.name] = {}
        all_metrics[folder.name] = {}
        for method in METHODS:
            score = aggregate(matrix, method)
            score_cache[folder.name][method] = score
            all_metrics[folder.name][method] = metrics(labels, score)
        print(f"computed {folder.name}", flush=True)

    targets = [folder.name for folder in folders]
    nested = {}
    selected_counts = {method: 0 for method in METHODS}
    for held in targets:
        seen = [target for target in targets if target != held]
        validation = {
            method: float(np.mean([all_metrics[target][method][args.selection_metric] for target in seen]))
            for method in METHODS
        }
        selected = max(METHODS, key=lambda method: (validation[method], -METHODS.index(method)))
        selected_counts[selected] += 1
        nested[held] = {
            "selected_method": selected,
            "seen_target_validation": validation,
            "held_metrics": all_metrics[held][selected],
        }

    names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.01"]
    nested_macro = {name: float(np.mean([nested[target]["held_metrics"][name] for target in targets])) for name in names}
    fixed_macro = {
        method: {name: float(np.mean([all_metrics[target][method][name] for target in targets])) for name in names}
        for method in METHODS
    }
    report = {
        "protocol": "Outer leave-one-target-out aggregation selection; held labels absent from method selection",
        "selection_metric": args.selection_metric,
        "n_targets": len(targets),
        "nested_macro": nested_macro,
        "fixed_method_macro": fixed_macro,
        "selected_method_counts": selected_counts,
        "folds": nested,
        "claim_boundary": "Aggregation baseline on official embeddings; no representation training or PAM efficacy claim.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
