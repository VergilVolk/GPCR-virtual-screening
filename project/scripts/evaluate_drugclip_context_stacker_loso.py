#!/usr/bin/env python3
"""Strict target-LOO context stacking for official DrugCLIP LIT-PCBA.

The model never sees assay labels from the held target.  It combines frozen
DrugCLIP similarities to the held target's pocket ensemble and co-crystal
reference-ligand ensemble.  Training is target-balanced and negative-sampled;
evaluation covers every molecule in the held target.
"""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
from scipy.stats import rankdata
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from evaluate_drugclip_full_litpcba_reference_fusion import pocket_map
from evaluate_drugclip_official_litpcba_embeddings import metrics


def row_stats(context: np.ndarray, molecules: np.ndarray) -> np.ndarray:
    sim = np.asarray(context, dtype=np.float32) @ molecules.T
    ordered = np.sort(sim, axis=0)
    top1 = ordered[-1]
    top2 = ordered[-min(2, len(ordered)) :].mean(axis=0)
    top3 = ordered[-min(3, len(ordered)) :].mean(axis=0)
    return np.column_stack([top1, top2, top3, sim.mean(axis=0), sim.std(axis=0)])


def percentile_columns(values: np.ndarray) -> np.ndarray:
    n = len(values)
    return np.column_stack([
        rankdata(values[:, index], method="average") / n for index in range(values.shape[1])
    ]).astype(np.float32)


def features(pockets: np.ndarray, references: np.ndarray, molecules: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pocket = row_stats(pockets, molecules)
    reference = row_stats(references, molecules)
    p_rank = percentile_columns(pocket)
    r_rank = percentile_columns(reference)
    cross = np.column_stack([
        np.minimum(p_rank[:, 0], r_rank[:, 0]),
        p_rank[:, 0] * r_rank[:, 0],
        np.abs(p_rank[:, 0] - r_rank[:, 0]),
        p_rank[:, 0] - r_rank[:, 0],
        (p_rank[:, 0] + r_rank[:, 0]) / 2,
    ])
    return np.column_stack([p_rank, r_rank, cross]).astype(np.float32), pocket[:, 0], reference[:, 0]


def balanced_indices(labels: np.ndarray, max_negatives: int, rng: np.random.Generator) -> np.ndarray:
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    count = min(len(negative), max(max_negatives, len(positive) * 20))
    chosen = rng.choice(negative, size=count, replace=False)
    return np.concatenate([positive, chosen])


def fit_model(x: np.ndarray, y: np.ndarray, target_ids: np.ndarray, args):
    # Equal total weight per training target and equal class mass within target.
    weight = np.zeros(len(y), dtype=np.float64)
    for target in np.unique(target_ids):
        mask = target_ids == target
        for label in (0, 1):
            subgroup = mask & (y == label)
            if subgroup.any():
                weight[subgroup] = 0.5 / (mask.sum() * subgroup.mean())
    if args.model == "linear":
        scaler = StandardScaler().fit(x)
        model = LogisticRegression(C=args.c, max_iter=2000, random_state=args.seed)
        model.fit(scaler.transform(x), y, sample_weight=weight)
        return lambda value: model.predict_proba(scaler.transform(value))[:, 1]
    model = HistGradientBoostingClassifier(
        learning_rate=args.learning_rate,
        max_iter=args.max_iter,
        max_leaf_nodes=args.max_leaf_nodes,
        l2_regularization=args.l2,
        min_samples_leaf=args.min_samples_leaf,
        random_state=args.seed,
    )
    model.fit(x, y, sample_weight=weight)
    return lambda value: model.predict_proba(value)[:, 1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pocket-root", type=Path, required=True)
    parser.add_argument("--pocket-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", choices=["linear", "histgb"], default="histgb")
    parser.add_argument("--max-negatives-per-target", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--c", type=float, default=1.0)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--max-iter", type=int, default=150)
    parser.add_argument("--max-leaf-nodes", type=int, default=7)
    parser.add_argument("--min-samples-leaf", type=int, default=40)
    parser.add_argument("--l2", type=float, default=1.0)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    pockets = pocket_map(args.pocket_root, args.pocket_archive)
    cache = {}
    for target in sorted(pockets):
        folder = args.root / target / "drugclip_emb"
        molecules, _, labels = pickle.load((folder / "mols.lmdb.pkl").open("rb"))
        references, _, _ = pickle.load((folder / "ligand.lmdb.pkl").open("rb"))
        x, pocket_score, reference_score = features(
            pockets[target], np.asarray(references, dtype=np.float32), np.asarray(molecules, dtype=np.float32)
        )
        cache[target] = {
            "x": x,
            "labels": np.asarray(labels, dtype=np.int64),
            "pocket": pocket_score,
            "reference": reference_score,
        }
        print(f"features {target}: {len(labels):,}", flush=True)

    rng = np.random.default_rng(args.seed)
    targets = sorted(cache)
    folds = {}
    score_dir = args.output.with_suffix("")
    score_dir.mkdir(parents=True, exist_ok=True)
    for held in targets:
        train_x, train_y, train_target = [], [], []
        for target_id, target in enumerate(targets):
            if target == held:
                continue
            idx = balanced_indices(cache[target]["labels"], args.max_negatives_per_target, rng)
            train_x.append(cache[target]["x"][idx])
            train_y.append(cache[target]["labels"][idx])
            train_target.append(np.full(len(idx), target_id, dtype=np.int64))
        predictor = fit_model(np.concatenate(train_x), np.concatenate(train_y), np.concatenate(train_target), args)
        score = predictor(cache[held]["x"])
        labels = cache[held]["labels"]
        folds[held] = {
            "n": int(len(labels)),
            "positives": int(labels.sum()),
            "pocket": metrics(labels, cache[held]["pocket"]),
            "reference": metrics(labels, cache[held]["reference"]),
            "stacked": metrics(labels, score),
        }
        np.savez_compressed(score_dir / f"{held}.npz", labels=labels, stacked=score,
                            pocket=cache[held]["pocket"], reference=cache[held]["reference"])
        print(f"completed {held}: BEDROC={folds[held]['stacked']['bedroc_alpha80_5']:.4f}", flush=True)

    names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    macro = {method: {name: float(np.mean([folds[t][method][name] for t in targets])) for name in names}
             for method in ("pocket", "reference", "stacked")}
    report = {
        "protocol": "Strict leave-one-entire-target-out target-balanced context stacker",
        "model": args.model,
        "n_targets": len(targets),
        "macro": macro,
        "folds": folds,
        "parameters": vars(args) | {"root": str(args.root), "pocket_root": str(args.pocket_root),
                                    "pocket_archive": str(args.pocket_archive), "output": str(args.output)},
        "claim_boundary": "Reference-ligand-assisted retrospective screening. Held-target labels are evaluation-only; held co-crystal reference ligands are target context, not assay actives. This does not predict PAM efficacy.",
    }
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(macro, indent=2))


if __name__ == "__main__":
    main()
