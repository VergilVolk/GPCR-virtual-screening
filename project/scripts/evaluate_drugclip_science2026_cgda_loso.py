#!/usr/bin/env python3
"""Strict target-LOO CGDA on full Science-2026 DrugCLIP LIT-PCBA embeddings."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from evaluate_drugclip_cgda_loso import evaluate_score, fit_fold
from evaluate_drugclip_official_litpcba_embeddings import metrics


def load_archives(root: Path) -> dict:
    data = {}
    for candidate_path in sorted(root.glob("*.npz")):
        if candidate_path.name.endswith(".references.npz"):
            continue
        target = candidate_path.stem
        reference_path = root / f"{target}.references.npz"
        metric_path = root / f"{target}.metrics.json"
        if not reference_path.exists() or not metric_path.exists():
            continue
        candidate = np.load(candidate_path, allow_pickle=False)
        reference = np.load(reference_path, allow_pickle=False)
        required = {"molecule_embeddings", "pocket_embeddings", "labels"}
        if not required.issubset(candidate.files):
            continue
        data[target] = {
            "molecules": candidate["molecule_embeddings"].astype(np.float32),
            "labels": candidate["labels"].astype(np.int64),
            "pockets": candidate["pocket_embeddings"].astype(np.float32),
            "references": reference["reference_embeddings"].astype(np.float32),
        }
        data[target]["reference_score"] = (
            data[target]["references"] @ data[target]["molecules"].T
        ).max(axis=0)
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--embedding-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minimum-targets", type=int, default=15)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--rank", type=int, default=2)
    parser.add_argument("--experts", type=int, default=2)
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--max-negatives-per-target", type=int, default=1000)
    parser.add_argument("--hard-negative-fraction", type=float, default=0.5)
    parser.add_argument("--ranking-weight", type=float, default=0.2)
    parser.add_argument("--ranking-margin", type=float, default=0.1)
    parser.add_argument("--lr", type=float, default=2e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.12)
    parser.add_argument("--retrieval-temperature", type=float, default=0.07)
    parser.add_argument("--retrieval-weight", type=float, default=0.2)
    parser.add_argument("--preserve-weight", type=float, default=0.5)
    parser.add_argument("--validation-fraction", type=float, default=0.2)
    parser.add_argument("--chunk-size", type=int, default=65536)
    args = parser.parse_args()

    data = load_archives(args.embedding_root)
    if len(data) < args.minimum_targets:
        missing = sorted(
            p.name.removesuffix(".references.npz")
            for p in args.embedding_root.glob("*.references.npz")
            if p.name.removesuffix(".references.npz") not in data
        )
        raise RuntimeError(
            f"Full-data gate failed: {len(data)}/{args.minimum_targets} complete targets; "
            f"missing candidate archives for {missing}"
        )

    targets = sorted(data)
    folds = {}
    score_dir = args.output.with_suffix("")
    score_dir.mkdir(parents=True, exist_ok=True)
    for fold_index, held in enumerate(targets):
        model, training = fit_fold(data, held, args, args.seed + fold_index * 1009)
        adapted, gate = evaluate_score(model, data[held], args.chunk_size)
        labels = data[held]["labels"]
        pocket = (data[held]["pockets"] @ data[held]["molecules"].T).max(axis=0)
        reference = data[held]["reference_score"]
        folds[held] = {
            "n": int(len(labels)),
            "positives": int(labels.sum()),
            "pocket": metrics(labels, pocket),
            "reference": metrics(labels, reference),
            "cgda": metrics(labels, adapted),
            "gate": gate.tolist(),
            "training": training,
        }
        np.savez_compressed(
            score_dir / f"{held}.npz", labels=labels, pocket=pocket,
            reference=reference, cgda=adapted,
        )
        print(f"completed={held}", flush=True)

    names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    macro = {
        method: {name: float(np.mean([folds[t][method][name] for t in targets])) for name in names}
        for method in ("pocket", "reference", "cgda")
    }
    report = {
        "method": "CGDA on verified Science-2026 DrugCLIP embeddings",
        "protocol": "Strict leave-one-entire-target-out; held labels evaluation-only",
        "targets": targets,
        "n_candidates": int(sum(len(data[t]["labels"]) for t in targets)),
        "macro": macro,
        "folds": folds,
        "parameters": vars(args) | {"embedding_root": str(args.embedding_root), "output": str(args.output)},
        "claim_boundary": "Reference-assisted retrospective binding screen; not PAM function or prospective proof.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(macro, indent=2))


if __name__ == "__main__":
    main()
