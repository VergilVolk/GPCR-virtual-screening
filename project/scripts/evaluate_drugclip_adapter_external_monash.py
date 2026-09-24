#!/usr/bin/env python3
"""Train frozen-DrugCLIP adapters on history and test on Monash allostery data."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from train_drugclip_m4_triplet_adapter import canonical, train_one


def load_features(path: Path, smiles: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    archive = np.load(path, allow_pickle=False)
    ids = list(map(str, archive["molecule_ids"]))
    lookup = {item: i for i, item in enumerate(ids)}
    canon = smiles.map(canonical)
    if missing := [item for item in canon if item not in lookup]:
        raise ValueError(f"Missing embeddings: {missing[:5]}")
    order = np.asarray([lookup[item] for item in canon])
    mol = archive["molecule_embeddings"].astype(np.float32)[order]
    scores = archive["scores"].astype(np.float32).T[order]
    return np.concatenate((mol, scores), axis=1), scores.max(1)


def bootstrap(y: np.ndarray, score: np.ndarray, repeats: int, seed: int) -> list[float]:
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    values = []
    for _ in range(repeats):
        idx = np.r_[rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)]
        values.append(float(roc_auc_score(y[idx], score[idx])))
    return list(map(float, np.quantile(values, [0.025, 0.975])))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-embeddings", type=Path, required=True)
    ap.add_argument("--train-benchmark", type=Path, required=True)
    ap.add_argument("--external-embeddings", type=Path, required=True)
    ap.add_argument("--external-benchmark", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--seeds", default="20260924,20260925,20260926")
    ap.add_argument("--bootstrap", type=int, default=5000)
    args = ap.parse_args()

    train = pd.read_csv(args.train_benchmark)
    external = pd.read_csv(args.external_benchmark)
    x_train, hardness = load_features(args.train_embeddings, train.canonical_smiles)
    x_external, external_hardness = load_features(args.external_embeddings, external.canonical_smiles)
    y_train = train.target.to_numpy(int)
    y_external = external.pam_label.to_numpy(int)
    source_names = sorted(train.source_component.unique())
    source_map = {name: i for i, name in enumerate(source_names)}
    source = train.source_component.map(source_map).to_numpy(int)
    seeds = [int(value) for value in args.seeds.split(",")]
    modes = ["bce", "triplet", "triplet_domain"]
    hp = SimpleNamespace(bottleneck=16, lr=1e-3, weight_decay=1e-3, epochs=250,
                         domain_weight=0.2, domain_loss_weight=0.2, margin=0.1,
                         triplet_weight=0.5, max_anchors=128, hard_negative_top_k=8)
    predictions = {mode: [] for mode in modes}
    for seed in seeds:
        for mode in modes:
            predictions[mode].append(train_one(x_train, y_train, source, hardness,
                                                x_external, mode, seed, hp))
    ensemble = {mode: np.mean(values, axis=0) for mode, values in predictions.items()}
    # The one exact historical overlap is excluded from the primary external estimate.
    primary = ~external.exact_history_overlap.astype(bool).to_numpy()
    report = {
        "evidence_level": "post_failure_exploratory_external_endpoint_transfer",
        "n_external": int(len(external)), "n_primary_nonoverlap": int(primary.sum()),
        "n_primary_positive": int(y_external[primary].sum()),
        "n_primary_negative": int((y_external[primary] == 0).sum()),
        "models": {},
        "claim_boundary": (
            "Independent publication but one narrow chemotype and only three inactive compounds; "
            "the set has been used previously in PACER analyses and is not an untouched confirmation set."
        ),
    }
    for offset, (mode, score) in enumerate(ensemble.items()):
        report["models"][mode] = {
            "roc_auc_nonoverlap": float(roc_auc_score(y_external[primary], score[primary])),
            "roc_auc_stratified_bootstrap_95ci": bootstrap(
                y_external[primary], score[primary], args.bootstrap, 20260924 + offset),
            "average_precision_nonoverlap": float(average_precision_score(y_external[primary], score[primary])),
        }
    # Unsupervised max-pocket DrugCLIP score is retained as a binding baseline.
    report["models"]["zero_shot_state_max"] = {
        "roc_auc_nonoverlap": float(roc_auc_score(y_external[primary], external_hardness[primary])),
        "roc_auc_stratified_bootstrap_95ci": bootstrap(
            y_external[primary], external_hardness[primary], args.bootstrap, 20260930),
        "average_precision_nonoverlap": float(average_precision_score(
            y_external[primary], external_hardness[primary])),
    }
    output = external[["compound_id", "pam_label", "exact_history_overlap"]].copy()
    for mode, score in ensemble.items(): output[mode] = score
    output["zero_shot_state_max"] = external_hardness
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    output.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
