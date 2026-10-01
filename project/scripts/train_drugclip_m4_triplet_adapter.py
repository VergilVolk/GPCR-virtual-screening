#!/usr/bin/env python3
"""Source-aware CPU adapter for frozen DrugCLIP representations.

The outer test folds are never used for tuning. Three preregistered ablations
are evaluated with repeated seeds: BCE only, cross-source functional triplets,
and triplets plus source-adversarial regularisation. This is a retrospective
functional-transfer experiment, not experimental PAM validation.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from sklearn.metrics import average_precision_score, balanced_accuracy_score, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler
from torch import nn
from torch.nn import functional as F


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)


class ReverseGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x: torch.Tensor, weight: float) -> torch.Tensor:
        ctx.weight = weight
        return x.view_as(x)

    @staticmethod
    def backward(ctx, gradient: torch.Tensor):
        return -ctx.weight * gradient, None


class ResidualAdapter(nn.Module):
    def __init__(self, n_input: int, bottleneck: int, n_sources: int):
        super().__init__()
        self.down = nn.Linear(n_input, bottleneck, bias=False)
        self.up = nn.Linear(bottleneck, n_input, bias=False)
        nn.init.zeros_(self.up.weight)
        self.classifier = nn.Linear(n_input, 1)
        self.source_classifier = nn.Linear(n_input, n_sources)

    def forward(self, x: torch.Tensor, adversarial_weight: float = 0.0):
        z = F.normalize(x + 0.1 * self.up(F.gelu(self.down(x))), dim=1)
        functional = self.classifier(z).squeeze(1)
        source = self.source_classifier(ReverseGradient.apply(z, adversarial_weight))
        return z, functional, source


def make_triplets(y: np.ndarray, source: np.ndarray, hardness: np.ndarray,
                  rng: np.random.Generator, max_anchors: int, top_k: int) -> np.ndarray:
    positives = np.flatnonzero(y == 1)
    negatives = np.flatnonzero(y == 0)
    if len(positives) > max_anchors:
        positives = rng.choice(positives, max_anchors, replace=False)
    triplets = []
    for anchor in positives:
        positive_pool = np.flatnonzero((y == 1) & (source != source[anchor]))
        if not len(positive_pool):
            continue
        same_source_negative = negatives[source[negatives] == source[anchor]]
        negative_pool = same_source_negative if len(same_source_negative) else negatives
        hard_order = negative_pool[np.argsort(hardness[negative_pool])[::-1]]
        hard_order = hard_order[:max(1, min(top_k, len(hard_order)))]
        triplets.append((anchor, int(rng.choice(positive_pool)), int(rng.choice(hard_order))))
    return np.asarray(triplets, dtype=np.int64)


def train_one(x_train: np.ndarray, y_train: np.ndarray, source_train: np.ndarray,
              hardness_train: np.ndarray, x_test: np.ndarray, mode: str,
              seed: int, args) -> np.ndarray:
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    scaler = StandardScaler().fit(x_train)
    xt = torch.tensor(scaler.transform(x_train), dtype=torch.float32)
    xv = torch.tensor(scaler.transform(x_test), dtype=torch.float32)
    yt = torch.tensor(y_train, dtype=torch.float32)
    st = torch.tensor(source_train, dtype=torch.long)
    model = ResidualAdapter(xt.shape[1], args.bottleneck, int(source_train.max()) + 1)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    pos_weight = torch.tensor(float((y_train == 0).sum() / max(1, (y_train == 1).sum())))
    rng = np.random.default_rng(seed)
    use_triplet = mode in {"triplet", "triplet_domain"}
    use_domain = mode == "triplet_domain"

    for _ in range(args.epochs):
        model.train(); optimizer.zero_grad()
        z, logit, source_logit = model(xt, args.domain_weight if use_domain else 0.0)
        loss = F.binary_cross_entropy_with_logits(logit, yt, pos_weight=pos_weight)
        if use_triplet:
            tri = make_triplets(y_train, source_train, hardness_train, rng,
                                args.max_anchors, args.hard_negative_top_k)
            if len(tri):
                a, p, n = map(torch.tensor, tri.T)
                dap = 1 - (z[a] * z[p]).sum(1)
                dan = 1 - (z[a] * z[n]).sum(1)
                loss = loss + args.triplet_weight * F.relu(args.margin + dap - dan).mean()
        if use_domain:
            loss = loss + args.domain_loss_weight * F.cross_entropy(source_logit, st)
        loss.backward(); optimizer.step()

    model.eval()
    with torch.no_grad():
        return torch.sigmoid(model(xv)[1]).numpy()


def metrics(y: np.ndarray, p: np.ndarray) -> dict:
    hard = p >= 0.5
    return {
        "roc_auc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
        "prevalence": float(y.mean()),
        "ap_lift_over_prevalence": float(average_precision_score(y, p) - y.mean()),
        "balanced_accuracy_0p5": float(balanced_accuracy_score(y, hard)),
        "mcc_0p5": float(matthews_corrcoef(y, hard)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings", type=Path, required=True)
    ap.add_argument("--benchmark", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--protocol", choices=["source", "series"], default="source")
    ap.add_argument("--seeds", default="20260924,20260925,20260926")
    ap.add_argument("--epochs", type=int, default=250)
    ap.add_argument("--bottleneck", type=int, default=16)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--weight-decay", type=float, default=1e-3)
    ap.add_argument("--margin", type=float, default=0.1)
    ap.add_argument("--triplet-weight", type=float, default=0.5)
    ap.add_argument("--domain-weight", type=float, default=0.2)
    ap.add_argument("--domain-loss-weight", type=float, default=0.2)
    ap.add_argument("--max-anchors", type=int, default=128)
    ap.add_argument("--hard-negative-top-k", type=int, default=8)
    args = ap.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    archive = np.load(args.embeddings, allow_pickle=False)
    embedded_ids = list(map(str, archive["molecule_ids"]))
    mol = archive["molecule_embeddings"].astype(np.float32)
    scores = archive["scores"].astype(np.float32).T
    table = pd.read_csv(args.benchmark)
    table["drugclip_smiles"] = table.canonical_smiles.map(canonical)
    lookup = {s: i for i, s in enumerate(embedded_ids)}
    if missing := table.loc[~table.drugclip_smiles.isin(lookup), "canonical_molecule_id"].tolist():
        raise ValueError(f"Missing embeddings: {missing[:5]}")
    order = np.asarray([lookup[s] for s in table.drugclip_smiles])
    # Static pocket scores are included as three low-dimensional binding priors.
    x = np.concatenate((mol[order], scores[order]), axis=1)
    y = table.target.to_numpy(np.int64)
    source_names = sorted(table.source_component.unique())
    source_map = {name: i for i, name in enumerate(source_names)}
    source = table.source_component.map(source_map).to_numpy(np.int64)
    hardness = scores[order].max(1)
    folds = table.source_fold.to_numpy(int) if args.protocol == "source" else table.series_holdout_fold.to_numpy(int)
    assigned = folds >= 0
    seeds = [int(v) for v in args.seeds.split(",")]
    modes = ["bce", "triplet", "triplet_domain"]
    seed_reports, prediction_rows = [], []

    for seed in seeds:
        mode_predictions = {mode: np.full(len(y), np.nan) for mode in modes}
        fold_audit = []
        for fold in sorted(np.unique(folds[assigned])):
            train = assigned & (folds != fold)
            test = assigned & (folds == fold)
            if len(np.unique(y[train])) != 2 or len(np.unique(y[test])) != 2:
                raise ValueError(f"Fold {fold} has a single class")
            fold_audit.append({"fold": int(fold), "n_train": int(train.sum()), "n_test": int(test.sum()),
                               "train_sources": int(len(np.unique(source[train]))),
                               "test_sources": int(len(np.unique(source[test])))})
            # Remap training source IDs to a dense range for the adversarial head.
            train_sources = sorted(np.unique(source[train]))
            dense = {old: i for i, old in enumerate(train_sources)}
            dense_source = np.asarray([dense[v] for v in source[train]], dtype=np.int64)
            for mode in modes:
                mode_predictions[mode][test] = train_one(
                    x[train], y[train], dense_source, hardness[train], x[test], mode,
                    seed + int(fold) * 100, args)
        result = {mode: metrics(y[assigned], p[assigned]) for mode, p in mode_predictions.items()}
        seed_reports.append({"seed": seed, "metrics": result, "fold_audit": fold_audit})
        for i in np.flatnonzero(assigned):
            prediction_rows.append({"canonical_molecule_id": table.iloc[i].canonical_molecule_id,
                                    "target": int(y[i]), "fold": int(folds[i]), "seed": seed,
                                    **{mode: float(mode_predictions[mode][i]) for mode in modes}})

    aggregate = {}
    for mode in modes:
        values = [entry["metrics"][mode]["roc_auc"] for entry in seed_reports]
        aggregate[mode] = {"roc_auc_mean_across_seeds": float(np.mean(values)),
                           "roc_auc_min": float(np.min(values)), "roc_auc_max": float(np.max(values)),
                           "per_seed": [entry["metrics"][mode] for entry in seed_reports]}
    report = {
        "evidence_level": "retrospective_frozen_drugclip_source_aware_adapter",
        "encoder_update": "none_frozen", "protocol": args.protocol,
        "n_assigned": int(assigned.sum()), "n_positive": int(y[assigned].sum()),
        "n_negative": int((y[assigned] == 0).sum()), "seeds": seeds,
        "hyperparameters": {k: v for k, v in vars(args).items()
                            if k not in {"embeddings", "benchmark", "output"}},
        "aggregate": aggregate, "seed_reports": seed_reports,
        "claim_boundary": "No test-fold tuning; retrospective labels; not experimental PAM validation.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(prediction_rows).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps({"protocol": args.protocol, "aggregate": aggregate}, indent=2))


if __name__ == "__main__":
    main()
