#!/usr/bin/env python3
"""Hybrid ECFP + frozen DrugCLIP molecule residual under strict target LOSO."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F

from evaluate_muscarinic_ecfp_pocket_baseline import fps


TARGETS = ["M1", "M2", "M3", "M4", "M5"]


class HybridMolecule(nn.Module):
    def __init__(self, seed: int, initial_drug_scale: float):
        super().__init__(); torch.manual_seed(seed)
        self.ecfp = nn.Linear(2048, 128, bias=False)
        nn.init.normal_(self.ecfp.weight, std=0.01)
        # softplus inverse; non-negative residual cannot arbitrarily invert the
        # frozen pretrained molecule geometry.
        raw = np.log(np.expm1(initial_drug_scale))
        self.raw_drug_scale = nn.Parameter(torch.tensor(float(raw)))

    def forward(self, x, frozen_drug):
        scale = F.softplus(self.raw_drug_scale)
        return F.normalize(self.ecfp(x) + scale * frozen_drug, dim=-1)


def fit_model(x, frozen_drug, pocket, frame, mi, pi, seed, args):
    model = HybridMolecule(seed, args.initial_drug_scale)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    mol = torch.as_tensor([mi[str(v)] for v in frame.canonical_smiles], dtype=torch.long)
    pos = torch.as_tensor([pi[str(v)] for v in frame.positive_subtype], dtype=torch.long)
    neg = torch.as_tensor([pi[str(v)] for v in frame.negative_subtype], dtype=torch.long)
    for _ in range(args.epochs):
        opt.zero_grad()
        score = model(x, frozen_drug) @ pocket.T
        loss = F.relu(args.margin - score[mol, pos] + score[mol, neg]).mean()
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); opt.step()
    return model


def pair_scores(model, x, frozen_drug, pocket, frame, mi, pi):
    with torch.inference_mode(): matrix = (model(x, frozen_drug) @ pocket.T).numpy()
    return np.asarray([matrix[mi[str(r.canonical_smiles)], pi[str(r.positive_subtype)]] -
                       matrix[mi[str(r.canonical_smiles)], pi[str(r.negative_subtype)]]
                       for r in frame.itertuples()])


def metrics(frame, score, held):
    preferred = frame.positive_subtype.to_numpy() == held
    disfavored = frame.negative_subtype.to_numpy() == held
    a = float((score[preferred] > 0).mean()); b = float((score[disfavored] > 0).mean())
    return {"accuracy": float((score > 0).mean()), "held_preferred_accuracy": a,
            "held_disfavored_accuracy": b, "direction_balanced_accuracy": (a + b) / 2}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--representations", type=Path, required=True)
    p.add_argument("--benchmark-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seeds", default="20260924,20260925,20260926")
    p.add_argument("--epochs", type=int, default=100)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight-decay", type=float, default=1e-2)
    p.add_argument("--margin", type=float, default=0.1)
    p.add_argument("--initial-drug-scale", type=float, default=0.1)
    args = p.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    arc = np.load(args.representations, allow_pickle=False)
    ids = list(map(str, arc["molecule_ids"])); mi = {v: i for i, v in enumerate(ids)}
    pids = [str(v).split("_")[0] for v in arc["pocket_ids"]]; pi = {v: i for i, v in enumerate(pids)}
    x = fps(ids)
    frozen_drug = F.normalize(torch.as_tensor(arc["molecule_embeddings"].astype(np.float32)), dim=-1)
    pocket = F.normalize(torch.as_tensor(arc["pocket_embeddings"].astype(np.float32)), dim=-1)
    reports, rows = {}, []
    for held in TARGETS:
        train_raw = pd.read_csv(args.benchmark_dir / f"{held}_train.csv")
        test_raw = pd.read_csv(args.benchmark_dir / f"{held}_test.csv")
        train = train_raw[train_raw.canonical_smiles.astype(str).isin(mi)].reset_index(drop=True)
        test = test_raw[test_raw.canonical_smiles.astype(str).isin(mi)].reset_index(drop=True)
        scores, scales = [], []
        for seed in map(int, args.seeds.split(",")):
            model = fit_model(x, frozen_drug, pocket, train, mi, pi, seed, args)
            value = pair_scores(model, x, frozen_drug, pocket, test, mi, pi)
            scores.append(value); scales.append(float(F.softplus(model.raw_drug_scale).detach()))
            for i, row in test.iterrows():
                rows.append({"held_subtype": held, "canonical_smiles": row.canonical_smiles,
                             "positive_subtype": row.positive_subtype, "negative_subtype": row.negative_subtype,
                             "seed": seed, "hybrid_delta": float(value[i])})
        mean = np.mean(scores, axis=0)
        reports[held] = {"train_pairs": int(len(train)), "test_pairs": int(len(test)),
                         "test_molecules": int(test.canonical_smiles.nunique()),
                         "hybrid": metrics(test, mean, held),
                         "per_seed_accuracy": [float((v > 0).mean()) for v in scores],
                         "learned_drug_scales": scales}
    macro = {key: float(np.mean([reports[h]["hybrid"][key] for h in TARGETS]))
             for key in ["accuracy", "direction_balanced_accuracy"]}
    report = {"method": "ECFP projection plus non-negative frozen DrugCLIP molecule residual; frozen DrugCLIP pockets",
              "split": "strict leave-one-muscarinic-target-out with molecule-disjoint folds",
              "folds": reports, "macro_average": macro,
              "hyperparameters": {"epochs": args.epochs, "lr": args.lr, "margin": args.margin,
                                      "initial_drug_scale": args.initial_drug_scale, "seeds": args.seeds},
              "claim_boundary": "Muscarinic subtype activity transfer; not PAM efficacy or binding affinity."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(rows).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
