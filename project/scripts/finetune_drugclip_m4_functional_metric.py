#!/usr/bin/env python3
"""Targeted internal DrugCLIP metric fine-tuning for M4 functional PAMs.

Unlike pocket-anchor tuning, this objective reshapes DrugCLIP's own molecule
projection so that functional PAMs from different literature/source components
are close while experimentally inactive molecules are separated.  Hard
negatives are selected with the *frozen* DrugCLIP geometry, preventing the
sampler from chasing its own changing scores.  An equal-size random-negative
arm is the causal control for the sampling strategy.
"""
from __future__ import annotations

import argparse
import copy
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from sklearn.metrics import average_precision_score, roc_auc_score
from torch import nn
from torch.nn import functional as F


def canonical(s):
    m = Chem.MolFromSmiles(str(s))
    if m is None: raise ValueError(s)
    return Chem.MolToSmiles(m, canonical=True)


class Projection(nn.Module):
    def __init__(self, state):
        super().__init__()
        self.linear1 = nn.Linear(state["linear1.weight"].shape[1], state["linear1.weight"].shape[0])
        self.linear2 = nn.Linear(state["linear2.weight"].shape[1], state["linear2.weight"].shape[0])
        self.load_state_dict(state)

    def forward(self, x):
        return F.normalize(self.linear2(F.relu(self.linear1(x))), dim=-1)


def metrics(y, p):
    return {"roc_auc": float(roc_auc_score(y, p)),
            "average_precision": float(average_precision_score(y, p)),
            "prevalence": float(np.mean(y))}


def neighbor_score(train_z, train_y, test_z, k):
    sim = test_z @ train_z.T
    pos = np.sort(sim[:, train_y == 1], axis=1)[:, -min(k, int((train_y == 1).sum())):]
    neg = np.sort(sim[:, train_y == 0], axis=1)[:, -min(k, int((train_y == 0).sum())):]
    return pos.mean(1) - neg.mean(1)


def fit(x, y, source, state, mode, seed, args):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    x = torch.as_tensor(x, dtype=torch.float32)
    y_t = torch.as_tensor(y, dtype=torch.long)
    source_t = torch.as_tensor(source, dtype=torch.long)
    model = Projection(copy.deepcopy(state))
    frozen = Projection(copy.deepcopy(state)).eval()
    for p in frozen.parameters(): p.requires_grad_(False)
    with torch.no_grad(): z0 = frozen(x)
    pos_idx = torch.nonzero(y_t == 1, as_tuple=False).reshape(-1)
    neg_idx = torch.nonzero(y_t == 0, as_tuple=False).reshape(-1)
    if not len(pos_idx) or not len(neg_idx): raise ValueError("Both classes required")
    # Fixed sampler: each PAM gets positives from another source and equally many
    # frozen-geometry hard or random inactive negatives.
    positive_sets, negative_sets = [], []
    generator = torch.Generator().manual_seed(seed)
    for anchor in pos_idx:
        candidates = pos_idx[source_t[pos_idx] != source_t[anchor]]
        if not len(candidates): candidates = pos_idx[pos_idx != anchor]
        original = z0[anchor] @ z0[candidates].T
        positive_sets.append(candidates[torch.topk(original, min(args.positive_k, len(candidates))).indices])
        if mode == "hard":
            order = torch.argsort(z0[anchor] @ z0[neg_idx].T, descending=True)
        else:
            order = torch.randperm(len(neg_idx), generator=generator)
        negative_sets.append(neg_idx[order[:min(args.negative_k, len(order))]])
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        opt.zero_grad(); z = model(x)
        losses = []
        for anchor, positives, negatives in zip(pos_idx, positive_sets, negative_sets):
            pos_sim = (z[anchor] @ z[positives].T).mean()
            neg_sim = torch.logsumexp((z[anchor] @ z[negatives].T) / args.temperature, 0) * args.temperature
            losses.append(F.relu(args.margin - pos_sim + neg_sim))
        triplet = torch.stack(losses).mean()
        preserve = (1 - (z * z0).sum(1)).mean()
        loss = triplet + args.preserve_weight * preserve
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); opt.step()
    return model


def bootstrap(y, a, b, n, seed=20260924):
    rng = np.random.default_rng(seed); out = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) == 2:
            out.append(roc_auc_score(y[idx], a[idx]) - roc_auc_score(y[idx], b[idx]))
    return list(map(float, np.quantile(out, [0.025, 0.5, 0.975])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--representations", type=Path, required=True)
    ap.add_argument("--projection", type=Path, required=True)
    ap.add_argument("--benchmark", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--protocol", choices=["source", "series"], required=True)
    ap.add_argument("--seeds", default="20260924,20260925,20260926")
    ap.add_argument("--epochs", type=int, default=120)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--weight-decay", type=float, default=1e-3)
    ap.add_argument("--margin", type=float, default=0.05)
    ap.add_argument("--temperature", type=float, default=0.05)
    ap.add_argument("--positive-k", type=int, default=4)
    ap.add_argument("--negative-k", type=int, default=8)
    ap.add_argument("--neighbor-k", type=int, default=5)
    ap.add_argument("--preserve-weight", type=float, default=1.0)
    ap.add_argument("--bootstrap", type=int, default=1000)
    args = ap.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    arc = np.load(args.representations, allow_pickle=False)
    initial = torch.load(args.projection, map_location="cpu")
    table = pd.read_csv(args.benchmark)
    ids = list(map(str, arc["molecule_ids"])); lookup = {v:i for i,v in enumerate(ids)}
    order = np.asarray([lookup[canonical(s)] for s in table.canonical_smiles])
    x = arc["molecule_representations"].astype(np.float32)[order]
    y = table.target.to_numpy(int)
    sources = pd.Categorical(table.source_component).codes.astype(int)
    fold_col = "source_fold" if args.protocol == "source" else "series_holdout_fold"
    folds = table[fold_col].to_numpy(int); assigned = folds >= 0
    modes = ["frozen", "random", "hard"]; rows = []
    for seed in map(int, args.seeds.split(",")):
        pred = {m: np.full(len(y), np.nan) for m in modes}
        for fold in sorted(np.unique(folds[assigned])):
            tr, te = assigned & (folds != fold), assigned & (folds == fold)
            frozen = Projection(copy.deepcopy(initial["mol_project"])).eval()
            with torch.no_grad():
                ztr0 = frozen(torch.as_tensor(x[tr])).numpy(); zte0 = frozen(torch.as_tensor(x[te])).numpy()
            pred["frozen"][te] = neighbor_score(ztr0, y[tr], zte0, args.neighbor_k)
            for mode in ["random", "hard"]:
                model = fit(x[tr], y[tr], sources[tr], initial["mol_project"], mode,
                            seed + 101 * int(fold), args)
                with torch.no_grad():
                    ztr = model(torch.as_tensor(x[tr])).numpy(); zte = model(torch.as_tensor(x[te])).numpy()
                pred[mode][te] = neighbor_score(ztr, y[tr], zte, args.neighbor_k)
        for i in np.flatnonzero(assigned):
            rows.append({"canonical_molecule_id": table.iloc[i].canonical_molecule_id,
                         "target": int(y[i]), "fold": int(folds[i]), "seed": seed,
                         **{m: float(pred[m][i]) for m in modes}})
    frame = pd.DataFrame(rows)
    ensemble = frame.groupby("canonical_molecule_id", sort=False)[modes].mean()
    truth = table.set_index("canonical_molecule_id").loc[ensemble.index, "target"].to_numpy(int)
    report = {"method": "DrugCLIP internal molecule projection cross-source functional metric tuning",
              "protocol": args.protocol, "n_assigned": int(assigned.sum()),
              "metrics": {m: metrics(truth, ensemble[m].to_numpy()) for m in modes},
              "hard_minus_random_auc_bootstrap_95ci": bootstrap(truth, ensemble.hard.to_numpy(),
                                                                 ensemble.random.to_numpy(), args.bootstrap),
              "hard_minus_frozen_auc_bootstrap_95ci": bootstrap(truth, ensemble.hard.to_numpy(),
                                                                 ensemble.frozen.to_numpy(), args.bootstrap),
              "claim_boundary": "Internal DrugCLIP projection fine-tune; retrospective grouped OOF, not wet-lab PAM proof."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    frame.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__": main()
