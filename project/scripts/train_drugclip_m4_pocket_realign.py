#!/usr/bin/env python3
"""M4-specific DrugCLIP pocket--molecule re-alignment on frozen embeddings.

DrugCLIP's dense-retrieval geometry is retained while small residual adapters
learn a multi-positive M4 allosteric pocket ensemble. Functional inactive
molecules with high pretrained pocket similarity are explicit hard negatives.
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
from torch import nn
from torch.nn import functional as F


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None: raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)


class ReverseGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, weight): ctx.weight = weight; return x.view_as(x)
    @staticmethod
    def backward(ctx, grad): return -ctx.weight * grad, None


class LowRankResidual(nn.Module):
    def __init__(self, dim: int, rank: int):
        super().__init__()
        self.down, self.up = nn.Linear(dim, rank, bias=False), nn.Linear(rank, dim, bias=False)
        nn.init.zeros_(self.up.weight)

    def forward(self, x): return F.normalize(x + 0.1 * self.up(F.gelu(self.down(x))), dim=-1)


class PocketRealign(nn.Module):
    def __init__(self, dim: int, rank: int, n_sources: int):
        super().__init__()
        self.molecule_adapter = LowRankResidual(dim, rank)
        self.pocket_adapter = LowRankResidual(dim, rank)
        self.log_scale = nn.Parameter(torch.tensor(2.3))
        self.bias = nn.Parameter(torch.tensor(0.0))
        self.source_head = nn.Linear(dim, n_sources)

    def forward(self, molecule, pocket, domain_weight=0.0):
        zm, zp = self.molecule_adapter(molecule), self.pocket_adapter(pocket)
        state_similarity = zm @ zp.T
        # Soft maximum keeps all three M4 states differentiable without forcing
        # every PAM to prefer the same receptor state.
        retrieval = 0.1 * torch.logsumexp(state_similarity / 0.1, dim=1)
        logit = F.softplus(self.log_scale) * retrieval + self.bias
        source_logit = self.source_head(ReverseGradient.apply(zm, domain_weight))
        return zm, zp, state_similarity, retrieval, logit, source_logit


def metrics(y, p):
    hard = p >= 0.5
    return {"roc_auc": float(roc_auc_score(y, p)),
            "average_precision": float(average_precision_score(y, p)),
            "prevalence": float(y.mean()),
            "ap_lift_over_prevalence": float(average_precision_score(y, p) - y.mean()),
            "balanced_accuracy_0p5": float(balanced_accuracy_score(y, hard)),
            "mcc_0p5": float(matthews_corrcoef(y, hard))}


def train_predict(x_train, pockets, y_train, source_train, x_test, mode, seed, args):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    x, p, y = map(lambda a: torch.tensor(a, dtype=torch.float32), (x_train, pockets, y_train))
    xv = torch.tensor(x_test, dtype=torch.float32)
    s = torch.tensor(source_train, dtype=torch.long)
    model = PocketRealign(x.shape[1], args.rank, int(source_train.max()) + 1)
    if args.family_checkpoint is not None:
        family = torch.load(args.family_checkpoint, map_location="cpu")
        state = family["model"]
        with torch.no_grad():
            for target, prefix in [(model.molecule_adapter, "mol"), (model.pocket_adapter, "pocket")]:
                target.down.weight.copy_(state[f"{prefix}.down.weight"])
                target.up.weight.copy_(state[f"{prefix}.up.weight"])
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    pos_weight = torch.tensor(float((y_train == 0).sum() / max(1, (y_train == 1).sum())))
    original_x, original_p = F.normalize(x, dim=1), F.normalize(p, dim=1)
    use_rank = mode in {"pocket_rank", "pocket_rank_domain"}
    use_domain = mode == "pocket_rank_domain"
    for _ in range(args.epochs):
        model.train(); opt.zero_grad()
        zm, zp, _, retrieval, logit, source_logit = model(
            x, p, args.domain_weight if use_domain else 0.0)
        loss = F.binary_cross_entropy_with_logits(logit, y, pos_weight=pos_weight)
        if use_rank:
            positive = retrieval[y.bool()]
            negative = retrieval[~y.bool()]
            hard_negative = torch.topk(negative, min(args.hard_negative_top_k, len(negative))).values
            if len(positive) > args.max_positive:
                positive = positive[torch.randperm(len(positive))[:args.max_positive]]
            rank = F.relu(args.margin - positive[:, None] + hard_negative[None, :]).mean()
            loss = loss + args.rank_weight * rank
        preserve = (1 - (zm * original_x).sum(1)).mean() + (1 - (zp * original_p).sum(1)).mean()
        loss = loss + args.preserve_weight * preserve
        if use_domain:
            loss = loss + args.domain_loss_weight * F.cross_entropy(source_logit, s)
        loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        _, _, state_sim, _, logit, _ = model(xv, p)
    return torch.sigmoid(logit).numpy(), state_sim.numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings", type=Path, required=True)
    ap.add_argument("--benchmark", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--protocol", choices=["source", "series"], required=True)
    ap.add_argument("--seeds", default="20260924,20260925,20260926")
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--rank", type=int, default=16)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--weight-decay", type=float, default=1e-3)
    ap.add_argument("--margin", type=float, default=0.05)
    ap.add_argument("--rank-weight", type=float, default=0.5)
    ap.add_argument("--preserve-weight", type=float, default=0.2)
    ap.add_argument("--domain-weight", type=float, default=0.2)
    ap.add_argument("--domain-loss-weight", type=float, default=0.2)
    ap.add_argument("--hard-negative-top-k", type=int, default=16)
    ap.add_argument("--max-positive", type=int, default=128)
    ap.add_argument("--family-checkpoint", type=Path)
    args = ap.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    archive = np.load(args.embeddings, allow_pickle=False)
    ids = list(map(str, archive["molecule_ids"])); lookup = {v: i for i, v in enumerate(ids)}
    table = pd.read_csv(args.benchmark); canon = table.canonical_smiles.map(canonical)
    if missing := [v for v in canon if v not in lookup]: raise ValueError(f"Missing: {missing[:5]}")
    order = np.asarray([lookup[v] for v in canon])
    x = archive["molecule_embeddings"].astype(np.float32)[order]
    pockets = archive["pocket_embeddings"].astype(np.float32)
    y = table.target.to_numpy(np.int64)
    source_names = sorted(table.source_component.unique()); smap = {v: i for i, v in enumerate(source_names)}
    source = table.source_component.map(smap).to_numpy(int)
    folds = (table.source_fold if args.protocol == "source" else table.series_holdout_fold).to_numpy(int)
    assigned = folds >= 0; seeds = [int(v) for v in args.seeds.split(",")]
    modes = ["pocket_bce", "pocket_rank", "pocket_rank_domain"]
    records, seed_reports = [], []
    for seed in seeds:
        pred = {m: np.full(len(y), np.nan) for m in modes}
        for fold in sorted(np.unique(folds[assigned])):
            tr, te = assigned & (folds != fold), assigned & (folds == fold)
            train_sources = sorted(np.unique(source[tr])); dense = {v: i for i, v in enumerate(train_sources)}
            ds = np.asarray([dense[v] for v in source[tr]])
            for mode in modes:
                probability, _ = train_predict(x[tr], pockets, y[tr], ds, x[te], mode,
                                                seed + int(fold) * 100, args)
                pred[mode][te] = probability
        result = {m: metrics(y[assigned], pred[m][assigned]) for m in modes}
        seed_reports.append({"seed": seed, "metrics": result})
        for i in np.flatnonzero(assigned):
            records.append({"canonical_molecule_id": table.iloc[i].canonical_molecule_id,
                            "target": int(y[i]), "fold": int(folds[i]), "seed": seed,
                            **{m: float(pred[m][i]) for m in modes}})
    aggregate = {}
    for mode in modes:
        auc = [item["metrics"][mode]["roc_auc"] for item in seed_reports]
        aggregate[mode] = {"roc_auc_mean": float(np.mean(auc)), "roc_auc_min": float(np.min(auc)),
                           "roc_auc_max": float(np.max(auc)),
                           "per_seed": [item["metrics"][mode] for item in seed_reports]}
    report = {"evidence_level": "retrospective_m4_specific_drugclip_realign",
              "protocol": args.protocol, "n_assigned": int(assigned.sum()),
              "pocket_ids": list(map(str, archive["pocket_ids"])), "aggregate": aggregate,
              "hyperparameters": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()
                                  if k not in {"embeddings", "benchmark", "output"}},
              "claim_boundary": "Frozen-embedding target adaptation; no wet-lab PAM confirmation."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(records).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__": main()
