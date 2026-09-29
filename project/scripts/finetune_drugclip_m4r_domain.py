#!/usr/bin/env python3
"""M4R-domain adaptation of DrugCLIP internal projection heads.

Uses broad M4R allosteric active/decoy membership as an intermediate task.
Evaluation is strictly Murcko-scaffold-held-out and reports early enrichment.
The downstream functional PAM set was excluded when the benchmark was built.
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
from sklearn.metrics import average_precision_score, roc_auc_score
from torch import nn
from torch.nn import functional as F


class Projection(nn.Module):
    def __init__(self, state):
        super().__init__()
        self.linear1 = nn.Linear(state["linear1.weight"].shape[1], state["linear1.weight"].shape[0])
        self.linear2 = nn.Linear(state["linear2.weight"].shape[1], state["linear2.weight"].shape[0])
        self.load_state_dict(state)
    def forward(self, x): return F.normalize(self.linear2(F.relu(self.linear1(x))), dim=-1)


class Model(nn.Module):
    def __init__(self, initial):
        super().__init__()
        self.mol_project = Projection(initial["mol_project"])
        self.pocket_project = Projection(initial["pocket_project"])
        self.log_scale = nn.Parameter(initial["logit_scale"].float().reshape(()).clone())
        self.bias = nn.Parameter(torch.tensor(0.0))
    def forward(self, x, p):
        zm, zp = self.mol_project(x), self.pocket_project(p)
        state = zm @ zp.T
        retrieval = 0.1 * torch.logsumexp(state / 0.1, dim=1)
        return zm, zp, retrieval, self.log_scale.exp().clamp(max=100) * retrieval + self.bias


def ef(y, score, fraction=0.01):
    n = max(1, int(np.ceil(len(y) * fraction)))
    top = y[np.argsort(-score)[:n]].mean(); base = y.mean()
    return float(top / base) if base else float("nan")


def metrics(y, score):
    return {"roc_auc": float(roc_auc_score(y, score)),
            "average_precision": float(average_precision_score(y, score)),
            "ef1pct": ef(y, score), "prevalence": float(y.mean()), "n": int(len(y))}


def train(x, p, y, initial, mode, seed, args):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    x = torch.as_tensor(x, dtype=torch.float32); p = torch.as_tensor(p, dtype=torch.float32)
    y = np.asarray(y, int); pos = np.flatnonzero(y == 1); neg = np.flatnonzero(y == 0)
    model = Model(copy.deepcopy(initial)); frozen = Model(copy.deepcopy(initial)).eval()
    for q in frozen.parameters(): q.requires_grad_(False)
    with torch.no_grad(): z0m, z0p, r0, _ = frozen(x, p)
    hard_pool = neg[np.argsort(-r0.numpy()[neg])[:min(args.hard_pool, len(neg))]]
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    rng = np.random.default_rng(seed)
    half = args.batch_size // 2
    for _ in range(args.epochs):
        for _ in range(args.steps_per_epoch):
            pi = rng.choice(pos, half, replace=len(pos) < half)
            pool = hard_pool if mode == "hard" else neg
            ni = rng.choice(pool, half, replace=len(pool) < half)
            idx = np.concatenate([pi, ni]); rng.shuffle(idx)
            idx_t = torch.as_tensor(idx, dtype=torch.long)
            target = torch.as_tensor(y[idx], dtype=torch.float32)
            opt.zero_grad(); zm, zp, retrieval, logits = model(x[idx_t], p)
            loss = F.binary_cross_entropy_with_logits(logits, target)
            if mode != "bce":
                rp, rn = retrieval[target > .5], retrieval[target < .5]
                loss = loss + args.rank_weight * F.relu(args.margin - rp[:, None] + rn[None, :]).mean()
            preserve = (1 - (zm * z0m[idx_t]).sum(1)).mean() + (1 - (zp * z0p).sum(1)).mean()
            loss = loss + args.preserve_weight * preserve
            loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); opt.step()
    return model


def score(model, x, p):
    model.eval()
    with torch.no_grad(): return model(torch.as_tensor(x), torch.as_tensor(p))[2].numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--representations", type=Path, required=True)
    ap.add_argument("--projection", type=Path, required=True)
    ap.add_argument("--benchmark", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=20260924)
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--steps-per-epoch", type=int, default=20)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--weight-decay", type=float, default=1e-3)
    ap.add_argument("--margin", type=float, default=0.05)
    ap.add_argument("--rank-weight", type=float, default=0.5)
    ap.add_argument("--preserve-weight", type=float, default=0.5)
    ap.add_argument("--hard-pool", type=int, default=2048)
    args = ap.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    arc = np.load(args.representations, allow_pickle=False); initial = torch.load(args.projection, map_location="cpu")
    table = pd.read_csv(args.benchmark); lookup = {str(v): i for i, v in enumerate(arc["molecule_ids"])}
    keep = table.canonical_smiles.astype(str).isin(lookup)
    omitted = table.loc[~keep, "canonical_molecule_id"].astype(str).tolist(); table = table[keep].reset_index(drop=True)
    order = np.asarray([lookup[str(v)] for v in table.canonical_smiles], dtype=np.int64)
    x = arc["molecule_representations"].astype(np.float32)[order]
    p = arc["pocket_representations"].astype(np.float32); y = table.target.to_numpy(int)
    folds = table.scaffold_fold.to_numpy(int); modes = ["frozen", "bce", "random", "hard"]
    predictions = {m: np.full(len(y), np.nan) for m in modes}
    frozen = Model(copy.deepcopy(initial)).eval()
    for fold in sorted(np.unique(folds)):
        tr, te = folds != fold, folds == fold
        predictions["frozen"][te] = score(frozen, x[te], p)
        for mode in modes[1:]:
            fitted = train(x[tr], p, y[tr], initial, mode, args.seed + 101*int(fold), args)
            predictions[mode][te] = score(fitted, x[te], p)
    report = {"method": "DrugCLIP internal projection M4R-domain adaptation",
              "split": "5-fold Murcko scaffold OOF", "n": len(table), "omitted_input_ids": omitted,
              "metrics": {m: metrics(y, predictions[m]) for m in modes},
              "claim_boundary": "Broad allosteric active/decoy retrieval; not functional PAM efficacy."}
    # Deployable checkpoint is trained after OOF evaluation with frozen settings.
    final = train(x, p, y, initial, "hard", args.seed, args)
    checkpoint = {"mol_project": final.mol_project.state_dict(),
                  "pocket_project": final.pocket_project.state_dict(),
                  "logit_scale": final.log_scale.detach().cpu().reshape(1),
                  "classification_bias": final.bias.detach().cpu(),
                  "source_checkpoint": initial.get("source_checkpoint"),
                  "metadata": {"stage": "M4R broad allosteric domain", "n": len(table)}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    ckpt_path = args.output.with_suffix(".projection.pt"); torch.save(checkpoint, ckpt_path)
    report["checkpoint"] = str(ckpt_path)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame({"canonical_molecule_id": table.canonical_molecule_id, "target": y,
                  "fold": folds, **predictions}).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__": main()
