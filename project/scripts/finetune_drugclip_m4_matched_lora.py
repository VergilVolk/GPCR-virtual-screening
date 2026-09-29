#!/usr/bin/env python3
"""Parameter-efficient M4 functional triplet adaptation of DrugCLIP.

The GPCR-specialized pocket space is frozen. A rank-r LoRA delta is attached
to the final DrugCLIP molecule projection and trained on source/chemotype-
matched PAM-vs-inactive pairs.
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
from torch import nn
from torch.nn import functional as F

from finetune_drugclip_m4_matched_triplet import Proj, bootstrap, make_pairs, metric


class MolLoRA(nn.Module):
    def __init__(self, state, rank, seed):
        super().__init__(); torch.manual_seed(seed)
        self.base = Proj(copy.deepcopy(state))
        for parameter in self.base.parameters(): parameter.requires_grad_(False)
        hidden = state["linear2.weight"].shape[1]; output = state["linear2.weight"].shape[0]
        self.down = nn.Linear(hidden, rank, bias=False); self.up = nn.Linear(rank, output, bias=False)
        nn.init.normal_(self.down.weight, std=0.02); nn.init.zeros_(self.up.weight)

    def forward(self, x):
        hidden = F.relu(self.base.linear1(x))
        return F.normalize(self.base.linear2(hidden) + self.up(self.down(hidden)), dim=-1)

    def materialized_state(self):
        state = copy.deepcopy(self.base.state_dict())
        state["linear2.weight"] = state["linear2.weight"] + self.up.weight.detach().cpu() @ self.down.weight.detach().cpu()
        return state


def retrieval(model, x, zp):
    zm = model(x); states = zm @ zp.T
    return zm, 0.1 * torch.logsumexp(states / 0.1, dim=1)


def fit(x, zp, frame, initial, mode, seed, args):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    model = MolLoRA(initial["mol_project"], args.rank, seed)
    with torch.no_grad(): z0 = model.base(x).detach().clone()
    pairs = make_pairs(frame, mode, args.positive_k, seed); pair = torch.as_tensor(pairs, dtype=torch.long)
    optimizer = torch.optim.AdamW([*model.down.parameters(), *model.up.parameters()],
                                  lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        optimizer.zero_grad(); zm, score = retrieval(model, x, zp)
        rank_loss = F.relu(args.margin - (score[pair[:, 0]] - score[pair[:, 1]])).mean()
        preserve = (1 - (zm * z0).sum(1)).mean()
        loss = rank_loss + args.preserve_weight * preserve
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
    return model, pairs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol", choices=["series", "source"], default="series")
    parser.add_argument("--seeds", default="20260924,20260925,20260926")
    parser.add_argument("--rank", type=int, default=4)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--lr", type=float, default=1e-2)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--margin", type=float, default=0.05)
    parser.add_argument("--preserve-weight", type=float, default=0.5)
    parser.add_argument("--positive-k", type=int, default=3)
    parser.add_argument("--bootstrap", type=int, default=3000)
    args = parser.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    arc = np.load(args.representations, allow_pickle=False); initial = torch.load(args.projection, map_location="cpu")
    table = pd.read_csv(args.benchmark); lookup = {str(v): i for i, v in enumerate(arc["molecule_ids"])}
    order = np.asarray([lookup[str(v)] for v in table.canonical_smiles])
    x = torch.as_tensor(arc["molecule_representations"].astype(np.float32)[order])
    pocket_rep = torch.as_tensor(arc["pocket_representations"].astype(np.float32))
    pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode(): zp = pocket(pocket_rep)
    zp = zp.detach().clone()
    fold_column = "series_holdout_fold" if args.protocol == "series" else "source_fold"
    folds = table[fold_column].to_numpy(int); assigned = folds >= 0; y = table.target.to_numpy(int)
    base_mol = Proj(copy.deepcopy(initial["mol_project"])).eval()
    with torch.inference_mode(): _, base_score = retrieval(base_mol, x, zp)
    rows = []
    for seed in map(int, args.seeds.split(",")):
        predictions = {mode: np.full(len(table), np.nan) for mode in ["frozen", "random", "matched"]}
        predictions["frozen"][assigned] = base_score.numpy()[assigned]
        for fold in sorted(np.unique(folds[assigned])):
            train_mask = assigned & (folds != fold); test_mask = assigned & (folds == fold)
            frame = table.loc[train_mask].reset_index(drop=True); train_x = x[train_mask]
            for mode in ["random", "matched"]:
                model, _ = fit(train_x, zp, frame, initial, mode, seed + 101 * int(fold), args)
                with torch.no_grad(): _, score = retrieval(model, x[test_mask], zp)
                predictions[mode][test_mask] = score.numpy()
        for i in np.flatnonzero(assigned):
            rows.append({"canonical_molecule_id": table.iloc[i].canonical_molecule_id,
                         "target": int(y[i]), "fold": int(folds[i]), "seed": seed,
                         **{mode: float(predictions[mode][i]) for mode in predictions}})
    frame = pd.DataFrame(rows); modes = ["frozen", "random", "matched"]
    ensemble = frame.groupby("canonical_molecule_id", sort=False)[modes].mean()
    truth = table.set_index("canonical_molecule_id").loc[ensemble.index, "target"].to_numpy(int)
    args.output.parent.mkdir(parents=True, exist_ok=True); checkpoints = []
    full_table = table.loc[assigned].reset_index(drop=True); full_x = x[assigned]
    for seed in map(int, args.seeds.split(",")):
        model, pairs = fit(full_x, zp, full_table, initial, "matched", seed, args)
        path = args.output.with_name(f"{args.output.stem}.seed{seed}.projection.pt")
        torch.save({"mol_project": model.materialized_state(),
                    "pocket_project": copy.deepcopy(initial["pocket_project"]),
                    "logit_scale": initial.get("logit_scale"), "seed": seed, "rank": args.rank,
                    "training_pairs": int(len(pairs)), "source_projection": str(args.projection)}, path)
        checkpoints.append(str(path))
    report = {
        "method": "GPCR-initialized DrugCLIP M4 matched functional rank-4 molecule LoRA",
        "protocol": args.protocol + " holdout", "n": int(len(truth)), "trainable_parameters": int(args.rank * 640),
        "metrics": {mode: metric(truth, ensemble[mode].to_numpy()) for mode in modes},
        "matched_minus_frozen_auc_bootstrap_95ci": bootstrap(
            truth, ensemble.matched.to_numpy(), ensemble.frozen.to_numpy(), args.bootstrap),
        "matched_minus_random_auc_bootstrap_95ci": bootstrap(
            truth, ensemble.matched.to_numpy(), ensemble.random.to_numpy(), args.bootstrap),
        "deployable_checkpoints": checkpoints,
        "claim_boundary": "Retrospective PAM-vs-inactive ranking; external labels not used in this fit.",
    }
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    frame.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
