#!/usr/bin/env python3
"""Low-rank DrugCLIP pocket-projection adaptation under target LOSO.

Only a rank-r delta on the final pocket projection is learned. The official
DrugCLIP molecule projection and both encoders remain frozen.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F

from finetune_drugclip_muscarinic_triplet import Proj, clustered_bootstrap


class PocketLoRA(nn.Module):
    def __init__(self, state: dict, rank: int, seed: int):
        super().__init__(); torch.manual_seed(seed)
        self.base = Proj(copy.deepcopy(state))
        for parameter in self.base.parameters(): parameter.requires_grad_(False)
        hidden = state["linear2.weight"].shape[1]; output = state["linear2.weight"].shape[0]
        self.down = nn.Linear(hidden, rank, bias=False); self.up = nn.Linear(rank, output, bias=False)
        nn.init.normal_(self.down.weight, std=0.02); nn.init.zeros_(self.up.weight)

    def forward(self, x):
        hidden = F.relu(self.base.linear1(x))
        return F.normalize(self.base.linear2(hidden) + self.up(self.down(hidden)), dim=-1)


def orient_indices(frame, mi, pi, seed, random_orientation):
    mol = torch.tensor([mi[s] for s in frame.canonical_smiles], dtype=torch.long)
    pos = torch.tensor([pi[s] for s in frame.positive_subtype], dtype=torch.long)
    neg = torch.tensor([pi[s] for s in frame.negative_subtype], dtype=torch.long)
    if random_orientation:
        generator = torch.Generator().manual_seed(seed); swap = torch.rand(len(frame), generator=generator) < 0.5
        original = pos.clone(); pos[swap] = neg[swap]; neg[swap] = original[swap]
    return mol, pos, neg


def train_lora(zm, pocket_rep, base_zp, train, mi, pi, state, seed, random_orientation, args):
    model = PocketLoRA(state, args.rank, seed)
    mol, pos, neg = orient_indices(train, mi, pi, seed, random_orientation)
    optimizer = torch.optim.AdamW([*model.down.parameters(), *model.up.parameters()],
                                  lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        optimizer.zero_grad(); zp = model(pocket_rep)
        delta = (zm[mol] * zp[pos]).sum(1) - (zm[mol] * zp[neg]).sum(1)
        rank_loss = F.relu(args.margin - delta).mean()
        preserve = (1 - (zp * base_zp).sum(1)).mean()
        loss = rank_loss + args.preserve_weight * preserve
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
    return model


def pair_delta(zm, zp, frame, mi, pi):
    scores = (zm @ zp.T).detach().numpy()
    return np.asarray([scores[mi[r.canonical_smiles], pi[r.positive_subtype]] -
                       scores[mi[r.canonical_smiles], pi[r.negative_subtype]]
                       for r in frame.itertuples()])


def direction(frame, score, held):
    preferred = frame.positive_subtype.to_numpy() == held; disfavored = ~preferred
    a = float((score[preferred] > 0).mean()); b = float((score[disfavored] > 0).mean())
    return {"accuracy": float((score > 0).mean()), "held_preferred_accuracy": a,
            "held_disfavored_accuracy": b, "direction_balanced_accuracy": (a + b) / 2}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", default="20260924,20260925,20260926")
    parser.add_argument("--rank", type=int, default=4)
    parser.add_argument("--epochs", type=int, default=500)
    parser.add_argument("--lr", type=float, default=1e-2)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--margin", type=float, default=0.1)
    parser.add_argument("--preserve-weight", type=float, default=0.2)
    args = parser.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    arc = np.load(args.representations, allow_pickle=False); initial = torch.load(args.projection, map_location="cpu")
    mids = list(map(str, arc["molecule_ids"])); pids = list(map(str, arc["pocket_ids"]))
    mi = {value: i for i, value in enumerate(mids)}; pi = {value.split("_")[0]: i for i, value in enumerate(pids)}
    molecule_rep = torch.as_tensor(arc["molecule_representations"], dtype=torch.float32)
    pocket_rep = torch.as_tensor(arc["pocket_representations"], dtype=torch.float32)
    mol_projection = Proj(copy.deepcopy(initial["mol_project"])).eval()
    pocket_projection = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode():
        zm = mol_projection(molecule_rep).clone()
        base_zp = pocket_projection(pocket_rep).clone()
    # Convert inference tensors to ordinary detached tensors before they
    # participate in an autograd graph through the LoRA branch.
    zm = zm.detach().clone(); base_zp = base_zp.detach().clone()
    reports, rows = {}, []
    for held in ["M1", "M2", "M3", "M4", "M5"]:
        train_raw = pd.read_csv(args.benchmark_dir / f"{held}_train.csv")
        test_raw = pd.read_csv(args.benchmark_dir / f"{held}_test.csv")
        train = train_raw[train_raw.canonical_smiles.isin(mi)].reset_index(drop=True)
        test = test_raw[test_raw.canonical_smiles.isin(mi)].reset_index(drop=True)
        base = pair_delta(zm, base_zp, test, mi, pi); targeted = []; controls = []
        for seed in map(int, args.seeds.split(",")):
            target_model = train_lora(zm, pocket_rep, base_zp, train, mi, pi,
                                      initial["pocket_project"], seed, False, args)
            random_model = train_lora(zm, pocket_rep, base_zp, train, mi, pi,
                                      initial["pocket_project"], seed, True, args)
            with torch.inference_mode():
                target = pair_delta(zm, target_model(pocket_rep), test, mi, pi)
                control = pair_delta(zm, random_model(pocket_rep), test, mi, pi)
            targeted.append(target); controls.append(control)
            for i, row in test.iterrows():
                rows.append({"held_subtype": held, "canonical_smiles": row.canonical_smiles,
                             "positive_subtype": row.positive_subtype, "negative_subtype": row.negative_subtype,
                             "seed": seed, "frozen_delta": float(base[i]),
                             "random_control_delta": float(control[i]), "targeted_delta": float(target[i])})
        target_mean = np.mean(targeted, axis=0); control_mean = np.mean(controls, axis=0)
        reports[held] = {
            "train_pairs": int(len(train)), "test_pairs": int(len(test)),
            "frozen": direction(test, base, held), "random_control": direction(test, control_mean, held),
            "targeted": direction(test, target_mean, held),
            "targeted_minus_frozen_molecule_cluster_95ci": clustered_bootstrap(test, base, target_mean),
            "targeted_minus_random_molecule_cluster_95ci": clustered_bootstrap(test, control_mean, target_mean),
            "per_seed_targeted_accuracy": [float((value > 0).mean()) for value in targeted],
        }
    macro = {method: {
        "accuracy": float(np.mean([reports[h][method]["accuracy"] for h in reports])),
        "direction_balanced_accuracy": float(np.mean([reports[h][method]["direction_balanced_accuracy"] for h in reports])),
    } for method in ["frozen", "random_control", "targeted"]}
    trainable = args.rank * (initial["pocket_project"]["linear2.weight"].shape[0] +
                             initial["pocket_project"]["linear2.weight"].shape[1])
    report = {"method": "DrugCLIP internal pocket projection rank-4 LoRA triplet",
              "trainable_parameters": int(trainable), "folds": reports, "macro_average": macro,
              "hyperparameters": vars(args) | {"representations": str(args.representations),
                                                "projection": str(args.projection),
                                                "benchmark_dir": str(args.benchmark_dir),
                                                "output": str(args.output)},
              "claim_boundary": "Posthoc target-LOSO development benchmark; not PAM efficacy."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(rows).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
