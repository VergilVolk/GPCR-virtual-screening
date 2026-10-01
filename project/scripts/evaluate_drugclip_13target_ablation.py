#!/usr/bin/env python3
"""Strict 13-target DrugCLIP adapter ablations and target-robust training.

Every reported prediction is leave-one-entire-target-out.  Scaffolds observed
for the held target are removed from all training targets, and the held pocket
is excluded from every optimization term.  This script tests which adapter
side and which loss are actually responsible for transfer.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from torch.nn import functional as F

from finetune_drugclip_gpcr_screening import FastProjectionLoRA, binary_metrics
from finetune_drugclip_muscarinic_triplet import Proj


MODES = {
    "mol_only": {"mol": True, "pocket": False, "retrieval": 0.25, "preserve_mol": 0.2, "preserve_pocket": 0.2, "dro": False},
    "pocket_only": {"mol": False, "pocket": True, "retrieval": 0.25, "preserve_mol": 0.2, "preserve_pocket": 0.2, "dro": False},
    "dual_bce": {"mol": True, "pocket": True, "retrieval": 0.0, "preserve_mol": 0.2, "preserve_pocket": 0.2, "dro": False},
    "dual_no_preserve": {"mol": True, "pocket": True, "retrieval": 0.25, "preserve_mol": 0.0, "preserve_pocket": 0.0, "dro": False},
    "dual_mol_anchor": {"mol": True, "pocket": True, "retrieval": 0.25, "preserve_mol": 0.2, "preserve_pocket": 0.0, "dro": False},
    "dual_pocket_anchor": {"mol": True, "pocket": True, "retrieval": 0.25, "preserve_mol": 0.0, "preserve_pocket": 0.2, "dro": False},
    "dual_weak_anchor": {"mol": True, "pocket": True, "retrieval": 0.25, "preserve_mol": 0.05, "preserve_pocket": 0.05, "dro": False},
    "dual_standard": {"mol": True, "pocket": True, "retrieval": 0.25, "preserve_mol": 0.2, "preserve_pocket": 0.2, "dro": False},
    "dual_target_dro": {"mol": True, "pocket": True, "retrieval": 0.25, "preserve_mol": 0.2, "preserve_pocket": 0.2, "dro": True},
    "full_last_no_preserve": {"mol": True, "pocket": True, "retrieval": 0.25, "preserve_mol": 0.0, "preserve_pocket": 0.0, "dro": False, "adapter": "full_last", "lr_scale": 0.1},
}


class FullLastProjection(torch.nn.Module):
    """Train the complete final projection matrix; keep the first layer frozen."""
    def __init__(self, state):
        super().__init__()
        self.register_buffer("w1", state["linear1.weight"].clone())
        self.register_buffer("b1", state["linear1.bias"].clone())
        self.w2 = torch.nn.Parameter(state["linear2.weight"].clone())
        self.b2 = torch.nn.Parameter(state["linear2.bias"].clone())
        self.register_buffer("base_w2", state["linear2.weight"].clone())
        self.register_buffer("base_b2", state["linear2.bias"].clone())

    def hidden(self, representation):
        return F.relu(F.linear(representation, self.w1, self.b1))

    def from_hidden(self, hidden):
        return F.normalize(F.linear(hidden, self.w2, self.b2), dim=-1)

    def base_from_hidden(self, hidden):
        return F.normalize(F.linear(hidden, self.base_w2, self.base_b2), dim=-1)

    def materialized_state(self):
        return {"linear1.weight": self.w1.detach().cpu(), "linear1.bias": self.b1.detach().cpu(),
                "linear2.weight": self.w2.detach().cpu(), "linear2.bias": self.b2.detach().cpu()}


class Adapter(torch.nn.Module):
    def __init__(self, state, rank, seed, train_mol=True, train_pocket=True, adapter="lora"):
        super().__init__()
        if adapter == "full_last":
            self.mol = FullLastProjection(state["mol_project"])
            self.pocket = FullLastProjection(state["pocket_project"])
        else:
            self.mol = FastProjectionLoRA(state["mol_project"], rank, seed)
            self.pocket = FastProjectionLoRA(state["pocket_project"], rank, seed + 17)
        if not train_mol:
            for parameter in self.mol.parameters():
                parameter.requires_grad_(False)
        if not train_pocket:
            for parameter in self.pocket.parameters():
                parameter.requires_grad_(False)


def balanced_weights(targets, labels, n_targets):
    out = np.zeros(len(labels), dtype=np.float32)
    present = sorted(set(map(int, targets)))
    for target in present:
        for label in (0, 1):
            keep = (targets == target) & (labels == label)
            if keep.any():
                out[keep] = 1.0 / (len(present) * 2 * keep.sum())
    return out * len(labels) / out.sum()


def fit(mol, pocket, pm, pt, labels, seen, state, seed, args, mode):
    config = MODES[mode]
    model = Adapter(state, args.rank, seed, config["mol"], config["pocket"], config.get("adapter", "lora"))
    with torch.no_grad():
        hm = model.mol.hidden(mol)
        hp = model.pocket.hidden(pocket)
        z0m = model.mol.base_from_hidden(hm)
        z0p = model.pocket.base_from_hidden(hp)

    pmt = torch.as_tensor(np.asarray(pm, dtype=np.int64).copy())
    ptt = torch.as_tensor(np.asarray(pt, dtype=np.int64).copy())
    y = torch.as_tensor(np.asarray(labels, dtype=np.float32).copy())
    wt = torch.as_tensor(balanced_weights(pt, labels, len(pocket)))
    train_mol = torch.unique(pmt)
    seen_t = torch.as_tensor(seen, dtype=torch.long)

    remap = {target: i for i, target in enumerate(seen)}
    active = {}
    for molecule, target, label in zip(pm, pt, labels):
        if label:
            active.setdefault(int(molecule), set()).add(int(target))
    retrieval = [(m, remap[next(iter(ts))]) for m, ts in active.items()
                 if len(ts) == 1 and next(iter(ts)) in remap]
    rm = torch.as_tensor([v[0] for v in retrieval], dtype=torch.long)
    rt = torch.as_tensor([v[1] for v in retrieval], dtype=torch.long)

    target_rows = {target: torch.as_tensor(np.flatnonzero(pt == target), dtype=torch.long)
                   for target in seen}
    dro_q = torch.ones(len(seen), dtype=torch.float32) / len(seen)
    optimizer = torch.optim.AdamW(
        [p for p in model.parameters() if p.requires_grad],
        lr=args.lr * config.get("lr_scale", 1.0),
        weight_decay=args.weight_decay,
    )
    for _ in range(args.epochs):
        optimizer.zero_grad()
        zm = model.mol.from_hidden(hm)
        zp = model.pocket.from_hidden(hp)
        score = zm @ zp.T
        pair_loss = F.binary_cross_entropy_with_logits(
            score[pmt, ptt] / args.temperature, y, reduction="none")
        if config["dro"]:
            group_loss = torch.stack([(pair_loss[rows] * wt[rows]).mean() for rows in target_rows.values()])
            with torch.no_grad():
                dro_q *= torch.exp(args.dro_eta * group_loss.detach())
                dro_q /= dro_q.sum()
            bce = (dro_q * group_loss).sum()
        else:
            bce = (pair_loss * wt).mean()
        retrieval_loss = (F.cross_entropy(score[rm][:, seen_t] / args.temperature, rt)
                          if len(rm) else torch.zeros((), dtype=bce.dtype))
        preserve_mol = (1 - (zm[train_mol] * z0m[train_mol]).sum(1)).mean()
        preserve_pocket = (1 - (zp[seen_t] * z0p[seen_t]).sum(1)).mean()
        loss = (bce + config["retrieval"] * retrieval_loss
                + config["preserve_mol"] * preserve_mol
                + config["preserve_pocket"] * preserve_pocket)
        loss.backward()
        torch.nn.utils.clip_grad_norm_([p for p in model.parameters() if p.requires_grad], 1.0)
        optimizer.step()
    return model, {str(target): float(dro_q[i]) for i, target in enumerate(seen)} if config["dro"] else None


def summarize(table, scores, targets):
    per_target = {}
    for target in targets:
        keep = table.target.eq(target).to_numpy()
        per_target[target] = binary_metrics(table.loc[keep, "label"].to_numpy(int), scores[keep])
    names = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]
    return {
        "macro": {name: float(np.mean([per_target[target][name] for target in targets])) for name in names},
        "per_target": per_target,
    }


def load_data(args):
    ga = np.load(args.gpcr_representations, allow_pickle=False)
    ea = np.load(args.external_representations, allow_pickle=False)
    gp = pd.read_csv(args.gpcr_pairs)
    ep = pd.read_csv(args.external_pairs)
    gids = list(map(str, ga["molecule_ids"]))
    eids = list(map(str, ea["molecule_ids"]))
    gmi = {value: i for i, value in enumerate(gids)}
    emi = {value: i for i, value in enumerate(eids)}
    gp = gp[gp.canonical_smiles.astype(str).isin(gmi)].copy()
    ep = ep[ep.canonical_smiles.astype(str).isin(emi)].copy()
    offset = len(gids)
    gp["mol_index"] = [gmi[v] for v in gp.canonical_smiles.astype(str)]
    ep["mol_index"] = [offset + emi[v] for v in ep.canonical_smiles.astype(str)]
    gp_targets = ["B2AR", "CCR2", "M2R", "M4R"]
    gp_ids = list(map(str, ga["pocket_ids"]))
    gp_order = [gp_ids.index(f"{target}_cluster0") for target in gp_targets]
    external_targets = list(map(str, ea["pocket_ids"]))
    targets = gp_targets + external_targets
    target_index = {value: i for i, value in enumerate(targets)}
    gp["target_index"] = [target_index[v] for v in gp.target]
    ep["target_index"] = [target_index[v] for v in ep.target]
    table = pd.concat([gp, ep], ignore_index=True, sort=False)
    table["murcko_scaffold"] = [
        MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(value)) or value
        for value in table.canonical_smiles.astype(str)
    ]
    mol = torch.as_tensor(np.concatenate([
        ga["molecule_representations"], ea["molecule_representations"]]).astype(np.float32))
    pocket = torch.as_tensor(np.concatenate([
        ga["pocket_representations"][gp_order], ea["pocket_representations"]]).astype(np.float32))
    return table, mol, pocket, targets


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpcr-representations", type=Path, required=True)
    parser.add_argument("--gpcr-pairs", type=Path, required=True)
    parser.add_argument("--external-representations", type=Path, required=True)
    parser.add_argument("--external-pairs", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--modes", default=",".join(MODES))
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--rank", type=int, default=8)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--lr", type=float, default=5e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--dro-eta", type=float, default=0.05)
    parser.add_argument("--full-only", action="store_true")
    parser.add_argument("--save-full-checkpoint", type=Path)
    args = parser.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    modes = [value.strip() for value in args.modes.split(",") if value.strip()]
    unknown = sorted(set(modes) - set(MODES))
    if unknown:
        raise ValueError(f"Unknown modes: {unknown}")

    table, mol, pocket, targets = load_data(args)
    state = torch.load(args.projection, map_location="cpu")
    pm = table.mol_index.to_numpy(int)
    pt = table.target_index.to_numpy(int)
    labels = table.label.to_numpy(int)
    scaffolds = table.murcko_scaffold.astype(str).to_numpy()
    frozen_m = Proj(copy.deepcopy(state["mol_project"])).eval()
    frozen_p = Proj(copy.deepcopy(state["pocket_project"])).eval()
    with torch.inference_mode():
        official_matrix = frozen_m(mol) @ frozen_p(pocket).T
        official = official_matrix.numpy()[pm, pt]

    if args.full_only:
        if len(modes) != 1 or args.save_full_checkpoint is None:
            raise ValueError("--full-only requires exactly one --modes value and --save-full-checkpoint")
        model, dro_weights = fit(
            mol, pocket, pm, pt, labels, list(range(len(targets))), state,
            args.seed, args, modes[0])
        args.save_full_checkpoint.parent.mkdir(parents=True, exist_ok=True)
        torch.save({
            "mol_project": model.mol.materialized_state(),
            "pocket_project": model.pocket.materialized_state(),
            "base_projection": str(args.projection),
            "mode": modes[0],
            "training_targets": targets,
            "n_pairs": len(table),
            "seed": args.seed,
            "rank": args.rank,
            "epochs": args.epochs,
            "dro_weights": dro_weights,
            "claim_boundary": "Full-data deployment checkpoint; performance claims come only from LOSO outputs.",
        }, args.save_full_checkpoint)
        print(json.dumps({
            "checkpoint": str(args.save_full_checkpoint),
            "mode": modes[0],
            "n_pairs": len(table),
            "n_targets": len(targets),
            "seed": args.seed,
        }, indent=2))
        return

    predictions = {mode: np.full(len(table), np.nan) for mode in modes}
    audits = {}
    for held, target in enumerate(targets):
        test = pt == held
        held_scaffolds = set(scaffolds[test])
        train = (pt != held) & ~np.asarray([value in held_scaffolds for value in scaffolds])
        seen = [i for i in range(len(targets)) if i != held]
        audits[target] = {
            "train_pairs": int(train.sum()),
            "test_pairs": int(test.sum()),
            "scaffold_overlap": 0,
            "held_target_in_loss": False,
            "modes": {},
        }
        for mode in modes:
            model, dro_weights = fit(
                mol, pocket, pm[train], pt[train], labels[train], seen, state,
                args.seed + held * 1009, args, mode)
            with torch.inference_mode():
                matrix = model.mol.from_hidden(model.mol.hidden(mol)) @ model.pocket.from_hidden(model.pocket.hidden(pocket)).T
                predictions[mode][test] = matrix.numpy()[pm[test], pt[test]]
            audits[target]["modes"][mode] = {
                "trainable_parameters": int(sum(p.numel() for p in model.parameters() if p.requires_grad)),
                "final_dro_weights": dro_weights,
            }
        print(f"completed held target {held + 1}/{len(targets)}: {target}", flush=True)

    metrics = {"official": summarize(table, official, targets)}
    metrics.update({mode: summarize(table, values, targets) for mode, values in predictions.items()})
    report = {
        "protocol": "13-target leave-one-entire-target-out with global held-target scaffold purge",
        "n_pairs": len(table),
        "n_targets": len(targets),
        "targets": targets,
        "metrics": metrics,
        "training_audit": audits,
        "hyperparameters": {
            "seed": args.seed,
            "rank": args.rank,
            "epochs": args.epochs,
            "lr": args.lr,
            "temperature": args.temperature,
            "dro_eta": args.dro_eta,
        },
        "claim_boundary": "Retrospective adapter ablation. Target-DRO is exploratory until multi-seed confidence intervals and external confirmation.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    out = table[["target", "label", "canonical_smiles", "murcko_scaffold"]].copy()
    out["official"] = official
    for mode, values in predictions.items():
        out[mode] = values
    out.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps({name: value["macro"] for name, value in metrics.items()}, indent=2))


if __name__ == "__main__":
    main()
