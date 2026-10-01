#!/usr/bin/env python3
"""Episodic target-domain generalization for DrugCLIP projection LoRA.

The outer protocol leaves one entire target out and purges its scaffolds.  In
each training epoch, one of the remaining targets is treated as a pseudo-unseen
query domain.  A differentiable virtual update is computed on all other target
domains, and the original adapter is optimized for performance on the query
target after that update (MLDG-style one-step meta objective).

The truly held target and its pocket never enter training or model selection.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.func import functional_call
from torch.nn import functional as F

from evaluate_drugclip_13target_ablation import Adapter, load_data, summarize
from finetune_drugclip_muscarinic_triplet import Proj


class EpisodicAdapter(Adapter):
    def forward(self, molecule_hidden, pocket_hidden):
        return self.mol.from_hidden(molecule_hidden) @ self.pocket.from_hidden(pocket_hidden).T


def target_label_weights(targets, labels):
    result = np.zeros(len(labels), dtype=np.float32)
    target_values = sorted(set(map(int, targets)))
    for target in target_values:
        for label in (0, 1):
            keep = (targets == target) & (labels == label)
            if keep.any():
                result[keep] = 1.0 / (len(target_values) * 2 * keep.sum())
    return result * len(labels) / result.sum()


def retrieval_indices(pm, pt, labels, allowed_targets):
    allowed = set(map(int, allowed_targets))
    remap = {target: i for i, target in enumerate(allowed_targets)}
    active = {}
    for molecule, target, label in zip(pm, pt, labels):
        if label and int(target) in allowed:
            active.setdefault(int(molecule), set()).add(int(target))
    pairs = [(molecule, remap[next(iter(targets))]) for molecule, targets in active.items()
             if len(targets) == 1]
    return (
        torch.as_tensor([value[0] for value in pairs], dtype=torch.long),
        torch.as_tensor([value[1] for value in pairs], dtype=torch.long),
        torch.as_tensor(list(allowed_targets), dtype=torch.long),
    )


def objective(score, pm, pt, labels, weights, retrieval, temperature, retrieval_weight):
    pair_score = score[pm, pt] / temperature
    bce = (F.binary_cross_entropy_with_logits(pair_score, labels, reduction="none") * weights).mean()
    rm, rt, allowed = retrieval
    if len(rm):
        ret = F.cross_entropy(score[rm][:, allowed] / temperature, rt)
    else:
        ret = torch.zeros((), dtype=bce.dtype)
    return bce + retrieval_weight * ret, {"bce": bce, "retrieval": ret}


def fit(mol, pocket, pm, pt, labels, seen, state, seed, args):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = EpisodicAdapter(state, args.rank, seed, True, True, "lora")
    with torch.no_grad():
        molecule_hidden = model.mol.hidden(mol)
        pocket_hidden = model.pocket.hidden(pocket)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    pm = np.asarray(pm, dtype=np.int64)
    pt = np.asarray(pt, dtype=np.int64)
    labels = np.asarray(labels, dtype=np.int64)
    history = []

    # Every block visits each pseudo-query target once in a random order.
    schedule = []
    while len(schedule) < args.epochs:
        schedule.extend(rng.permutation(seen).tolist())
    schedule = schedule[:args.epochs]

    for epoch, pseudo_query in enumerate(schedule):
        support_mask = pt != pseudo_query
        query_mask = pt == pseudo_query
        support_targets = [target for target in seen if target != pseudo_query]

        support_pm = torch.as_tensor(pm[support_mask], dtype=torch.long)
        support_pt = torch.as_tensor(pt[support_mask], dtype=torch.long)
        support_y = torch.as_tensor(labels[support_mask], dtype=torch.float32)
        support_w = torch.as_tensor(target_label_weights(pt[support_mask], labels[support_mask]))
        support_retrieval = retrieval_indices(
            pm[support_mask], pt[support_mask], labels[support_mask], support_targets)

        query_pm = torch.as_tensor(pm[query_mask], dtype=torch.long)
        query_pt = torch.as_tensor(pt[query_mask], dtype=torch.long)
        query_y = torch.as_tensor(labels[query_mask], dtype=torch.float32)
        query_w = torch.as_tensor(target_label_weights(pt[query_mask], labels[query_mask]))

        optimizer.zero_grad()
        score = model(molecule_hidden, pocket_hidden)
        support_loss, support_parts = objective(
            score, support_pm, support_pt, support_y, support_w, support_retrieval,
            args.temperature, args.retrieval_weight)

        named_parameters = dict(model.named_parameters())
        gradients = torch.autograd.grad(
            support_loss, tuple(named_parameters.values()),
            create_graph=not args.first_order, retain_graph=True)
        virtual_parameters = {
            name: parameter - args.inner_lr * gradient
            for (name, parameter), gradient in zip(named_parameters.items(), gradients)
        }
        query_score = functional_call(
            model, virtual_parameters, (molecule_hidden, pocket_hidden), strict=False)
        query_loss = (F.binary_cross_entropy_with_logits(
            query_score[query_pm, query_pt] / args.temperature,
            query_y,
            reduction="none") * query_w).mean()
        total_loss = support_loss + args.meta_weight * query_loss
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        history.append({
            "epoch": epoch,
            "pseudo_query": int(pseudo_query),
            "support_loss": float(support_loss.detach()),
            "support_bce": float(support_parts["bce"].detach()),
            "support_retrieval": float(support_parts["retrieval"].detach()),
            "query_after_virtual_update": float(query_loss.detach()),
        })
    return model, history


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpcr-representations", type=Path, required=True)
    parser.add_argument("--gpcr-pairs", type=Path, required=True)
    parser.add_argument("--external-representations", type=Path, required=True)
    parser.add_argument("--external-pairs", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--rank", type=int, default=8)
    parser.add_argument("--epochs", type=int, default=48)
    parser.add_argument("--lr", type=float, default=5e-3)
    parser.add_argument("--inner-lr", type=float, default=0.1)
    parser.add_argument("--meta-weight", type=float, default=0.5)
    parser.add_argument("--retrieval-weight", type=float, default=0.25)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--first-order", action="store_true")
    args = parser.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

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

    episodic = np.full(len(table), np.nan)
    audit = {}
    for held, target in enumerate(targets):
        test = pt == held
        held_scaffolds = set(scaffolds[test])
        train = (pt != held) & ~np.asarray([value in held_scaffolds for value in scaffolds])
        seen = [index for index in range(len(targets)) if index != held]
        model, history = fit(
            mol, pocket, pm[train], pt[train], labels[train], seen, state,
            args.seed + held * 1009, args)
        with torch.inference_mode():
            score = model(model.mol.hidden(mol), model.pocket.hidden(pocket)).numpy()
            episodic[test] = score[pm[test], pt[test]]
        audit[target] = {
            "train_pairs": int(train.sum()),
            "test_pairs": int(test.sum()),
            "scaffold_overlap": 0,
            "held_target_in_loss": False,
            "pseudo_query_targets": sorted(set(item["pseudo_query"] for item in history)),
            "last_epoch": history[-1],
            "trainable_parameters": int(sum(p.numel() for p in model.parameters() if p.requires_grad)),
        }
        print(f"completed held target {held + 1}/{len(targets)}: {target}", flush=True)

    report = {
        "method": "DrugCLIP dual projection LoRA with MLDG-style pseudo-unseen-target episodic objective",
        "protocol": "13-target leave-one-entire-target-out; global held-target scaffold purge; held pocket absent from training",
        "n_pairs": len(table),
        "n_targets": len(targets),
        "targets": targets,
        "metrics": {
            "official": summarize(table, official, targets),
            "episodic": summarize(table, episodic, targets),
        },
        "training_audit": audit,
        "hyperparameters": {key: getattr(args, key) for key in [
            "seed", "rank", "epochs", "lr", "inner_lr", "meta_weight",
            "retrieval_weight", "temperature", "first_order"]},
        "claim_boundary": "Retrospective target-domain generalization pilot; no PAM efficacy, prospective, or SOTA claim.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    predictions = table[["target", "label", "canonical_smiles", "murcko_scaffold"]].copy()
    predictions["official"] = official
    predictions["episodic"] = episodic
    predictions.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps({key: value["macro"] for key, value in report["metrics"].items()}, indent=2))


if __name__ == "__main__":
    main()
