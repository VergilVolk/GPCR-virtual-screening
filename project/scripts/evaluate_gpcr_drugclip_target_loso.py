#!/usr/bin/env python3
"""Leave-one-GPCR-target-out transfer audit for DrugCLIP projection tuning.

For each held target, all molecules from that target are excluded.  Training
uses only the other three targets and excludes the held pocket from the CE
denominator, so it is neither treated as a positive nor an artificial negative.
The adapted shared projections are then asked to retrieve the unseen held
pocket among all four pockets.  This is a target-transfer audit, not efficacy
or prospective-binding validation.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, roc_auc_score
from torch.nn import functional as F

from finetune_drugclip_gpcr_retrieval import TARGETS, RetrievalLoRA
from finetune_drugclip_muscarinic_triplet import Proj


def fit_seen_targets(molecule_rep, pocket_rep, labels, seen, initial, seed, args, random_labels=False):
    model = RetrievalLoRA(initial, args.rank, seed)
    frozen_mol = Proj(copy.deepcopy(initial["mol_project"])).eval()
    frozen_pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.no_grad():
        z0m = frozen_mol(molecule_rep).detach()
        z0p = frozen_pocket(pocket_rep).detach()
    remap = {original: local for local, original in enumerate(seen)}
    local_labels = torch.as_tensor([remap[int(value)] for value in labels], dtype=torch.long)
    if random_labels:
        generator = torch.Generator().manual_seed(seed)
        local_labels = local_labels[torch.randperm(len(local_labels), generator=generator)]
    counts = torch.bincount(local_labels, minlength=len(seen)).float()
    weights = counts.sum() / (len(seen) * counts.clamp_min(1))
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    seen_tensor = torch.as_tensor(seen, dtype=torch.long)
    for _ in range(args.epochs):
        optimizer.zero_grad()
        zm, zp, all_scores = model(molecule_rep, pocket_rep)
        scores = all_scores[:, seen_tensor]
        ce = F.cross_entropy(scores / args.temperature, local_labels, weight=weights)
        preserve = (1 - (zm * z0m).sum(1)).mean() + (1 - (zp * z0p).sum(1)).mean()
        loss = ce + args.preserve_weight * preserve
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
    return model


def metrics(labels, scores):
    order = np.argsort(-scores, axis=1)
    ranks = np.asarray([int(np.flatnonzero(order[i] == y)[0]) + 1 for i, y in enumerate(labels)])
    per_target = {}
    for j, target in enumerate(TARGETS):
        keep = labels == j
        per_target[target] = {
            "n": int(keep.sum()),
            "recall_at_1": float((ranks[keep] == 1).mean()),
            "recall_at_2": float((ranks[keep] <= 2).mean()),
            "mrr": float((1.0 / ranks[keep]).mean()),
        }
    one_hot = np.eye(len(TARGETS), dtype=int)[labels]
    return {
        "n": int(len(labels)),
        "recall_at_1": float((ranks == 1).mean()),
        "mrr": float((1.0 / ranks).mean()),
        "macro_recall_at_1": float(np.mean([v["recall_at_1"] for v in per_target.values()])),
        "macro_mrr": float(np.mean([v["mrr"] for v in per_target.values()])),
        "pair_roc_auc": float(roc_auc_score(one_hot.ravel(), scores.ravel())),
        "pair_pr_auc": float(average_precision_score(one_hot.ravel(), scores.ravel())),
        "per_target": per_target,
    }


def scaffold_bootstrap_delta(table, labels, tuned, official, repeats, seed):
    groups = table.murcko_scaffold.astype(str).unique()
    scaffold_values = table.murcko_scaffold.astype(str).to_numpy()
    group_indices = {group: np.flatnonzero(scaffold_values == group) for group in groups}
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(repeats):
        sampled = rng.choice(groups, len(groups), replace=True)
        idx = np.concatenate([group_indices[group] for group in sampled])
        if len(np.unique(labels[idx])) != len(TARGETS):
            continue
        values.append(metrics(labels[idx], tuned[idx])["macro_recall_at_1"] -
                      metrics(labels[idx], official[idx])["macro_recall_at_1"])
    return list(map(float, np.quantile(values, [0.025, 0.5, 0.975])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", default="20260925,20260926,20260927")
    parser.add_argument("--rank", type=int, default=4)
    parser.add_argument("--epochs", type=int, default=250)
    parser.add_argument("--lr", type=float, default=5e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--preserve-weight", type=float, default=0.2)
    parser.add_argument("--bootstrap", type=int, default=1000)
    args = parser.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    archive = np.load(args.representations, allow_pickle=False)
    initial = torch.load(args.projection, map_location="cpu")
    table_raw = pd.read_csv(args.benchmark)
    index = {value: i for i, value in enumerate(map(str, archive["molecule_ids"]))}
    table = table_raw[table_raw.canonical_smiles.astype(str).isin(index)].copy().reset_index(drop=True)
    order = np.asarray([index[value] for value in table.canonical_smiles.astype(str)])
    molecule_rep = torch.as_tensor(archive["molecule_representations"].astype(np.float32)[order])
    pocket_rep = torch.as_tensor(archive["pocket_representations"].astype(np.float32))
    labels = np.asarray([TARGETS.index(value) for value in table.target], dtype=int)
    frozen_mol = Proj(copy.deepcopy(initial["mol_project"])).eval()
    frozen_pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode():
        official = (frozen_mol(molecule_rep) @ frozen_pocket(pocket_rep).T).numpy()

    seed_values = list(map(int, args.seeds.split(",")))
    tuned_by_seed = [np.full_like(official, np.nan) for _ in seed_values]
    random_by_seed = [np.full_like(official, np.nan) for _ in seed_values]
    training_audit = {}
    scaffolds = table.murcko_scaffold.astype(str).to_numpy()
    for held in range(len(TARGETS)):
        test = labels == held
        held_scaffolds = set(scaffolds[test])
        train = (labels != held) & ~np.asarray([value in held_scaffolds for value in scaffolds])
        seen = [j for j in range(len(TARGETS)) if j != held]
        training_audit[TARGETS[held]] = {
            "test_n": int(test.sum()),
            "train_n_after_held_scaffold_purge": int(train.sum()),
            "purged_cross_target_rows": int(((labels != held) & ~train).sum()),
            "seen_targets": [TARGETS[j] for j in seen],
        }
        train_labels = torch.as_tensor(labels[train], dtype=torch.long)
        for seed_index, seed in enumerate(seed_values):
            model = fit_seen_targets(molecule_rep[train], pocket_rep, train_labels, seen, initial,
                                     seed + 1009 * held, args, random_labels=False)
            random_model = fit_seen_targets(molecule_rep[train], pocket_rep, train_labels, seen, initial,
                                            seed + 1009 * held, args, random_labels=True)
            with torch.inference_mode():
                tuned_by_seed[seed_index][test] = model(molecule_rep[test], pocket_rep)[2].numpy()
                random_by_seed[seed_index][test] = random_model(molecule_rep[test], pocket_rep)[2].numpy()

    tuned = np.mean(tuned_by_seed, axis=0)
    random = np.mean(random_by_seed, axis=0)
    report = {
        "protocol": "leave-one-entire-GPCR-target-out; held-target scaffolds purged from the other targets; held pocket absent from training loss",
        "training_audit": training_audit,
        "metrics": {
            "official_drugclip": metrics(labels, official),
            "three_target_balanced_ce_transfer": metrics(labels, tuned),
            "three_target_random_label_control": metrics(labels, random),
        },
        "tuned_minus_official_macro_recall1_scaffold_bootstrap_95ci": scaffold_bootstrap_delta(
            table, labels, tuned, official, args.bootstrap, 20260925),
        "claim_boundary": "Strict retrospective target-transfer audit across four GPCRs. No potency, efficacy, prospective-binding, or broad GPCR-generalization claim.",
        "hyperparameters": {"rank": args.rank, "epochs": args.epochs, "lr": args.lr,
                            "temperature": args.temperature, "preserve_weight": args.preserve_weight,
                            "seeds": seed_values},
    }
    rows = table[["canonical_molecule_id", "target", "murcko_scaffold"]].copy()
    for name, values in [("official", official), ("tuned", tuned), ("random", random)]:
        for j, target in enumerate(TARGETS):
            rows[f"{name}_{target}"] = values[:, j]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    rows.to_csv(args.output.with_suffix(".csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
