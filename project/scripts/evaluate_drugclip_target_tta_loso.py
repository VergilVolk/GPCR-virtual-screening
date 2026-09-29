#!/usr/bin/env python3
"""Label-free target-level test-time adaptation of a DrugCLIP query.

For each target, the DrugCLIP encoders stay frozen.  A rank-r residual is
optimized from co-crystal reference ligands (positives), the pocket ensemble
(structural anchor), and unlabeled screening candidates (negatives).  Assay
labels are used only after adaptation for evaluation.

Hyperparameter variants are selected for every outer held target using only
the other targets' retrospective metrics (nested target LOO).
"""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import torch
from torch.nn import functional as F

from evaluate_drugclip_full_litpcba_reference_fusion import pocket_map
from evaluate_drugclip_official_litpcba_embeddings import metrics


VARIANTS = {
    "ref": {"pocket_weight": 0.0, "hard_fraction": 0.0, "anchor": 1e9},
    "ref_pocket25": {"pocket_weight": 0.25, "hard_fraction": 0.0, "anchor": 1e9},
    "tta_random": {"pocket_weight": 0.25, "hard_fraction": 0.0, "anchor": 1.0},
    "tta_mixed": {"pocket_weight": 0.25, "hard_fraction": 0.5, "anchor": 1.0},
    "tta_hard": {"pocket_weight": 0.25, "hard_fraction": 1.0, "anchor": 1.0},
    "tta_mixed_strong_anchor": {"pocket_weight": 0.25, "hard_fraction": 0.5, "anchor": 5.0},
}


def unit(value):
    return F.normalize(value, dim=-1)


def adapt_query(references, pockets, candidates, cfg, args, seed):
    torch.manual_seed(seed)
    references = unit(torch.as_tensor(references, dtype=torch.float32))
    pockets = unit(torch.as_tensor(pockets, dtype=torch.float32))
    candidates = unit(torch.as_tensor(candidates, dtype=torch.float32))
    pocket_center = unit(pockets.mean(dim=0, keepdim=True))[0]
    positive_sets = [references] if args.reference_mode == "centroid" else [value[None] for value in references]
    queries = []
    for query_index, positive in enumerate(positive_sets):
        torch.manual_seed(seed + query_index)
        ref_center = unit(positive.mean(dim=0, keepdim=True))[0]
        initial = unit(((1 - cfg["pocket_weight"]) * ref_center + cfg["pocket_weight"] * pocket_center)[None])[0]
        if cfg["anchor"] > 1e8:
            queries.append(initial.numpy())
            continue

        raw = candidates @ initial
        n_hard = int(args.n_negatives * cfg["hard_fraction"])
        n_random = args.n_negatives - n_hard
        chosen = []
        if n_hard:
            # Exclude the extreme top tail to reduce the chance of treating true
            # analogues/actives as negatives, then sample from the difficult band.
            order = torch.argsort(raw, descending=True)
            start = max(1, int(len(order) * args.hard_exclude_top))
            pool = order[start : min(len(order), start + max(n_hard * 20, n_hard))]
            chosen.append(pool[torch.randperm(len(pool))[:n_hard]])
        if n_random:
            chosen.append(torch.randperm(len(candidates))[:n_random])
        negative = candidates[torch.cat(chosen)]

        delta = torch.nn.Parameter(torch.zeros_like(initial))
        optimizer = torch.optim.AdamW([delta], lr=args.lr, weight_decay=0.0)
        for _ in range(args.steps):
            optimizer.zero_grad()
            query = unit((initial + delta)[None])[0]
            positive_score = positive @ query
            negative_score = negative @ query
            ranking = F.softplus(args.margin - positive_score[:, None] + negative_score[None, :]).mean()
            structural = 1 - (query * pocket_center).sum()
            preserve = 1 - (query * initial).sum()
            loss = ranking + args.pocket_anchor * structural + cfg["anchor"] * preserve
            loss.backward()
            optimizer.step()
        queries.append(unit((initial + delta).detach()[None])[0].numpy())
    return np.stack(queries)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pocket-root", type=Path, required=True)
    parser.add_argument("--pocket-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--n-negatives", type=int, default=2048)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--lr", type=float, default=0.03)
    parser.add_argument("--margin", type=float, default=0.15)
    parser.add_argument("--pocket-anchor", type=float, default=0.1)
    parser.add_argument("--hard-exclude-top", type=float, default=0.001)
    parser.add_argument("--reference-mode", choices=["centroid", "individual"], default="individual")
    parser.add_argument("--selection-metric", default="bedroc_alpha80_5")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    pockets = pocket_map(args.pocket_root, args.pocket_archive)
    per_target = {}
    score_dir = args.output.with_suffix("")
    score_dir.mkdir(parents=True, exist_ok=True)
    for target_index, target in enumerate(sorted(pockets)):
        folder = args.root / target / "drugclip_emb"
        molecules, _, labels = pickle.load((folder / "mols.lmdb.pkl").open("rb"))
        references, _, _ = pickle.load((folder / "ligand.lmdb.pkl").open("rb"))
        molecules = np.asarray(molecules, dtype=np.float32)
        references = np.asarray(references, dtype=np.float32)
        labels = np.asarray(labels, dtype=np.int64)
        result, scores = {}, {}
        for variant_index, (name, cfg) in enumerate(VARIANTS.items()):
            query = adapt_query(references, pockets[target], molecules, cfg, args,
                                args.seed + target_index * 100 + variant_index)
            score = (molecules @ query.T).max(axis=1)
            scores[name] = score
            result[name] = metrics(labels, score)
        per_target[target] = result
        np.savez_compressed(score_dir / f"{target}.npz", labels=labels, **scores)
        print(f"completed {target}", flush=True)

    targets = sorted(per_target)
    folds, selection_counts = {}, {name: 0 for name in VARIANTS}
    for held in targets:
        seen = [target for target in targets if target != held]
        validation = {name: float(np.mean([per_target[target][name][args.selection_metric] for target in seen]))
                      for name in VARIANTS}
        selected = max(VARIANTS, key=lambda name: validation[name])
        selection_counts[selected] += 1
        folds[held] = {"selected": selected, "seen_validation": validation,
                       "held_metrics": per_target[held][selected]}

    names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    fixed_macro = {variant: {name: float(np.mean([per_target[t][variant][name] for t in targets])) for name in names}
                   for variant in VARIANTS}
    nested_macro = {name: float(np.mean([folds[t]["held_metrics"][name] for t in targets])) for name in names}
    report = {
        "protocol": "Label-free per-target TTA with outer target-LOO hyperparameter selection",
        "n_targets": len(targets),
        "variants": VARIANTS,
        "fixed_macro": fixed_macro,
        "nested_macro": nested_macro,
        "selection_counts": selection_counts,
        "folds": folds,
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "claim_boundary": "Reference-ligand-assisted transductive screening; candidate labels are evaluation-only. Candidate embeddings are visible without labels during TTA. Not PAM efficacy prediction.",
    }
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"fixed_macro": fixed_macro, "nested_macro": nested_macro,
                      "selection_counts": selection_counts}, indent=2))


if __name__ == "__main__":
    main()
