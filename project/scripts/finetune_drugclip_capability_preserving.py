#!/usr/bin/env python3
"""Capability-preserving DrugCLIP adaptation for GPCR and muscarinic pockets.

The adapter combines target-balanced four-GPCR retrieval with same-molecule
cross-muscarinic-subtype triplets.  Frozen-score/embedding preservation limits
catastrophic forgetting.  Every reported retrieval prediction is global
Murcko-scaffold OOF; auxiliary triplets sharing a held-out scaffold are removed.
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

from finetune_drugclip_gpcr_retrieval import (
    TARGETS, RetrievalLoRA, retrieval_metrics,
)
from finetune_drugclip_muscarinic_triplet import Proj


def scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    core = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    return core or Chem.MolToSmiles(mol, canonical=True)


def load_task(archive_path: Path, projection_path: Path, benchmark_path: Path,
              triplet_archive_path: Path, triplet_path: Path):
    arc = np.load(archive_path, allow_pickle=False)
    aux = np.load(triplet_archive_path, allow_pickle=False)
    initial = torch.load(projection_path, map_location="cpu")
    table = pd.read_csv(benchmark_path)
    index = {str(v): i for i, v in enumerate(arc["molecule_ids"])}
    table = table[table.canonical_smiles.astype(str).isin(index)].copy().reset_index(drop=True)
    order = [index[str(v)] for v in table.canonical_smiles]
    main_mol = torch.as_tensor(arc["molecule_representations"][order].astype(np.float32))
    main_pocket = torch.as_tensor(arc["pocket_representations"].astype(np.float32))
    pocket_ids = [str(v).split("_")[0] for v in arc["pocket_ids"]]
    if pocket_ids != TARGETS:
        raise ValueError(f"Unexpected GPCR pocket order: {pocket_ids}")

    aux_index = {str(v): i for i, v in enumerate(aux["molecule_ids"])}
    aux_pocket_ids = [str(v).split("_")[0] for v in aux["pocket_ids"]]
    aux_pocket_index = {v: i for i, v in enumerate(aux_pocket_ids)}
    triplets = pd.read_csv(triplet_path)
    keep = (triplets.canonical_smiles.astype(str).isin(aux_index)
            & triplets.positive_subtype.astype(str).isin(aux_pocket_index)
            & triplets.negative_subtype.astype(str).isin(aux_pocket_index))
    triplets = triplets[keep].copy().reset_index(drop=True)
    triplets["mol_index"] = [aux_index[str(v)] for v in triplets.canonical_smiles]
    triplets["positive_index"] = [aux_pocket_index[str(v)] for v in triplets.positive_subtype]
    triplets["negative_index"] = [aux_pocket_index[str(v)] for v in triplets.negative_subtype]
    triplets["murcko_scaffold"] = triplets.canonical_smiles.astype(str).map(scaffold)
    return (initial, table, main_mol, main_pocket,
            torch.as_tensor(aux["molecule_representations"].astype(np.float32)),
            torch.as_tensor(aux["pocket_representations"].astype(np.float32)), triplets)


def fit(main_mol, main_pocket, labels, aux_mol, aux_pocket, triplets,
        initial, seed, mode, args):
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    model = RetrievalLoRA(initial, args.rank, seed)
    frozen_mol = Proj(copy.deepcopy(initial["mol_project"])).eval()
    frozen_pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode():
        main_zm0 = frozen_mol(main_mol); main_zp0 = frozen_pocket(main_pocket)
        aux_zm0 = frozen_mol(aux_mol); aux_zp0 = frozen_pocket(aux_pocket)
        official_main_scores = main_zm0 @ main_zp0.T
        official_aux_scores = aux_zm0 @ aux_zp0.T
    main_zm0 = main_zm0.clone(); main_zp0 = main_zp0.clone()
    aux_zm0 = aux_zm0.clone(); aux_zp0 = aux_zp0.clone()
    official_main_scores = official_main_scores.clone()
    official_aux_scores = official_aux_scores.clone()
    n_classes = main_pocket.shape[0]
    counts = torch.bincount(labels, minlength=n_classes).float()
    class_weights = counts.sum() / (n_classes * counts.clamp_min(1))
    mi = torch.as_tensor(triplets.mol_index.to_numpy(), dtype=torch.long)
    pi = torch.as_tensor(triplets.positive_index.to_numpy(), dtype=torch.long)
    ni = torch.as_tensor(triplets.negative_index.to_numpy(), dtype=torch.long)
    if mode == "random_aux":
        swap = torch.as_tensor(rng.random(len(triplets)) < 0.5)
        pi, ni = torch.where(swap, ni, pi), torch.where(swap, pi, ni)
    used_aux = torch.unique(mi)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        optimizer.zero_grad()
        zm, zp, scores = model(main_mol, main_pocket)
        azm, azp, aux_scores = model(aux_mol, aux_pocket)
        ce = F.cross_entropy(scores / args.temperature, labels, weight=class_weights)
        triplet = F.relu(args.margin - aux_scores[mi, pi] + aux_scores[mi, ni]).mean()
        preserve = ((1 - (zm * main_zm0).sum(1)).mean()
                    + (1 - (zp * main_zp0).sum(1)).mean()
                    + (1 - (azm[used_aux] * aux_zm0[used_aux]).sum(1)).mean()
                    + (1 - (azp * aux_zp0).sum(1)).mean())
        distill = (F.mse_loss(scores, official_main_scores)
                   + F.mse_loss(aux_scores[used_aux], official_aux_scores[used_aux]))
        auxiliary_weight = 0.0 if mode == "ce_only" else args.triplet_weight
        loss = (ce + auxiliary_weight * triplet
                + args.preserve_weight * preserve + args.distill_weight * distill)
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
    return model


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-representations", type=Path, required=True)
    p.add_argument("--projection", type=Path, required=True)
    p.add_argument("--retrieval-benchmark", type=Path, required=True)
    p.add_argument("--triplet-representations", type=Path, required=True)
    p.add_argument("--triplets", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--zero-shot-test", type=Path)
    p.add_argument("--seed", type=int, default=20260925)
    p.add_argument("--rank", type=int, default=4)
    p.add_argument("--epochs", type=int, default=100)
    p.add_argument("--lr", type=float, default=5e-3)
    p.add_argument("--weight-decay", type=float, default=1e-3)
    p.add_argument("--temperature", type=float, default=0.07)
    p.add_argument("--margin", type=float, default=0.1)
    p.add_argument("--triplet-weight", type=float, default=0.5)
    p.add_argument("--preserve-weight", type=float, default=0.2)
    p.add_argument("--distill-weight", type=float, default=0.2)
    args = p.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    initial, table, main_mol, main_pocket, aux_mol, aux_pocket, triplets = load_task(
        args.retrieval_representations, args.projection, args.retrieval_benchmark,
        args.triplet_representations, args.triplets)
    labels = torch.as_tensor([TARGETS.index(v) for v in table.target], dtype=torch.long)
    folds = table.scaffold_fold.to_numpy(int)
    frozen_mol = Proj(copy.deepcopy(initial["mol_project"])).eval()
    frozen_pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode(): official = (frozen_mol(main_mol) @ frozen_pocket(main_pocket).T).numpy()
    modes = ["ce_only", "multitask", "random_aux"]
    predictions = {m: np.full((len(table), len(TARGETS)), np.nan) for m in modes}
    fold_audit = []
    for fold in sorted(np.unique(folds)):
        train = folds != fold; test = folds == fold
        held_scaffolds = set(table.loc[test, "murcko_scaffold"].astype(str))
        aux_train = triplets[~triplets.murcko_scaffold.astype(str).isin(held_scaffolds)].copy()
        fold_audit.append({"fold": int(fold), "main_train": int(train.sum()),
                           "main_test": int(test.sum()), "aux_triplets": int(len(aux_train)),
                           "aux_scaffold_overlap": 0})
        for mode in modes:
            model = fit(main_mol[train], main_pocket, labels[train], aux_mol, aux_pocket,
                        aux_train, initial, args.seed + 101 * int(fold), mode, args)
            with torch.inference_mode():
                predictions[mode][test] = model(main_mol[test], main_pocket)[2].numpy()
    report = {
        "method": "Capability-preserving DrugCLIP GPCR CE plus muscarinic same-molecule triplets",
        "split": "Five-fold global Murcko-scaffold OOF; overlapping auxiliary scaffolds removed per fold",
        "metrics": {"official": retrieval_metrics(labels.numpy(), official),
                    **{m: retrieval_metrics(labels.numpy(), v) for m, v in predictions.items()}},
        "fold_audit": fold_audit,
        "hyperparameters": {k: getattr(args, k) for k in ["rank", "epochs", "lr", "temperature",
            "margin", "triplet_weight", "preserve_weight", "distill_weight", "seed"]},
        "claim_boundary": "Retrospective GPCR target retrieval; no binding, PAM identity, or efficacy claim.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    rows = table[["canonical_smiles", "target", "murcko_scaffold", "scaffold_fold"]].copy()
    for name, values in {"official": official, **predictions}.items():
        for i, target in enumerate(TARGETS): rows[f"{name}_{target}"] = values[:, i]
    rows.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    if args.zero_shot_test:
        zero = pd.read_csv(args.zero_shot_test)
        zero_scaffolds = set(zero.canonical_smiles.astype(str).map(scaffold))
        safe_main = ((table.target != "M4R")
                     & ~table.murcko_scaffold.astype(str).isin(zero_scaffolds)).to_numpy()
        safe_aux = triplets[~triplets.murcko_scaffold.astype(str).isin(zero_scaffolds)].copy()
        zero_model = fit(main_mol[safe_main], main_pocket[:3], labels[safe_main],
                         aux_mol, aux_pocket, safe_aux, initial, args.seed,
                         "multitask", args)
        checkpoint = args.output.with_name(args.output.stem + ".zero_shot_safe.projection.pt")
        torch.save({"mol_project": zero_model.mol.materialized_state(),
                    "pocket_project": zero_model.pocket.materialized_state(),
                    "logit_scale": initial.get("logit_scale"),
                    "training_task": "B2AR/CCR2/M2R CE plus M1/M2/M3/M5 triplet",
                    "zero_shot_exclusions": {"M4R_target": True,
                                             "test_scaffolds": len(zero_scaffolds),
                                             "main_rows": int(safe_main.sum()),
                                             "aux_triplets": int(len(safe_aux))}}, checkpoint)
        report["zero_shot_safe_checkpoint"] = str(checkpoint)
        report["zero_shot_safe_training"] = {
            "main_rows": int(safe_main.sum()), "aux_triplets": int(len(safe_aux)),
            "excluded_test_scaffolds": int(len(zero_scaffolds))}
        args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
