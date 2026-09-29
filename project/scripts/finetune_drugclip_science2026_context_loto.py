#!/usr/bin/env python3
"""Context-gated target-LOTO adaptation of frozen Science-2026 DrugCLIP.

This model changes only the similarity relation.  A pocket-conditioned gate
mixes low-rank bilinear experts on top of the frozen cosine score.  Evaluation
holds out an entire target and purges its Murcko scaffolds from all training
rows.  No held-target labels, pocket vectors, or hyperparameter selection are
used during fitting.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from torch.nn import functional as F

from finetune_drugclip_gpcr_screening import binary_metrics


def scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or Chem.MolToSmiles(mol, canonical=True)


class ContextMetric(torch.nn.Module):
    def __init__(self, dim: int, rank: int, experts: int, seed: int):
        super().__init__()
        torch.manual_seed(seed)
        self.mol = torch.nn.Parameter(torch.empty(experts, dim, rank))
        self.pocket = torch.nn.Parameter(torch.empty(experts, dim, rank))
        self.gate = torch.nn.Linear(dim, experts)
        self.scale = torch.nn.Parameter(torch.tensor(-2.0))
        torch.nn.init.normal_(self.mol, std=0.02)
        torch.nn.init.normal_(self.pocket, std=0.02)
        torch.nn.init.zeros_(self.gate.weight)
        torch.nn.init.zeros_(self.gate.bias)
        self.rank = rank

    def score_pairs(self, molecules: torch.Tensor, pockets: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        base = (molecules * pockets).sum(dim=-1)
        m = torch.einsum("nd,edr->ner", molecules, self.mol)
        p = torch.einsum("nd,edr->ner", pockets, self.pocket)
        expert = (m * p).sum(dim=-1) / math.sqrt(self.rank)
        gate = torch.softmax(self.gate(pockets), dim=-1)
        residual = torch.tanh(self.scale) * (gate * expert).sum(dim=-1)
        return base + residual, residual


def balanced_weights(target: np.ndarray, label: np.ndarray, targets: np.ndarray) -> np.ndarray:
    out = np.zeros(len(label), dtype=np.float32)
    for t in targets:
        for y in (0, 1):
            keep = (target == t) & (label == y)
            if keep.any():
                out[keep] = 1.0 / (2 * len(targets) * keep.sum())
    out *= len(out) / out.sum()
    return out


def train(
    molecule_embedding: torch.Tensor,
    pocket_embedding: torch.Tensor,
    pair_mol: np.ndarray,
    pair_target: np.ndarray,
    labels: np.ndarray,
    train_rows: np.ndarray,
    seen_targets: np.ndarray,
    args: argparse.Namespace,
    seed: int,
    shuffle: bool,
) -> ContextMetric:
    rng = np.random.default_rng(seed)
    fit_label = labels[train_rows].copy()
    fit_target = pair_target[train_rows]
    if shuffle:
        for target in seen_targets:
            local = np.flatnonzero(fit_target == target)
            fit_label[local] = rng.permutation(fit_label[local])

    mol = molecule_embedding[torch.as_tensor(pair_mol[train_rows])]
    pocket = pocket_embedding[torch.as_tensor(fit_target)]
    y = torch.as_tensor(fit_label, dtype=torch.float32)
    weight = torch.as_tensor(balanced_weights(fit_target, fit_label, seen_targets))
    model = ContextMetric(mol.shape[1], args.rank, args.experts, seed)

    with torch.no_grad():
        base = (mol * pocket).sum(dim=-1).cpu().numpy()
    rank_groups = []
    for target in seen_targets:
        local = np.flatnonzero(fit_target == target)
        positive = local[fit_label[local] == 1]
        negative = local[fit_label[local] == 0]
        if len(positive) and len(negative):
            hard = negative[np.argsort(base[negative])[::-1][: args.hard_negatives]]
            rank_groups.append((torch.as_tensor(positive), torch.as_tensor(hard)))

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        optimizer.zero_grad()
        score, residual = model.score_pairs(mol, pocket)
        bce = (
            F.binary_cross_entropy_with_logits(score / args.temperature, y, reduction="none") * weight
        ).mean()
        ranking = torch.stack([
            F.softplus(args.margin - score[pos, None] + score[neg][None, :]).mean()
            for pos, neg in rank_groups
        ]).mean()
        # The residual remains a bounded correction to the frozen model.
        preservation = residual.square().mean()
        loss = bce + args.ranking_weight * ranking + args.preserve_weight * preservation
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
    return model.eval()


def metrics(frame: pd.DataFrame, score: np.ndarray, targets: list[str]) -> dict:
    per = {}
    for target in targets:
        keep = frame.target.eq(target).to_numpy()
        per[target] = binary_metrics(frame.loc[keep, "label"].to_numpy(int), score[keep])
    names = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]
    return {
        "macro": {name: float(np.mean([per[t][name] for t in targets])) for name in names},
        "per_target": per,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", default="20260925,20260926,20260927")
    parser.add_argument("--rank", type=int, default=8)
    parser.add_argument("--experts", type=int, default=3)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--lr", type=float, default=3e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--margin", type=float, default=0.10)
    parser.add_argument("--hard-negatives", type=int, default=256)
    parser.add_argument("--ranking-weight", type=float, default=0.50)
    parser.add_argument("--preserve-weight", type=float, default=0.25)
    parser.add_argument("--bootstrap", type=int, default=10000)
    parser.add_argument(
        "--deployable-dir",
        type=Path,
        default=None,
        help="Optionally fit one final adapter per seed on all development targets and save checkpoints.",
    )
    args = parser.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    archive = np.load(args.representations, allow_pickle=False)
    molecule_ids = list(map(str, archive["molecule_ids"]))
    pocket_ids = list(map(str, archive["pocket_ids"]))
    frame = pd.read_csv(args.pairs)
    targets = sorted(set(frame.target.astype(str)) & set(pocket_ids))
    mol_index = {value: i for i, value in enumerate(molecule_ids)}
    pocket_index = {value: pocket_ids.index(value) for value in targets}
    frame = frame[
        frame.target.astype(str).isin(targets)
        & frame.canonical_smiles.astype(str).isin(mol_index)
    ].copy().reset_index(drop=True)
    frame["scaffold"] = [scaffold(value) for value in frame.canonical_smiles.astype(str)]
    pair_mol = np.asarray([mol_index[value] for value in frame.canonical_smiles.astype(str)])
    pair_target = np.asarray([pocket_index[value] for value in frame.target.astype(str)])
    labels = frame.label.to_numpy(int)
    molecule_embedding = torch.as_tensor(archive["molecule_embeddings"].astype(np.float32))
    pocket_embedding = torch.as_tensor(archive["pocket_embeddings"].astype(np.float32))
    raw = (
        molecule_embedding[torch.as_tensor(pair_mol)]
        * pocket_embedding[torch.as_tensor(pair_target)]
    ).sum(dim=-1).numpy()

    seeds = [int(value) for value in args.seeds.split(",")]
    true_predictions = {seed: np.full(len(frame), np.nan, dtype=np.float32) for seed in seeds}
    # One shuffled seed is sufficient as a falsification control; it is not a competing model.
    shuffled_prediction = np.full(len(frame), np.nan, dtype=np.float32)
    audit = []
    for held in targets:
        test = frame.target.eq(held).to_numpy()
        held_scaffolds = set(frame.loc[test, "scaffold"])
        train_mask = (~test) & (~frame.scaffold.isin(held_scaffolds).to_numpy())
        train_rows = np.flatnonzero(train_mask)
        seen_targets = np.unique(pair_target[train_rows])
        audit.append({
            "held_target": held,
            "train_pairs": int(train_mask.sum()),
            "test_pairs": int(test.sum()),
            "purged_pairs": int((~test).sum() - train_mask.sum()),
            "scaffold_overlap": 0,
            "held_labels_used": False,
            "held_pocket_used": False,
        })
        for seed in seeds:
            model = train(molecule_embedding, pocket_embedding, pair_mol, pair_target, labels,
                          train_rows, seen_targets, args, seed + 101 * pocket_index[held], False)
            with torch.inference_mode():
                score, _ = model.score_pairs(
                    molecule_embedding[torch.as_tensor(pair_mol[test])],
                    pocket_embedding[torch.as_tensor(pair_target[test])],
                )
            true_predictions[seed][test] = score.numpy()
        control = train(molecule_embedding, pocket_embedding, pair_mol, pair_target, labels,
                        train_rows, seen_targets, args, seeds[0] + 101 * pocket_index[held], True)
        with torch.inference_mode():
            score, _ = control.score_pairs(
                molecule_embedding[torch.as_tensor(pair_mol[test])],
                pocket_embedding[torch.as_tensor(pair_target[test])],
            )
        shuffled_prediction[test] = score.numpy()
        print(f"held={held} train={train_mask.sum()} test={test.sum()} purged={audit[-1]['purged_pairs']}", flush=True)

    candidate_arrays = [*true_predictions.values(), shuffled_prediction]
    if any(not np.isfinite(value).all() for value in candidate_arrays):
        raise RuntimeError("Missing or non-finite held-target predictions")
    adapted = np.mean(np.stack(list(true_predictions.values())), axis=0)
    raw_report = metrics(frame, raw, targets)
    adapted_report = metrics(frame, adapted, targets)
    shuffled_report = metrics(frame, shuffled_prediction, targets)
    per_seed = {str(seed): metrics(frame, prediction, targets)["macro"]
                for seed, prediction in true_predictions.items()}

    metric_names = list(raw_report["macro"])
    per_target_delta = {
        metric: np.asarray([
            adapted_report["per_target"][target][metric] - raw_report["per_target"][target][metric]
            for target in targets
        ]) for metric in metric_names
    }
    rng = np.random.default_rng(20260926)
    draws = rng.integers(0, len(targets), size=(args.bootstrap, len(targets)))
    ci = {
        metric: list(map(float, np.quantile(delta[draws].mean(axis=1), [0.025, 0.5, 0.975])))
        for metric, delta in per_target_delta.items()
    }
    deployable_checkpoints = []
    if args.deployable_dir is not None:
        args.deployable_dir.mkdir(parents=True, exist_ok=True)
        all_rows = np.arange(len(frame), dtype=np.int64)
        all_targets = np.unique(pair_target)
        for seed in seeds:
            deployable = train(
                molecule_embedding,
                pocket_embedding,
                pair_mol,
                pair_target,
                labels,
                all_rows,
                all_targets,
                args,
                seed,
                False,
            )
            checkpoint_path = args.deployable_dir / f"pacer_cgm_triplet_seed{seed}.pt"
            torch.save(
                {
                    "format_version": 1,
                    "method": "PACER-CGM pocket-anchored hard-triplet adapter",
                    "backbone": "Science-2026 litpcba_identity_90.pt frozen embeddings",
                    "state_dict": deployable.state_dict(),
                    "dim": int(molecule_embedding.shape[1]),
                    "rank": int(args.rank),
                    "experts": int(args.experts),
                    "seed": int(seed),
                    "targets": targets,
                    "objective": {
                        "classification": "target-balanced BCE",
                        "triplet": "anchor=pocket; positive=active; negative=top frozen-score inactive",
                        "margin": float(args.margin),
                        "hard_negatives_per_target": int(args.hard_negatives),
                        "ranking_weight": float(args.ranking_weight),
                        "preservation_weight": float(args.preserve_weight),
                    },
                    "claim_boundary": (
                        "Binding-ranking adapter only; not affinity, intrinsic agonism, "
                        "cooperativity, or PAM efficacy."
                    ),
                },
                checkpoint_path,
            )
            deployable_checkpoints.append(str(checkpoint_path))

    report = {
        "method": "PACER-CGM: pocket-context-gated low-rank bilinear metric residual",
        "backbone": "verified Science-2026 litpcba_identity_90.pt embeddings",
        "training_objective": (
            "target-balanced BCE + pocket-anchored hard-triplet ranking + bounded residual preservation"
        ),
        "triplet_definition": {
            "anchor": "target pocket embedding",
            "positive": "experimentally active molecule for that target",
            "negative": "highest frozen-DrugCLIP-scoring inactive molecule for that target",
        },
        "protocol": "leave-one-target-out with held-target Murcko scaffold purge; three-seed ensemble",
        "metrics": {
            "science2026_raw": raw_report,
            "pacer_cgm": adapted_report,
            "shuffled_label_control": shuffled_report,
            "per_seed": per_seed,
        },
        "delta_pacer_vs_raw_target_bootstrap_95ci": ci,
        "fold_audit": audit,
        "deployable_checkpoints": deployable_checkpoints,
        "hyperparameters": {name: getattr(args, name) for name in [
            "rank", "experts", "epochs", "lr", "temperature", "margin",
            "hard_negatives", "ranking_weight", "preserve_weight", "seeds"
        ]},
        "claim_boundary": "Checkpoint-matched target-transfer development pilot; not official full LIT-PCBA, SOTA, affinity, or PAM efficacy.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    prediction_table = frame[["pair_id", "target", "label", "scaffold"]].copy()
    prediction_table["science2026_raw"] = raw
    prediction_table["pacer_cgm"] = adapted
    prediction_table["shuffled_label_control"] = shuffled_prediction
    for seed, prediction in true_predictions.items():
        prediction_table[f"pacer_seed_{seed}"] = prediction
    prediction_table.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
