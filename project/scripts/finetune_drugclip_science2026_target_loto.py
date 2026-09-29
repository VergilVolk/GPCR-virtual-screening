#!/usr/bin/env python3
"""Target-held-out, scaffold-purged adaptation of Science-2026 DrugCLIP.

The encoder is frozen.  Shared rank-r residuals on the molecule and pocket
projection heads are trained on all but one target.  Any training row whose
Murcko scaffold occurs in the held target is removed.  The held target pocket
and labels are never used for fitting or model selection.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from torch.nn import functional as F

from finetune_drugclip_gpcr_screening import FastProjectionLoRA, binary_metrics
from finetune_drugclip_muscarinic_triplet import Proj


class Adapter(torch.nn.Module):
    def __init__(self, projection: dict, rank: int, seed: int,
                 pocket_experts: int = 1, pocket_expert_rank: int = 2):
        super().__init__()
        self.mol = FastProjectionLoRA(projection["mol_project"], rank, seed)
        self.pocket = (
            PocketConditionedMoE(projection["pocket_project"], pocket_experts,
                                 pocket_expert_rank, seed + 17)
            if pocket_experts > 1
            else FastProjectionLoRA(projection["pocket_project"], rank, seed + 17)
        )


class PocketConditionedMoE(torch.nn.Module):
    """Low-rank pocket experts selected from the frozen pocket context."""

    def __init__(self, state: dict, experts: int, rank: int, seed: int):
        super().__init__()
        torch.manual_seed(seed)
        self.experts = experts
        self.register_buffer("w1", state["linear1.weight"].clone())
        self.register_buffer("b1", state["linear1.bias"].clone())
        self.register_buffer("w2", state["linear2.weight"].clone())
        self.register_buffer("b2", state["linear2.bias"].clone())
        self.down = torch.nn.Parameter(torch.empty(experts, rank, self.w2.shape[1]))
        self.up = torch.nn.Parameter(torch.zeros(experts, self.w2.shape[0], rank))
        torch.nn.init.normal_(self.down, std=0.02)
        self.gate = torch.nn.Sequential(
            torch.nn.Linear(self.w2.shape[0], 32),
            torch.nn.GELU(),
            torch.nn.Linear(32, experts),
        )
        torch.nn.init.zeros_(self.gate[-1].weight)
        torch.nn.init.zeros_(self.gate[-1].bias)

    def hidden(self, representation: torch.Tensor) -> torch.Tensor:
        return F.relu(F.linear(representation, self.w1, self.b1))

    def base_from_hidden(self, hidden: torch.Tensor) -> torch.Tensor:
        return F.normalize(F.linear(hidden, self.w2, self.b2), dim=-1)

    def gate_weights(self, hidden: torch.Tensor) -> torch.Tensor:
        return torch.softmax(self.gate(self.base_from_hidden(hidden)), dim=-1)

    def from_hidden(self, hidden: torch.Tensor) -> torch.Tensor:
        base = F.linear(hidden, self.w2, self.b2)
        latent = torch.einsum("erh,nh->enr", self.down, hidden)
        delta = torch.einsum("edr,enr->end", self.up, latent)
        weights = self.gate_weights(hidden)
        return F.normalize(base + torch.einsum("ne,end->nd", weights, delta), dim=-1)


def scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or Chem.MolToSmiles(mol, canonical=True)


def balanced_weights(target: np.ndarray, label: np.ndarray, target_ids: np.ndarray) -> np.ndarray:
    out = np.zeros(len(label), dtype=np.float32)
    for t in target_ids:
        for y in (0, 1):
            keep = (target == t) & (label == y)
            if keep.any():
                out[keep] = 1.0 / (2 * len(target_ids) * keep.sum())
    out *= len(out) / out.sum()
    return out


def macro_metrics(frame: pd.DataFrame, score: np.ndarray, targets: list[str]) -> dict:
    per = {}
    for target in targets:
        keep = frame.target.eq(target).to_numpy()
        per[target] = binary_metrics(frame.loc[keep, "label"].to_numpy(int), score[keep])
    names = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]
    return {
        "macro": {name: float(np.mean([per[t][name] for t in targets])) for name in names},
        "per_target": per,
    }


def train_fold(
    mol_rep: torch.Tensor,
    pocket_rep: torch.Tensor,
    pair_mol: np.ndarray,
    pair_target: np.ndarray,
    labels: np.ndarray,
    train_rows: np.ndarray,
    seen_targets: np.ndarray,
    projection: dict,
    seed: int,
    args: argparse.Namespace,
    pair_fingerprints: list | None = None,
    shuffle_labels: bool = False,
) -> Adapter:
    rng = np.random.default_rng(seed)
    fit_labels = labels[train_rows].copy()
    fit_targets = pair_target[train_rows]
    if shuffle_labels:
        for target in seen_targets:
            local = np.flatnonzero(fit_targets == target)
            fit_labels[local] = rng.permutation(fit_labels[local])

    model = Adapter(projection, args.rank, seed, args.pocket_experts, args.pocket_expert_rank)
    hm = model.mol.hidden(mol_rep).detach()
    hp = model.pocket.hidden(pocket_rep).detach()
    with torch.no_grad():
        base_mol = model.mol.base_from_hidden(hm)
        base_pocket = model.pocket.base_from_hidden(hp)
        base_score = (base_mol @ base_pocket.T).cpu().numpy()

    rows_t = torch.as_tensor(train_rows, dtype=torch.long)
    mol_t = torch.as_tensor(pair_mol[train_rows], dtype=torch.long)
    target_t = torch.as_tensor(fit_targets, dtype=torch.long)
    label_t = torch.as_tensor(fit_labels, dtype=torch.float32)
    weight_t = torch.as_tensor(balanced_weights(fit_targets, fit_labels, seen_targets))
    preserve_mol = torch.as_tensor(np.unique(pair_mol[train_rows]), dtype=torch.long)
    preserve_target = torch.as_tensor(seen_targets, dtype=torch.long)

    rank_pairs = []
    matched_pairs = []
    for target in seen_targets:
        local = np.flatnonzero(fit_targets == target)
        pos = local[fit_labels[local] == 1]
        neg = local[fit_labels[local] == 0]
        if len(pos) and len(neg):
            if args.matched_hard_negatives or args.hybrid_matched_triplet:
                if pair_fingerprints is None:
                    raise ValueError("pair_fingerprints are required for matched hard negatives")
                paired_pos = []
                paired_neg = []
                global_neg = train_rows[neg]
                neg_fps = [pair_fingerprints[i] for i in global_neg]
                for positive in pos:
                    similarities = np.asarray(
                        DataStructs.BulkTanimotoSimilarity(
                            pair_fingerprints[train_rows[positive]], neg_fps
                        ),
                        dtype=np.float32,
                    )
                    order = np.argsort(similarities)[::-1][
                        : min(args.matched_negative_k, len(neg))
                    ]
                    paired_pos.extend([int(positive)] * len(order))
                    paired_neg.extend(neg[order].tolist())
                matched_pairs.append((
                    torch.as_tensor(paired_pos, dtype=torch.long),
                    torch.as_tensor(paired_neg, dtype=torch.long),
                ))
                if args.matched_hard_negatives:
                    continue
            if not args.dynamic_hard_negatives:
                hard_order = np.argsort(base_score[pair_mol[train_rows][neg], target])[::-1]
                neg = neg[hard_order[: min(args.hard_negatives, len(neg))]]
            rank_pairs.append((torch.as_tensor(pos), torch.as_tensor(neg), False))

    unique_active = {}
    for m, t, y in zip(pair_mol[train_rows], fit_targets, fit_labels):
        if y:
            unique_active.setdefault(int(m), set()).add(int(t))
    retrieval = [(m, next(iter(ts))) for m, ts in unique_active.items() if len(ts) == 1]
    seen_map = {int(t): i for i, t in enumerate(seen_targets)}
    retrieval = [(m, seen_map[t]) for m, t in retrieval if t in seen_map]
    retrieval_mol = torch.as_tensor([v[0] for v in retrieval], dtype=torch.long)
    retrieval_y = torch.as_tensor([v[1] for v in retrieval], dtype=torch.long)
    seen_t = torch.as_tensor(seen_targets, dtype=torch.long)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        optimizer.zero_grad()
        zm = model.mol.from_hidden(hm)
        zp = model.pocket.from_hidden(hp)
        matrix = zm @ zp.T
        pair_score = matrix[mol_t, target_t]
        bce = (
            F.binary_cross_entropy_with_logits(pair_score / args.temperature, label_t, reduction="none")
            * weight_t
        ).mean()
        ranking_terms = []
        listwise_terms = []
        for pos, negative_pool, paired in rank_pairs:
            if paired:
                ranking_terms.append(
                    F.softplus(args.margin - pair_score[pos] + pair_score[negative_pool]).mean()
                )
                continue
            if args.dynamic_hard_negatives and len(negative_pool) > args.hard_negatives:
                local_score = pair_score[negative_pool].detach()
                hard_local = torch.topk(local_score, args.hard_negatives, largest=True).indices
                neg = negative_pool[hard_local]
            else:
                neg = negative_pool
            ranking_terms.append(
                F.softplus(args.margin - pair_score[pos, None] + pair_score[neg][None, :]).mean()
            )
            if args.listwise_weight > 0:
                tau = args.listwise_temperature
                positive_peak = tau * (
                    torch.logsumexp(pair_score[pos] / tau, dim=0) - math.log(len(pos))
                )
                negative_peak = tau * (
                    torch.logsumexp(pair_score[neg] / tau, dim=0) - math.log(len(neg))
                )
                listwise_terms.append(F.softplus((negative_peak - positive_peak) / tau))
        ranking = torch.stack(ranking_terms).mean() if ranking_terms else pair_score.sum() * 0
        matched_terms = [
            F.softplus(args.margin - pair_score[pos] + pair_score[neg]).mean()
            for pos, neg in matched_pairs
        ]
        matched_ranking = (
            torch.stack(matched_terms).mean() if matched_terms else pair_score.sum() * 0
        )
        listwise = torch.stack(listwise_terms).mean() if listwise_terms else pair_score.sum() * 0
        retrieval_loss = (
            F.cross_entropy(matrix[retrieval_mol][:, seen_t] / args.retrieval_temperature, retrieval_y)
            if len(retrieval_mol) else pair_score.sum() * 0
        )
        preserve = (
            (1 - (zm[preserve_mol] * base_mol[preserve_mol]).sum(1)).mean()
            + (1 - (zp[preserve_target] * base_pocket[preserve_target]).sum(1)).mean()
        )
        if args.pocket_experts > 1:
            gate_mean = model.pocket.gate_weights(hp[seen_t]).mean(dim=0)
            gate_balance = ((gate_mean - 1.0 / args.pocket_experts) ** 2).sum()
        else:
            gate_balance = pair_score.sum() * 0
        loss = (
            bce
            + args.ranking_weight * ranking
            + args.matched_ranking_weight * matched_ranking
            + args.listwise_weight * listwise
            + args.retrieval_weight * retrieval_loss
            + args.preserve_weight * preserve
            + args.gate_balance_weight * gate_balance
        )
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
    return model


def predict_target(model: Adapter, mol_rep: torch.Tensor, pocket_rep: torch.Tensor,
                   molecule_indices: np.ndarray, target_index: int) -> np.ndarray:
    with torch.inference_mode():
        zm = model.mol.from_hidden(model.mol.hidden(mol_rep))
        zp = model.pocket.from_hidden(model.pocket.hidden(pocket_rep))
        return (zm[torch.as_tensor(molecule_indices)] * zp[target_index]).sum(1).cpu().numpy()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", default="20260925,20260926,20260927")
    parser.add_argument("--rank", type=int, default=8)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--lr", type=float, default=3e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--retrieval-temperature", type=float, default=0.10)
    parser.add_argument("--margin", type=float, default=0.10)
    parser.add_argument("--hard-negatives", type=int, default=256)
    parser.add_argument("--ranking-weight", type=float, default=0.25)
    parser.add_argument("--listwise-weight", type=float, default=0.0)
    parser.add_argument("--listwise-temperature", type=float, default=0.05)
    parser.add_argument("--dynamic-hard-negatives", action="store_true")
    parser.add_argument(
        "--matched-hard-negatives",
        action="store_true",
        help="Pair every active with the most ECFP-similar inactive inside the training fold.",
    )
    parser.add_argument("--matched-negative-k", type=int, default=1)
    parser.add_argument(
        "--hybrid-matched-triplet",
        action="store_true",
        help="Keep global score-hard triplets and add ECFP-matched counterfactual triplets.",
    )
    parser.add_argument("--matched-ranking-weight", type=float, default=0.25)
    parser.add_argument("--pocket-experts", type=int, default=1)
    parser.add_argument("--pocket-expert-rank", type=int, default=2)
    parser.add_argument("--gate-balance-weight", type=float, default=0.05)
    parser.add_argument("--retrieval-weight", type=float, default=0.20)
    parser.add_argument("--preserve-weight", type=float, default=0.25)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument(
        "--deployable-dir",
        type=Path,
        default=None,
        help="Fit final all-development-target adapters and save materialized projection checkpoints.",
    )
    args = parser.parse_args()
    if args.matched_hard_negatives and args.dynamic_hard_negatives:
        parser.error("--matched-hard-negatives and --dynamic-hard-negatives are mutually exclusive")
    if args.matched_hard_negatives and args.hybrid_matched_triplet:
        parser.error("choose either matched-only or hybrid matched triplets")
    if args.matched_negative_k < 1:
        parser.error("--matched-negative-k must be positive")
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    archive = np.load(args.representations, allow_pickle=False)
    projection = torch.load(args.projection, map_location="cpu")
    molecule_ids = list(map(str, archive["molecule_ids"]))
    pocket_ids = list(map(str, archive["pocket_ids"]))
    targets = sorted(set(pd.read_csv(args.pairs).target.astype(str)) & set(pocket_ids))
    target_index = {target: pocket_ids.index(target) for target in targets}
    molecule_index = {smiles: i for i, smiles in enumerate(molecule_ids)}
    frame = pd.read_csv(args.pairs)
    frame = frame[
        frame.target.astype(str).isin(targets)
        & frame.canonical_smiles.astype(str).isin(molecule_index)
    ].copy().reset_index(drop=True)
    frame["scaffold"] = [scaffold(s) for s in frame.canonical_smiles.astype(str)]
    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    pair_fingerprints = [
        morgan.GetFingerprint(Chem.MolFromSmiles(smiles))
        for smiles in frame.canonical_smiles.astype(str)
    ]
    pair_mol = np.asarray([molecule_index[s] for s in frame.canonical_smiles.astype(str)])
    pair_target = np.asarray([target_index[t] for t in frame.target.astype(str)])
    labels = frame.label.to_numpy(int)
    mol_rep = torch.as_tensor(archive["molecule_representations"].astype(np.float32))
    pocket_rep = torch.as_tensor(archive["pocket_representations"].astype(np.float32))

    frozen_mol = Proj(copy.deepcopy(projection["mol_project"])).eval()
    frozen_pocket = Proj(copy.deepcopy(projection["pocket_project"])).eval()
    with torch.inference_mode():
        raw_matrix = frozen_mol(mol_rep) @ frozen_pocket(pocket_rep).T
        raw = raw_matrix[torch.as_tensor(pair_mol), torch.as_tensor(pair_target)].cpu().numpy()

    seeds = [int(value) for value in args.seeds.split(",")]
    predictions = {seed: np.full(len(frame), np.nan, dtype=np.float32) for seed in seeds}
    shuffled = {seed: np.full(len(frame), np.nan, dtype=np.float32) for seed in seeds}
    fold_audit = []
    for held in targets:
        test = frame.target.eq(held).to_numpy()
        held_scaffolds = set(frame.loc[test, "scaffold"])
        train = (~test) & (~frame.scaffold.isin(held_scaffolds).to_numpy())
        train_rows = np.flatnonzero(train)
        seen_targets = np.unique(pair_target[train_rows])
        fold_audit.append({
            "held_target": held,
            "train_pairs": int(train.sum()),
            "test_pairs": int(test.sum()),
            "purged_pairs": int((~test).sum() - train.sum()),
            "held_labels_used": False,
            "held_pocket_used_in_loss": False,
            "scaffold_overlap": 0,
        })
        for seed in seeds:
            model = train_fold(mol_rep, pocket_rep, pair_mol, pair_target, labels, train_rows,
                               seen_targets, projection, seed + 101 * target_index[held], args,
                               pair_fingerprints=pair_fingerprints)
            predictions[seed][test] = predict_target(
                model, mol_rep, pocket_rep, pair_mol[test], target_index[held]
            )
            control = train_fold(mol_rep, pocket_rep, pair_mol, pair_target, labels, train_rows,
                                 seen_targets, projection, seed + 101 * target_index[held], args,
                                 pair_fingerprints=pair_fingerprints,
                                 shuffle_labels=True)
            shuffled[seed][test] = predict_target(
                control, mol_rep, pocket_rep, pair_mol[test], target_index[held]
            )
        print(f"held={held} train={train.sum()} test={test.sum()} purged={fold_audit[-1]['purged_pairs']}", flush=True)

    if any(not np.isfinite(value).all() for value in [*predictions.values(), *shuffled.values()]):
        raise RuntimeError("Non-finite or missing out-of-target predictions")
    adapted = np.mean(np.stack(list(predictions.values())), axis=0)
    label_control = np.mean(np.stack(list(shuffled.values())), axis=0)
    reports = {
        "science2026_raw": macro_metrics(frame, raw, targets),
        "pacer_loto": macro_metrics(frame, adapted, targets),
        "shuffled_label_control": macro_metrics(frame, label_control, targets),
        "per_seed": {str(seed): macro_metrics(frame, score, targets)["macro"]
                     for seed, score in predictions.items()},
    }

    metric_names = list(reports["science2026_raw"]["macro"])
    rng = np.random.default_rng(20260926)
    # Per-target metrics are already frozen above. Bootstrap those nine paired
    # deltas directly instead of recomputing ranking metrics hundreds of
    # thousands of times inside the resampling loop.
    per_target_delta = {
        name: np.asarray([
            reports["pacer_loto"]["per_target"][target][name]
            - reports["science2026_raw"]["per_target"][target][name]
            for target in targets
        ], dtype=np.float64)
        for name in metric_names
    }
    draws = rng.integers(0, len(targets), size=(args.bootstrap, len(targets)))
    ci = {
        name: list(map(float, np.quantile(delta[draws].mean(axis=1), [0.025, 0.5, 0.975])))
        for name, delta in per_target_delta.items()
    }
    deployable_checkpoints = []
    if args.deployable_dir is not None:
        args.deployable_dir.mkdir(parents=True, exist_ok=True)
        all_rows = np.arange(len(frame), dtype=np.int64)
        all_targets = np.unique(pair_target)
        for seed in seeds:
            deployable = train_fold(
                mol_rep,
                pocket_rep,
                pair_mol,
                pair_target,
                labels,
                all_rows,
                all_targets,
                projection,
                seed,
                args,
                pair_fingerprints=pair_fingerprints,
                shuffle_labels=False,
            )
            checkpoint_path = args.deployable_dir / f"pacer_projection_triplet_seed{seed}.pt"
            torch.save(
                {
                    "format_version": 1,
                    "method": "Science-2026 DrugCLIP dual-projection hard-triplet LoRA",
                    "mol_project": deployable.mol.materialized_state(),
                    "pocket_project": deployable.pocket.materialized_state(),
                    "base_projection": str(args.projection),
                    "seed": int(seed),
                    "training_targets": targets,
                    "n_pairs": int(len(frame)),
                    "objective": {
                        "classification": "target-balanced BCE",
                        "triplet": (
                            "anchor=pocket; positive=active; negative=ECFP-nearest inactive"
                            if args.matched_hard_negatives
                            else "anchor=pocket; positive=active; negative=top frozen-score inactive"
                        ),
                        "target_retrieval_weight": float(args.retrieval_weight),
                        "ranking_weight": float(args.ranking_weight),
                        "listwise_weight": float(args.listwise_weight),
                        "listwise_temperature": float(args.listwise_temperature),
                        "dynamic_hard_negatives": bool(args.dynamic_hard_negatives),
                        "matched_hard_negatives": bool(args.matched_hard_negatives),
                        "matched_negative_k": int(args.matched_negative_k),
                        "hybrid_matched_triplet": bool(args.hybrid_matched_triplet),
                        "matched_ranking_weight": float(args.matched_ranking_weight),
                        "preservation_weight": float(args.preserve_weight),
                    },
                    "claim_boundary": (
                        "Binding-ranking projection only; not affinity, cooperativity, "
                        "intrinsic agonism, or PAM efficacy."
                    ),
                },
                checkpoint_path,
            )
            deployable_checkpoints.append(str(checkpoint_path))

    report = {
        "method": (
            "Science-2026 frozen encoder + global molecule LoRA + pocket-conditioned "
            "low-rank experts + dual global/chemotype triplets"
            if args.pocket_experts > 1 and args.hybrid_matched_triplet
            else
            "Science-2026 frozen encoder + shared dual-projection rank-8 LoRA + "
            "global-hard and chemotype-matched dual triplet"
            if args.hybrid_matched_triplet
            else
            "Science-2026 frozen encoder + shared dual-projection rank-8 LoRA + "
            "chemotype-matched hard triplet"
            if args.matched_hard_negatives
            else
            "Science-2026 frozen encoder + shared dual-projection rank-8 LoRA + "
            "dynamic hard-negative top-heavy triplet"
            if args.listwise_weight > 0 or args.dynamic_hard_negatives
            else "Science-2026 frozen encoder + shared dual-projection rank-8 LoRA"
        ),
        "protocol": "Leave-one-target-out; all held-target Murcko scaffolds purged from training; 3-seed mean",
        "targets": targets,
        "n_pairs": len(frame),
        "n_molecules": len(np.unique(pair_mol)),
        "metrics": reports,
        "delta_pacer_vs_raw_target_bootstrap_95ci": ci,
        "fold_audit": fold_audit,
        "deployable_checkpoints": deployable_checkpoints,
        "hyperparameters": {key: getattr(args, key) for key in [
            "rank", "epochs", "lr", "temperature", "retrieval_temperature", "margin",
            "hard_negatives", "ranking_weight", "listwise_weight", "listwise_temperature",
            "dynamic_hard_negatives", "retrieval_weight", "preserve_weight", "seeds"
            , "matched_hard_negatives", "matched_negative_k", "hybrid_matched_triplet",
            "matched_ranking_weight", "pocket_experts", "pocket_expert_rank",
            "gate_balance_weight"
        ]},
        "claim_boundary": "Checkpoint-matched target-transfer pilot; not official full LIT-PCBA, SOTA, affinity, or PAM efficacy.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pred = frame[["pair_id", "target", "label", "scaffold"]].copy()
    pred["science2026_raw"] = raw
    pred["pacer_loto"] = adapted
    pred["shuffled_label_control"] = label_control
    for seed, score in predictions.items():
        pred[f"pacer_seed_{seed}"] = score
    pred.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
