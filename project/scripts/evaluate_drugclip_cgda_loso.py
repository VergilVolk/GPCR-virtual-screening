#!/usr/bin/env python3
"""Context-Gated DrugCLIP Adapter (CGDA), strict target-LOO evaluation.

The frozen DrugCLIP molecular vectors are scored against target-specific
reference-ligand queries.  A small bank of shared low-rank experts adapts each
query; mixture weights are generated from the target pocket and reference
ligand context.  The held target contributes no activity labels to fitting.
"""
from __future__ import annotations

import argparse
import copy
import json
import pickle
from pathlib import Path

import numpy as np
import torch
from torch.nn import functional as F

from evaluate_drugclip_full_litpcba_reference_fusion import pocket_map
from evaluate_drugclip_official_litpcba_embeddings import metrics


def unit(x: torch.Tensor) -> torch.Tensor:
    return F.normalize(x, dim=-1)


class CGDA(torch.nn.Module):
    def __init__(self, dimension: int, rank: int, experts: int, seed: int):
        super().__init__()
        torch.manual_seed(seed)
        self.down = torch.nn.Parameter(torch.empty(experts, rank, dimension))
        self.up = torch.nn.Parameter(torch.zeros(experts, dimension, rank))
        torch.nn.init.normal_(self.down, std=0.02)
        self.gate = torch.nn.Sequential(
            torch.nn.Linear(dimension * 4, 32),
            torch.nn.GELU(),
            torch.nn.Linear(32, experts),
        )
        # Uniform gate and identity adapter at initialization.
        torch.nn.init.zeros_(self.gate[-1].weight)
        torch.nn.init.zeros_(self.gate[-1].bias)

    def context(self, pockets: torch.Tensor, references: torch.Tensor) -> torch.Tensor:
        p = unit(pockets.mean(dim=0, keepdim=True))[0]
        r = unit(references.mean(dim=0, keepdim=True))[0]
        return torch.cat([p, r, torch.abs(p - r), p * r])

    def adapted_references(self, pockets: torch.Tensor, references: torch.Tensor):
        references = unit(references)
        weights = torch.softmax(self.gate(self.context(pockets, references)), dim=-1)
        # E x R x rank -> E x R x D, followed by gated expert sum.
        latent = torch.einsum("ekd,nd->enk", self.down, references)
        delta = torch.einsum("edk,enk->end", self.up, latent)
        adapted = unit(references + torch.einsum("e,end->nd", weights, delta))
        return adapted, weights

    def score(self, molecules: torch.Tensor, pockets: torch.Tensor, references: torch.Tensor):
        queries, weights = self.adapted_references(pockets, references)
        return (unit(molecules) @ queries.T).max(dim=1).values, weights


def balanced_sample(entry, max_negatives: int, hard_fraction: float, rng: np.random.Generator) -> np.ndarray:
    labels = entry["labels"]
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    n_negative = min(len(negative), max(max_negatives, len(positive) * 20))
    n_hard = min(n_negative, int(round(n_negative * hard_fraction)))
    hard = np.empty(0, dtype=np.int64)
    if n_hard:
        order = np.argsort(entry["reference_score"][negative])[::-1]
        hard = negative[order[:n_hard]]
    remaining = np.setdiff1d(negative, hard, assume_unique=False)
    random = rng.choice(remaining, n_negative - n_hard, replace=False)
    return np.concatenate([positive, hard, random])


def load_data(root: Path, pockets: dict[str, np.ndarray]):
    data = {}
    for target in sorted(pockets):
        folder = root / target / "drugclip_emb"
        molecules, _, labels = pickle.load((folder / "mols.lmdb.pkl").open("rb"))
        references, _, _ = pickle.load((folder / "ligand.lmdb.pkl").open("rb"))
        data[target] = {
            "molecules": np.asarray(molecules, dtype=np.float32),
            "labels": np.asarray(labels, dtype=np.int64),
            "references": np.asarray(references, dtype=np.float32),
            "pockets": np.asarray(pockets[target], dtype=np.float32),
        }
        data[target]["reference_score"] = (
            data[target]["references"] @ data[target]["molecules"].T
        ).max(axis=0)
    return data


def fit_fold(data, held: str, args, fold_seed: int):
    rng = np.random.default_rng(fold_seed)
    seen = [target for target in sorted(data) if target != held]
    sampled = {}
    for target in seen:
        idx = balanced_sample(data[target], args.max_negatives_per_target,
                              args.hard_negative_fraction, rng)
        # Stratified internal validation is used only for early stopping.  The
        # outer held target remains completely invisible.
        train_parts, validation_parts = [], []
        for label in (0, 1):
            group = idx[data[target]["labels"][idx] == label]
            group = group[rng.permutation(len(group))]
            n_validation = max(1, int(round(len(group) * args.validation_fraction)))
            validation_parts.append(group[:n_validation])
            train_parts.append(group[n_validation:])
        sampled[target] = {
            "train": np.concatenate(train_parts),
            "validation": np.concatenate(validation_parts),
        }

    model = CGDA(128, args.rank, args.experts, fold_seed)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    bce = torch.nn.BCEWithLogitsLoss()
    def validation_score():
        values = []
        model.eval()
        with torch.no_grad():
            for target in seen:
                idx = sampled[target]["validation"]
                score, _ = model.score(
                    torch.as_tensor(data[target]["molecules"][idx]),
                    torch.as_tensor(data[target]["pockets"]),
                    torch.as_tensor(data[target]["references"]),
                )
                values.append(metrics(data[target]["labels"][idx], score.numpy())["bedroc_alpha80_5"])
        model.train()
        return float(np.mean(values))

    best = {"epoch": -1, "validation_bedroc": validation_score(),
            "state": copy.deepcopy(model.state_dict())}
    last = {}
    for epoch in range(args.epochs):
        epoch_bce, epoch_retrieval = [], []
        for target_index, target in enumerate(seen):
            idx = sampled[target]["train"]
            start = (epoch * args.batch_size) % len(idx)
            take = np.concatenate([idx[start:], idx[:start]])[: min(args.batch_size, len(idx))]
            molecules = torch.as_tensor(data[target]["molecules"][take])
            labels = torch.as_tensor(data[target]["labels"][take], dtype=torch.float32)
            pockets = torch.as_tensor(data[target]["pockets"])
            references = torch.as_tensor(data[target]["references"])

            optimizer.zero_grad()
            score, _ = model.score(molecules, pockets, references)
            loss_bce = bce(score / args.temperature, labels)
            positive_score = score[labels > 0.5]
            negative_score = score[labels < 0.5]
            if len(positive_score) and len(negative_score):
                loss_ranking = F.softplus(
                    args.ranking_margin - positive_score[:, None] + negative_score[None, :]
                ).mean()
            else:
                loss_ranking = score.sum() * 0

            active = molecules[labels > 0.5]
            if len(active):
                target_logits = []
                for other in seen:
                    other_score, _ = model.score(
                        active,
                        torch.as_tensor(data[other]["pockets"]),
                        torch.as_tensor(data[other]["references"]),
                    )
                    target_logits.append(other_score)
                logits = torch.stack(target_logits, dim=1) / args.retrieval_temperature
                truth = torch.full((len(active),), target_index, dtype=torch.long)
                loss_retrieval = F.cross_entropy(logits, truth)
            else:
                loss_retrieval = score.sum() * 0
            adapted, _ = model.adapted_references(pockets, references)
            preserve = (1 - (adapted * unit(references)).sum(dim=1)).mean()
            loss = (loss_bce + args.ranking_weight * loss_ranking +
                    args.retrieval_weight * loss_retrieval + args.preserve_weight * preserve)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            epoch_bce.append(float(loss_bce.detach()))
            epoch_retrieval.append(float(loss_retrieval.detach()))
        last = {"epoch": epoch, "bce": float(np.mean(epoch_bce)),
                "retrieval": float(np.mean(epoch_retrieval))}
        current_validation = validation_score()
        if current_validation > best["validation_bedroc"]:
            best = {"epoch": epoch, "validation_bedroc": current_validation,
                    "state": copy.deepcopy(model.state_dict())}
    model.load_state_dict(best.pop("state"))
    return model, {"final": last, "selected": best}


def evaluate_score(model, entry, chunk_size: int):
    pocket = torch.as_tensor(entry["pockets"])
    reference = torch.as_tensor(entry["references"])
    adapted, gate = model.adapted_references(pocket, reference)
    output = []
    with torch.no_grad():
        for start in range(0, len(entry["molecules"]), chunk_size):
            molecules = unit(torch.as_tensor(entry["molecules"][start : start + chunk_size]))
            output.append((molecules @ adapted.T).max(dim=1).values.numpy())
    return np.concatenate(output), gate.detach().numpy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pocket-root", type=Path, required=True)
    parser.add_argument("--pocket-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--rank", type=int, default=4)
    parser.add_argument("--experts", type=int, default=4)
    parser.add_argument("--epochs", type=int, default=120)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--max-negatives-per-target", type=int, default=3000)
    parser.add_argument("--hard-negative-fraction", type=float, default=0.0)
    parser.add_argument("--ranking-weight", type=float, default=0.0)
    parser.add_argument("--ranking-margin", type=float, default=0.1)
    parser.add_argument("--lr", type=float, default=2e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.12)
    parser.add_argument("--retrieval-temperature", type=float, default=0.07)
    parser.add_argument("--retrieval-weight", type=float, default=0.2)
    parser.add_argument("--preserve-weight", type=float, default=0.5)
    parser.add_argument("--validation-fraction", type=float, default=0.2)
    parser.add_argument("--chunk-size", type=int, default=65536)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    pockets = pocket_map(args.pocket_root, args.pocket_archive)
    data = load_data(args.root, pockets)
    targets = sorted(data)
    folds = {}
    score_dir = args.output.with_suffix("")
    score_dir.mkdir(parents=True, exist_ok=True)
    for fold_index, held in enumerate(targets):
        model, training = fit_fold(data, held, args, args.seed + fold_index * 1009)
        adapted_score, gate = evaluate_score(model, data[held], args.chunk_size)
        molecules = data[held]["molecules"]
        labels = data[held]["labels"]
        pocket_score = (data[held]["pockets"] @ molecules.T).max(axis=0)
        reference_score = (data[held]["references"] @ molecules.T).max(axis=0)
        folds[held] = {
            "n": int(len(labels)), "positives": int(labels.sum()),
            "pocket": metrics(labels, pocket_score),
            "reference": metrics(labels, reference_score),
            "cgda": metrics(labels, adapted_score),
            "gate": gate.tolist(), "training_final": training,
        }
        np.savez_compressed(score_dir / f"{held}.npz", labels=labels, pocket=pocket_score,
                            reference=reference_score, cgda=adapted_score)
        print(f"completed {fold_index + 1}/{len(targets)} {held}: "
              f"BEDROC={folds[held]['cgda']['bedroc_alpha80_5']:.4f}", flush=True)

    names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    macro = {method: {name: float(np.mean([folds[t][method][name] for t in targets])) for name in names}
             for method in ("pocket", "reference", "cgda")}
    report = {
        "method": "Context-Gated DrugCLIP Adapter (CGDA)",
        "protocol": "strict leave-one-entire-target-out; held labels evaluation-only",
        "n_targets": len(targets), "macro": macro, "folds": folds,
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "claim_boundary": "Reference-ligand-assisted retrospective binding screen; not PAM identity, efficacy, or prospective validation.",
    }
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(macro, indent=2))


if __name__ == "__main__":
    main()
