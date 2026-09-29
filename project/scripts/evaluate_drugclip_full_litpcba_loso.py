#!/usr/bin/env python3
"""Full 15-target LIT-PCBA leave-one-target-out low-rank DrugCLIP adaptation.

Uses the official 128-dimensional precomputed DrugCLIP embeddings.  For each
fold, one complete target is excluded, exact held-target SMILES are purged from
the other targets, and a shared molecule/pocket low-rank residual is trained on
the remaining targets.  Evaluation uses every molecule of the held target.
"""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import lmdb
import torch
from rdkit.ML.Scoring.Scoring import CalcBEDROC
from sklearn.metrics import average_precision_score, roc_auc_score
from torch.nn import functional as F


def split_smiles(names):
    return np.asarray([str(value).rsplit(" ", 1)[0] for value in names], dtype=object)


def load_target(folder, pocket_map):
    mol, names, labels = pickle.load((folder / "drugclip_emb" / "mols.lmdb.pkl").open("rb"))
    return (np.asarray(mol, dtype=np.float32), split_smiles(names),
            np.asarray(labels, dtype=np.int64), pocket_map[folder.name])


def load_pocket_map(pocket_root, pocket_archive):
    archive = np.load(pocket_archive, allow_pickle=False)
    embeddings = np.asarray(archive["pocket_embeddings"], dtype=np.float32)
    files = sorted(pocket_root.glob("*.pockets.lmdb"))
    counts = []
    for path in files:
        env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
        with env.begin() as txn:
            counts.append(txn.stat()["entries"])
        env.close()
    if sum(counts) != len(embeddings):
        raise ValueError(f"Pocket count mismatch: {sum(counts)} != {len(embeddings)}")
    offsets = np.cumsum([0] + counts)
    return {path.name.removesuffix(".pockets.lmdb"): embeddings[start:end]
            for path, start, end in zip(files, offsets[:-1], offsets[1:])}


class LowRankResidual(torch.nn.Module):
    def __init__(self, dimension, rank, seed):
        super().__init__()
        torch.manual_seed(seed)
        self.down = torch.nn.Linear(dimension, rank, bias=False)
        self.up = torch.nn.Linear(rank, dimension, bias=False)
        torch.nn.init.normal_(self.down.weight, std=0.02)
        torch.nn.init.zeros_(self.up.weight)

    def forward(self, values):
        return F.normalize(values + self.up(self.down(values)), dim=-1)


class Model(torch.nn.Module):
    def __init__(self, dimension, rank, seed, adapter_side):
        super().__init__()
        self.molecule = LowRankResidual(dimension, rank, seed)
        self.pocket = LowRankResidual(dimension, rank, seed + 17)
        if adapter_side == "molecule":
            for parameter in self.pocket.parameters():
                parameter.requires_grad_(False)
        elif adapter_side == "pocket":
            for parameter in self.molecule.parameters():
                parameter.requires_grad_(False)


def target_scores(molecules, transformed_pockets):
    return torch.stack([(molecules @ pocket.T).max(dim=1).values for pocket in transformed_pockets], dim=1)


def metrics(labels, scores):
    ranked = [[float(score), int(label)] for score, label in sorted(
        zip(scores, labels), key=lambda item: item[0], reverse=True)]
    out = {
        "n": int(len(labels)), "actives": int(labels.sum()),
        "roc_auc": float(roc_auc_score(labels, scores)),
        "pr_auc": float(average_precision_score(labels, scores)),
        "bedroc_alpha80_5": float(CalcBEDROC(ranked, 1, 80.5)),
    }
    for fraction in (0.005, 0.01, 0.02, 0.05):
        n = max(1, int(np.ceil(len(labels) * fraction)))
        order = np.argsort(-scores, kind="mergesort")[:n]
        out[f"ef{fraction:g}"] = float(labels[order].mean() / labels.mean())
    return out


def build_training_pools(root, targets, held_smiles, seed, negative_pool, pocket_map):
    rng = np.random.default_rng(seed)
    pools = []
    audit = {}
    for target in targets:
        molecules, smiles, labels, pockets = load_target(root / target, pocket_map)
        keep = ~np.isin(smiles, held_smiles)
        purged = int((~keep).sum())
        molecules, labels = molecules[keep], labels[keep]
        positive = np.flatnonzero(labels == 1)
        negative = np.flatnonzero(labels == 0)
        if len(negative) > negative_pool:
            negative = rng.choice(negative, size=negative_pool, replace=False)
        pools.append({
            "positive": torch.as_tensor(molecules[positive]),
            "negative": torch.as_tensor(molecules[negative]),
            "pocket": torch.as_tensor(pockets),
        })
        audit[target] = {
            "positive_pool": int(len(positive)),
            "negative_pool": int(len(negative)),
            "exact_smiles_purged": purged,
        }
    return pools, audit


def sample_rows(values, n, generator):
    index = torch.randint(len(values), (n,), generator=generator)
    return values[index]


def fit(pools, dimension, seed, args):
    model = Model(dimension, args.rank, seed, args.adapter_side)
    optimizer = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=args.lr, weight_decay=args.weight_decay)
    generator = torch.Generator().manual_seed(seed)
    last = None
    for epoch in range(args.epochs):
        optimizer.zero_grad()
        transformed_pockets = [model.pocket(pool["pocket"]) for pool in pools]
        losses = []
        ranking_losses = []
        preserve_terms = []
        active_batches = []
        active_targets = []
        for target_index, pool in enumerate(pools):
            positive = sample_rows(pool["positive"], args.batch_per_class, generator)
            hard_count = int(round(args.batch_per_class * args.hard_negative_fraction))
            random_count = args.batch_per_class - hard_count
            negative_parts = []
            if hard_count:
                with torch.no_grad():
                    negative_embedding = model.molecule(pool["negative"])
                    negative_score = (negative_embedding @ transformed_pockets[target_index].T).max(dim=1).values
                    hard_index = torch.topk(negative_score, k=min(hard_count, len(negative_score))).indices
                    negative_parts.append(pool["negative"][hard_index])
            if random_count:
                negative_parts.append(sample_rows(pool["negative"], random_count, generator))
            negative = torch.cat(negative_parts)
            batch = torch.cat([positive, negative])
            labels = torch.cat([torch.ones(len(positive)), torch.zeros(len(negative))])
            molecule = model.molecule(batch)
            preserve_terms.append((1 - (molecule * batch).sum(dim=1)).mean())
            score = (molecule @ transformed_pockets[target_index].T).max(dim=1).values
            losses.append(F.binary_cross_entropy_with_logits(score / args.temperature, labels))
            positive_score = score[:len(positive)]
            negative_score = score[len(positive):]
            ranking_losses.append(F.relu(args.ranking_margin - positive_score + negative_score).mean())
            active_batches.append(model.molecule(positive))
            active_targets.extend([target_index] * len(positive))
        bce = torch.stack(losses).mean()
        ranking = torch.stack(ranking_losses).mean()
        pocket_preserve = torch.stack([
            (1 - (transformed * pool["pocket"]).sum(dim=1)).mean()
            for transformed, pool in zip(transformed_pockets, pools)]).mean()
        preserve = torch.stack(preserve_terms).mean() + pocket_preserve
        active = torch.cat(active_batches)
        retrieval_score = target_scores(active, transformed_pockets) / args.temperature
        retrieval = F.cross_entropy(retrieval_score, torch.as_tensor(active_targets, dtype=torch.long))
        loss = (bce + args.retrieval_weight * retrieval + args.ranking_weight * ranking
                + args.preserve_weight * preserve)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last = {"epoch": epoch, "bce": float(bce.detach()), "retrieval": float(retrieval.detach()),
                "ranking": float(ranking.detach()),
                "preserve": float(preserve.detach()), "loss": float(loss.detach())}
    return model, last


def score_full(model, molecules, pockets, chunk_size):
    pocket = model.pocket(torch.as_tensor(pockets))
    scores = []
    with torch.inference_mode():
        for start in range(0, len(molecules), chunk_size):
            batch = model.molecule(torch.as_tensor(molecules[start:start + chunk_size]))
            scores.append((batch @ pocket.T).max(dim=1).values.numpy())
    return np.concatenate(scores)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pocket-root", type=Path, required=True)
    parser.add_argument("--pocket-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--rank", type=int, default=8)
    parser.add_argument("--epochs", type=int, default=60)
    parser.add_argument("--batch-per-class", type=int, default=128)
    parser.add_argument("--negative-pool", type=int, default=10000)
    parser.add_argument("--lr", type=float, default=5e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--retrieval-weight", type=float, default=0.25)
    parser.add_argument("--preserve-weight", type=float, default=0.0)
    parser.add_argument("--adapter-side", choices=["dual", "molecule", "pocket"], default="dual")
    parser.add_argument("--hard-negative-fraction", type=float, default=0.0)
    parser.add_argument("--ranking-weight", type=float, default=0.0)
    parser.add_argument("--ranking-margin", type=float, default=0.1)
    parser.add_argument("--chunk-size", type=int, default=32768)
    args = parser.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    targets = sorted(path.name for path in args.root.iterdir()
                     if path.is_dir() and (path / "drugclip_emb" / "mols.lmdb.pkl").exists())
    pocket_map = load_pocket_map(args.pocket_root, args.pocket_archive)
    if set(targets) != set(pocket_map):
        raise ValueError("Target mismatch between molecule and pocket archives")
    official_metrics, adapted_metrics, fold_audit = {}, {}, {}
    prediction_files = {}
    for held_index, held in enumerate(targets):
        held_mol, held_smiles, held_labels, held_pocket = load_target(args.root / held, pocket_map)
        official_score = (held_pocket @ held_mol.T).max(axis=0)
        seen = [target for target in targets if target != held]
        pools, pool_audit = build_training_pools(
            args.root, seen, set(held_smiles.tolist()),
            args.seed + held_index * 1009, args.negative_pool, pocket_map)
        model, last = fit(pools, held_mol.shape[1], args.seed + held_index * 1009, args)
        adapted_score = score_full(model, held_mol, held_pocket, args.chunk_size)
        official_metrics[held] = metrics(held_labels, official_score)
        adapted_metrics[held] = metrics(held_labels, adapted_score)
        prediction_path = args.output.with_name(f"{args.output.stem}.{held}.scores.npz")
        np.savez_compressed(prediction_path, labels=held_labels, official=official_score, adapted=adapted_score)
        prediction_files[held] = str(prediction_path)
        fold_audit[held] = {
            "held_target_in_loss": False,
            "held_pairs": int(len(held_labels)),
            "seen_targets": seen,
            "training_pools": pool_audit,
            "last_epoch": last,
            "trainable_parameters": int(sum(p.numel() for p in model.parameters() if p.requires_grad)),
        }
        print(f"completed {held_index + 1}/{len(targets)} {held}", flush=True)

    names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    macro = lambda block: {name: float(np.mean([block[target][name] for target in targets])) for name in names}
    report = {
        "method": "Official 128D DrugCLIP molecule and pocket embedding low-rank residual with BCE, target retrieval, and optional anchoring",
        "protocol": "Full 15-target LIT-PCBA leave-one-entire-target-out; exact held-target SMILES purged; all held candidates evaluated",
        "n_targets": len(targets),
        "n_pairs": int(sum(value["n"] for value in official_metrics.values())),
        "metrics": {
            "official": {"macro": macro(official_metrics), "per_target": official_metrics},
            "adapted": {"macro": macro(adapted_metrics), "per_target": adapted_metrics},
        },
        "fold_audit": fold_audit,
        "prediction_files": prediction_files,
        "hyperparameters": {key: getattr(args, key) for key in [
            "seed", "rank", "epochs", "batch_per_class", "negative_pool", "lr",
            "temperature", "retrieval_weight", "preserve_weight", "adapter_side",
            "hard_negative_fraction", "ranking_weight", "ranking_margin"]},
        "claim_boundary": "Retrospective full-benchmark target transfer; exact-SMILES rather than scaffold purge; no PAM efficacy or prospective claim.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({key: value["macro"] for key, value in report["metrics"].items()}, indent=2))


if __name__ == "__main__":
    main()
