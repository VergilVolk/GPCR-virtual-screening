#!/usr/bin/env python3
"""Target-LOO pocket-to-co-crystal-ligand alignment for full LIT-PCBA.

Training supervision consists only of pocket/reference-ligand set membership on
seen targets.  The held target contributes neither reference ligands nor assay
labels.  A rank-r residual maps pocket embeddings toward the molecular shared
space with a target-balanced multi-positive contrastive objective.
"""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import lmdb
import numpy as np
import torch
from torch.nn import functional as F

from evaluate_drugclip_official_litpcba_embeddings import metrics


class PocketAligner(torch.nn.Module):
    def __init__(self, dimension, rank, seed):
        super().__init__()
        torch.manual_seed(seed)
        self.down = torch.nn.Linear(dimension, rank, bias=False)
        self.up = torch.nn.Linear(rank, dimension, bias=False)
        torch.nn.init.normal_(self.down.weight, std=0.02)
        torch.nn.init.zeros_(self.up.weight)

    def forward(self, pocket):
        return F.normalize(pocket + self.up(self.down(pocket)), dim=-1)


def load_pockets(root, archive_path):
    archive = np.load(archive_path, allow_pickle=False)
    embeddings = np.asarray(archive["pocket_embeddings"], dtype=np.float32)
    files = sorted(root.glob("*.pockets.lmdb"))
    counts = []
    for path in files:
        env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
        with env.begin() as txn:
            counts.append(txn.stat()["entries"])
        env.close()
    offsets = np.cumsum([0] + counts)
    return {path.name.removesuffix(".pockets.lmdb"): embeddings[start:end]
            for path, start, end in zip(files, offsets[:-1], offsets[1:])}


def target_balanced_multi_positive(logits, pocket_targets, reference_targets, seen):
    losses = []
    for target in seen:
        query = pocket_targets == target
        positive = reference_targets == target
        target_logits = logits[query]
        numerator = torch.logsumexp(target_logits[:, positive], dim=1)
        denominator = torch.logsumexp(target_logits, dim=1)
        losses.append(-(numerator - denominator).mean())
    return torch.stack(losses).mean()


def fit(pockets, references, pocket_targets, reference_targets, seen, seed, args):
    model = PocketAligner(pockets.shape[1], args.rank, seed)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    pocket = torch.as_tensor(pockets)
    reference = F.normalize(torch.as_tensor(references), dim=-1)
    pt = torch.as_tensor(pocket_targets, dtype=torch.long)
    rt = torch.as_tensor(reference_targets, dtype=torch.long)
    seen_pocket = torch.as_tensor(np.isin(pocket_targets, seen))
    seen_reference = torch.as_tensor(np.isin(reference_targets, seen))
    p = pocket[seen_pocket]
    r = reference[seen_reference]
    p_target = pt[seen_pocket]
    r_target = rt[seen_reference]
    remap = {target: index for index, target in enumerate(seen)}
    p_target = torch.as_tensor([remap[int(value)] for value in p_target], dtype=torch.long)
    r_target = torch.as_tensor([remap[int(value)] for value in r_target], dtype=torch.long)
    local_seen = list(range(len(seen)))
    last = None
    for epoch in range(args.epochs):
        optimizer.zero_grad()
        mapped = model(p)
        logits = mapped @ r.T / args.temperature
        contrastive = target_balanced_multi_positive(logits, p_target, r_target, local_seen)
        preserve = (1 - (mapped * p).sum(dim=1)).mean()
        loss = contrastive + args.preserve_weight * preserve
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last = {"epoch": epoch, "contrastive": float(contrastive.detach()),
                "preserve": float(preserve.detach()), "loss": float(loss.detach())}
    return model, last


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pocket-root", type=Path, required=True)
    parser.add_argument("--pocket-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--rank", type=int, default=8)
    parser.add_argument("--epochs", type=int, default=500)
    parser.add_argument("--lr", type=float, default=1e-2)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--preserve-weight", type=float, default=0.05)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    pocket_map = load_pockets(args.pocket_root, args.pocket_archive)
    targets = sorted(pocket_map)
    pocket_rows, pocket_targets, reference_rows, reference_targets = [], [], [], []
    for target_index, target in enumerate(targets):
        references, _, _ = pickle.load(
            (args.root / target / "drugclip_emb" / "ligand.lmdb.pkl").open("rb"))
        references = np.asarray(references, dtype=np.float32)
        pocket_rows.append(pocket_map[target])
        reference_rows.append(references)
        pocket_targets.extend([target_index] * len(pocket_map[target]))
        reference_targets.extend([target_index] * len(references))
    pockets = np.concatenate(pocket_rows)
    references = np.concatenate(reference_rows)
    pocket_targets = np.asarray(pocket_targets)
    reference_targets = np.asarray(reference_targets)

    official_metrics, aligned_metrics, reference_metrics, audit = {}, {}, {}, {}
    for held, target in enumerate(targets):
        seen = [index for index in range(len(targets)) if index != held]
        model, last = fit(
            pockets, references, pocket_targets, reference_targets, seen,
            args.seed + held * 1009, args)
        molecules, _, labels = pickle.load(
            (args.root / target / "drugclip_emb" / "mols.lmdb.pkl").open("rb"))
        target_references, _, _ = pickle.load(
            (args.root / target / "drugclip_emb" / "ligand.lmdb.pkl").open("rb"))
        molecules = np.asarray(molecules, dtype=np.float32)
        labels = np.asarray(labels, dtype=np.int64)
        held_pockets = torch.as_tensor(pocket_map[target])
        with torch.inference_mode():
            aligned_pockets = model(held_pockets).numpy()
        official_score = (pocket_map[target] @ molecules.T).max(axis=0)
        aligned_score = (aligned_pockets @ molecules.T).max(axis=0)
        reference_score = (np.asarray(target_references, dtype=np.float32) @ molecules.T).max(axis=0)
        official_metrics[target] = metrics(labels, official_score)
        aligned_metrics[target] = metrics(labels, aligned_score)
        reference_metrics[target] = metrics(labels, reference_score)
        score_path = args.output.with_name(f"{args.output.stem}.{target}.scores.npz")
        np.savez_compressed(score_path, labels=labels, official=official_score,
                            aligned=aligned_score, reference=reference_score)
        audit[target] = {
            "held_reference_in_training": False,
            "held_assay_labels_in_training": False,
            "seen_targets": [targets[index] for index in seen],
            "last_epoch": last,
            "trainable_parameters": int(sum(p.numel() for p in model.parameters() if p.requires_grad)),
        }
        print(f"completed {held + 1}/{len(targets)} {target}", flush=True)

    names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    macro = lambda block: {name: float(np.mean([block[target][name] for target in targets])) for name in names}
    report = {
        "method": "Rank-8 pocket-to-reference-ligand multi-positive contrastive alignment",
        "protocol": "Full LIT-PCBA leave-one-entire-target-out; held references and assay labels absent from training",
        "n_targets": len(targets),
        "n_pairs": int(sum(value["n"] for value in official_metrics.values())),
        "metrics": {
            "official_pocket": {"macro": macro(official_metrics), "per_target": official_metrics},
            "aligned_pocket": {"macro": macro(aligned_metrics), "per_target": aligned_metrics},
            "reference_ligand_upper_context": {"macro": macro(reference_metrics), "per_target": reference_metrics},
        },
        "audit": audit,
        "hyperparameters": {key: getattr(args, key) for key in [
            "seed", "rank", "epochs", "lr", "temperature", "preserve_weight"]},
        "claim_boundary": "Retrospective cross-modal target transfer; reference-free on held target but trained with co-crystal ligands from other targets; no PAM efficacy claim.",
    }
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({key: value["macro"] for key, value in report["metrics"].items()}, indent=2))


if __name__ == "__main__":
    main()
