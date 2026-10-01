#!/usr/bin/env python3
"""Fine-tune DrugCLIP's own projection heads for M4 functional PAM ranking.

This is deliberately different from a frozen-embedding classifier/adapter:
the trainable modules are initialized from and exported back to the official
``mol_project`` and ``pocket_project`` keys.  The large Uni-Mol encoders stay
frozen in this CPU gate.  A later GPU gate may unfreeze their final blocks only
after this low-capacity experiment passes grouped out-of-fold evaluation.

The targeted loss mirrors the local metric-learning logic used in DreaMS:
M4 pocket states are anchors, functional PAMs are positives, and experimental
inactive molecules that the original DrugCLIP ranks highly are hard negatives.
A random-negative arm with the same number of pairs is an explicit control.
"""
from __future__ import annotations

import argparse
import copy
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from sklearn.metrics import average_precision_score, roc_auc_score
from torch import nn
from torch.nn import functional as F


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)


class DrugCLIPProjection(nn.Module):
    """Exact architecture of DrugCLIP's internal NonLinearHead."""

    def __init__(self, state: dict[str, torch.Tensor]):
        super().__init__()
        in_dim = state["linear1.weight"].shape[1]
        hidden = state["linear1.weight"].shape[0]
        out_dim = state["linear2.weight"].shape[0]
        self.linear1 = nn.Linear(in_dim, hidden)
        self.linear2 = nn.Linear(hidden, out_dim)
        self.load_state_dict(state)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return F.normalize(self.linear2(F.relu(self.linear1(x))), dim=-1)


class M4ProjectionTune(nn.Module):
    def __init__(self, initial: dict):
        super().__init__()
        self.mol_project = DrugCLIPProjection(initial["mol_project"])
        self.pocket_project = DrugCLIPProjection(initial["pocket_project"])
        self.log_scale = nn.Parameter(initial["logit_scale"].float().reshape(()).clone())
        self.bias = nn.Parameter(torch.tensor(0.0))

    def encode(self, molecule_rep: torch.Tensor, pocket_rep: torch.Tensor):
        return self.mol_project(molecule_rep), self.pocket_project(pocket_rep)

    def scores(self, molecule_rep: torch.Tensor, pocket_rep: torch.Tensor):
        zm, zp = self.encode(molecule_rep, pocket_rep)
        state_scores = zm @ zp.T
        # Multi-instance state aggregation: a PAM need not prefer every state.
        retrieval = 0.1 * torch.logsumexp(state_scores / 0.1, dim=1)
        logits = self.log_scale.exp().clamp(max=100.0) * retrieval + self.bias
        return zm, zp, state_scores, retrieval, logits


def metric(y: np.ndarray, p: np.ndarray) -> dict:
    return {
        "roc_auc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
        "prevalence": float(np.mean(y)),
    }


def paired_bootstrap_auc(y, a, b, seed=20260924, n=2000):
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if np.unique(y[idx]).size < 2:
            continue
        values.append(roc_auc_score(y[idx], a[idx]) - roc_auc_score(y[idx], b[idx]))
    return [float(x) for x in np.quantile(values, [0.025, 0.5, 0.975])]


def train_one(x, pockets, y, mode, seed, initial, args):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    x = torch.as_tensor(x, dtype=torch.float32)
    pockets = torch.as_tensor(pockets, dtype=torch.float32)
    y = torch.as_tensor(y, dtype=torch.float32)
    model = M4ProjectionTune(copy.deepcopy(initial))
    frozen = M4ProjectionTune(copy.deepcopy(initial)).eval()
    for parameter in frozen.parameters():
        parameter.requires_grad_(False)
    with torch.no_grad():
        z0m, z0p = frozen.encode(x, pockets)
        _, _, _, original_retrieval, _ = frozen.scores(x, pockets)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    # DrugCLIP's frozen environment uses an older PyTorch without flatnonzero.
    positive_idx = torch.nonzero(y > 0.5, as_tuple=False).reshape(-1)
    negative_idx = torch.nonzero(y < 0.5, as_tuple=False).reshape(-1)
    if not len(positive_idx) or not len(negative_idx):
        raise ValueError("Training fold must contain both functional PAM and inactive molecules")
    if mode == "hard_triplet":
        negative_order = negative_idx[torch.argsort(original_retrieval[negative_idx], descending=True)]
    else:
        generator = torch.Generator().manual_seed(seed)
        negative_order = negative_idx[torch.randperm(len(negative_idx), generator=generator)]
    selected_negative = negative_order[: min(args.negative_count, len(negative_order))]
    pos_weight = torch.tensor(float(len(negative_idx) / len(positive_idx)))
    for _ in range(args.epochs):
        model.train(); opt.zero_grad()
        zm, zp, _, retrieval, logits = model.scores(x, pockets)
        classification = F.binary_cross_entropy_with_logits(logits, y, pos_weight=pos_weight)
        loss = classification
        if mode != "bce":
            positives = retrieval[positive_idx]
            if len(positives) > args.positive_count:
                choose = torch.randperm(len(positives))[: args.positive_count]
                positives = positives[choose]
            negatives = retrieval[selected_negative]
            triplet = F.relu(args.margin - positives[:, None] + negatives[None, :]).mean()
            loss = loss + args.triplet_weight * triplet
        preserve = (1 - (zm * z0m).sum(1)).mean() + (1 - (zp * z0p).sum(1)).mean()
        loss = loss + args.preserve_weight * preserve
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
    return model


def predict(model, x, pockets):
    model.eval()
    with torch.no_grad():
        *_, logits = model.scores(torch.as_tensor(x, dtype=torch.float32),
                                  torch.as_tensor(pockets, dtype=torch.float32))
    return torch.sigmoid(logits).cpu().numpy()


def exported_state(model, initial, metadata):
    return {
        "mol_project": model.mol_project.state_dict(),
        "pocket_project": model.pocket_project.state_dict(),
        "logit_scale": model.log_scale.detach().cpu().reshape(1),
        "classification_bias": model.bias.detach().cpu(),
        "source_checkpoint": initial.get("source_checkpoint"),
        "metadata": metadata,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--representations", type=Path, required=True)
    ap.add_argument("--projection", type=Path, required=True)
    ap.add_argument("--benchmark", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--protocol", choices=["source", "series"], required=True)
    ap.add_argument("--seeds", default="20260924,20260925,20260926")
    ap.add_argument("--epochs", type=int, default=250)
    ap.add_argument("--lr", type=float, default=2e-5)
    ap.add_argument("--weight-decay", type=float, default=1e-3)
    ap.add_argument("--margin", type=float, default=0.05)
    ap.add_argument("--triplet-weight", type=float, default=0.5)
    ap.add_argument("--preserve-weight", type=float, default=0.5)
    ap.add_argument("--negative-count", type=int, default=32)
    ap.add_argument("--positive-count", type=int, default=128)
    ap.add_argument("--bootstrap", type=int, default=2000)
    args = ap.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    archive = np.load(args.representations, allow_pickle=False)
    if "molecule_representations" not in archive or "pocket_representations" not in archive:
        raise ValueError("Archive lacks pre-projection representations; rerun extraction script")
    initial = torch.load(args.projection, map_location="cpu")
    ids = list(map(str, archive["molecule_ids"])); lookup = {v: i for i, v in enumerate(ids)}
    table = pd.read_csv(args.benchmark)
    smiles = table.canonical_smiles.map(canonical)
    missing = [v for v in smiles if v not in lookup]
    if missing:
        raise ValueError(f"Missing molecules: {missing[:5]}")
    order = np.asarray([lookup[v] for v in smiles])
    x = archive["molecule_representations"].astype(np.float32)[order]
    pockets = archive["pocket_representations"].astype(np.float32)
    y = table.target.to_numpy(np.int64)
    fold_col = "source_fold" if args.protocol == "source" else "series_holdout_fold"
    folds = table[fold_col].to_numpy(int)
    assigned = folds >= 0
    seeds = [int(s) for s in args.seeds.split(",")]
    modes = ["bce", "random_triplet", "hard_triplet"]
    records, per_seed = [], []
    for seed in seeds:
        predictions = {mode: np.full(len(y), np.nan) for mode in modes}
        for fold in sorted(np.unique(folds[assigned])):
            train = assigned & (folds != fold); test = assigned & (folds == fold)
            for mode in modes:
                model = train_one(x[train], pockets, y[train], mode, seed + int(fold) * 101,
                                  initial, args)
                predictions[mode][test] = predict(model, x[test], pockets)
        seed_metrics = {mode: metric(y[assigned], predictions[mode][assigned]) for mode in modes}
        per_seed.append({"seed": seed, "metrics": seed_metrics})
        for index in np.flatnonzero(assigned):
            records.append({
                "canonical_molecule_id": table.iloc[index].canonical_molecule_id,
                "target": int(y[index]), "fold": int(folds[index]), "seed": seed,
                **{mode: float(predictions[mode][index]) for mode in modes},
            })
    pred_table = pd.DataFrame(records)
    ensemble = pred_table.groupby("canonical_molecule_id", sort=False)[modes].mean()
    truth = table.set_index("canonical_molecule_id").loc[ensemble.index, "target"].to_numpy(int)
    aggregate = {mode: metric(truth, ensemble[mode].to_numpy()) for mode in modes}
    aggregate["hard_minus_random_auc_bootstrap_95ci"] = paired_bootstrap_auc(
        truth, ensemble.hard_triplet.to_numpy(), ensemble.random_triplet.to_numpy(), n=args.bootstrap)
    aggregate["hard_minus_bce_auc_bootstrap_95ci"] = paired_bootstrap_auc(
        truth, ensemble.hard_triplet.to_numpy(), ensemble.bce.to_numpy(), n=args.bootstrap)

    # Fit a deployable projection checkpoint only after all OOF predictions are frozen.
    final_model = train_one(x[assigned], pockets, y[assigned], "hard_triplet", seeds[0], initial, args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_path = args.output.with_suffix(".projection_finetuned.pt")
    metadata = {"target": "CHRM4 functional PAM", "protocol": args.protocol,
                "training_molecules": int(assigned.sum()), "mode": "hard_triplet",
                "evidence_level": "retrospective_internal_projection_finetune"}
    torch.save(exported_state(final_model, initial, metadata), checkpoint_path)
    report = {
        "method": "DrugCLIP internal projection-head fine-tuning",
        "protocol": args.protocol,
        "n_assigned": int(assigned.sum()),
        "class_counts": {"pam": int(y[assigned].sum()), "inactive": int((y[assigned] == 0).sum())},
        "pocket_ids": list(map(str, archive["pocket_ids"])),
        "aggregate_ensemble": aggregate,
        "per_seed": per_seed,
        "checkpoint": str(checkpoint_path),
        "trainable_scope": ["mol_project", "pocket_project", "logit_scale", "classification_bias"],
        "claim_boundary": (
            "True DrugCLIP projection-space fine-tuning, but frozen Uni-Mol encoders and retrospective "
            "labels. Performance is not wet-lab PAM validation and not evidence of full-backbone gain."
        ),
    }
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pred_table.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
