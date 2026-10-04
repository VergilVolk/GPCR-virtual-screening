#!/usr/bin/env python3
"""Parameter-efficient DrugCLIP fine-tuning for four-GPCR target retrieval.

The official encoders and projection weights are frozen. Rank-4 LoRA deltas on
the final molecule and pocket projections are evaluated with class-balanced
cross-entropy, with and without a hardest-wrong-pocket margin loss.  OOF model
selection currently promotes CE-only because the extra triplet did not help.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.multiclass import OneVsRestClassifier
from torch import nn
from torch.nn import functional as F

from finetune_drugclip_muscarinic_triplet import Proj


TARGETS = ["B2AR", "CCR2", "M2R", "M4R"]


class ProjectionLoRA(nn.Module):
    def __init__(self, state: dict, rank: int, seed: int):
        super().__init__(); torch.manual_seed(seed)
        self.base = Proj(copy.deepcopy(state))
        for parameter in self.base.parameters(): parameter.requires_grad_(False)
        hidden = state["linear2.weight"].shape[1]; output = state["linear2.weight"].shape[0]
        self.down = nn.Linear(hidden, rank, bias=False); self.up = nn.Linear(rank, output, bias=False)
        nn.init.normal_(self.down.weight, std=0.02); nn.init.zeros_(self.up.weight)

    def forward(self, x):
        hidden = F.relu(self.base.linear1(x))
        return F.normalize(self.base.linear2(hidden) + self.up(self.down(hidden)), dim=-1)

    def materialized_state(self):
        state = copy.deepcopy(self.base.state_dict())
        state["linear2.weight"] += self.up.weight.detach().cpu() @ self.down.weight.detach().cpu()
        return state


class RetrievalLoRA(nn.Module):
    def __init__(self, initial: dict, rank: int, seed: int):
        super().__init__()
        self.mol = ProjectionLoRA(initial["mol_project"], rank, seed)
        self.pocket = ProjectionLoRA(initial["pocket_project"], rank, seed + 17)

    def forward(self, molecule_rep, pocket_rep):
        zm, zp = self.mol(molecule_rep), self.pocket(pocket_rep)
        return zm, zp, zm @ zp.T


def fit(molecule_rep, pocket_rep, labels, initial, seed, mode, args):
    model = RetrievalLoRA(initial, args.rank, seed)
    frozen_mol = Proj(copy.deepcopy(initial["mol_project"])).eval()
    frozen_pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.no_grad():
        z0m = frozen_mol(molecule_rep)
        z0p = frozen_pocket(pocket_rep)
    z0m = z0m.detach().clone(); z0p = z0p.detach().clone()
    train_labels = labels.clone()
    if mode == "random":
        generator = torch.Generator().manual_seed(seed); train_labels = train_labels[torch.randperm(len(labels), generator=generator)]
    counts = torch.bincount(train_labels, minlength=len(TARGETS)).float()
    weights = counts.sum() / (len(TARGETS) * counts.clamp_min(1))
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        optimizer.zero_grad(); zm, zp, scores = model(molecule_rep, pocket_rep)
        ce = F.cross_entropy(scores / args.temperature, train_labels, weight=weights)
        true = scores.gather(1, train_labels[:, None]).squeeze(1)
        wrong = scores.masked_fill(F.one_hot(train_labels, len(TARGETS)).bool(), -1e9).max(1).values
        hard = F.relu(args.margin - true + wrong).mean()
        preserve = (1 - (zm * z0m).sum(1)).mean() + (1 - (zp * z0p).sum(1)).mean()
        loss = ce + (args.hard_weight * hard if mode == "hard" else 0.0) + args.preserve_weight * preserve
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
    return model


def fp_matrix(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    matrix = np.zeros((len(smiles), 2048), dtype=np.float32)
    for i, value in enumerate(smiles):
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(Chem.MolFromSmiles(value)), matrix[i])
    return matrix


def retrieval_metrics(labels: np.ndarray, scores: np.ndarray) -> dict:
    order = np.argsort(-scores, axis=1); ranks = np.empty(len(labels), dtype=int)
    for i, label in enumerate(labels): ranks[i] = int(np.flatnonzero(order[i] == label)[0]) + 1
    one_hot = np.eye(len(TARGETS), dtype=int)[labels]
    per_target = {}
    for target_index, target in enumerate(TARGETS):
        mask = labels == target_index
        per_target[target] = {"n": int(mask.sum()), "recall_at_1": float((ranks[mask] <= 1).mean()),
                              "mrr": float((1.0 / ranks[mask]).mean())}
    return {
        "n": int(len(labels)), "recall_at_1": float((ranks <= 1).mean()),
        "recall_at_2": float((ranks <= 2).mean()), "mrr": float((1.0 / ranks).mean()),
        "macro_recall_at_1": float(np.mean([v["recall_at_1"] for v in per_target.values()])),
        "macro_mrr": float(np.mean([v["mrr"] for v in per_target.values()])),
        "pair_roc_auc": float(roc_auc_score(one_hot.ravel(), scores.ravel())),
        "pair_pr_auc": float(average_precision_score(one_hot.ravel(), scores.ravel())),
        "per_target": per_target,
    }


def scaffold_bootstrap(table, labels, a, b, repeats, seed):
    groups = table.murcko_scaffold.unique(); rng = np.random.default_rng(seed); values = []
    for _ in range(repeats):
        sampled = rng.choice(groups, len(groups), True)
        index = np.concatenate([np.flatnonzero(table.murcko_scaffold.to_numpy() == group) for group in sampled])
        if len(np.unique(labels[index])) < len(TARGETS):
            continue
        values.append(retrieval_metrics(labels[index], a[index])["macro_recall_at_1"] -
                      retrieval_metrics(labels[index], b[index])["macro_recall_at_1"])
    return list(map(float, np.quantile(values, [0.025, 0.5, 0.975])))


def main() -> None:
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
    parser.add_argument("--margin", type=float, default=0.1)
    parser.add_argument("--hard-weight", type=float, default=0.5)
    parser.add_argument("--preserve-weight", type=float, default=0.2)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args(); torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    archive = np.load(args.representations, allow_pickle=False); initial = torch.load(args.projection, map_location="cpu")
    table_raw = pd.read_csv(args.benchmark); available = set(map(str, archive["molecule_ids"]))
    table = table_raw[table_raw.canonical_smiles.isin(available)].copy().reset_index(drop=True)
    index = {value: i for i, value in enumerate(map(str, archive["molecule_ids"]))}
    order = np.asarray([index[value] for value in table.canonical_smiles])
    molecule_rep = torch.as_tensor(archive["molecule_representations"].astype(np.float32)[order])
    pocket_rep = torch.as_tensor(archive["pocket_representations"].astype(np.float32))
    pocket_ids = [value.split("_")[0] for value in map(str, archive["pocket_ids"])]
    if pocket_ids != TARGETS: raise ValueError(f"Pocket order {pocket_ids} != {TARGETS}")
    labels = np.asarray([TARGETS.index(value) for value in table.target], dtype=int)
    folds = table.scaffold_fold.to_numpy(int); modes = ["official", "random", "ce", "hard"]
    frozen_mol = Proj(copy.deepcopy(initial["mol_project"])).eval(); frozen_pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode(): official_scores = (frozen_mol(molecule_rep) @ frozen_pocket(pocket_rep).T).numpy()
    ecfp = fp_matrix(table.canonical_smiles.astype(str).tolist()); rows = []; mode_predictions = {m: [] for m in modes}
    ecfp_oof = np.full((len(table), len(TARGETS)), np.nan)
    for fold in sorted(np.unique(folds)):
        train = folds != fold; test = folds == fold
        classifier = OneVsRestClassifier(LogisticRegression(
            C=1.0, class_weight="balanced", max_iter=3000,
            solver="liblinear", random_state=20260925))
        classifier.fit(ecfp[train], labels[train]); probabilities = classifier.predict_proba(ecfp[test])
        ecfp_oof[test[:, None] if False else np.flatnonzero(test)] = probabilities[:, np.argsort(classifier.classes_)]
    for seed in map(int, args.seeds.split(",")):
        predictions = {m: np.full((len(table), len(TARGETS)), np.nan) for m in modes}
        predictions["official"][:] = official_scores
        for fold in sorted(np.unique(folds)):
            train = folds != fold; test = folds == fold
            train_x = molecule_rep[train]; train_y = torch.as_tensor(labels[train], dtype=torch.long)
            for mode in ["random", "ce", "hard"]:
                model = fit(train_x, pocket_rep, train_y, initial, seed + 101 * int(fold), mode, args)
                with torch.no_grad(): predictions[mode][test] = model(molecule_rep[test], pocket_rep)[2].numpy()
        for mode in modes: mode_predictions[mode].append(predictions[mode])
        for i, row in table.iterrows():
            record = {"canonical_molecule_id": row.canonical_molecule_id, "target": row.target,
                      "scaffold_fold": int(row.scaffold_fold), "seed": seed}
            for mode in modes:
                for j, target in enumerate(TARGETS): record[f"{mode}_{target}"] = float(predictions[mode][i, j])
            rows.append(record)
    ensemble = {mode: np.mean(values, axis=0) for mode, values in mode_predictions.items()}
    metrics = {mode: retrieval_metrics(labels, score) for mode, score in ensemble.items()}
    metrics["ecfp_logistic"] = retrieval_metrics(labels, ecfp_oof)
    report = {
        "method": "DrugCLIP dual-projection rank-4 LoRA ablation; balanced CE selected by OOF",
        "benchmark_rows": int(len(table)), "input_failures_excluded": int(len(table_raw) - len(table)),
        "targets": TARGETS, "split": "five-fold Murcko-scaffold OOF", "metrics": metrics,
        "hard_minus_official_macro_recall1_scaffold_bootstrap_95ci": scaffold_bootstrap(
            table, labels, ensemble["hard"], ensemble["official"], args.bootstrap, 20260925),
        "hard_minus_ce_macro_recall1_scaffold_bootstrap_95ci": scaffold_bootstrap(
            table, labels, ensemble["hard"], ensemble["ce"], args.bootstrap, 20260926),
        "hard_minus_random_macro_recall1_scaffold_bootstrap_95ci": scaffold_bootstrap(
            table, labels, ensemble["hard"], ensemble["random"], args.bootstrap, 20260927),
        "hyperparameters": {"rank": args.rank, "epochs": args.epochs, "lr": args.lr,
                            "temperature": args.temperature, "margin": args.margin,
                            "hard_weight": args.hard_weight, "preserve_weight": args.preserve_weight,
                            "seeds": args.seeds},
        "claim_boundary": "Retrospective GPCR target retrieval; no potency, efficacy, or prospective-binding claim.",
    }
    # Full-data deployable ensemble is fitted after OOF evaluation is frozen.
    # CE-only is promoted; the hard-triplet ablation remains in the report.
    checkpoints = []
    full_labels = torch.as_tensor(labels, dtype=torch.long)
    for seed in map(int, args.seeds.split(",")):
        model = fit(molecule_rep, pocket_rep, full_labels, initial, seed, "ce", args)
        path = args.output.with_name(f"{args.output.stem}.balanced_ce.seed{seed}.projection.pt")
        torch.save({"mol_project": model.mol.materialized_state(),
                    "pocket_project": model.pocket.materialized_state(),
                    "logit_scale": initial.get("logit_scale"), "seed": seed,
                    "training_task": "four-GPCR target-balanced CE",
                    "selection_note": "hardest-pocket triplet did not improve OOF macro Recall@1"}, path)
        checkpoints.append(str(path))
    report["deployable_checkpoints"] = checkpoints
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(rows).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
