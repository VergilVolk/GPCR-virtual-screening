#!/usr/bin/env python3
"""Parameter-efficient DrugCLIP fine-tuning for GPCR active/decoy screening.

Official encoders and projection weights remain frozen.  Rank-r residuals on
the final molecule and pocket projections are trained with target- and
label-balanced BCE under global Murcko-scaffold OOF.  A shuffled-label adapter
and ECFP4 logistic models are evaluated under the identical split.
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
from rdkit.ML.Scoring.Scoring import CalcBEDROC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from torch import nn
from torch.nn import functional as F

from finetune_drugclip_muscarinic_triplet import Proj


TARGETS = ["B2AR", "CCR2", "M2R", "M4R"]


class FastProjectionLoRA(nn.Module):
    def __init__(self, state: dict, rank: int, seed: int):
        super().__init__()
        torch.manual_seed(seed)
        self.register_buffer("w1", state["linear1.weight"].clone())
        self.register_buffer("b1", state["linear1.bias"].clone())
        self.register_buffer("w2", state["linear2.weight"].clone())
        self.register_buffer("b2", state["linear2.bias"].clone())
        self.down = nn.Linear(self.w2.shape[1], rank, bias=False)
        self.up = nn.Linear(rank, self.w2.shape[0], bias=False)
        nn.init.normal_(self.down.weight, std=0.02)
        nn.init.zeros_(self.up.weight)

    def hidden(self, representation):
        return F.relu(F.linear(representation, self.w1, self.b1))

    def from_hidden(self, hidden):
        return F.normalize(F.linear(hidden, self.w2, self.b2) + self.up(self.down(hidden)), dim=-1)

    def base_from_hidden(self, hidden):
        return F.normalize(F.linear(hidden, self.w2, self.b2), dim=-1)

    def materialized_state(self):
        return {"linear1.weight": self.w1.detach().cpu(), "linear1.bias": self.b1.detach().cpu(),
                "linear2.weight": (self.w2 + self.up.weight @ self.down.weight).detach().cpu(),
                "linear2.bias": self.b2.detach().cpu()}


class ScreeningAdapter(nn.Module):
    def __init__(self, initial, rank, seed):
        super().__init__()
        self.molecule = FastProjectionLoRA(initial["mol_project"], rank, seed)
        self.pocket = FastProjectionLoRA(initial["pocket_project"], rank, seed + 17)


def ef(labels, scores, fraction):
    n = max(1, int(np.ceil(len(labels) * fraction)))
    order = np.argsort(-scores, kind="mergesort")[:n]
    return float(labels[order].mean() / labels.mean())


def binary_metrics(labels, scores):
    ranked = [[float(score), int(label)] for score, label in sorted(
        zip(scores, labels), key=lambda pair: pair[0], reverse=True)]
    return {"n": int(len(labels)), "actives": int(labels.sum()),
            "roc_auc": float(roc_auc_score(labels, scores)),
            "pr_auc": float(average_precision_score(labels, scores)),
            "bedroc_alpha20": float(CalcBEDROC(ranked, 1, 20.0)),
            "ef1pct": ef(labels, scores, 0.01), "ef5pct": ef(labels, scores, 0.05)}


def screening_metrics(table, scores):
    per_target = {target: binary_metrics(table.loc[table.target == target, "label"].to_numpy(int),
                                         scores[table.target.to_numpy() == target]) for target in TARGETS}
    names = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]
    return {"macro": {name: float(np.mean([per_target[target][name] for target in TARGETS])) for name in names},
            "per_target": per_target}


def fp_matrix(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    out = np.zeros((len(smiles), 2048), dtype=np.float32)
    for i, value in enumerate(smiles):
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(Chem.MolFromSmiles(value)), out[i])
    return out


def balanced_weights(targets, labels):
    weights = np.zeros(len(labels), dtype=np.float32)
    for target in range(len(TARGETS)):
        for label in [0, 1]:
            keep = (targets == target) & (labels == label)
            if keep.any():
                weights[keep] = 1.0 / (len(TARGETS) * 2 * keep.sum())
    weights *= len(labels) / weights.sum()
    return weights


def fit(train_pair_mol, train_pair_target, train_labels, molecule_rep, pocket_rep,
        initial, seed, args, random_labels=False):
    model = ScreeningAdapter(initial, args.rank, seed)
    with torch.no_grad():
        hm = model.molecule.hidden(molecule_rep)
        hp = model.pocket.hidden(pocket_rep)
        z0m = model.molecule.base_from_hidden(hm)
        z0p = model.pocket.base_from_hidden(hp)
    labels = train_labels.copy()
    if random_labels:
        rng = np.random.default_rng(seed)
        for target in range(len(TARGETS)):
            keep = np.flatnonzero(train_pair_target == target)
            labels[keep] = rng.permutation(labels[keep])
    weights = torch.as_tensor(balanced_weights(train_pair_target, labels))
    pair_mol = torch.as_tensor(train_pair_mol, dtype=torch.long)
    pair_target = torch.as_tensor(train_pair_target, dtype=torch.long)
    # The preservation regularizer must not inspect held-out molecules.  Using
    # every molecule representation here would make the otherwise scaffold-OOF
    # experiment transductive, even though no held-out labels are used.
    preserve_mol = torch.as_tensor(np.unique(train_pair_mol), dtype=torch.long)
    preserve_target = torch.as_tensor(np.unique(train_pair_target), dtype=torch.long)
    y = torch.as_tensor(labels, dtype=torch.float32)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        optimizer.zero_grad()
        zm, zp = model.molecule.from_hidden(hm), model.pocket.from_hidden(hp)
        scores = (zm[pair_mol] * zp[pair_target]).sum(1) / args.temperature
        bce = (F.binary_cross_entropy_with_logits(scores, y, reduction="none") * weights).mean()
        preserve = ((1 - (zm[preserve_mol] * z0m[preserve_mol]).sum(1)).mean()
                    + (1 - (zp[preserve_target] * z0p[preserve_target]).sum(1)).mean())
        loss = bce + args.preserve_weight * preserve
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
    return model


def predict(model, molecule_rep, pocket_rep, pair_mol, pair_target):
    with torch.inference_mode():
        hm = model.molecule.hidden(molecule_rep); hp = model.pocket.hidden(pocket_rep)
        zm = model.molecule.from_hidden(hm); zp = model.pocket.from_hidden(hp)
        return (zm[torch.as_tensor(pair_mol)] * zp[torch.as_tensor(pair_target)]).sum(1).numpy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", default="20260925")
    parser.add_argument("--rank", type=int, default=4)
    parser.add_argument("--epochs", type=int, default=150)
    parser.add_argument("--lr", type=float, default=5e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--preserve-weight", type=float, default=0.2)
    args = parser.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    archive = np.load(args.representations, allow_pickle=False)
    initial = torch.load(args.projection, map_location="cpu")
    pairs_raw = pd.read_csv(args.pairs)
    mol_index = {value: i for i, value in enumerate(map(str, archive["molecule_ids"]))}
    pairs = pairs_raw[pairs_raw.canonical_smiles.astype(str).isin(mol_index)].copy().reset_index(drop=True)
    pair_mol = np.asarray([mol_index[value] for value in pairs.canonical_smiles.astype(str)], dtype=int)
    pair_target = np.asarray([TARGETS.index(value) for value in pairs.target], dtype=int)
    labels = pairs.label.to_numpy(int); folds = pairs.scaffold_fold.to_numpy(int)
    molecule_rep = torch.as_tensor(archive["molecule_representations"].astype(np.float32))
    all_pocket_rep = archive["pocket_representations"].astype(np.float32)
    raw_pocket_ids = list(map(str, archive["pocket_ids"]))
    pocket_ids = [value.split("_")[0] for value in raw_pocket_ids]
    if pocket_ids == TARGETS:
        pocket_rep = torch.as_tensor(all_pocket_rep)
    elif set(raw_pocket_ids) == {
            f"{target}_cluster{cluster}" for target in TARGETS for cluster in range(10)}:
        cluster0_ids = [f"{target}_cluster0" for target in TARGETS]
        pocket_rep = torch.as_tensor(
            all_pocket_rep[[raw_pocket_ids.index(value) for value in cluster0_ids]]
        )
    else:
        raise ValueError(f"Unsupported pocket order: {raw_pocket_ids[:12]}")

    frozen_mol = Proj(copy.deepcopy(initial["mol_project"])).eval()
    frozen_pocket = Proj(copy.deepcopy(initial["pocket_project"])).eval()
    with torch.inference_mode():
        all_scores = (frozen_mol(molecule_rep) @ frozen_pocket(pocket_rep).T).numpy()
    official = all_scores[pair_mol, pair_target]

    # Strong ligand-only baseline under the identical global scaffold folds.
    smiles_by_index = list(map(str, archive["molecule_ids"]))
    fp = fp_matrix(smiles_by_index)
    ecfp = np.full(len(pairs), np.nan)
    for fold in sorted(np.unique(folds)):
        for target_index, target in enumerate(TARGETS):
            train = (folds != fold) & (pair_target == target_index)
            test = (folds == fold) & (pair_target == target_index)
            classifier = LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000,
                                            solver="liblinear", random_state=20260925)
            classifier.fit(fp[pair_mol[train]], labels[train])
            ecfp[test] = classifier.predict_proba(fp[pair_mol[test]])[:, 1]

    seed_values = list(map(int, args.seeds.split(",")))
    tuned_seeds, random_seeds = [], []
    prediction_rows = []
    for seed in seed_values:
        tuned = np.full(len(pairs), np.nan); random = np.full(len(pairs), np.nan)
        for fold in sorted(np.unique(folds)):
            train = folds != fold; test = folds == fold
            model = fit(pair_mol[train], pair_target[train], labels[train], molecule_rep, pocket_rep,
                        initial, seed + 101 * int(fold), args, False)
            random_model = fit(pair_mol[train], pair_target[train], labels[train], molecule_rep, pocket_rep,
                               initial, seed + 101 * int(fold), args, True)
            tuned[test] = predict(model, molecule_rep, pocket_rep, pair_mol[test], pair_target[test])
            random[test] = predict(random_model, molecule_rep, pocket_rep, pair_mol[test], pair_target[test])
        tuned_seeds.append(tuned); random_seeds.append(random)
        for i, row in pairs.iterrows():
            prediction_rows.append({"pair_id": row.pair_id, "seed": seed, "target": row.target,
                                    "label": int(row.label), "scaffold_fold": int(row.scaffold_fold),
                                    "official": float(official[i]), "ecfp": float(ecfp[i]),
                                    "balanced_bce": float(tuned[i]), "random_label": float(random[i])})
    tuned = np.mean(tuned_seeds, axis=0); random = np.mean(random_seeds, axis=0)
    report = {
        "method": "DrugCLIP dual-projection rank-4 LoRA with target/label-balanced BCE",
        "protocol": "Five-fold global Murcko-scaffold OOF; identical folds for every baseline.",
        "pairs_evaluated": int(len(pairs)), "input_pairs_excluded": int(len(pairs_raw) - len(pairs)),
        "metrics": {"official_drugclip": screening_metrics(pairs, official),
                    "ecfp4_logistic": screening_metrics(pairs, ecfp),
                    "random_label_lora": screening_metrics(pairs, random),
                    "balanced_bce_lora": screening_metrics(pairs, tuned)},
        "hyperparameters": {"rank": args.rank, "epochs": args.epochs, "lr": args.lr,
                            "temperature": args.temperature, "preserve_weight": args.preserve_weight,
                            "seeds": seed_values},
        "claim_boundary": "Retrospective active/decoy screening; not affinity, efficacy, or prospective validation.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(prediction_rows).to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
