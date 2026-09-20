# -*- coding: utf-8 -*-
"""Source-domain generalization prototype for strict M4 PAM classification.

Compares identical neural networks trained with:
  - ERM: class-balanced empirical risk minimization
  - GroupDRO: dynamically upweights high-loss source components

All preprocessing is fitted inside each outer training fold. Test folds are the
frozen PACER-M4 scaffold/source folds. This is a prototype gate: GroupDRO is kept
only if it improves source-held-out performance without damaging scaffold-held-
out performance.
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from sklearn.metrics import (
    average_precision_score, balanced_accuracy_score, brier_score_loss,
    matthews_corrcoef, roc_auc_score,
)

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "pam_vs_inactive.csv"
OUT = PROJECT / "results" / "pacer_groupdro_v1"
TOP_FP_BITS = 512

DESC = [
    Descriptors.MolWt, Descriptors.MolLogP, Descriptors.NumHDonors,
    Descriptors.NumHAcceptors, Descriptors.TPSA,
    Descriptors.NumRotatableBonds, Descriptors.NumAromaticRings,
    Descriptors.FractionCSP3, Descriptors.HeavyAtomCount,
]


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(min(14, torch.get_num_threads()))


def raw_features(smiles: list[str]) -> tuple[np.ndarray, np.ndarray]:
    fps, descs = [], []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smi}")
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros(2048, dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fps.append(arr)
        descs.append([fn(mol) for fn in DESC])
    return np.asarray(fps, dtype=np.float32), np.asarray(descs, dtype=np.float32)


def preprocess(train: np.ndarray, test: np.ndarray, train_desc: np.ndarray,
               test_desc: np.ndarray) -> tuple[np.ndarray, np.ndarray, list[int]]:
    variance = train.var(axis=0)
    selected = np.argsort(variance)[-TOP_FP_BITS:]
    mean = train_desc.mean(axis=0)
    std = train_desc.std(axis=0)
    std[std < 1e-8] = 1.0
    train_x = np.hstack([train[:, selected], (train_desc - mean) / std]).astype(np.float32)
    test_x = np.hstack([test[:, selected], (test_desc - mean) / std]).astype(np.float32)
    return train_x, test_x, selected.tolist()


class Classifier(torch.nn.Module):
    def __init__(self, n_features: int):
        super().__init__()
        self.network = torch.nn.Sequential(
            torch.nn.Linear(n_features, 64),
            torch.nn.LayerNorm(64),
            torch.nn.GELU(),
            torch.nn.Dropout(0.20),
            torch.nn.Linear(64, 32),
            torch.nn.GELU(),
            torch.nn.Dropout(0.10),
            torch.nn.Linear(32, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x).squeeze(-1)


def class_weights(y: np.ndarray) -> torch.Tensor:
    counts = np.bincount(y.astype(int), minlength=2).astype(float)
    weights = np.asarray([len(y) / (2.0 * counts[int(label)]) for label in y], dtype=np.float32)
    return torch.from_numpy(weights)


def hard_negative_pairs(fp: np.ndarray, y: np.ndarray, k: int = 5) -> tuple[torch.Tensor, torch.Tensor]:
    """For each inactive, select the k most structurally similar PAMs."""
    negatives = np.where(y == 0)[0]
    positives = np.where(y == 1)[0]
    pairs_neg, pairs_pos = [], []
    pos_fp = fp[positives]
    for neg in negatives:
        intersection = np.logical_and(pos_fp > 0, fp[neg] > 0).sum(axis=1)
        union = np.logical_or(pos_fp > 0, fp[neg] > 0).sum(axis=1)
        similarity = intersection / np.maximum(union, 1)
        chosen = np.argsort(similarity)[-min(k, len(positives)):]
        pairs_neg.extend([int(neg)] * len(chosen))
        pairs_pos.extend(positives[chosen].astype(int).tolist())
    return torch.tensor(pairs_neg, dtype=torch.long), torch.tensor(pairs_pos, dtype=torch.long)


def train_model(x: np.ndarray, raw_fp: np.ndarray, y: np.ndarray, groups: np.ndarray,
                method: str, seed: int, epochs: int) -> Classifier:
    seed_everything(seed)
    xt = torch.from_numpy(x)
    yt = torch.from_numpy(y.astype(np.float32))
    sample_weights = class_weights(y)
    group_names = sorted(set(map(str, groups)))
    group_to_idx = {g: i for i, g in enumerate(group_names)}
    gt = torch.tensor([group_to_idx[str(g)] for g in groups], dtype=torch.long)
    q = torch.ones(len(group_names), dtype=torch.float32) / len(group_names)
    hard_neg, hard_pos = hard_negative_pairs(raw_fp, y)

    model = Classifier(x.shape[1])
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-3)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    for _ in range(epochs):
        model.train()
        optimizer.zero_grad()
        logits = model(xt)
        losses = torch.nn.functional.binary_cross_entropy_with_logits(logits, yt, reduction="none")
        losses = losses * sample_weights
        use_groupdro = method in {"GroupDRO", "GroupDROHardNeg"}
        use_hardneg = method in {"HardNeg", "GroupDROHardNeg"}
        if not use_groupdro:
            objective = losses.mean()
        else:
            group_losses = torch.stack([losses[gt == idx].mean() for idx in range(len(group_names))])
            with torch.no_grad():
                q *= torch.exp(0.05 * group_losses.detach())
                q /= q.sum()
            objective = torch.sum(q * group_losses)
        if use_hardneg and len(hard_neg):
            # PAM logit should exceed the similar inactive logit by a margin.
            rank_loss = torch.nn.functional.softplus(
                0.5 - (logits[hard_pos] - logits[hard_neg])
            ).mean()
            objective = objective + 0.5 * rank_loss
        objective.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step()
        scheduler.step()
    return model


def expected_calibration_error(y: np.ndarray, prob: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0, 1, bins + 1)
    ece = 0.0
    for low, high in zip(edges[:-1], edges[1:]):
        mask = (prob >= low) & (prob < high if high < 1 else prob <= high)
        if mask.any():
            ece += mask.mean() * abs(prob[mask].mean() - y[mask].mean())
    return float(ece)


def metrics(y: np.ndarray, prob: np.ndarray) -> dict[str, float]:
    pred = (prob >= 0.5).astype(int)
    return {
        "ROC_AUC": float(roc_auc_score(y, prob)),
        "PR_AUC": float(average_precision_score(y, prob)),
        "MCC": float(matthews_corrcoef(y, pred)),
        "BalancedAcc": float(balanced_accuracy_score(y, pred)),
        "Brier": float(brier_score_loss(y, prob)),
        "ECE10": expected_calibration_error(y, prob),
        "positive_prediction_rate": float(pred.mean()),
    }


def run(split: str, method: str, df: pd.DataFrame, fps: np.ndarray,
        descs: np.ndarray, seeds: list[int], epochs: int) -> tuple[dict, pd.DataFrame]:
    fold_col = "series_holdout_fold" if split == "series" else f"{split}_fold"
    evaluation_mask = df[fold_col].to_numpy() >= 0
    seed_predictions = []
    fold_reports = []
    for seed in seeds:
        oof = np.full(len(df), np.nan, dtype=float)
        for fold in sorted(x for x in df[fold_col].unique() if x >= 0):
            test = df[fold_col].to_numpy() == fold
            train = ~test
            train_x, test_x, _ = preprocess(fps[train], fps[test], descs[train], descs[test])
            model = train_model(
                train_x, fps[train], df.loc[train, "target"].to_numpy(int),
                df.loc[train, "source_component"].to_numpy(str), method, seed, epochs,
            )
            model.eval()
            with torch.no_grad():
                prob = torch.sigmoid(model(torch.from_numpy(test_x))).numpy()
            oof[test] = prob
        if np.isnan(oof[evaluation_mask]).any():
            raise RuntimeError(f"Missing OOF predictions for {split}/{method}/{seed}")
        seed_predictions.append(oof[evaluation_mask])

    matrix = np.vstack(seed_predictions)
    mean_prob = matrix.mean(axis=0)
    eval_df = df[evaluation_mask].reset_index(drop=True)
    y = eval_df["target"].to_numpy(int)
    for fold in sorted(eval_df[fold_col].unique()):
        mask = eval_df[fold_col].to_numpy() == fold
        fold_reports.append({"fold": int(fold), "n": int(mask.sum()), **metrics(y[mask], mean_prob[mask])})

    group_reports = []
    for group, idx in eval_df.groupby("source_component").groups.items():
        idx = np.asarray(list(idx), dtype=int)
        item = {"source_component": group, "n": int(len(idx)),
                "positive_rate": float(y[idx].mean()),
                "predicted_positive_rate": float((mean_prob[idx] >= 0.5).mean()),
                "Brier": float(brier_score_loss(y[idx], mean_prob[idx]))}
        if len(np.unique(y[idx])) == 2:
            item["ROC_AUC"] = float(roc_auc_score(y[idx], mean_prob[idx]))
            item["BalancedAcc"] = float(balanced_accuracy_score(y[idx], mean_prob[idx] >= 0.5))
        group_reports.append(item)

    report = {
        "aggregate_oof": metrics(y, mean_prob),
        "seeds": seeds,
        "epochs": epochs,
        "folds": fold_reports,
        "source_components": group_reports,
    }
    output = eval_df[["canonical_molecule_id", "canonical_smiles", "target", fold_col,
                 "source_component"]].copy()
    output["probability"] = mean_prob
    output["seed_std"] = matrix.std(axis=0)
    return report, output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=120)
    parser.add_argument("--seeds", default="42", help="Comma-separated integer seeds")
    parser.add_argument("--methods", default="ERM,GroupDRO,HardNeg,GroupDROHardNeg")
    parser.add_argument("--splits", default="scaffold,source")
    args = parser.parse_args()
    seeds = [int(x) for x in args.seeds.split(",")]
    methods = [x.strip() for x in args.methods.split(",") if x.strip()]
    splits = [x.strip() for x in args.splits.split(",") if x.strip()]
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA).reset_index(drop=True)
    fps, descs = raw_features(df["canonical_smiles"].tolist())
    metrics_path = OUT / "groupdro_metrics.json"
    if metrics_path.exists():
        with open(metrics_path, encoding="utf-8") as handle:
            results = json.load(handle)
    else:
        results = {}
    for split in splits:
        for method in methods:
            report, predictions = run(split, method, df, fps, descs, seeds, args.epochs)
            results[f"{split}|{method}"] = report
            predictions.to_csv(OUT / f"{split}_{method}_oof.csv", index=False)
            print(split, method, report["aggregate_oof"], flush=True)
    with open(metrics_path, "w", encoding="utf-8") as handle:
        json.dump(results, handle, ensure_ascii=False, indent=2)
    print(f"Saved source-domain generalization results to {OUT}")


if __name__ == "__main__":
    main()
