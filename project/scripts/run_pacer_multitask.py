# -*- coding: utf-8 -*-
"""PACER-M4 v0.1: joint functional classification, potency and delta-SAR learning.

Tasks share one molecular encoder:
  - strict PAM vs experimentally inactive classification (529 molecules)
  - exact calcium/ACh pEC50 regression (430 molecules, masked loss)
  - optional local delta-pEC50 loss on near-neighbor pairs

Every outer split is inherited from M4-PAM-Benchmark v1.0. Regression labels
and pair labels are used only when their molecules are in the outer training
partition.
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
    matthews_corrcoef, mean_absolute_error, mean_squared_error, r2_score,
    roc_auc_score,
)

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(__file__).resolve().parents[1]
BENCH = PROJECT / "data" / "benchmarks" / "m4_pam_v1"
OUT = PROJECT / "results" / "pacer_multitask_v01"
TOP_FP_BITS = 512

DESC = [
    Descriptors.MolWt, Descriptors.MolLogP, Descriptors.NumHDonors,
    Descriptors.NumHAcceptors, Descriptors.TPSA,
    Descriptors.NumRotatableBonds, Descriptors.NumAromaticRings,
    Descriptors.FractionCSP3, Descriptors.HeavyAtomCount,
]


def set_seed(seed: int) -> None:
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    torch.set_num_threads(6)


def raw_features(smiles: list[str]) -> tuple[np.ndarray, np.ndarray]:
    fps, desc = [], []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros(2048, dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fps.append(arr); desc.append([fn(mol) for fn in DESC])
    return np.asarray(fps, np.float32), np.asarray(desc, np.float32)


def preprocess(train_fp, test_fp, train_desc, test_desc):
    selected = np.argsort(train_fp.var(axis=0))[-TOP_FP_BITS:]
    mean, std = train_desc.mean(0), train_desc.std(0)
    std[std < 1e-8] = 1.0
    train = np.hstack([train_fp[:, selected], (train_desc - mean) / std]).astype(np.float32)
    test = np.hstack([test_fp[:, selected], (test_desc - mean) / std]).astype(np.float32)
    return train, test


class PACERMultiTask(torch.nn.Module):
    def __init__(self, n_features: int):
        super().__init__()
        self.encoder = torch.nn.Sequential(
            torch.nn.Linear(n_features, 96), torch.nn.LayerNorm(96),
            torch.nn.GELU(), torch.nn.Dropout(0.20),
            torch.nn.Linear(96, 48), torch.nn.GELU(),
        )
        self.classifier = torch.nn.Linear(48, 1)
        self.regressor = torch.nn.Linear(48, 1)

    def forward(self, x):
        z = self.encoder(x)
        return self.classifier(z).squeeze(-1), self.regressor(z).squeeze(-1)


def class_sample_weights(y: np.ndarray) -> torch.Tensor:
    counts = np.bincount(y.astype(int), minlength=2)
    return torch.tensor([len(y) / (2 * counts[int(v)]) for v in y], dtype=torch.float32)


def train_model(x: np.ndarray, y_cls: np.ndarray, y_reg: np.ndarray,
                train_ids: list[str], pair_table: pd.DataFrame, method: str,
                seed: int, epochs: int, reg_weight: float = 0.5,
                pair_weight_value: float = 0.25) -> tuple[PACERMultiTask, float, float]:
    set_seed(seed)
    xt = torch.from_numpy(x)
    cls_t = torch.from_numpy(y_cls.astype(np.float32))
    cls_w = class_sample_weights(y_cls)
    reg_mask_np = np.isfinite(y_reg)
    reg_mean = float(np.nanmean(y_reg)); reg_std = float(np.nanstd(y_reg))
    reg_std = max(reg_std, 1e-6)
    reg_scaled = np.nan_to_num((y_reg - reg_mean) / reg_std, nan=0.0).astype(np.float32)
    reg_t = torch.from_numpy(reg_scaled)
    reg_mask = torch.from_numpy(reg_mask_np)

    id_to_local = {mid: idx for idx, mid in enumerate(train_ids)}
    valid_pairs = pair_table[
        pair_table["mol_a"].isin(id_to_local) & pair_table["mol_b"].isin(id_to_local)
    ]
    pair_a = torch.tensor([id_to_local[x] for x in valid_pairs["mol_a"]], dtype=torch.long)
    pair_b = torch.tensor([id_to_local[x] for x in valid_pairs["mol_b"]], dtype=torch.long)
    pair_delta = torch.tensor(
        valid_pairs["delta_pEC50_b_minus_a"].to_numpy(np.float32) / reg_std,
        dtype=torch.float32,
    )
    pair_weight = torch.tensor(
        np.where(valid_pairs["is_cliff"].to_numpy(int) == 1, 2.0, 1.0), dtype=torch.float32
    )

    model = PACERMultiTask(x.shape[1])
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-3)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)
    for _ in range(epochs):
        model.train(); optimizer.zero_grad()
        logits, reg = model(xt)
        cls_loss = (
            torch.nn.functional.binary_cross_entropy_with_logits(logits, cls_t, reduction="none")
            * cls_w
        ).mean()
        objective = cls_loss
        if method in {"MultiTask", "MultiTaskPair"}:
            reg_loss = torch.nn.functional.smooth_l1_loss(reg[reg_mask], reg_t[reg_mask])
            objective = objective + reg_weight * reg_loss
        if method == "MultiTaskPair" and len(pair_a):
            predicted_delta = reg[pair_b] - reg[pair_a]
            pair_loss = (
                torch.nn.functional.smooth_l1_loss(predicted_delta, pair_delta, reduction="none")
                * pair_weight
            ).mean()
            objective = objective + pair_weight_value * pair_loss
        objective.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step(); scheduler.step()
    return model, reg_mean, reg_std


def ece(y, p, bins=10):
    edges = np.linspace(0, 1, bins + 1); result = 0.0
    for low, high in zip(edges[:-1], edges[1:]):
        mask = (p >= low) & (p < high if high < 1 else p <= high)
        if mask.any(): result += mask.mean() * abs(p[mask].mean() - y[mask].mean())
    return float(result)


def cls_metrics(y, p):
    pred = p >= 0.5
    return {
        "ROC_AUC": float(roc_auc_score(y, p)),
        "PR_AUC": float(average_precision_score(y, p)),
        "MCC": float(matthews_corrcoef(y, pred)),
        "BalancedAcc": float(balanced_accuracy_score(y, pred)),
        "Brier": float(brier_score_loss(y, p)), "ECE10": ece(y, p),
        "positive_prediction_rate": float(pred.mean()),
    }


def reg_metrics(y, p):
    return {
        "n": int(len(y)), "MAE": float(mean_absolute_error(y, p)),
        "RMSE": float(np.sqrt(mean_squared_error(y, p))),
        "R2": float(r2_score(y, p)),
        "Spearman": float(pd.Series(y).corr(pd.Series(p), method="spearman")),
    }


def run(split, method, df, fps, desc, pair_table, seeds, epochs):
    fold_col = "series_holdout_fold" if split == "series" else f"{split}_fold"
    eval_mask = df[fold_col].to_numpy() >= 0
    all_cls, all_reg = [], []
    for seed in seeds:
        cls_oof = np.full(len(df), np.nan); reg_oof = np.full(len(df), np.nan)
        for fold in sorted(x for x in df[fold_col].unique() if x >= 0):
            test = df[fold_col].to_numpy() == fold; train = ~test
            train_x, test_x = preprocess(fps[train], fps[test], desc[train], desc[test])
            model, mean, std = train_model(
                train_x, df.loc[train, "target"].to_numpy(int),
                df.loc[train, "pEC50"].to_numpy(float),
                df.loc[train, "canonical_molecule_id"].tolist(), pair_table,
                method, seed, epochs,
            )
            model.eval()
            with torch.no_grad():
                logits, reg = model(torch.from_numpy(test_x))
                cls_oof[test] = torch.sigmoid(logits).numpy()
                reg_oof[test] = reg.numpy() * std + mean
        all_cls.append(cls_oof[eval_mask]); all_reg.append(reg_oof[eval_mask])
    eval_df = df[eval_mask].reset_index(drop=True)
    cls_pred = np.vstack(all_cls).mean(0); reg_pred = np.vstack(all_reg).mean(0)
    reg_mask = eval_df["pEC50"].notna().to_numpy()
    report = {"classification": cls_metrics(eval_df["target"].to_numpy(int), cls_pred)}
    if reg_mask.any():
        report["potency"] = reg_metrics(eval_df.loc[reg_mask, "pEC50"].to_numpy(float), reg_pred[reg_mask])
    output = eval_df[["canonical_molecule_id", "canonical_smiles", "target", "pEC50", fold_col]].copy()
    output["pam_probability"] = cls_pred; output["pEC50_prediction"] = reg_pred
    return report, output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=120)
    parser.add_argument("--seeds", default="42")
    parser.add_argument("--splits", default="scaffold,source,series")
    parser.add_argument("--methods", default="ClassifierOnly,MultiTask,MultiTaskPair")
    args = parser.parse_args()
    seeds = [int(x) for x in args.seeds.split(",")]
    splits = [x for x in args.splits.split(",") if x]
    methods = [x for x in args.methods.split(",") if x]
    OUT.mkdir(parents=True, exist_ok=True)

    strict = pd.read_csv(BENCH / "pam_vs_inactive.csv")
    potency = pd.read_csv(BENCH / "potency_molecules.csv")[["canonical_molecule_id", "pEC50"]]
    df = strict.merge(potency, on="canonical_molecule_id", how="left").reset_index(drop=True)
    pairs = pd.read_csv(BENCH / "activity_cliffs" / "near_neighbor_pairs.csv")
    fps, desc = raw_features(df["canonical_smiles"].tolist())
    metrics_path = OUT / "multitask_metrics.json"
    results = json.load(open(metrics_path, encoding="utf-8")) if metrics_path.exists() else {}
    for split in splits:
        for method in methods:
            report, predictions = run(split, method, df, fps, desc, pairs, seeds, args.epochs)
            results[f"{split}|{method}"] = {"seeds": seeds, "epochs": args.epochs, **report}
            predictions.to_csv(OUT / f"{split}_{method}_oof.csv", index=False)
            print(split, method, report, flush=True)
    with open(metrics_path, "w", encoding="utf-8") as handle:
        json.dump(results, handle, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
