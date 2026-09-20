# -*- coding: utf-8 -*-
"""Summarize PACER-M4 domain methods from saved OOF predictions.

This comparison is intentionally prediction-file based so partial experiments
cannot erase previous methods. It also computes paired bootstrap delta-AUC
against ERM on the exact same molecules.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score, balanced_accuracy_score, brier_score_loss,
    matthews_corrcoef, roc_auc_score,
)

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "results" / "pacer_groupdro_v1"
METHODS = ["ERM", "GroupDRO", "HardNeg", "GroupDROHardNeg"]
SEED = 42


def ece(y: np.ndarray, p: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0, 1, bins + 1)
    value = 0.0
    for low, high in zip(edges[:-1], edges[1:]):
        mask = (p >= low) & (p < high if high < 1 else p <= high)
        if mask.any():
            value += mask.mean() * abs(p[mask].mean() - y[mask].mean())
    return float(value)


def metrics(y: np.ndarray, p: np.ndarray) -> dict:
    pred = p >= 0.5
    return {
        "ROC_AUC": float(roc_auc_score(y, p)),
        "PR_AUC": float(average_precision_score(y, p)),
        "MCC": float(matthews_corrcoef(y, pred)),
        "BalancedAcc": float(balanced_accuracy_score(y, pred)),
        "Brier": float(brier_score_loss(y, p)),
        "ECE10": ece(y, p),
        "positive_prediction_rate": float(pred.mean()),
    }


def paired_auc_delta(y: np.ndarray, candidate: np.ndarray, reference: np.ndarray,
                     n_boot: int = 2000) -> dict:
    rng = np.random.RandomState(SEED)
    values = []
    for _ in range(n_boot):
        idx = rng.randint(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        values.append(roc_auc_score(y[idx], candidate[idx]) - roc_auc_score(y[idx], reference[idx]))
    return {
        "estimate": float(roc_auc_score(y, candidate) - roc_auc_score(y, reference)),
        "ci95": [float(x) for x in np.percentile(values, [2.5, 97.5])],
        "bootstrap_probability_gt_zero": float(np.mean(np.asarray(values) > 0)),
    }


def main() -> None:
    result = {}
    for split in ["scaffold", "source", "series"]:
        frames = {}
        for method in METHODS:
            path = ROOT / f"{split}_{method}_oof.csv"
            if path.exists():
                frames[method] = pd.read_csv(path).sort_values("canonical_molecule_id").reset_index(drop=True)
        reference = frames["ERM"]
        y = reference["target"].to_numpy(int)
        ref_prob = reference["probability"].to_numpy(float)
        for method, frame in frames.items():
            if not np.array_equal(reference["canonical_molecule_id"], frame["canonical_molecule_id"]):
                raise ValueError(f"Molecule order mismatch: {split}/{method}")
            prob = frame["probability"].to_numpy(float)
            item = {"metrics": metrics(y, prob), "mean_seed_std": float(frame["seed_std"].mean())}
            if method != "ERM":
                item["paired_delta_auc_vs_ERM"] = paired_auc_delta(y, prob, ref_prob)
            result[f"{split}|{method}"] = item
    with open(ROOT / "domain_method_comparison.json", "w", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    for key, value in result.items():
        print(key, value)


if __name__ == "__main__":
    main()
