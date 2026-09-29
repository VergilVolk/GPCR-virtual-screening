#!/usr/bin/env python3
"""Evaluate published static and GaMD-ensemble docking on a frozen subset.

Missing docking rows are assigned a score below the worst observed score for
that target/method.  This avoids the optimistic bias of complete-case-only
evaluation while preserving an explicit coverage audit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit.ML.Scoring.Scoring import CalcBEDROC
from sklearn.metrics import average_precision_score, roc_auc_score


TARGETS = ["B2AR", "CCR2", "M2R", "M4R"]
METHODS = {
    "PDB_Glide": ("PDB/{t}_PDB_Glide_scores.csv", "glide_gscore", "lower"),
    "PDB_Vina": ("PDB/{t}_PDB_Vina_scores.csv", "vina_score", "lower"),
    "Ensemble_Glide_BEavg": ("Ensemble/{t}_Ensemble_Glide_BEavg_ranked.csv", "BE_avg", "lower"),
    "Ensemble_Glide_BEmin": ("Ensemble/{t}_Ensemble_Glide_BEmin_ranked.csv", "BE_min", "lower"),
    "Ensemble_Vina_BEavg": ("Ensemble/{t}_Ensemble_Vina_BEavg_ranked.csv", "BE_avg", "lower"),
    "Ensemble_Vina_BEmin": ("Ensemble/{t}_Ensemble_Vina_BEmin_ranked.csv", "BE_min", "lower"),
}


def ef(labels: np.ndarray, scores: np.ndarray, fraction: float) -> float:
    n = max(1, int(np.ceil(len(labels) * fraction)))
    order = np.argsort(-scores, kind="mergesort")[:n]
    prevalence = labels.mean()
    return float(labels[order].mean() / prevalence) if prevalence else float("nan")


def metrics(labels: np.ndarray, scores: np.ndarray) -> dict:
    ranked = [[float(score), int(label)] for score, label in sorted(
        zip(scores, labels), key=lambda pair: pair[0], reverse=True)]
    return {
        "n": int(len(labels)), "actives": int(labels.sum()),
        "roc_auc": float(roc_auc_score(labels, scores)),
        "pr_auc": float(average_precision_score(labels, scores)),
        "bedroc_alpha20": float(CalcBEDROC(ranked, 1, 20.0)),
        "ef1pct": ef(labels, scores, 0.01),
        "ef5pct": ef(labels, scores, 0.05),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--scores-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pairs = pd.read_csv(args.pairs)
    report = {"protocol": "Frozen benchmark subset; missing docking scores rank below every observed score.",
              "methods": {}}
    long_rows = []
    for method, (pattern, score_column, direction) in METHODS.items():
        per_target, coverage = {}, {}
        for target in TARGETS:
            target_pairs = pairs[pairs.target == target].copy()
            path = args.scores_root / target / pattern.format(t=target)
            scores = pd.read_csv(path)
            if score_column not in scores.columns:
                raise ValueError(f"{path}: expected {score_column}; got {list(scores.columns)}")
            values = scores[["ligand_id", score_column]].drop_duplicates("ligand_id")
            merged = target_pairs.merge(values, left_on="source_id", right_on="ligand_id", how="left")
            observed = merged[score_column].notna().to_numpy()
            raw = merged[score_column].to_numpy(float)
            transformed = -raw if direction == "lower" else raw
            floor = float(np.nanmin(transformed) - max(1e-6, np.nanstd(transformed) * 1e-6))
            transformed[~observed] = floor
            labels = merged.label.to_numpy(int)
            per_target[target] = metrics(labels, transformed)
            coverage[target] = {
                "n": int(len(merged)), "observed": int(observed.sum()),
                "coverage": float(observed.mean()),
                "missing_actives": int(((labels == 1) & ~observed).sum()),
                "missing_decoys": int(((labels == 0) & ~observed).sum()),
            }
            for pair_id, label, score, is_observed in zip(merged.pair_id, labels, transformed, observed):
                long_rows.append({"pair_id": pair_id, "method": method, "label": int(label),
                                  "score": float(score), "observed": bool(is_observed)})
        metric_names = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]
        report["methods"][method] = {
            "macro": {name: float(np.mean([per_target[t][name] for t in TARGETS])) for name in metric_names},
            "per_target": per_target, "coverage": coverage,
        }
    report["claim_boundary"] = "Published retrospective docking baselines only; no prospective or PAM-efficacy claim."
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(long_rows).to_csv(args.output.with_suffix(".csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
