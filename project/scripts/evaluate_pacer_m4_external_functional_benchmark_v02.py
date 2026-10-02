#!/usr/bin/env python3
"""Evaluate frozen scores on PACER-M4 external functional benchmark v02."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score


METHODS = [
    "drugclip_official", "drugclip_gpcr_triplet", "drugclip_m4_matched",
    "ecfp_rf_task", "max_train_similarity",
]


def concordance(y: np.ndarray, score: np.ndarray) -> float:
    dy = y[:, None] - y[None, :]
    ds = score[:, None] - score[None, :]
    mask = np.triu(dy != 0, 1)
    if not mask.any():
        return float("nan")
    product = dy[mask] * ds[mask]
    return float(np.mean(product > 0) + 0.5 * np.mean(product == 0))


def primary_metric(y: np.ndarray, score: np.ndarray, task: str) -> tuple[str, float, float]:
    if task == "binary":
        value = float(roc_auc_score(y, score))
        return "roc_auc", value, 2.0 * value - 1.0
    value = concordance(y, score)
    return "pairwise_concordance", value, 2.0 * value - 1.0


def bootstrap(y: np.ndarray, score: np.ndarray, task: str, seed: int,
              n_boot: int) -> list[float]:
    rng = np.random.default_rng(seed)
    values = []
    pos = np.flatnonzero(y == 1) if task == "binary" else None
    neg = np.flatnonzero(y == 0) if task == "binary" else None
    for _ in range(n_boot):
        if task == "binary":
            if len(pos) == 0 or len(neg) == 0:
                continue
            idx = np.concatenate([rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True)])
        else:
            idx = rng.integers(0, len(y), len(y))
        _, value, _ = primary_metric(y[idx], score[idx], task)
        if np.isfinite(value):
            values.append(value)
    return list(map(float, np.quantile(values, [0.025, 0.5, 0.975]))) if values else [float("nan")] * 3


def endpoint_result(frame: pd.DataFrame, method: str, n_boot: int, seed: int) -> dict:
    y = frame.value.to_numpy(float)
    score = frame[method].to_numpy(float)
    task = str(frame.task_type.iloc[0])
    metric, value, skill = primary_metric(y, score, task)
    out = {
        metric: value, "directional_skill": skill,
        f"{metric}_bootstrap_95ci": bootstrap(y, score, task, seed, n_boot),
    }
    if task == "binary":
        out["average_precision"] = float(average_precision_score(y, score))
    else:
        out["spearman"] = float(spearmanr(y, score).statistic)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--frozen-predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=3000)
    args = parser.parse_args()

    records = pd.read_csv(args.records)
    endpoint_summary = pd.read_csv(args.summary)
    predictions = pd.read_csv(args.frozen_predictions)
    score_cols = ["drugclip_official", "drugclip_gpcr_triplet", "drugclip_m4_matched",
                  "ecfp_rf_pam", "ecfp_rf_potency"]
    frozen = predictions.groupby("external_molecule_id", as_index=False)[score_cols].mean()
    table = records.merge(frozen, on="external_molecule_id", how="left", validate="many_to_one")
    if table[score_cols].isna().any().any():
        missing = table.loc[table[score_cols].isna().any(axis=1), "external_molecule_id"].unique()
        raise ValueError(f"Missing frozen scores for {len(missing)} molecules: {missing[:5]}")
    table["max_train_similarity"] = table.max_train_tanimoto_ecfp4
    table["ecfp_rf_task"] = np.where(table.task_type.eq("binary"),
                                      table.ecfp_rf_pam, table.ecfp_rf_potency)
    table = table[table.eligible_zero_shot.astype(bool)].copy()
    table = table.sort_values(["dataset", "endpoint", "external_molecule_id"]).drop_duplicates(
        ["dataset", "endpoint", "external_molecule_id"], keep="first")

    endpoint_reports = []
    skill_rows = []
    for i, ((dataset, endpoint), frame) in enumerate(table.groupby(["dataset", "endpoint"], sort=False)):
        meta = endpoint_summary[(endpoint_summary.dataset == dataset) &
                                (endpoint_summary.endpoint == endpoint)].iloc[0]
        y = frame.value.to_numpy(float)
        task = str(frame.task_type.iloc[0])
        report = {
            "dataset": dataset, "endpoint": endpoint, "task_type": task,
            "endpoint_role": str(frame.endpoint_role.iloc[0]), "n": int(len(frame)),
            "n_scaffolds": int(frame.murcko_scaffold.nunique()),
            "confirmatory_macro_eligible": bool(meta.confirmatory_macro_eligible),
            "methods": {},
        }
        if task == "binary":
            report["positive"] = int((y == 1).sum())
            report["negative"] = int((y == 0).sum())
            if min(report["positive"], report["negative"]) < 10:
                report["evidence_warning"] = "descriptive_only_fewer_than_10_per_class"
        for j, method in enumerate(METHODS):
            result = endpoint_result(frame, method, args.bootstrap, 20261002 + i * 17 + j)
            report["methods"][method] = result
            if report["confirmatory_macro_eligible"]:
                skill_rows.append({"dataset": dataset, "endpoint": endpoint,
                                   "method": method, "directional_skill": result["directional_skill"]})
        endpoint_reports.append(report)

    skills = pd.DataFrame(skill_rows)
    macro = {}
    for method, frame in skills.groupby("method"):
        vals = frame.directional_skill.to_numpy(float)
        macro[method] = {
            "macro_directional_skill": float(np.mean(vals)),
            "n_endpoints": int(len(vals)),
            "endpoint_min": float(np.min(vals)),
            "endpoint_max": float(np.max(vals)),
        }
    report = {
        "benchmark": "PACER-M4 external functional benchmark v02",
        "n_unique_molecules": int(table.external_molecule_id.nunique()),
        "n_endpoint_units": int(len(table)),
        "n_endpoints": int(table.endpoint.nunique()),
        "confirmatory_macro_definition": (
            "Unweighted mean of endpoint directional skill. Binary skill=2*AUC-1; "
            "continuous/ordinal skill=2*pairwise-concordance-1. Only primary functional "
            "endpoints with n>=15 and adequate class balance enter the macro score."
        ),
        "macro": macro,
        "endpoints": endpoint_reports,
        "claim_boundary": (
            "Retrospective multi-source functional transfer benchmark. It tests ranking and "
            "classification of known public compounds; it does not validate generated candidates "
            "or establish prospective/wet-lab PAM activity."
        ),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "benchmark_result.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    table.to_csv(args.output_dir / "frozen_predictions_by_endpoint.csv", index=False)
    pd.DataFrame(skill_rows).to_csv(args.output_dir / "confirmatory_endpoint_skills.csv", index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
