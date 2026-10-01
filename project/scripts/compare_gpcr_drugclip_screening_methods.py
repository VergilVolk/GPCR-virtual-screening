#!/usr/bin/env python3
"""Freeze common-pair comparisons and scaffold-bootstrap uncertainty."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from finetune_drugclip_gpcr_screening import TARGETS, screening_metrics


def valid_for_metrics(table):
    for target in TARGETS:
        labels = table.loc[table.target == target, "label"]
        if labels.nunique() < 2:
            return False
    return True


def bootstrap(table, scores, proposed, comparator, repeats, seed):
    scaffold_values = table.murcko_scaffold.astype(str).to_numpy()
    groups = np.unique(scaffold_values)
    group_indices = {group: np.flatnonzero(scaffold_values == group) for group in groups}
    rng = np.random.default_rng(seed); rows = []
    for _ in range(repeats):
        sampled = rng.choice(groups, len(groups), replace=True)
        idx = np.concatenate([group_indices[group] for group in sampled])
        boot = table.iloc[idx].reset_index(drop=True)
        if not valid_for_metrics(boot):
            continue
        a = screening_metrics(boot, scores[proposed][idx])["macro"]
        b = screening_metrics(boot, scores[comparator][idx])["macro"]
        rows.append({metric: a[metric] - b[metric] for metric in ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]})
    return {metric: {"q025": float(np.quantile([row[metric] for row in rows], 0.025)),
                     "median": float(np.median([row[metric] for row in rows])),
                     "q975": float(np.quantile([row[metric] for row in rows], 0.975)),
                     "p_one_sided_le_zero": float((1 + sum(row[metric] <= 0 for row in rows)) / (1 + len(rows)))}
            for metric in rows[0]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--single", type=Path, required=True)
    parser.add_argument("--ensemble", type=Path, required=True)
    parser.add_argument("--docking", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=500)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs)
    single = pd.read_csv(args.single).groupby("pair_id", as_index=False).mean(numeric_only=True)
    ensemble = pd.read_csv(args.ensemble).groupby("pair_id", as_index=False).mean(numeric_only=True)
    docking = pd.read_csv(args.docking)
    common = set(single.pair_id) & set(ensemble.pair_id)
    table = pairs[pairs.pair_id.isin(common)].copy().sort_values("pair_id").reset_index(drop=True)
    single = single.set_index("pair_id").loc[table.pair_id]
    ensemble = ensemble.set_index("pair_id").loc[table.pair_id]
    scores = {
        "official_cluster0": ensemble.official_cluster0.to_numpy(float),
        "official_population_mean": ensemble.official_population_mean.to_numpy(float),
        "official_max": ensemble.official_max.to_numpy(float),
        "official_population_lse": ensemble.official_population_lse.to_numpy(float),
        "finetuned_cluster0": ensemble.finetuned_cluster0.to_numpy(float),
        "finetuned_population_lse": ensemble.finetuned_population_lse.to_numpy(float),
        "random_population_lse": ensemble.random_population_lse.to_numpy(float),
        "ecfp4_logistic": single.ecfp.to_numpy(float),
        "single_balanced_bce": single.balanced_bce.to_numpy(float),
    }
    for method, frame in docking[docking.pair_id.isin(common)].groupby("method"):
        values = frame.drop_duplicates("pair_id").set_index("pair_id").reindex(table.pair_id)
        if values.score.isna().any():
            raise ValueError(f"Docking method {method} lacks common-pair rows")
        scores[f"docking_{method}"] = values.score.to_numpy(float)
    metrics = {name: screening_metrics(table, values) for name, values in scores.items()}
    proposed = "finetuned_population_lse"
    comparisons = {}
    for comparator in ["finetuned_cluster0", "official_population_lse", "random_population_lse",
                       "ecfp4_logistic", "docking_Ensemble_Glide_BEavg", "docking_Ensemble_Glide_BEmin"]:
        comparisons[f"{proposed}_minus_{comparator}"] = bootstrap(
            table, scores, proposed, comparator, args.bootstrap, 20260925 + len(comparisons))
    report = {
        "common_pairs": int(len(table)), "excluded_for_missing_ai_input": int(len(pairs) - len(table)),
        "metrics": metrics, "scaffold_bootstrap_differences": comparisons,
        "claim_boundary": "Retrospective common-pair active/decoy screening; no PAM-efficacy or prospective claim.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    summary = []
    for method, value in metrics.items():
        summary.append({"method": method, **value["macro"]})
    pd.DataFrame(summary).sort_values("bedroc_alpha20", ascending=False).to_csv(args.output.with_suffix(".csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
