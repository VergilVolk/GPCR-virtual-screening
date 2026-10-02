#!/usr/bin/env python3
"""Unified retrospective bake-off for DrugCLIP generations and baselines.

This script deliberately distinguishes a *development/selection benchmark*
from a held-out lockbox.  It never calls the bundled four-GPCR benchmark
blind because its labels have already been used during model development.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp

from finetune_drugclip_gpcr_screening import screening_metrics


TARGETS = ["B2AR", "CCR2", "M2R", "M4R"]
METRICS = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_archive_scores(
    scores: dict[str, np.ndarray], prefix: str, archive: np.lib.npyio.NpzFile,
    pairs: pd.DataFrame, pocket_meta: pd.DataFrame, tau: float = 0.1,
) -> None:
    mol_index = {str(v): i for i, v in enumerate(archive["molecule_ids"])}
    pocket_ids = list(map(str, archive["pocket_ids"]))
    values = archive["scores"].astype(np.float64)
    pair_mol = np.asarray([mol_index[str(v)] for v in pairs.canonical_smiles], dtype=int)
    cluster = np.empty((len(pairs), 10), dtype=np.float64)
    weights = np.empty((len(pairs), 10), dtype=np.float64)
    for target in TARGETS:
        mask = pairs.target.to_numpy() == target
        ids = [f"{target}_cluster{i}" for i in range(10)]
        rows = [pocket_ids.index(v) for v in ids]
        cluster[mask] = values[np.ix_(rows, pair_mol[mask])].T
        p = pocket_meta.set_index("pocket_id").loc[ids, "population"].to_numpy(float).copy()
        p /= p.sum()
        weights[mask] = p
    scores[f"{prefix}_cluster0"] = cluster[:, 0]
    scores[f"{prefix}_population_mean"] = (cluster * weights).sum(axis=1)
    scores[f"{prefix}_population_lse"] = tau * logsumexp(
        np.log(np.maximum(weights, 1e-12)) + cluster / tau, axis=1
    )
    scores[f"{prefix}_max"] = cluster.max(axis=1)


def merge_prediction(
    table: pd.DataFrame, files: list[Path], score_column: str, output_name: str
) -> np.ndarray:
    keys = ["target", "canonical_smiles"]
    merged = table[keys].copy()
    values = []
    for i, path in enumerate(files):
        frame = pd.read_csv(path)
        if score_column not in frame:
            raise ValueError(f"{path} lacks score column {score_column}")
        # One score per target--molecule pair is required.
        frame = frame[keys + [score_column]].drop_duplicates(keys)
        col = f"seed_{i}"
        merged = merged.merge(frame.rename(columns={score_column: col}), on=keys, how="left")
        values.append(col)
    if merged[values].isna().any().any():
        missing = int(merged[values].isna().any(axis=1).sum())
        raise ValueError(f"{output_name}: {missing} benchmark rows lack predictions")
    return merged[values].mean(axis=1).to_numpy(float)


def method_metrics(table: pd.DataFrame, scores: dict[str, np.ndarray]) -> tuple[pd.DataFrame, pd.DataFrame]:
    macro_rows, target_rows = [], []
    for name, value in scores.items():
        result = screening_metrics(table, value)
        macro_rows.append({"method": name, **result["macro"]})
        for target, metrics in result["per_target"].items():
            target_rows.append({"method": name, "target": target, **metrics})
    return pd.DataFrame(macro_rows), pd.DataFrame(target_rows)


def target_percentile_rank(table: pd.DataFrame, values: np.ndarray) -> np.ndarray:
    frame = pd.DataFrame({"target": table.target.to_numpy(), "score": values})
    return frame.groupby("target", sort=False).score.rank(method="average", pct=True).to_numpy(float)


def matched_docking_metrics(
    table: pd.DataFrame, docking: pd.DataFrame, scores: dict[str, np.ndarray]
) -> pd.DataFrame:
    """Compare docking and AI only on rows observed for each docking method."""
    rows = []
    ai_methods = ["drugclip2023_raw_population_lse", "drugclip2023_gpcr_loto",
                  "drugclip2026_13t_famaug", "exploratory_rank_fusion_2023loto_2026famaug",
                  "frozen_m4_safe_router"]
    pair_index = pd.Series(np.arange(len(table)), index=table.pair_id)
    for docking_name, frame in docking[docking.observed.astype(bool)].groupby("method"):
        frame = frame.drop_duplicates("pair_id")
        frame = frame[frame.pair_id.isin(pair_index.index)].copy()
        idx = pair_index.loc[frame.pair_id].to_numpy(int)
        subset = table.iloc[idx].reset_index(drop=True)
        if set(subset.target) != set(TARGETS):
            continue
        methods = {f"docking_{docking_name}": frame.score.to_numpy(float)}
        methods.update({name: scores[name][idx] for name in ai_methods})
        for name, value in methods.items():
            try:
                result = screening_metrics(subset, value)["macro"]
            except ValueError:
                continue
            rows.append({"docking_cohort": docking_name, "n_pairs": len(subset),
                         "method": name, **result})
    return pd.DataFrame(rows)


def bootstrap_ci(
    table: pd.DataFrame, scores: dict[str, np.ndarray], reference: str,
    repeats: int, seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    target_values = table.target.to_numpy()
    scaffolds = table.murcko_scaffold.fillna("").astype(str).to_numpy()
    groups: dict[str, list[np.ndarray]] = {}
    for target in TARGETS:
        target_idx = np.flatnonzero(target_values == target)
        target_scaffolds = scaffolds[target_idx]
        groups[target] = [target_idx[target_scaffolds == s] for s in np.unique(target_scaffolds)]
    samples = {name: {metric: [] for metric in METRICS} for name in scores}
    deltas = {name: {metric: [] for metric in METRICS} for name in scores}
    for _ in range(repeats):
        idx_parts = []
        for target in TARGETS:
            target_groups = groups[target]
            selected = rng.integers(0, len(target_groups), size=len(target_groups))
            idx_parts.extend(target_groups[i] for i in selected)
        idx = np.concatenate(idx_parts)
        sample = table.iloc[idx].reset_index(drop=True)
        results = {}
        try:
            for name, value in scores.items():
                results[name] = screening_metrics(sample, value[idx])["macro"]
        except ValueError:
            continue
        for name in scores:
            for metric in METRICS:
                samples[name][metric].append(results[name][metric])
                deltas[name][metric].append(results[name][metric] - results[reference][metric])
    rows = []
    for name in scores:
        for metric in METRICS:
            lo, med, hi = np.quantile(samples[name][metric], [0.025, 0.5, 0.975])
            dlo, dmed, dhi = np.quantile(deltas[name][metric], [0.025, 0.5, 0.975])
            rows.append({"method": name, "metric": metric, "ci_low": lo, "median": med,
                         "ci_high": hi, "delta_reference": reference,
                         "delta_ci_low": dlo, "delta_median": dmed, "delta_ci_high": dhi})
    return pd.DataFrame(rows)


def plot_summary(macro: pd.DataFrame, per_target: pd.DataFrame, output_dir: Path) -> None:
    preferred = [
        "drugclip2023_raw_population_lse", "drugclip2026_raw_population_lse",
        "drugclip2023_gpcr_loto", "drugclip2026_gpcr_only_ep80",
        "drugclip2026_13t", "drugclip2026_13t_famaug", "drugclip2026_20t",
        "exploratory_rank_fusion_2023loto_2026famaug",
        "frozen_m4_safe_router",
        "ecfp4_target_loto", "iid_random_score", "random_label_control",
    ]
    chosen = macro[macro.method.isin(preferred)].copy()
    chosen["order"] = chosen.method.map({v: i for i, v in enumerate(preferred)})
    chosen = chosen.sort_values("order")
    display = {
        "drugclip2023_raw_population_lse": "2023 raw (ensemble LSE)",
        "drugclip2026_raw_population_lse": "2026 raw (ensemble LSE)",
        "drugclip2023_gpcr_loto": "2023 GPCR LOTO",
        "drugclip2026_gpcr_only_ep80": "2026 GPCR-only",
        "drugclip2026_13t": "2026 13-target",
        "drugclip2026_13t_famaug": "2026 13T family-aug",
        "drugclip2026_20t": "2026 20-target",
        "exploratory_rank_fusion_2023loto_2026famaug": "Fixed rank fusion",
        "frozen_m4_safe_router": "M4-safe routed fusion",
        "ecfp4_target_loto": "ECFP4 target-LOTO",
        "iid_random_score": "IID random",
        "random_label_control": "Random-label adapter",
    }
    labels = [display.get(v, v) for v in chosen.method]
    colors = ["#d73027" if v == "frozen_m4_safe_router" else
              "#fc8d59" if "fusion" in v else
              "#4575b4" if "drugclip" in v else
              "#1a9850" if "ecfp" in v else "#bdbdbd" for v in chosen.method]
    fig, axes = plt.subplots(1, 3, figsize=(17, 6))
    for ax, metric, title in zip(axes, ["roc_auc", "pr_auc", "bedroc_alpha20"],
                                  ["ROC-AUC", "PR-AUC", "BEDROC (alpha=20)"]):
        bars = ax.barh(np.arange(len(chosen)), chosen[metric], color=colors)
        ax.set_yticks(np.arange(len(chosen)), labels if ax is axes[0] else [""] * len(chosen))
        ax.invert_yaxis(); ax.set_title(title); ax.set_xlim(0, 1)
        for bar, value in zip(bars, chosen[metric]):
            ax.text(min(value + .01, .96), bar.get_y() + bar.get_height()/2, f"{value:.3f}", va="center", fontsize=8)
    fig.suptitle("Four-GPCR retrospective selection benchmark (not a blind lockbox)")
    fig.tight_layout(); fig.savefig(output_dir / "macro_metric_bars.png", dpi=200); plt.close(fig)

    heat_methods = [v for v in preferred if v in set(per_target.method)]
    pivot = per_target[per_target.method.isin(heat_methods)].pivot(index="method", columns="target", values="roc_auc").loc[heat_methods]
    fig, ax = plt.subplots(figsize=(8, max(4, .48 * len(pivot))))
    image = ax.imshow(pivot.to_numpy(), vmin=.35, vmax=1.0, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(pivot.columns)), pivot.columns)
    ax.set_yticks(range(len(pivot.index)), [v.replace("drugclip", "DC") for v in pivot.index])
    for i in range(len(pivot)):
        for j in range(len(pivot.columns)):
            ax.text(j, i, f"{pivot.iloc[i,j]:.3f}", ha="center", va="center",
                    color="white" if pivot.iloc[i,j] < .72 else "black", fontsize=8)
    fig.colorbar(image, ax=ax, label="ROC-AUC"); ax.set_title("Per-target ROC-AUC")
    fig.tight_layout(); fig.savefig(output_dir / "per_target_roc_heatmap.png", dpi=200); plt.close(fig)

    raw = macro.set_index("method")
    pairs = [("drugclip2023_raw_population_lse", "drugclip2026_raw_population_lse", "Raw"),
             ("drugclip2023_gpcr_loto", "drugclip2026_gpcr_only_ep80", "GPCR adapted")]
    fig, ax = plt.subplots(figsize=(8, 4))
    for y, (old, new, label) in enumerate(pairs):
        if old not in raw.index or new not in raw.index:
            continue
        x1, x2 = raw.loc[old, "roc_auc"], raw.loc[new, "roc_auc"]
        ax.plot([x1, x2], [y, y], color="#777", lw=2); ax.scatter(x1, y, s=80, label="2023" if y == 0 else None)
        ax.scatter(x2, y, s=80, label="2026" if y == 0 else None)
        ax.text(x1, y+.12, f"{x1:.3f}", ha="center"); ax.text(x2, y-.20, f"{x2:.3f}", ha="center")
    ax.set_yticks(range(len(pairs)), [p[2] for p in pairs]); ax.set_xlim(.45, .8); ax.set_xlabel("Macro ROC-AUC")
    ax.set_title("Checkpoint generation is not a universal ordering"); ax.legend(); fig.tight_layout()
    fig.savefig(output_dir / "generation_dumbbell.png", dpi=200); plt.close(fig)


def plot_delta_ci(ci: pd.DataFrame, output_dir: Path) -> None:
    methods = ["drugclip2026_13t_famaug", "exploratory_rank_fusion_2023loto_2026famaug",
               "frozen_m4_safe_router"]
    metrics = ["roc_auc", "pr_auc", "bedroc_alpha20"]
    names = {"drugclip2026_13t_famaug": "2026 13T family-aug",
             "exploratory_rank_fusion_2023loto_2026famaug": "Fixed rank fusion",
             "frozen_m4_safe_router": "M4-safe routed fusion"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for ax, metric in zip(axes, metrics):
        frame = ci[(ci.method.isin(methods)) & (ci.metric == metric)].set_index("method").loc[methods]
        y = np.arange(len(methods))
        center = frame.delta_median.to_numpy(float)
        lower = center - frame.delta_ci_low.to_numpy(float)
        upper = frame.delta_ci_high.to_numpy(float) - center
        ax.errorbar(center, y, xerr=[lower, upper], fmt="o", color="#2166ac", capsize=4)
        ax.axvline(0, color="#555", ls="--", lw=1)
        ax.set_title({"roc_auc": "Delta ROC-AUC", "pr_auc": "Delta PR-AUC",
                      "bedroc_alpha20": "Delta BEDROC20"}[metric])
        ax.set_yticks(y); ax.grid(axis="x", alpha=.2)
    axes[0].set_yticks(np.arange(len(methods)), [names[v] for v in methods])
    axes[0].invert_yaxis()
    fig.suptitle("Paired scaffold-bootstrap delta vs 2023 GPCR LOTO (95% CI)")
    fig.tight_layout(); fig.savefig(output_dir / "paired_delta_forest.png", dpi=200); plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--old-archive", type=Path, required=True)
    p.add_argument("--new-archive", type=Path, required=True)
    p.add_argument("--pairs", type=Path, required=True)
    p.add_argument("--pocket-metadata", type=Path, required=True)
    p.add_argument("--old-loto", type=Path, required=True)
    p.add_argument("--new-gpcr", type=Path, required=True)
    p.add_argument("--new-13t", type=Path, nargs=3, required=True)
    p.add_argument("--new-13t-famaug", type=Path, nargs=3, required=True)
    p.add_argument("--new-20t", type=Path, nargs=3, required=True)
    p.add_argument("--docking", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=300)
    p.add_argument("--seed", type=int, default=20261002)
    a = p.parse_args()

    old = np.load(a.old_archive, allow_pickle=False); new = np.load(a.new_archive, allow_pickle=False)
    if not np.array_equal(old["molecule_ids"], new["molecule_ids"]):
        raise ValueError("Old/new molecule order differs")
    if not np.array_equal(old["pocket_ids"], new["pocket_ids"]):
        raise ValueError("Old/new pocket order differs")
    mols = set(map(str, old["molecule_ids"]))
    pairs0 = pd.read_csv(a.pairs)
    pairs = pairs0[pairs0.canonical_smiles.astype(str).isin(mols)].copy().reset_index(drop=True)
    pocket_meta = pd.read_csv(a.pocket_metadata)
    scores: dict[str, np.ndarray] = {}
    add_archive_scores(scores, "drugclip2023_raw", old, pairs, pocket_meta)
    add_archive_scores(scores, "drugclip2026_raw", new, pairs, pocket_meta)

    old_loto = pd.read_csv(a.old_loto).set_index("pair_id").reindex(pairs.pair_id)
    if old_loto[["tuned", "ecfp", "random"]].isna().any().any():
        raise ValueError("Old LOTO predictions do not cover the common benchmark")
    scores["drugclip2023_gpcr_loto"] = old_loto.tuned.to_numpy(float)
    scores["ecfp4_target_loto"] = old_loto.ecfp.to_numpy(float)
    scores["random_label_control"] = old_loto.random.to_numpy(float)
    scores["drugclip2026_gpcr_only_ep80"] = merge_prediction(pairs, [a.new_gpcr], "tuned", "new GPCR-only")
    scores["drugclip2026_13t"] = merge_prediction(pairs, list(a.new_13t), "tuned", "new 13-target")
    scores["drugclip2026_13t_famaug"] = merge_prediction(pairs, list(a.new_13t_famaug), "tuned", "new 13-target family aug")
    scores["drugclip2026_20t"] = merge_prediction(pairs, list(a.new_20t), "tuned", "new 20-target")

    # Prespecified 50:50 rank fusion.  It is exploratory on this development
    # set and must be frozen before evaluation on a new lockbox.
    old_rank = target_percentile_rank(pairs, scores["drugclip2023_gpcr_loto"])
    fam_rank = target_percentile_rank(pairs, scores["drugclip2026_13t_famaug"])
    scores["exploratory_rank_fusion_2023loto_2026famaug"] = 0.5 * old_rank + 0.5 * fam_rank
    # The M4-safe route was fixed from the earlier M4 degradation audit:
    # retain the 2023 GPCR adapter for M4, use the complementary fusion elsewhere.
    scores["frozen_m4_safe_router"] = np.where(
        pairs.target.to_numpy() == "M4R",
        scores["drugclip2023_gpcr_loto"],
        scores["exploratory_rank_fusion_2023loto_2026famaug"],
    )
    scores["iid_random_score"] = np.random.default_rng(a.seed).normal(size=len(pairs))

    docking = pd.read_csv(a.docking)

    macro, per_target = method_metrics(pairs, scores)
    # The strongest deployable 2023 GPCR transfer model is the primary comparator.
    reference = "drugclip2023_gpcr_loto"
    ci = bootstrap_ci(pairs, scores, reference, a.bootstrap, a.seed)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    macro.sort_values("bedroc_alpha20", ascending=False).to_csv(a.output_dir / "macro_metrics.csv", index=False)
    per_target.to_csv(a.output_dir / "per_target_metrics.csv", index=False)
    ci.to_csv(a.output_dir / "scaffold_bootstrap_ci.csv", index=False)
    matched_docking_metrics(pairs, docking, scores).to_csv(
        a.output_dir / "matched_docking_comparison.csv", index=False
    )
    pred = pairs[["pair_id", "target", "label", "canonical_smiles", "murcko_scaffold"]].copy()
    for name, value in scores.items(): pred[name] = value
    pred.to_csv(a.output_dir / "predictions.csv", index=False)
    plot_summary(macro, per_target, a.output_dir)
    plot_delta_ci(ci, a.output_dir)
    report = {
        "benchmark_role": "retrospective_selection_benchmark_not_lockbox",
        "n_input_pairs": int(len(pairs0)), "n_common_pairs": int(len(pairs)),
        "targets": TARGETS, "bootstrap": {"unit": "target-stratified Murcko scaffold", "repeats": a.bootstrap, "seed": a.seed},
        "primary_metrics": ["bedroc_alpha20", "ef1pct", "ef5pct"],
        "secondary_metrics": ["pr_auc", "roc_auc"],
        "input_sha256": {str(path): sha256(path) for path in [a.old_archive, a.new_archive, a.pairs, a.pocket_metadata, a.old_loto, a.new_gpcr, *a.new_13t, *a.new_13t_famaug, *a.new_20t, a.docking]},
        "claim_boundary": "Same-input retrospective comparison. Labels and this benchmark were previously seen during development; results select deployment routes but cannot establish prospective superiority, affinity, or PAM efficacy.",
    }
    (a.output_dir / "audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(macro.sort_values("bedroc_alpha20", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
