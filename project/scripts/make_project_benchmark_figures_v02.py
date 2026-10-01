#!/usr/bin/env python3
"""Publication-style v02 figures using only protocol-compatible comparisons."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

METRICS = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]
LABEL = {"roc_auc": "ROC-AUC", "pr_auc": "PR-AUC", "bedroc_alpha20": "BEDROC20", "ef1pct": "EF1%", "ef5pct": "EF5%"}


def save(fig, path):
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = p.parse_args()
    root = args.repo_root
    out = root / "project/results/project_wide_integration_benchmark_v02"
    figdir = out / "figures"; figdir.mkdir(parents=True, exist_ok=True)
    bind = json.loads((out / "binding_recomputed_v02.json").read_text())
    master = json.loads((root / "project/results/project_wide_integration_benchmark_v01/project_wide_benchmark_master_v01.json").read_text())
    registry = pd.read_csv(out / "benchmark_method_registry_v02.csv")
    candidates = pd.read_csv(out / "candidate_handoff_all200_v02.csv")
    audit = json.loads((out / "candidate_handoff_audit_v02.json").read_text())
    plt.rcParams.update({"font.size": 8, "axes.titleweight": "bold"})

    # 1. Fair bars: protocols are deliberately separate panels.
    methods = ["official_drugclip", "ecfp4_logistic", "ep80_score_ensemble", "famaug_score_ensemble"]
    names = ["Official DrugCLIP", "ECFP4 logistic", "ep80 ensemble", "family-aug ensemble"]
    colors = ["#777777", "#4C72B0", "#DD8452", "#55A868"]
    fig, axes = plt.subplots(2, 5, figsize=(15, 7.4), constrained_layout=True)
    for row, protocol in enumerate(("13T", "20T")):
        for col, metric_name in enumerate(METRICS):
            vals = [bind["methods"][f"{protocol}/{m}"]["macro"][metric_name] for m in methods]
            bars = axes[row, col].bar(range(4), vals, color=colors)
            axes[row, col].set_xticks(range(4), names, rotation=35, ha="right", fontsize=6)
            axes[row, col].set_title(f"{protocol} {LABEL[metric_name]}")
            for bar, value in zip(bars, vals):
                axes[row, col].text(bar.get_x()+bar.get_width()/2, value, f"{value:.3f}", ha="center", va="bottom", fontsize=6)
    fig.suptitle("Binding benchmark v02 — identical metrics, separate 13T and 20T protocols")
    save(fig, figdir / "fig_v02_01_binding_fair_bars.png")

    # 2. Delta + CI heatmap. Dot marks intervals excluding zero.
    comparisons = bind["paired_target_bootstrap"]
    mat = np.array([[comparisons[c][m]["delta"] for m in METRICS] for c in comparisons])
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    lim = max(abs(mat.min()), abs(mat.max()))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(5), [LABEL[m] for m in METRICS])
    ax.set_yticks(range(len(comparisons)), [x.replace("/", "\n") for x in comparisons])
    for i, comp in enumerate(comparisons):
        for j, metric_name in enumerate(METRICS):
            node = comparisons[comp][metric_name]; lo, hi = node["target_bootstrap_95ci"]
            mark = "*" if lo > 0 or hi < 0 else ""
            ax.text(j, i, f"{node['delta']:+.3f}{mark}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, label="candidate − baseline")
    ax.set_title("Paired target-bootstrap deltas (* 95% CI excludes zero)")
    save(fig, figdir / "fig_v02_02_binding_delta_ci_heatmap.png")

    # 3. 13T per-target heatmap.
    targets = sorted(bind["methods"]["13T/official_drugclip"]["per_target"])
    mat = np.array([[bind["methods"][f"13T/{m}"]["per_target"][t]["roc_auc"] for t in targets] for m in methods])
    fig, ax = plt.subplots(figsize=(10, 3.6))
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0.3, vmax=0.85, aspect="auto")
    ax.set_xticks(range(len(targets)), targets, rotation=45, ha="right")
    ax.set_yticks(range(4), names)
    for i in range(4):
        for j in range(len(targets)): ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, label="ROC-AUC")
    ax.set_title("13-target strict LOSO: per-target ROC-AUC")
    save(fig, figdir / "fig_v02_03_binding_target_heatmap.png")

    # 4. Functional three-class direction heatmap (same frozen endpoint).
    regions = ["compound110_extension", "cooperativity_mutagenesis", "intracellular_microswitches",
               "orthosteric_activation_core", "orthosteric_contact_union", "pam_contact_consensus", "pam_contact_union"]
    cm = {r["region"]: r for r in master["functional"]["cm00734_region_rows"]}
    fmat = np.array([[cm[r][name] for r in regions] for name in ("compound110_R1R3", "LY2119620_R1R3", "CM00734_R1R3")])
    fig, ax = plt.subplots(figsize=(9, 3.2))
    im = ax.imshow(fmat, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
    ax.set_xticks(range(len(regions)), [r.replace("_", "\n") for r in regions], fontsize=6)
    ax.set_yticks(range(3), ["compound110", "LY2119620", "CM00734 inactive"])
    for i in range(3):
        for j in range(len(regions)): ax.text(j, i, f"{fmat[i,j]:+.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, label="R1/R3 DeltaINT direction cosine")
    ax.set_title("Functional line: frozen STATE_MOTION / DeltaINT endpoint")
    save(fig, figdir / "fig_v02_04_functional_three_class_heatmap.png")

    # 5. Baseline coverage grid; status is not performance.
    families = ["2D/QSAR", "static structure", "binding", "MD representation", "functional MD"]
    status_order = ["RECOMPUTED_RAW", "COMMITTED_RESULT", "COMMITTED_NEGATIVE", "RUN_AS_VinaOnly",
                    "RUN_SUPERSEDED_BY_C1_BS256", "SOFTWARE_SMOKE_ONLY", "PILOT_ONLY_NOT_COMMON_PROTOCOL",
                    "NOT_RUN_NO_MATCHED_BINARY_GPU", "NOT_RUN_LEGACY_CUDA_STACK", "NOT_RUN_NO_COMMON_PREDICTIONS", "NOT_RUN_DATA_GATE"]
    methods_short = registry.method.tolist()
    grid = np.zeros((len(methods_short), len(families)))
    for i, row in registry.iterrows():
        j = families.index(row.module)
        if row.status in ("RECOMPUTED_RAW", "COMMITTED_RESULT", "RUN_AS_VinaOnly"):
            grid[i, j] = 2
        elif row.status in ("RUN_SUPERSEDED_BY_C1_BS256", "SOFTWARE_SMOKE_ONLY", "PILOT_ONLY_NOT_COMMON_PROTOCOL"):
            grid[i, j] = 1
        else:
            grid[i, j] = -1
    fig, ax = plt.subplots(figsize=(8.5, max(8, len(methods_short)*0.22)))
    from matplotlib.colors import ListedColormap, BoundaryNorm
    im = ax.imshow(grid, cmap=ListedColormap(["#d95f5f", "#eeeeee", "#e6ab55", "#55a868"]), norm=BoundaryNorm([-1.5,-0.5,0.5,1.5,2.5],4), aspect="auto")
    ax.set_xticks(range(len(families)), families, rotation=25, ha="right")
    ax.set_yticks(range(len(methods_short)), methods_short, fontsize=5.5)
    ax.set_title("Benchmark coverage grid: green=full, orange=partial/superseded, red=negative/not run")
    save(fig, figdir / "fig_v02_05_baseline_coverage_grid.png")

    # 6. Candidate handoff: this is the correct same-200 comparison.
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    dual = (candidates.drugclip_rank <= 50) & (candidates.new_rank <= 25)
    ax.scatter(candidates.drugclip_rank, candidates.new_rank, s=16, alpha=.5, color="#777")
    ax.scatter(candidates.loc[dual, "drugclip_rank"], candidates.loc[dual, "new_rank"], s=45, color="#C44E52", label=f"dual gate n={dual.sum()}")
    offsets = [(2, -2), (2, 2), (2, -4), (2, 4), (2, 0)]
    for ((_, row), (dx, dy)) in zip(candidates.loc[dual].sort_values("new_rank").iterrows(), offsets):
        ax.annotate(row.candidate_id, (row.drugclip_rank, row.new_rank), xytext=(dx, dy),
                    textcoords="offset points", fontsize=6, ha="left", va="center")
    ax.axvline(50, ls="--", color="#4C72B0"); ax.axhline(25, ls="--", color="#55A868")
    ax.invert_xaxis(); ax.invert_yaxis(); ax.set_xlabel("old DrugCLIP rank (lower better)"); ax.set_ylabel("family-augmented rank (lower better)")
    ax.set_title(f"PACER-200 handoff, Spearman={audit['spearman_old_drugclip_vs_famaug_rank']:+.3f}")
    ax.legend()
    save(fig, figdir / "fig_v02_06_candidate_handoff_scatter.png")
    print(json.dumps({"figures": 6, "output": str(figdir)}, indent=2))


if __name__ == "__main__":
    main()
