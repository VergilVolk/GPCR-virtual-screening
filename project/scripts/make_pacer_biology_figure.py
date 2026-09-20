# -*- coding: utf-8 -*-
"""Create the PACER functional-method and biological-mechanism summary figure."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_biology_v01"
OUT.mkdir(parents=True, exist_ok=True)


def bootstrap_ci(values, seed=830, n=50000):
    values = np.asarray(values, float); rng = np.random.default_rng(seed)
    draws = values[rng.integers(0, len(values), size=(n, len(values)))].mean(1)
    return np.quantile(draws, [.025, .975])


def main():
    models = pd.read_csv(P / "results" / "pacer_meta_metric_v01" / "series_results.csv")
    residues = pd.read_csv(P / "results" / "m4_residue_potency_association_v01" / "residue_associations.csv")
    candidates = pd.read_csv(P / "results" / "pacer_candidate_7trq_mechanism_v01" / "candidate_mechanism_annotations.csv")

    order = ["Absolute", "AnchorResidual", "PACER-MetaMetric", "PACER-FS", "PACER-Hybrid"]
    labels = ["Absolute", "Anchor\nresidual", "MetaMetric", "PACER-FS", "Hybrid"]
    colors = ["#9aa0a6", "#74a9cf", "#41b6c4", "#1769aa", "#ef8a62"]
    means, lows, highs = [], [], []
    for method in order:
        vals = models.loc[models.method == method, "Spearman"].to_numpy()
        ci = bootstrap_ci(vals)
        means.append(vals.mean()); lows.append(vals.mean() - ci[0]); highs.append(ci[1] - vals.mean())

    key_features = ["7TRQ|res89|min_A", "7TRQ|res432|min_A", "7TRQ|res184|pairs_per_heavy"]
    key = residues.set_index("feature").loc[key_features].copy()
    key_labels = ["Y89 proximity", "D432 proximity", "Q184 contacts"]

    heat_cols = ["Y89_proximity_favorable_percentile", "Q184_contact_favorable_percentile",
                 "D432_proximity_favorable_percentile"]
    candidates = candidates.sort_values(
        ["mechanism_checks_above_median", "Y89_proximity_favorable_percentile"],
        ascending=[False, False],
    ).head(12)

    plt.rcParams.update({"font.size": 10, "axes.titleweight": "bold"})
    fig = plt.figure(figsize=(15, 5.4), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.05, 1.35])

    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(order))
    ax.bar(x, means, yerr=np.vstack([lows, highs]), capsize=4, color=colors, edgecolor="white")
    for i, method in enumerate(order):
        vals = models.loc[models.method == method, "Spearman"].to_numpy()
        ax.scatter(np.full(len(vals), i) + np.linspace(-.12, .12, len(vals)), vals,
                   s=14, color="#263238", alpha=.55, zorder=3)
    ax.axhline(0, color="#555", lw=.8); ax.axhline(.263, color="#1769aa", ls="--", lw=1)
    ax.set_xticks(x, labels); ax.set_ylabel("Series-held-out Spearman")
    ax.set_title("A  Functional few-shot benchmark")
    ax.text(3.05, .275, "frozen PACER-FS", color="#1769aa", fontsize=8)

    ax = fig.add_subplot(gs[0, 1])
    yy = np.arange(3)
    ax.barh(yy + .17, key.pooled_centered_Spearman, height=.30, color="#2b8cbe", label="series-centered")
    ax.barh(yy - .17, key.physchem_adjusted_Spearman, height=.30, color="#fdae61", label="+ physchem adjusted")
    ax.axvline(0, color="#555", lw=.8)
    ax.set_yticks(yy, key_labels); ax.invert_yaxis(); ax.set_xlim(-.23, .23)
    ax.set_xlabel("Spearman with pEC50")
    ax.set_title("B  7TRQ residue hypotheses")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.text(.02, .02, "Primary BH q=0.021 for all three;\nnone survives adjusted global FDR",
            transform=ax.transAxes, fontsize=8, va="bottom",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": .88, "pad": 2})

    ax = fig.add_subplot(gs[0, 2])
    matrix = candidates[heat_cols].to_numpy(float)
    im = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap="RdYlBu_r")
    ax.set_yticks(np.arange(len(candidates)), candidates.candidate_id)
    ax.set_xticks(np.arange(3), ["Y89\nproximity", "Q184\ncontacts", "D432\nproximity"])
    for i in range(matrix.shape[0]):
        for j in range(3):
            ax.text(j, i, f"{matrix[i,j]:.2f}", ha="center", va="center",
                    fontsize=7, color="black" if .2 < matrix[i,j] < .8 else "white")
    ax.set_title("C  Candidate geometry percentiles")
    fig.colorbar(im, ax=ax, fraction=.045, pad=.03, label="favorable percentile")

    fig.suptitle("PACER-M4: few-shot potency, residue coupling, and candidate falsifiability",
                 fontsize=14, fontweight="bold")
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"PACER_M4_BIOLOGY_SUMMARY.{ext}", dpi=240, bbox_inches="tight")
    print(OUT / "PACER_M4_BIOLOGY_SUMMARY.png")


if __name__ == "__main__":
    main()
