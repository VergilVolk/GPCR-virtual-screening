# -*- coding: utf-8 -*-
"""Publication figure for PACER-XR Cascade benchmark results."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import roc_curve


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_xr_v01"


def main():
    pred = pd.read_csv(OUT / "m4_blind_test_predictions.csv")
    met = pd.read_csv(OUT / "m4_blind_test_metrics.csv").set_index("method")
    methods = {
        "PACER-XR Cascade": ("PACER_XR_Cascade1pct", "#7a1fa2", 2.8),
        "PACER-XR": ("PACER_XR", "#d17c00", 2.2),
        "Glide BEmin": ("Glide_BEmin", "#276fbf", 1.8),
        "Vina BEmin": ("Vina_BEmin", "#2a9d5b", 1.8),
    }
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)

    y = pred.target.to_numpy()
    for label, (col, color, width) in methods.items():
        fpr, tpr, _ = roc_curve(y, pred[col])
        auc = met.loc[col, "ROC_AUC"]
        axes[0].plot(fpr, tpr, color=color, lw=width, label=f"{label} ({auc:.3f})")
    axes[0].plot([0, 1], [0, 1], "--", color="#888888", lw=1)
    axes[0].set(xlabel="False positive rate", ylabel="True positive rate", title="A  Global discrimination")
    axes[0].legend(frameon=False, loc="lower right")
    axes[0].grid(alpha=.18)

    cols = ["EF_0.5pct", "EF_1pct", "EF_2pct", "EF_5pct"]
    labels = ["0.5%", "1%", "2%", "5%"]
    x = range(len(cols)); width = .19
    for j, (label, (col, color, _)) in enumerate(methods.items()):
        axes[1].bar([i + (j - 1.5) * width for i in x], met.loc[col, cols], width,
                    color=color, label=label)
    axes[1].set_xticks(list(x), labels)
    axes[1].set(xlabel="Library sampled", ylabel="Enrichment factor", title="B  Early enrichment")
    axes[1].grid(axis="y", alpha=.18)
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("PACER-XR Cascade on the official M4 allosteric-modulator benchmark",
                 fontsize=14, fontweight="bold")
    fig.savefig(OUT / "PACER_XR_CASCADE_FIGURE.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "PACER_XR_CASCADE_FIGURE.svg", bbox_inches="tight")
    print(OUT / "PACER_XR_CASCADE_FIGURE.png")


if __name__ == "__main__":
    main()
