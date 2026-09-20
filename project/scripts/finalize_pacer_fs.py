# -*- coding: utf-8 -*-
"""Freeze PACER-FS evidence tables and publication-style summary figure."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_fs_final"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260829)


def boot(delta, n=100_000):
    delta = np.asarray(delta, float)
    ix = RNG.integers(0, len(delta), size=(n, len(delta)))
    b = delta[ix].mean(1)
    return {"estimate": float(delta.mean()), "ci95": [float(x) for x in np.quantile(b, [.025, .975])],
            "probability_gt_zero": float((b > 0).mean())}


def main():
    series = pd.read_csv(P / "results" / "pacer_series_normalized_v01" / "outer_episodes.csv")
    dual = pd.read_csv(P / "results" / "pacer_dual_channel_v01" / "outer_episodes.csv")
    robust = pd.read_csv(P / "results" / "pacer_centering_robustness_v01" / "summary.csv")
    domain = pd.read_csv(P / "results" / "pacer_centered_domain_audit_v01" / "summary.csv")
    temporal = pd.read_csv(P / "results" / "pacer_temporal_fewshot_v01" / "metrics.csv")

    sw = series.pivot(index="group", columns="method", values=["Spearman", "MAE"])
    dw = dual.pivot(index="group", columns="method", values=["Spearman", "MAE"])
    evidence = {
        "status": "strongest_internal_series_generalization_method_not_external_SOTA",
        "primary_protocol": "11 leave-one-source-component-out series, n>=12, 3 label-free predicted-span anchors",
        "centered_vs_absolute_spearman": boot(sw["Spearman"]["Center"] - sw["Spearman"]["Absolute"]),
        "centered_vs_absolute_mae": boot(-(sw["MAE"]["Center"] - sw["MAE"]["Absolute"])),
        "dual_vs_centered_spearman": boot(dw["Spearman"]["DualChannel-Centered"] - dw["Spearman"]["CenteredCalibrated"]),
        "hyperparameter_robustness": {
            "positive_configs": int((robust.delta_macro > 0).sum()), "n_configs": len(robust),
            "median_delta": float(robust.delta_macro.median()),
        },
        "temporal_single_series": temporal.to_dict("records"),
        "limitations": [
            "The temporal result contains one future series and is exploratory.",
            "The exact method was developed on the same 11 LOSO series; independent external validation is absent.",
            "Potency data are almost entirely ACh/calcium and cannot validate probe dependence.",
            "Potency ranking does not establish that a generated molecule is a PAM.",
        ],
    }
    (OUT / "evidence.json").write_text(json.dumps(evidence, indent=2), encoding="utf-8")

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.2))
    methods = ["Absolute", "Center", "DualChannel-Centered"]
    vals = [sw["Spearman"]["Absolute"].mean(), sw["Spearman"]["Center"].mean(),
            dw["Spearman"]["DualChannel-Centered"].mean()]
    ax[0, 0].bar(["Absolute\nQSAR", "Series-centered\nPACER-FS", "Dual-channel\n(exploratory)"], vals,
                 color=["#9aa0a6", "#1b9e77", "#7570b3"])
    ax[0, 0].axhline(0, color="black", lw=.8); ax[0, 0].set_ylabel("Macro series Spearman")
    ax[0, 0].set_title("A  New-series ranking")
    for i, v in enumerate(vals): ax[0, 0].text(i, v + .008, f"{v:.3f}", ha="center")

    a = sw["Spearman"]["Absolute"]; c = sw["Spearman"]["Center"]
    ax[0, 1].scatter(a, c, color="#1b9e77", s=42)
    lo = min(a.min(), c.min()) - .04; hi = max(a.max(), c.max()) + .04
    ax[0, 1].plot([lo, hi], [lo, hi], "--", color="#777777", lw=1)
    for group in a.index: ax[0, 1].annotate(group[-2:], (a[group], c[group]), fontsize=7, xytext=(3, 2), textcoords="offset points")
    ax[0, 1].set(xlabel="Absolute QSAR Spearman", ylabel="Centered PACER-FS Spearman",
                 title="B  Paired held-out series", xlim=(lo, hi), ylim=(lo, hi))

    ax[1, 0].bar(np.arange(len(robust)), robust.delta_macro, color=np.where(robust.delta_macro > 0, "#1b9e77", "#d95f02"))
    ax[1, 0].axhline(0, color="black", lw=.8); ax[1, 0].set_xlabel("Fixed LightGBM configuration")
    ax[1, 0].set_ylabel("Centered - absolute macro Spearman")
    ax[1, 0].set_title("C  Hyperparameter-region robustness (11/12 positive)")

    labels = {"all_query": "All", "novel_murcko_scaffold": "Novel scaffold",
              "tanimoto_lt_0.70": "Tanimoto <0.70",
              "novel_scaffold_and_tanimoto_lt_0.70": "Both"}
    domain = domain.set_index("domain").loc[list(labels)]
    x = np.arange(len(domain)); w = .36
    ax[1, 1].bar(x - w/2, domain.absolute_macro, w, label="Absolute", color="#9aa0a6")
    ax[1, 1].bar(x + w/2, domain.centered_macro, w, label="Centered", color="#1b9e77")
    ax[1, 1].set_xticks(x, [labels[z] for z in domain.index], rotation=18, ha="right")
    ax[1, 1].set_ylabel("Macro series Spearman"); ax[1, 1].set_title("D  Chemical-domain stress test")
    ax[1, 1].legend(frameon=False)
    fig.suptitle("PACER-FS: factorizing series shift from within-series M4 PAM SAR", fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / "pacer_fs_summary.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT / "pacer_fs_summary.svg", bbox_inches="tight")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
