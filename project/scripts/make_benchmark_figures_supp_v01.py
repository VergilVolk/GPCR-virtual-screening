#!/usr/bin/env python
"""Supplementary figures 12-13: bootstrap CIs and LY-vs-c110 direct agreement."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FUNC = Path(r"D:\CLC_geom2vec_pilot\project\results")
OUT = FUNC / "project_wide_integration_benchmark_v01"
FIG = OUT / "figures"
M = json.loads((OUT / "project_wide_benchmark_master_v01.json").read_text(encoding="utf-8"))

plt.rcParams.update({"figure.dpi": 160, "font.size": 9, "axes.titlesize": 10,
                     "axes.titleweight": "bold", "savefig.bbox": "tight"})

# fig12: delta bootstrap CIs (ep40 vs official) + headline deltas
s = json.loads(Path(r"D:\CLC\project\results\drugclip_science2026\migration_13target_v01\summary_preserve0_3seed.json").read_text())
ci = s["bootstrap_95ci"]["tuned_minus_official"]
e20 = json.loads(Path(r"D:\CLC\project\results\drugclip_science2026\extended20_finetune_v01\ecfp20_loso.json").read_text())

fig, ax = plt.subplots(figsize=(8, 3.6))
items = [("ROC-AUC", *ci["roc_auc"]), ("PR-AUC", *ci["pr_auc"]),
         ("BEDROC20", *ci["bedroc_alpha20"]), ("EF1%", *ci["ef1pct"])]
y = np.arange(len(items))
for i, (name, lo, mid, hi) in enumerate(items):
    color = "#55A868" if lo > 0 else ("#C44E52" if hi < 0 else "#DD8452")
    ax.plot([lo, hi], [i, i], color=color, lw=3, solid_capstyle="round")
    ax.plot(mid, i, "o", color=color, ms=7)
    ax.text(hi + 0.004, i, f"[{lo:+.3f}, {hi:+.3f}]", va="center", fontsize=7.5)
ax.axvline(0, color="k", lw=1)
ax.set_yticks(y, [i[0] for i in items])
ax.set_xlabel("finetuned(ep40 ens) - official raw, macro delta, target-stratified scaffold bootstrap 95% CI")
ax.set_title("Binding line significance: what the CI does and does not support\n(green=positive CI, red=negative CI, orange=crosses zero)")
fig.text(0.01, -0.08, "source: summary_preserve0_3seed.json bootstrap_95ci | EF1% negative CI is the known boundary", fontsize=6.5, color="gray")
fig.savefig(FIG / "fig12_binding_delta_ci.png"); plt.close(fig)
print("saved fig12")

# fig13: LY vs compound110 direct direction agreement (matched 20ns)
cmp_ = M["functional"]["ly_vs_c110"]["STATE_MOTION"]
fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
reps = ["R1", "R2", "R3"]
for ax, contrast in zip(axes, ["Delta_PAM", "Delta_AGO", "Delta_INT"]):
    regions = list(cmp_[contrast])
    vals = {rep: [cmp_[contrast][reg]["replicas"][rep]["direction_cosine_LY_vs_compound110"] for reg in regions] for rep in reps}
    positions = np.arange(len(regions))
    w = 0.25
    for k, rep in enumerate(reps):
        ax.bar(positions + (k - 1) * w, vals[rep], w, label=rep, color=["#4C72B0", "#DD8452", "#55A868"][k])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(positions, [r.replace("_", "\n") for r in regions], rotation=60, ha="right", fontsize=6)
    ax.set_title(contrast, fontsize=9)
axes[0].set_ylabel("cos(LY d-vector, c110 d-vector)")
axes[0].legend(fontsize=7)
fig.suptitle("Do the two functional molecules move the receptor the same way?\nSTATE_MOTION, LY2119620 vs compound110 direct direction cosine per region", y=1.06, fontsize=10)
fig.savefig(FIG / "fig13_ly_vs_c110_direction.png"); plt.close(fig)
print("saved fig13")

