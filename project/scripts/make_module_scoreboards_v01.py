#!/usr/bin/env python
"""Module scoreboards: per-module horizontal bar charts + cascade funnel.

Reads only committed numbers: method registry v02 (binding), master v01
(static/QSAR/functional), handoff audit v02 (funnel counts).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(r"D:\CLC_integration_benchmark\project\results\project_wide_integration_benchmark_v02")
MASTER = json.loads(Path(r"D:\CLC_geom2vec_pilot\project\results\project_wide_integration_benchmark_v01\project_wide_benchmark_master_v01.json").read_text(encoding="utf-8"))
FIG = HERE / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 160, "font.size": 9, "axes.titlesize": 10,
                     "axes.titleweight": "bold", "savefig.bbox": "tight"})

# ---- module A: binding (from teammate's registry, same protocol 13T) ----
reg = list(csv.DictReader((HERE / "benchmark_method_registry_v02.csv").open(encoding="utf-8")))
a_rows = [("random-label control", 0.49621, "#9aa0a6", "leakage floor"),
          ("official DrugCLIP", 0.54422, "#888d94", "starting point"),
          ("ECFP4 logistic (2D)", 0.56807, "#4C72B0", "strong ligand-only"),
          ("ep80 score ensemble", 0.60823, "#DD8452", "ours"),
          ("famaug score ensemble", 0.64136, "#C44E52", "ours, current headline")]

# ---- module B: static gating (macro Spearman, same task) ----
agg = MASTER["qsar_static"]["structure_loso_aggregate"]
b_rows = sorted([(k, v["macro_Spearman"]) for k, v in agg.items()], key=lambda x: x[1])

# ---- module C: dynamic discrimination (direction cosine) ----
c_rows = [("CM00734 (exp. inactive)", -0.327, "#C44E52"),
          ("compound110 (ago-PAM)", 0.584, "#DD8452"),
          ("LY2119620 (known PAM)", 0.634, "#55A868")]

# ---- module D: QSAR collapse (ROC, same models two splits) ----
base = MASTER["qsar_static"]["baselines"]["pam_vs_inactive"]
models = ["RandomForest", "ExtraTrees"]
d_scaffold = [base[f"scaffold|{m}"]["aggregate_oof"]["ROC_AUC"] for m in models]
d_source = [base[f"source|{m}"]["aggregate_oof"]["ROC_AUC"] for m in models]

fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.45)

# A
ax = fig.add_subplot(gs[0, 0])
names = [r[0] for r in a_rows][::-1]; vals = [r[1] for r in a_rows][::-1]
cols = [r[2] for r in a_rows][::-1]
ax.barh(names, vals, color=cols)
for i, v in enumerate(vals):
    ax.text(v + 0.004, i, f"{v:.3f}", va="center", fontsize=8)
ax.axvline(0.5, color="k", ls="--", lw=1)
ax.set_xlim(0.45, 0.68)
ax.set_xlabel("macro ROC-AUC (13-target strict LOSO)")
ax.set_title("A. Binding module | same-protocol leaderboard\n(ruler: retrieval quality)")

# B
ax = fig.add_subplot(gs[0, 1])
ax.barh([r[0] for r in b_rows], [r[1] for r in b_rows], color="#4C72B0")
for i, (n, v) in enumerate(b_rows):
    ax.text(v + 0.004, i, f"{v:.3f}", va="center", fontsize=7.5)
ax.set_xlabel("macro within-series Spearman")
ax.set_title("B. Static structure module | ALL methods near zero\n(-> used only as gate, never as score)")

# C
ax = fig.add_subplot(gs[1, 0])
ax.barh([r[0] for r in c_rows], [r[1] for r in c_rows], color=[r[2] for r in c_rows])
for i, (n, v, c) in enumerate(c_rows):
    ax.text(v + (0.02 if v > 0 else -0.02), i, f"{v:+.3f}", va="center",
            ha="left" if v > 0 else "right", fontsize=9, fontweight="bold")
ax.axvline(0, color="k", lw=1.2)
ax.axvline(0.5, color="gray", ls=":", lw=1)
ax.text(0.505, 2.35, "0.50 frozen descriptive line", fontsize=6.5, color="gray")
ax.set_xlim(-0.55, 0.8)
ax.set_xlabel("R1/R3 direction cosine (Delta_INT, compound110_extension)")
ax.set_title("C. Dynamic module | our endpoint, first of its kind\n(ruler: cross-replica direction reproducibility)")

# D
ax = fig.add_subplot(gs[1, 1])
y = np.arange(len(models)); h = 0.32
ax.barh(y + h/2, d_scaffold, h, label="scaffold holdout (in-distribution)", color="#8fbf9f")
ax.barh(y - h/2, d_source, h, label="source holdout (real generalization)", color="#C44E52")
for yy, v in zip(y + h/2, d_scaffold): ax.text(v + 0.005, yy, f"{v:.3f}", va="center", fontsize=8)
for yy, v in zip(y - h/2, d_source): ax.text(v + 0.005, yy, f"{v:.3f}", va="center", fontsize=8)
ax.set_yticks(y, models); ax.legend(fontsize=7, loc="lower right")
ax.set_xlabel("ROC-AUC (PAM vs inactive)")
ax.set_title("D. Why QSAR module was closed | the collapse that\nmotivated the whole dynamic line")

fig.suptitle("Module scoreboards: each module has its own ruler - numbers are NOT comparable across panels", fontsize=12, y=0.99)
fig.savefig(FIG / "fig_v02_07_module_scoreboards.png")
plt.close(fig)
print("saved fig_v02_07")

# ---- funnel: the only legitimate 'overall' comparison ----
fig, ax = plt.subplots(figsize=(9.5, 4.6))
stages = ["generated\n(new molecules)", "passed medchem +\nDrugCLIP dual-model", "passed docking/IFP\nstructure gate", "final candidates\n(frozen 4-ctx MD + FKG)"]
counts = [200, 24, 5, None]
survival = ["", "12%", "0.4% of 200\n(2 of top-100 dual)", "target 1-3\nPENDING"]
colors = ["#4C72B0", "#DD8452", "#55A868", "#bbbbbb"]
x = np.arange(4)
vals = [200, 24, 5, 1.5]
bars = ax.bar(x, vals, color=colors, width=0.6)
ax.set_yscale("log")
ax.set_ylim(0.8, 400)
ax.set_xticks(x, stages, fontsize=8)
for i, (b, v) in enumerate(zip(bars, vals)):
    label = f"{counts[i]}" if counts[i] else "1-3\n(pending)"
    ax.text(b.get_x() + b.get_width()/2, v * 1.15, label, ha="center", fontsize=10, fontweight="bold")
    if survival[i]:
        ax.text(b.get_x() + b.get_width()/2, v * 0.45, survival[i], ha="center", fontsize=7, color="white")
ax.set_ylabel("molecules remaining (log)")
ax.set_title("The overall 'score' is this funnel, not a single number\n(each stage removes molecules; final stage is the pending work)")
ax.annotate("", xy=(2.62, 3), xytext=(2.05, 8),
            arrowprops=dict(arrowstyle="->", color="#C44E52", lw=2))
ax.text(2.62, 1.6, "NEXT 72h", color="#C44E52", fontsize=9, fontweight="bold", ha="center")
fig.text(0.01, -0.04, "source: candidate_handoff_audit_v02.json (200/24/5 verified) + mainline doc v01", fontsize=6.5, color="gray")
fig.savefig(FIG / "fig_v02_08_cascade_funnel.png")
plt.close(fig)
print("saved fig_v02_08")
