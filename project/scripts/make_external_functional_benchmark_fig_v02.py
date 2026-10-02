#!/usr/bin/env python
"""External functional benchmark v02: two-layer summary figure."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R = Path(r"D:\CLC_integration_benchmark\project\results")
OUT = R / "pacer_m4_external_functional_v02"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 160, "font.size": 9, "axes.titlesize": 10,
                     "axes.titleweight": "bold", "savefig.bbox": "tight"})

b = json.loads((OUT / "benchmark_result.json").read_text(encoding="utf-8"))
p = json.loads((R / "pubchem_aid624126_m4_pam_v01" / "baseline_result.json").read_text(encoding="utf-8"))

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))

# Layer 1: confirmatory directional skill
ax = axes[0]
methods = ["drugclip_gpcr_triplet", "drugclip_m4_matched", "drugclip_official",
           "ecfp_rf_task", "max_train_similarity"]
labels = ["DrugCLIP\nGPCR-triplet", "DrugCLIP\nM4-matched", "DrugCLIP\nofficial",
          "ECFP-RF\n(2D)", "nearest-training\nsimilarity"]
vals = [b["macro"][m]["macro_directional_skill"] for m in methods]
lo = [b["macro"][m].get("macro_ci", [v, v])[0] if "macro_ci" in b["macro"][m] else v for m, v in zip(methods, vals)]
colors = ["#4C72B0"] * 3 + ["#DD8452", "#9aa0a6"]
y = np.arange(len(methods))
ax.barh(y, vals, color=colors, xerr=None)
for i, v in enumerate(vals):
    ax.text(v + (0.01 if v >= 0 else -0.01), i, f"{v:+.3f}", va="center",
            ha="left" if v >= 0 else "right", fontsize=8)
ax.axvline(0, color="k", lw=1.2)
ax.set_yticks(y, labels, fontsize=8)
ax.set_xlabel("macro directional skill (0 = random)")
ax.set_title(f"Layer 1: confirmatory functional endpoints\n{b['n_unique_molecules']} molecules / {b['n_endpoint_units']} records / {b['macro']['drugclip_official']['n_endpoints']} endpoints in macro")

# Layer 2: PubChem primary screen
ax = axes[1]
det = p["full_deterministic_subset"]
names = ["ECFP-Ridge", "ECFP-RF", "nearest-sim"]
keys = ["ecfp_ridge", "ecfp_rf", "max_train_similarity"]
metrics = [("roc_auc", "ROC-AUC"), ("EF@0.5pct", "EF@0.5%"), ("EF@1pct", "EF@1pct")]
x = np.arange(len(metrics))
w = 0.25
for k, (name, key) in enumerate(zip(names, keys)):
    vals = [det[key][m[0]] for m in metrics]
    bars = ax.bar(x + (k - 1) * w, vals, w, label=name,
                  color=["#4C72B0", "#DD8452", "#9aa0a6"][k])
    for bb, v in zip(bars, vals):
        ax.text(bb.get_x() + bb.get_width() / 2, v, f"{v:.2f}" if v < 5 else f"{v:.2f}",
                ha="center", va="bottom", fontsize=7)
ax.axhline(0.5, color="gray", ls=":", lw=1)
ax.set_xticks(x, [m[1] for m in metrics])
ax.set_ylabel("score / enrichment factor")
ax.set_title(f"Layer 2: PubChem AID 624126 M4 PAM primary screen\n21,445 ext. (1,445 act / 20,000 inact), 15,421 scaffolds")
ax.legend(fontsize=7.5, loc="upper left")

fig.suptitle("External functional benchmark v02: two layers, no mixed AUC\n(confirmatory truth vs large-scale primary screen - reported separately)", fontsize=11, y=1.03)
fig.text(0.01, -0.06, "layer 1: pacer_m4_external_functional_v02 | layer 2: pubchem_aid624126_m4_pam_v01 | actives' scaffold overlap with training: 0.76%, median max-sim 0.234", fontsize=6.5, color="gray")
fig.savefig(FIG / "fig_extfunc_v02_two_layers.png")
print("saved")
