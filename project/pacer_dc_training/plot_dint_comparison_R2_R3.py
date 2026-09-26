
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(
    "project/results/pacer_dc_four_context_v01/compound110"
    "/interaction_audit_R2_R3_v01"
)
INPUT = ROOT / "interaction_comparison_R2_R3_v01.json"

with INPUT.open(encoding="utf-8") as f:
    report = json.load(f)

assert report["validation_pass"] is True
records = report["window_metrics"]
means = report["mean_vector_metrics"]

names = ["dPAM", "dAGO", "dINT"]
colors = {
    "dPAM": "#277DA8",
    "dAGO": "#F8961E",
    "dINT": "#9B5DE5",
}
markers = {"dPAM": "o", "dAGO": "s", "dINT": "^"}

assert len(records) == 15
assert len(means) == 3
assert len(report["input_sha256"]) == 40

by_key = {(int(r["window"]), r["metric"]): r
          for r in records}
assert len(by_key) == 15

fig, axes = plt.subplots(
    2, 1, figsize=(9, 8), sharex=True,
    constrained_layout=True
)

windows = np.arange(5)

for name in names:
    rows = [by_key[(w, name)] for w in windows]

    cosines = [r["cosine"] for r in rows]
    distances = [r["l2_distance"] for r in rows]

    axes[0].plot(
        windows, cosines, marker=markers[name],
        color=colors[name], linewidth=2,
        markersize=7, label=name
    )
    axes[1].plot(
        windows, distances, marker=markers[name],
        color=colors[name], linewidth=2,
        markersize=7, label=name
    )

axes[0].axhline(0, color="gray",
                linestyle="--", linewidth=1)
axes[0].set_ylim(-0.65, 0.85)
axes[0].set_ylabel("Cosine similarity")
axes[0].set_title(
    "A. Cross-replica direction consistency"
)
axes[0].legend(ncol=3, loc="upper right")
axes[0].grid(alpha=0.2)

axes[1].set_ylabel("L2 distance")
axes[1].set_xlabel("Matched window")
axes[1].set_title(
    "B. Cross-replica difference magnitude"
)
axes[1].set_xticks(windows)
axes[1].set_xticklabels([f"W{w}" for w in windows])
axes[1].grid(alpha=0.2)

fig.suptitle(
    "Compound110: R2 vs R3 interaction contrasts",
    fontsize=15, fontweight="bold"
)

png = ROOT / "dint_comparison_R2_R3_v01.png"
svg = ROOT / "dint_comparison_R2_R3_v01.svg"

fig.savefig(png, dpi=220, bbox_inches="tight")
fig.savefig(svg, bbox_inches="tight")
plt.close(fig)

print("G0-5 FIGURE")
print("PNG:", png)
print("SVG:", svg)
print("INPUT_RECORDS:", len(records))
print("INPUT_HASHES:", len(report["input_sha256"]))
print(
    "G0_5_PASS:",
    png.is_file()
    and svg.is_file()
    and png.stat().st_size > 10000
    and svg.stat().st_size > 10000
)
