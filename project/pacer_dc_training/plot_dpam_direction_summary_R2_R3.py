from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("project/results/pacer_dc_four_context_v01/compound110")
IN_CSV = ROOT / "direction_audit_R2_R3_v01/dpam_direction_decomposition_R2_R3_v01.csv"
OUT_DIR = ROOT / "direction_audit_R2_R3_v01"

if not IN_CSV.is_file():
    raise SystemExit(f"MISSING_INPUT: {IN_CSV}")

rows = []
with IN_CSV.open("r", encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        parsed = {}
        for k, v in row.items():
            if k == "window":
                parsed[k] = int(v)
            else:
                parsed[k] = float(v)
        rows.append(parsed)

if len(rows) != 5:
    raise SystemExit(f"EXPECTED_5_ROWS: got {len(rows)}")

rows = sorted(rows, key=lambda x: x["window"])

w = np.array([r["window"] for r in rows], dtype=int)
labels = [f"W{i}" for i in w]

ca_cos = np.array([r["ca_cross_replica_cosine"] for r in rows])
a_cos = np.array([r["a_cross_replica_cosine"] for r in rows])
d_cos = np.array([r["dpam_cross_replica_cosine"] for r in rows])

delta_ca_a_cos = np.array([r["delta_ca_a_cosine"] for r in rows])
cross_term = np.array([r["cross_term"] for r in rows])

d_angle = np.array([r["dpam_cross_replica_angle_deg"] for r in rows])
delta_d_norm = np.array([r["delta_dpam_norm"] for r in rows])

png_path = OUT_DIR / "dpam_direction_summary_R2_R3_v01.png"
svg_path = OUT_DIR / "dpam_direction_summary_R2_R3_v01.svg"
json_path = OUT_DIR / "dpam_direction_summary_R2_R3_v01.json"

fig, axes = plt.subplots(3, 1, figsize=(10, 12), constrained_layout=True)

# Panel A
ax = axes[0]
ax.plot(w, ca_cos, marker="o", linewidth=2, label="candidate + ACh")
ax.plot(w, a_cos, marker="s", linewidth=2, label="ACh-only")
ax.plot(w, d_cos, marker="^", linewidth=2, label="dPAM")
ax.axhline(0.0, linewidth=1, linestyle="--")
ax.set_xticks(w, labels)
ax.set_ylim(-1.05, 1.05)
ax.set_ylabel("Cross-replica cosine")
ax.set_title("A. Cross-replica directional consistency of raw embeddings and dPAM")
ax.legend()
ax.grid(alpha=0.2)

# Panel B
ax = axes[1]
ax.bar(w, cross_term, width=0.55, alpha=0.75, label="Cross term: -2ΔzCA·ΔzA")
ax.plot(w, delta_ca_a_cos, marker="o", linewidth=2, label="cos(ΔzCA, ΔzA)")
ax.axhline(0.0, linewidth=1, linestyle="--")
ax.set_xticks(w, labels)
ax.set_ylabel("Geometric term")
ax.set_title("B. Relative direction of raw changes and cross term")
ax.legend()
ax.grid(alpha=0.2)

# Panel C
ax = axes[2]
ax.plot(w, d_angle, marker="o", linewidth=2, label="dPAM angle (deg)")
ax.set_xticks(w, labels)
ax.set_ylabel("Angle (deg)")
ax.set_title("C. Cross-replica directional difference and magnitude of dPAM")
ax.grid(alpha=0.2)

ax2 = ax.twinx()
ax2.plot(w, delta_d_norm, marker="s", linewidth=2, linestyle="--", label="||ΔdPAM||")
ax2.set_ylabel("L2 norm")

lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

fig.suptitle(
    "compound110: R2/R3 five-window dPAM direction sensitivity summary",
    fontsize=16
)

fig.savefig(png_path, dpi=220)
fig.savefig(svg_path)
plt.close(fig)

summary = {
    "candidate": "compound110",
    "replicas": [2, 3],
    "windows": labels,
    "figure_png": str(png_path),
    "figure_svg": str(svg_path),
    "notes": [
        "Panel A: Cross-replica cosine for candidate + ACh, ACh-only, and dPAM.",
        "Panel B: Cosine similarity between ΔzCA and ΔzA, together with the cross term -2ΔzCA·ΔzA.",
        "Panel C: Cross-replica dPAM angle and ||ΔdPAM||."
    ]
}

json_path.write_text(
    json.dumps(summary, indent=2, ensure_ascii=False),
    encoding="utf-8"
)

print("SUMMARY_FIGURE_REDRAW_PASS: True")
print("PNG:", png_path)
print("SVG:", svg_path)
print("JSON:", json_path)
