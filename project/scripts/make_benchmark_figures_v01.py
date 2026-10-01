#!/usr/bin/env python
"""Generate all benchmark figures (PNG) from the assembled master JSON.

Every figure reads only committed numbers; doc-sourced values carry a
provenance note inside the figure itself.
"""
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
FIG.mkdir(parents=True, exist_ok=True)

M = json.loads((OUT / "project_wide_benchmark_master_v01.json").read_text(encoding="utf-8"))

REGIONS = [
    "compound110_extension", "cooperativity_mutagenesis", "distal_control",
    "intracellular_microswitches", "orthosteric_activation_core",
    "orthosteric_contact_union", "pam_contact_consensus", "pam_contact_union",
    "stable_core_control",
]
RSHORT = {
    "compound110_extension": "c110_ext", "cooperativity_mutagenesis": "coop_mut",
    "distal_control": "distal", "intracellular_microswitches": "ic_switch",
    "orthosteric_activation_core": "ortho_core", "orthosteric_contact_union": "ortho_contact",
    "pam_contact_consensus": "pam_cons", "pam_contact_union": "pam_union",
    "stable_core_control": "stable_core",
}
plt.rcParams.update({"figure.dpi": 160, "font.size": 9, "axes.titlesize": 10,
                     "axes.titleweight": "bold", "savefig.bbox": "tight"})


def save(fig, name):
    fig.savefig(FIG / name)
    plt.close(fig)
    print("saved", name)


# ---------- fig01: binding macro bars ----------
def fig01():
    m13 = M["binding"]["m13_summary"]
    methods = [("official raw", m13["official"]["macro"]),
               ("ECFP4 logistic", m13["ecfp4_logistic"]["macro"]),
               ("random-label ctrl", m13["random"]["macro"]),
               ("ep40 tuned ens", m13["tuned"]["macro"])]
    seeds = sorted(M["binding"]["ep80_seeds"])
    ep80_vals = [M["binding"]["ep80_seeds"][s]["macro"]["roc_auc"] for s in seeds]
    methods.append((f"ep80 seeds ({np.mean(ep80_vals):.4f})", {"roc_auc": np.mean(ep80_vals)}))
    fa_vals = [M["binding"]["famaug_seeds"][s]["macro"]["roc_auc"] for s in seeds]
    methods.append((f"famaug seeds ({np.mean(fa_vals):.4f})", {"roc_auc": np.mean(fa_vals)}))
    methods.append(("famaug ensemble", M["binding"]["famaug_ensemble"]["ensemble_macro"]))

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    for ax, metric, title in zip(axes, ["roc_auc", "pr_auc", "ef1pct"],
                                 ["macro ROC-AUC (higher better)", "macro PR-AUC", "macro EF1%"]):
        names = [m[0] for m in methods]
        vals = [m[1].get(metric, np.nan) for m in methods]
        colors = ["#888"] * 3 + ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]
        bars = ax.bar(range(len(vals)), vals, color=colors[:len(vals)])
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(names, rotation=40, ha="right", fontsize=7)
        ax.set_title(title)
        for b, v in zip(bars, vals):
            if not np.isnan(v):
                ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.3f}", ha="center", va="bottom", fontsize=6)
        ax.axhline(0.5 if metric == "roc_auc" else 0, color="k", lw=0.6, ls="--", alpha=0.5)
    fig.suptitle("DrugCLIP binding line - 13-target strict LOSO (Science-2026 protocol)", y=1.06, fontsize=11)
    fig.text(0.01, -0.14, "source: migration_13target_v01 + family_aug_v01 committed JSONs", fontsize=6.5, color="gray")
    save(fig, "fig01_binding_macro_bars.png")


# ---------- fig02: per-target x method heatmap ----------
def fig02():
    m13 = M["binding"]["m13_summary"]
    per_t = m13["official"]["per_target"]
    targets = sorted(per_t)
    cols = {
        "official raw": m13["official"]["per_target"],
        "ECFP4 logistic": m13["ecfp4_logistic"]["per_target"],
        "ep40 ens": m13["tuned"]["per_target"],
        "ep80 seed25": M["binding"]["ep80_seeds"]["20260925"]["per_target"],
        "ep80 seed26": M["binding"]["ep80_seeds"]["20260926"]["per_target"],
        "ep80 seed27": M["binding"]["ep80_seeds"]["20260927"]["per_target"],
        "famaug seed25": M["binding"]["famaug_seeds"]["20260925"]["per_target"],
        "famaug ensemble*": None,
    }
    fam_pt = {}
    # famaug ensemble per-target approximated as mean of seeds (ensemble = mean scores)
    for s in M["binding"]["famaug_seeds"]:
        for t, v in M["binding"]["famaug_seeds"][s]["per_target"].items():
            fam_pt.setdefault(t, []).append(v["roc_auc"])
    cols["famaug ensemble*"] = {t: {"roc_auc": np.mean(v)} for t, v in fam_pt.items()}

    mat = np.array([[cols[c].get(t, {}).get("roc_auc", np.nan) for t in targets] for c in cols])
    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0.35, vmax=0.8, aspect="auto")
    ax.set_xticks(range(len(targets)), targets, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(cols)), list(cols), fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, shrink=0.8, label="ROC-AUC")
    ax.set_title("Per-target ROC-AUC: 13-target LOSO across methods")
    ax.axvline(3.5, color="k", lw=1.2)  # after M4R
    fig.text(0.5, 0.01, "* famaug ensemble per-target = mean of 3 per-seed ROC (score-level ensemble macro differs)", ha="center", fontsize=6.5, color="gray")
    save(fig, "fig02_binding_target_heatmap.png")


# ---------- fig03: seed variance ----------
def fig03():
    seeds = ["20260925", "20260926", "20260927"]
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    for ax, key, name in zip(axes, ["ep80_seeds", "famaug_seeds"], ["ep80 tuned", "family-augmented ep80"]):
        for metric, color in [("roc_auc", "#4C72B0"), ("pr_auc", "#DD8452")]:
            vals = [M["binding"][key][s]["macro"][metric] for s in seeds]
            ax.plot([1, 2, 3], vals, "o-", color=color, label=metric, ms=5)
            ax.axhline(np.mean(vals), color=color, ls=":", lw=1)
        ax.set_xticks([1, 2, 3], seeds, fontsize=8)
        ax.set_title(name)
        ax.legend(fontsize=7)
        ax.set_ylabel("macro")
    fig.suptitle("Seed stability - three training seeds, 13-target LOSO", y=1.03, fontsize=11)
    save(fig, "fig03_seed_variance.png")


# ---------- fig04: three-class cosine heatmap (money figure) ----------
def fig04():
    cm = {r["region"]: r for r in M["functional"]["cm00734_region_rows"]}
    sm = {}
    for row in M["functional"]["track_b_rows"]:
        if row["branch"] == "STATE_MOTION" and row["contrast"] == "Delta_INT":
            sm[(row["panel"], row["region"])] = row["R1R3_cosine"]
    mat = np.array([[sm[("compound110", r)], sm[("LY2119620", r)], cm[r]["CM00734_R1R3"] if r in cm else np.nan] for r in REGIONS])
    fig, ax = plt.subplots(figsize=(5.6, 5.8))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
    ax.set_xticks([0, 1, 2], ["compound110\n(ago/interacting)", "LY2119620\n(known PAM)", "CM00734\n(inactive binder)"], fontsize=8)
    ax.set_yticks(range(len(REGIONS)), [RSHORT[r] for r in REGIONS], fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > 0.55 else "black", fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.7, label="R1/R3 direction cosine")
    ax.set_title("STATE_MOTION / Delta_INT\nR1/R3 direction reproducibility (matched 20 ns)")
    fig.text(0.01, -0.04, "compound110 & LY: TRACK_B_MATCHED_20NS_RESULTS_v01.json | CM00734: synthesis doc table", fontsize=6.5, color="gray")
    save(fig, "fig04_three_class_cosine_heatmap.png")


# ---------- fig05: magnitude grid ----------
def fig05():
    sm = {}
    for row in M["functional"]["track_b_rows"]:
        if row["branch"] == "STATE_MOTION" and row["contrast"] == "Delta_INT":
            sm[(row["panel"], row["region"])] = row
    cs = M["functional"]["cm00734_status"]["primary_result"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    # left: compound110_extension magnitudes
    region = "compound110_extension"
    panels = ["compound110", "LY2119620", "CM00734"]
    reps = ["R1", "R2", "R3"]
    mag = np.array([[sm[("compound110", region)][f"{r}_norm"] for r in reps],
                    [sm[("LY2119620", region)][f"{r}_norm"] for r in reps],
                    [np.nan] * 3])
    # CM00734 magnitudes from synthesis table (R1/R2/R3 given in section 6)
    mag[2] = [0.067244, 0.068286, 0.055084]
    ax = axes[0]
    im = ax.imshow(mag, cmap="viridis", aspect="auto")
    ax.set_xticks(range(3), reps)
    ax.set_yticks(range(3), panels, fontsize=8)
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{mag[i, j]:.3f}", ha="center", va="center", fontsize=8,
                    color="white" if mag[i, j] < 0.09 else "black")
    fig.colorbar(im, ax=ax, shrink=0.8, label="mean ||Delta_INT vector||")
    ax.set_title(f"Delta_INT magnitude - {region}")
    # right: CM/LY ratio bar
    ax = axes[1]
    ratios = cs["CM_over_LY_magnitude_ratio"]
    ax.barh(range(3), [ratios["R1"], ratios["R2"], ratios["R3"]], color="#55A868")
    ax.axvline(1.0, color="k", ls="--", lw=1)
    ax.set_yticks(range(3), ["R1", "R2", "R3"])
    ax.set_xlabel("CM00734 / LY2119620 magnitude ratio")
    ax.set_title("Negative control is NOT magnitude-zero\n(direction, not size, separates classes)")
    for i, (k, v) in enumerate(ratios.items()):
        ax.text(v + 0.02, i, f"{v:.2f}", va="center", fontsize=8)
    fig.text(0.01, -0.05, "CM00734 magnitudes: PACER_DC_STAGE_A_B_FINAL_SYNTHESIS_v01.md section 6", fontsize=6.5, color="gray")
    save(fig, "fig05_magnitude_grid.png")


# ---------- fig06: branch x contrast grid ----------
def fig06():
    cm = {r["region"]: r for r in M["functional"]["cm00734_region_rows"]}
    branches = ["STATE_MOTION", "SIGNED_DRIFT"]
    contrasts = ["Delta_PAM", "Delta_AGO", "Delta_INT"]
    fig, axes = plt.subplots(2, 3, figsize=(11, 5.6), sharey=True)
    for i, br in enumerate(branches):
        for j, ct in enumerate(contrasts):
            ax = axes[i, j]
            lookup = {}
            for row in M["functional"]["track_b_rows"]:
                if row["branch"] == br and row["contrast"] == ct:
                    lookup[(row["panel"], row["region"])] = row["R1R3_cosine"]
            vec_c110 = [lookup.get(("compound110", r), np.nan) for r in REGIONS]
            vec_ly = [lookup.get(("LY2119620", r), np.nan) for r in REGIONS]
            vec_cm = [cm[r]["CM00734_R1R3"] if (br == "STATE_MOTION" and ct == "Delta_INT" and r in cm) else np.nan for r in REGIONS]
            mat = np.array([vec_c110, vec_ly, vec_cm])
            im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
            ax.set_xticks(range(len(REGIONS)), [RSHORT[r] for r in REGIONS], rotation=60, ha="right", fontsize=6)
            if j == 0:
                ax.set_yticks([0, 1, 2], ["c110", "LY", "CM00734"], fontsize=8)
            else:
                ax.set_yticks([])
            for a in range(3):
                for b in range(len(REGIONS)):
                    v = mat[a, b]
                    if not np.isnan(v):
                        ax.text(b, a, f"{v:+.2f}", ha="center", va="center", fontsize=5.5,
                                color="white" if abs(v) > 0.55 else "black")
            ax.set_title(f"{br} / {ct}", fontsize=9)
    fig.suptitle("R1/R3 direction cosine across branch x contrast x panel (CM00734 row only available for STATE_MOTION/Delta_INT)", fontsize=9, y=1.0)
    save(fig, "fig06_branch_contrast_grid.png")


# ---------- fig07: FKG v01 heatmap ----------
def fig07():
    rows = M["functional"]["fkg_v01_region_axis"]
    axes_names = ["dAGO", "dPAM", "dINT"]
    mat = np.full((3, len(REGIONS)), np.nan)
    for r in rows:
        if r["axis"] in axes_names and r["region"] in RSHORT:
            mat[axes_names.index(r["axis"]), list(REGIONS).index(r["region"])] = np.mean([float(r["R1_vs_R2_cosine"]), float(r["R1_vs_R3_cosine"]), float(r["R2_vs_R3_cosine"])])
    fig, ax = plt.subplots(figsize=(8.6, 3.2))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(len(REGIONS)), [RSHORT[r] for r in REGIONS], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(3), axes_names)
    for i in range(3):
        for j in range(len(REGIONS)):
            v = mat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(v) > 0.45 else "black")
    fig.colorbar(im, ax=ax, shrink=0.9, label="mean pooled pairwise cosine (R1R2/R1R3/R2R3)")
    ax.set_title("FKG v01 common-kernel, 600 ns compound110 - three axes x nine regions\n(historical 128D representation line)")
    save(fig, "fig07_fkg_v01_axis_heatmap.png")


# ---------- fig08: encoder/representation evolution ----------
def fig08():
    items = [
        ("OneProt-MD final head\ndPAM (5ns, doc)", 0.34, "PACER_DC_PROGRESS_20260923"),
        ("OneProt-MD dINT 5-window\n(L0-L4 all fail, doc)", 0.25, "progress doc"),
        ("Geom2Vec C0-M128 prod\n(common-kernel dINT best region)", 0.467, "FINAL_FREEZE mean pairwise"),
        ("C1-BS256 + FKG v02\nSTATE_MOTION dINT R1R3 (c110, 600ns)", 0.6197, "v02 audit report"),
        ("C1-BS256 + FKG v02\nmatched-20ns c110 dINT R1R3", 0.5837, "TRACK_B JSON"),
        ("C1-BS256 + FKG v02\nLY2119620 dINT R1R3 (new molecule)", 0.6343, "TRACK_B JSON"),
    ]
    fig, ax = plt.subplots(figsize=(9.5, 4))
    names = [i[0] for i in items][::-1]
    vals = [i[1] for i in items][::-1]
    srcs = [i[2] for i in items][::-1]
    colors = ["#999", "#999", "#4C72B0", "#DD8452", "#55A868", "#C44E52"][::-1]
    bars = ax.barh(range(len(vals)), vals, color=colors)
    ax.set_yticks(range(len(vals)), names, fontsize=7.5)
    for i, (b, v) in enumerate(zip(bars, vals)):
        ax.text(v + 0.01, i, f"{v:.3f}  [{srcs[i]}]", va="center", fontsize=6.5)
    ax.axvline(0.5, color="k", ls="--", lw=1)
    ax.text(0.502, -0.6, "0.50 descriptive reference (project convention)", fontsize=6.5, color="gray")
    ax.set_xlabel("best achieved cross-replica direction cosine")
    ax.set_title("Representation evolution: from failed OneProt pooling to frozen C1-BS256 three-class closure\n(metrics differ by generation; see report - not a single-protocol comparison)")
    save(fig, "fig08_representation_evolution.png")


# ---------- fig09: QSAR split collapse ----------
def fig09():
    b = M["qsar_static"]["baselines"]["pam_vs_inactive"]
    models = ["TanimotoKNN", "RandomForest", "ExtraTrees"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for metric, ax in [("ROC_AUC", axes[0]), ("MCC", axes[1])]:
        splits = ["scaffold", "source", "series"]
        x = np.arange(len(splits))
        w = 0.25
        for k, model in enumerate(models):
            vals = [b[f"{sp}|{model}"]["aggregate_oof"][metric] for sp in splits]
            bars = ax.bar(x + (k - 1) * w, vals, w, label=model)
            for bb, v in zip(bars, vals):
                ax.text(bb.get_x() + bb.get_width() / 2, v, f"{v:.2f}", ha="center", va="bottom", fontsize=6.5)
        ax.set_xticks(x, [f"{s} holdout" for s in splits])
        ax.set_title(f"PAM vs inactive - {metric}")
        ax.legend(fontsize=7)
    fig.suptitle("QSAR generalization collapse: scaffold-split performance dies at source/series holdout\n(the Aug bottleneck that motivated the dynamic line)", y=1.05, fontsize=10)
    fig.text(0.01, -0.08, "source: pacer_baselines_v1/baseline_metrics.json", fontsize=6.5, color="gray")
    save(fig, "fig09_qsar_split_collapse.png")


# ---------- fig10: static structure line ----------
def fig10():
    agg = M["qsar_static"]["structure_loso_aggregate"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax = axes[0]
    names = list(agg)
    vals = [agg[n]["macro_Spearman"] for n in names]
    colors = ["#C44E52" if "Vina" in n or "IFP" in n or "Ensemble" in n else "#4C72B0" for n in names]
    bars = ax.bar(range(len(vals)), vals, color=colors)
    ax.set_xticks(range(len(vals)), names, rotation=30, ha="right", fontsize=7.5)
    for bb, v in zip(bars, vals):
        ax.text(bb.get_x() + bb.get_width() / 2, v + 0.005, f"{v:.3f}", ha="center", fontsize=7)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title("Static structure line - macro within-series Spearman\n(430 molecules, leave-one-source-out)")
    ax = axes[1]
    miao = M["qsar_static"]["miao_external_macro"]
    mn = list(miao)
    mv = [miao[n]["ROC_AUC"] for n in mn]
    bars = ax.bar(range(len(mv)), mv, color="#55A868")
    ax.set_xticks(range(len(mv)), mn, rotation=30, ha="right", fontsize=7.5)
    for bb, v in zip(bars, mv):
        ax.text(bb.get_x() + bb.get_width() / 2, v + 0.005, f"{v:.3f}", ha="center", fontsize=7)
    ef = [miao[n].get("EF1%", np.nan) for n in mn]
    ax2 = ax.twinx()
    ax2.plot(range(len(ef)), ef, "ko--", ms=4, lw=1, label="EF1%")
    ax2.axhline(1.0, color="gray", ls=":", lw=1)
    ax2.set_ylabel("EF1% (1.0 = random)", fontsize=8)
    ax2.legend(fontsize=7)
    ax.set_title("Miao-2026 external screen (116-118k molecules)\nROC-AUC bars + EF1% line: high AUC, random early enrichment")
    save(fig, "fig10_static_and_external.png")


# ---------- fig11: project overview matrix ----------
def fig11():
    cells = [
        # (line, method, value, vmax, fmt)
        ("binding 13T", "official raw ROC", 0.544, 0.65),
        ("binding 13T", "ep40 ens ROC", 0.612, 0.65),
        ("binding 13T", "ep80 mean ROC", 0.624, 0.65),
        ("binding 13T", "famaug ens ROC", 0.641, 0.65),
        ("binding 20T", "official raw ROC", 0.582, 0.65),
        ("binding 20T", "ep80 ens ROC", 0.628, 0.65),
        ("functional", "c110 dINT R1R3", 0.584, 0.7),
        ("functional", "LY dINT R1R3", 0.634, 0.7),
        ("functional", "CM00734 dINT R1R3", -0.327, 0.7),
        ("QSAR series", "best model ROC (source holdout)", 0.630, 0.95),
        ("QSAR series", "best model ROC (scaffold holdout)", 0.919, 0.95),
        ("static", "VinaOnly macro Spearman", 0.043, 0.25),
        ("static", "best 2D macro Spearman", 0.207, 0.25),
        ("static", "PACER-FS (series-centered)", 0.263, 0.25),
    ]
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    rows = [c[1] for c in cells]
    vals = [c[2] for c in cells]
    lines = [c[0] for c in cells]
    vmaxs = [c[3] for c in cells]
    norm = [v / vm if v > 0 else 0 for v, vm in zip(vals, vmaxs)]
    colors = plt.cm.RdYlGn(np.clip(norm, 0, 1))
    ax.barh(range(len(vals)), [1] * len(vals), color=colors, height=0.7)
    ax.set_yticks(range(len(vals)), [f"[{l}]  {r}" for l, r in zip(lines, rows)], fontsize=7.5)
    for i, v in enumerate(vals):
        ax.text(1.01, i, f"{v:+.3f}", va="center", fontsize=8)
    ax.set_xticks([])
    ax.set_xlim(0, 1.15)
    ax.set_title("Project-wide benchmark overview (14 headline numbers, color = within-line normalized)")
    fig.text(0.01, -0.03, "values from committed result JSONs + synthesis docs; see master JSON for provenance", fontsize=6.5, color="gray")
    save(fig, "fig11_project_overview.png")


def main():
    fig01(); fig02(); fig03(); fig04(); fig05(); fig06()
    fig07(); fig08(); fig09(); fig10(); fig11()
    print("all figures done ->", FIG)


if __name__ == "__main__":
    main()



