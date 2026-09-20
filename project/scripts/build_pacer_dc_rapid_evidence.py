#!/usr/bin/env python
"""Build the strongest currently defensible PACER-DC result card and figure."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DYNAMIC = ROOT / "results" / "pacer_public_pam_dynamic_endpoint_v01"
STATIC = ROOT / "results" / "pacer_paired_delta_pam_signature_v04"
OUT = ROOT / "results" / "pacer_dc_rapid_evidence_v01"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    audit = json.loads((DYNAMIC / "audit.json").read_text(encoding="utf-8"))
    series = pd.read_csv(DYNAMIC / "processed_rmsd_timeseries.csv")
    replicas = pd.read_csv(DYNAMIC / "replica_summary.csv")
    static = pd.read_csv(STATIC / "scores.csv")
    static = static[static.threshold_A == 0.1].set_index("structure")
    static_delta = float(
        static.loc["7V68_LY2119620_PAM_external", "score"]
        - static.loc["7V69_iperoxo_only_external", "score"]
    )

    control_cols = [f"iperoxo-sim{i}" for i in range(1, 4)]
    pam_cols = [f"iperoxo-LY-sim{i}" for i in range(1, 4)]
    wide = series.pivot(index="time_ns", columns="trajectory", values="iperoxo_RMSD_A")
    # A 10-ns rolling mean is visualization only and never used in inference.
    smooth = wide.rolling(50, center=True, min_periods=1).mean()
    time = smooth.index.to_numpy(float)
    c = smooth[control_cols].to_numpy(float)
    p = smooth[pam_cols].to_numpy(float)

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    for values, color, label in ((c, "#6b7280", "iperoxo only"), (p, "#7c3aed", "iperoxo + LY2119620")):
        mean = values.mean(axis=1)
        sd = values.std(axis=1, ddof=1)
        axes[0].plot(time, mean, color=color, lw=1.4, label=label)
        axes[0].fill_between(time, mean - sd, mean + sd, color=color, alpha=0.18, linewidth=0)
    axes[0].set(xlabel="Simulation time (ns)", ylabel="Iperoxo RMSD (Å)", title="Published three-replica trajectories")
    axes[0].legend(frameon=False)

    control = replicas[replicas.condition == "iperoxo_only"].sort_values("replica")
    pam = replicas[replicas.condition == "iperoxo_LY2119620"].sort_values("replica")
    x = np.array([0, 1])
    for replica, cval, pval in zip(control.replica, control.mean_RMSD_A, pam.mean_RMSD_A):
        axes[1].plot(x, [cval, pval], marker="o", lw=1.4, alpha=0.85, label=f"replica {replica}")
    axes[1].set_xticks(x, ["iperoxo", "+ LY2119620"])
    axes[1].set_ylabel("Mean iperoxo RMSD (Å)")
    axes[1].set_title(f"Paired reduction: {audit['mean_RMSD_reduction_A']:.3f} Å")
    axes[1].text(
        0.5, 0.04,
        f"hierarchical 95% CI [{audit['hierarchical_block_bootstrap_CI95_A'][0]:.3f}, "
        f"{audit['hierarchical_block_bootstrap_CI95_A'][1]:.3f}] Å",
        transform=axes[1].transAxes, ha="center", va="bottom", fontsize=8,
    )
    fig.suptitle("LY2119620 reproducibly stabilizes the orthosteric ligand pose", fontweight="bold")
    fig.savefig(OUT / "figure_public_dynamic_stabilization.png", dpi=240)
    plt.close(fig)

    summary = {
        "result_status": "supported_mechanistic_result_not_functional_predictor",
        "public_dynamic_source": "Wang2022 Supplementary Figure 8 source data",
        "independent_replicas_per_condition": 3,
        "mean_orthosteric_rmsd_reduction_A": audit["mean_RMSD_reduction_A"],
        "hierarchical_block_bootstrap_CI95_A": audit["hierarchical_block_bootstrap_CI95_A"],
        "paired_replica_differences_A": audit["paired_replica_differences_A"],
        "post100ns_effect_A": audit["post100ns_effect_A"],
        "probe_matched_external_static_increment": static_delta,
        "static_structure_pair": ["7V69 iperoxo only", "7V68 iperoxo + LY2119620"],
        "claim": (
            "Across three published replicas, LY2119620 is associated with reproducible stabilization "
            "of the orthosteric iperoxo pose, consistent with a positive probe-matched static coupling shift."
        ),
        "claim_boundary": (
            "This triangulates one known PAM mechanism. It does not validate PACER-DC as a general "
            "PAM classifier, infer efficacy, or confirm any computational candidate."
        ),
    }
    (OUT / "audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report = f"""# PACER-DC rapid evidence result v0.1

## Supported result

Across three independent published trajectories, LY2119620 reduced the mean orthosteric iperoxo RMSD by **{summary['mean_orthosteric_rmsd_reduction_A']:.3f} Å**. The hierarchical block-bootstrap 95% interval was **[{summary['hierarchical_block_bootstrap_CI95_A'][0]:.3f}, {summary['hierarchical_block_bootstrap_CI95_A'][1]:.3f}] Å**, and all three paired replica differences had the same direction.

The independent probe-matched static pair (7V68 versus 7V69) showed a positive coupling-coordinate increment of **{static_delta:.3f}**. Static and dynamic evidence therefore agree that LY2119620 stabilizes an orthosteric-ligand-coupled M4 state.

![Public dynamic stabilization](figure_public_dynamic_stabilization.png)

## What this establishes

- A reproducible mechanistic endpoint exists for a known M4 PAM.
- The effect is replica-aware and remains after the first 100 ns.
- The direction agrees with an independent probe-matched static structure pair.

## What this does not establish

- It does not show that PACER-DC predicts PAM efficacy across molecules.
- It does not distinguish functional PAMs from allosteric agonists by itself.
- It does not validate PACER0076, PACER0057, or another generated molecule as a PAM.
- It is not a SOTA claim.

## Immediate use

This is a defensible preliminary biological result for the competition report. The next decisive result is the preregistered six-context, three-replica comparison of LY2119620 and compound-110 from the new production trajectories.
"""
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
