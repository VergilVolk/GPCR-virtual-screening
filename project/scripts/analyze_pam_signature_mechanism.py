"""Interpret residue-pair drivers and time blocks of the frozen PAM signature."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import test_consensus_pam_signature as sigmod


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "pacer_consensus_pam_signature_v01"


def main():
    q, p, s, pairs = sigmod.load_static()
    x = np.load(OUT / "ca_pair_distance_features_210.npy")
    frame = pd.read_csv(OUT / "frame_signature_scores.csv")
    dq, dp = q - s, p - s
    mask = (dq * dp > 0) & (np.abs(dq) >= .10) & (np.abs(dp) >= .10)
    agreement = np.minimum(np.abs(dq), np.abs(dp)) / np.maximum(np.abs(dq), np.abs(dp))
    v = ((dq + dp) / 2)[mask] * np.sqrt(agreement[mask])
    denom = float(np.dot(v, v))
    contributions = ((x[:, mask] - s[mask]) * np.sqrt(agreement[mask])) * v / denom
    selected_pairs = np.asarray(pairs, dtype=int)[mask]

    rep2 = frame.replica.to_numpy() == 2
    diff = contributions[rep2].mean(0) - contributions[~rep2].mean(0)
    rows = []
    for j in np.argsort(diff)[:20]:
        rows.append({"residue_i": int(selected_pairs[j, 0]), "residue_j": int(selected_pairs[j, 1]),
                     "replica2_minus_other_contribution": float(diff[j]),
                     "replica2_mean_contribution": float(contributions[rep2, j].mean()),
                     "other_mean_contribution": float(contributions[~rep2, j].mean()),
                     "7TRQ_minus_7TRS_A": float(dq[mask][j]),
                     "7TRP_minus_7TRS_A": float(dp[mask][j])})
    drivers = pd.DataFrame(rows)
    drivers.to_csv(OUT / "replica2_negative_signature_drivers.csv", index=False)

    blocks = frame.copy()
    blocks["block_50ns"] = blocks.frame_ns // 50
    block_summary = blocks.groupby(["replica", "block_50ns"]).signature_projection.agg(
        ["mean", "median", lambda z: float(np.mean(z > 0))]).reset_index()
    block_summary.columns = ["replica", "block_50ns", "mean", "median", "fraction_positive"]
    block_summary.to_csv(OUT / "signature_50ns_blocks.csv", index=False)

    report = "# Consensus-PAM Signature Mechanism Audit\n\n"
    report += "## 50 ns blocks\n\n" + block_summary.to_markdown(index=False, floatfmt=".4f")
    report += "\n\n## Replica 2 strongest negative residue-pair drivers\n\n"
    report += drivers.head(12).to_markdown(index=False, floatfmt=".4f")
    report += ("\n\nThese are contributors to a static-anchor structural projection, not causal "
               "residues and not efficacy determinants.\n")
    (OUT / "MECHANISM_AUDIT.md").write_text(report, encoding="utf-8")
    print(block_summary[block_summary.replica == 2].to_string(index=False))
    print(drivers.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
