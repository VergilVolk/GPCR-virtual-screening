#!/usr/bin/env python
"""Assemble only observed public/static evidence into a PACER-DC smoke-test table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data" / "benchmarks" / "pacer_dc_reference_evidence.csv",
    )
    args = parser.parse_args()

    rmsd = pd.read_csv(
        ROOT / "results" / "pacer_public_pam_dynamic_endpoint_v01" / "replica_summary.csv"
    )
    static = pd.read_csv(
        ROOT / "results" / "pacer_paired_delta_pam_signature_v04" / "scores.csv"
    )
    static = static[static.threshold_A == 0.1].set_index("structure")

    rows: list[dict[str, object]] = []
    for row in rmsd.itertuples(index=False):
        context = "probe_only" if row.condition == "iperoxo_only" else "candidate_probe"
        rows.append(
            {
                "candidate_id": "LY2119620",
                "context": context,
                "replicate_id": f"Wang2022_MD_{int(row.replica)}",
                "metric": "orthosteric_pose_rmsd_A",
                "value": row.mean_RMSD_A,
                "evidence_level": "trajectory",
                "source_id": "Wang2022_SuppFig8_SourceData",
            }
        )

    for context, structure in [
        ("probe_only", "7V69_iperoxo_only_external"),
        ("candidate_probe", "7V68_LY2119620_PAM_external"),
    ]:
        rows.append(
            {
                "candidate_id": "LY2119620",
                "context": context,
                "replicate_id": structure.split("_")[0],
                "metric": "coupling_coordinate",
                "value": float(static.loc[structure, "score"]),
                "evidence_level": "static_structure",
                "source_id": structure.split("_")[0],
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(f"Wrote {len(rows)} observed evidence rows to {args.output}")
    print("No candidate_no_probe, apo, or binding rows were invented.")


if __name__ == "__main__":
    main()
