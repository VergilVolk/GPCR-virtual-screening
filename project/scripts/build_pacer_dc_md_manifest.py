#!/usr/bin/env python
"""Freeze a matched pilot MD manifest for PACER-DC before trajectories are run."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PILOT_CANDIDATES = ["PACER0076", "PACER0057"]
CALIBRATION = [
    {
        "candidate_id": "LY2119620",
        "canonical_smiles": "Cc1c(Cl)c(OCC(=O)N2CCN(C)CC2)nc2sc(C(=O)NC3CC3)c(N)c12",
        "role": "known_PAM_positive_control",
        "ligand_pose_source": "7V68",
    },
    {
        "candidate_id": "compound110",
        "canonical_smiles": "COC(=O)N1CCC(N2CCC(n3c(=O)n(C)c4ccccc43)CC2)CC1",
        "role": "allosteric_agonist_specificity_control",
        "ligand_pose_source": "7V6A",
    },
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "config" / "pacer_dc_pilot_md_manifest.csv",
    )
    parser.add_argument("--replicas", type=int, default=3)
    parser.add_argument("--production-ns", type=int, default=100)
    args = parser.parse_args()

    portfolio = pd.read_csv(
        ROOT
        / "results"
        / "m4_gamd_ensemble"
        / "candidate_ensemble_evidence.csv"
    ).set_index("candidate_id")
    compounds = list(CALIBRATION)
    for candidate_id in PILOT_CANDIDATES:
        compounds.append(
            {
                "candidate_id": candidate_id,
                "canonical_smiles": portfolio.loc[candidate_id, "canonical_smiles"],
                "role": "prospective_computational_hypothesis",
                "ligand_pose_source": "7TRS_common_pocket_pose_proposal",
            }
        )

    rows: list[dict[str, object]] = []
    # Shared controls are run once per paired seed, not duplicated for each ligand.
    conditions = [
        ("SHARED_CONTROL", "probe_only", "7TRS", "ACh", "none"),
        ("SHARED_CONTROL", "apo", "7TRS", "none", "none"),
    ]
    for compound in compounds:
        conditions.extend(
            [
                (
                    compound["candidate_id"],
                    "candidate_probe",
                    "7TRS",
                    "ACh",
                    compound["candidate_id"],
                ),
                (
                    compound["candidate_id"],
                    "candidate_no_probe",
                    "7TRS",
                    "none",
                    compound["candidate_id"],
                ),
            ]
        )

    info = {c["candidate_id"]: c for c in compounds}
    for compound_id, context, template, probe, allosteric in conditions:
        compound = info.get(compound_id, {})
        for replica in range(1, args.replicas + 1):
            rows.append(
                {
                    "run_id": f"{compound_id}__{context}__r{replica}",
                    "candidate_id": compound_id,
                    "context": context,
                    "replicate_id": f"paired_seed_{replica}",
                    "role": compound.get("role", "shared_context_control"),
                    "canonical_smiles": compound.get("canonical_smiles", ""),
                    "receptor_template": template,
                    "orthosteric_probe": probe,
                    "allosteric_ligand": allosteric,
                    "ligand_pose_source": compound.get("ligand_pose_source", "not_applicable"),
                    "production_ns_initial": args.production_ns,
                    "extension_policy": "extend_to_500ns_if_endpoint_CI_or_state_ESS_fails_frozen_QC",
                    "paired_seed_group": replica,
                    "status": "not_started",
                }
            )

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(
        f"Wrote {len(out)} runs: {len(compounds)} ligands x 2 contexts x "
        f"{args.replicas} replicas + 2 shared controls x {args.replicas} replicas"
    )


if __name__ == "__main__":
    main()
