#!/usr/bin/env python
"""Add a same-source PAM/inactive matched pair to the frozen v1 MD ledger.

The original manifest is never overwritten.  CM00734 is an experimentally
inactive hard negative, not a confirmed binder; binding compatibility remains
an outcome of the structural/MD gate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


ROOT = Path(__file__).resolve().parents[1]
PAIR = {"positive": "CM00717", "inactive": "CM00734"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=ROOT / "config" / "pacer_dc_pilot_md_manifest.csv")
    parser.add_argument("--benchmark", type=Path, default=ROOT / "data" / "benchmarks" / "m4_pam_v1" / "pam_vs_inactive.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "config" / "pacer_dc_pilot_md_manifest_v2.csv")
    parser.add_argument("--audit", type=Path, default=ROOT / "config" / "pacer_dc_pilot_md_manifest_v2_audit.json")
    args = parser.parse_args()

    base = pd.read_csv(args.base)
    benchmark = pd.read_csv(args.benchmark).set_index("canonical_molecule_id")
    positive = benchmark.loc[PAIR["positive"]]
    inactive = benchmark.loc[PAIR["inactive"]]
    if int(positive.target) != 1 or int(inactive.target) != 0:
        raise ValueError("Frozen PAM/inactive labels no longer match the benchmark")
    if str(positive.primary_source_id) != str(inactive.primary_source_id):
        raise ValueError("The matched pair must come from the same experimental source")

    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [fpgen.GetFingerprint(Chem.MolFromSmiles(s)) for s in (positive.canonical_smiles, inactive.canonical_smiles)]
    similarity = float(DataStructs.TanimotoSimilarity(*fps))
    additions = []
    specifications = [
        (PAIR["positive"], str(positive.canonical_smiles), "same_source_PAM_positive_control"),
        (PAIR["inactive"], str(inactive.canonical_smiles), "strict_experimental_inactive_binding_to_be_tested"),
    ]
    for candidate_id, smiles, role in specifications:
        for context in ("candidate_probe", "candidate_no_probe"):
            for replica_number in (1, 2, 3):
                probe = "ACh" if context == "candidate_probe" else "none"
                additions.append({
                    "run_id": f"{candidate_id}__{context}__r{replica_number}",
                    "candidate_id": candidate_id,
                    "context": context,
                    "replicate_id": f"paired_seed_{replica_number}",
                    "role": role,
                    "canonical_smiles": smiles,
                    "receptor_template": "7TRS",
                    "orthosteric_probe": probe,
                    "allosteric_ligand": candidate_id,
                    "ligand_pose_source": "7TRS_common_pocket_pose_proposal",
                    "production_ns_initial": 100,
                    "extension_policy": "extend_to_500ns_if_endpoint_CI_or_state_ESS_fails_frozen_QC",
                    "paired_seed_group": replica_number,
                    "status": "not_started",
                })
    result = pd.concat([base, pd.DataFrame(additions)], ignore_index=True)
    if result.run_id.duplicated().any():
        raise ValueError("Generated v2 manifest contains duplicate run_id")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    audit = {
        "base_manifest": str(args.base),
        "output_manifest": str(args.output),
        "added_pair": PAIR,
        "same_primary_source_id": str(positive.primary_source_id),
        "ecfp4_tanimoto": similarity,
        "positive_target": int(positive.target),
        "inactive_target": int(inactive.target),
        "row_count_before": int(len(base)),
        "row_count_after": int(len(result)),
        "claim_boundary": "CM00734 is experimentally inactive but is not called a binder until the binding/MD gate passes.",
    }
    args.audit.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
