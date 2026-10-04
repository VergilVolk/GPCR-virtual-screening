from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
PRIMARY_BRANCH = "STATE_MOTION"
PRIMARY_CONTRAST = "Delta_INT"
PRIMARY_REGION = "compound110_extension"


def _prospective_rows(data: dict) -> list[dict]:
    if data.get("status") != "STAGE4_FULL_FROZEN_EVALUATION_COMPLETE":
        raise ValueError("Stage 4 frozen evaluation is not complete")
    rows = []
    for legacy_id, candidate in data["results"].items():
        regions = candidate["branches"][PRIMARY_BRANCH][PRIMARY_CONTRAST]["regions"]
        primary = regions[PRIMARY_REGION]
        rows.append({
            "legacy_id": legacy_id,
            "cluster": int(candidate["cluster"]),
            "four_context_status": "complete_unresolved",
            "delta_int_primary_region": PRIMARY_REGION,
            "delta_int_r1_r3_direction_cosine": float(primary["R1_R3_direction_cosine"]),
            "delta_int_pam_contact_consensus_cosine": float(
                regions["pam_contact_consensus"]["R1_R3_direction_cosine"]),
            "delta_int_pam_contact_union_cosine": float(
                regions["pam_contact_union"]["R1_R3_direction_cosine"]),
            "functional_prediction": "unresolved",
            "experimental_pam_validation": "not_performed",
            "stage4_claim_boundary": data["claim_boundary"],
        })
    return rows


def _legacy_rows(data: dict) -> list[dict]:
    if data.get("source_status") != "STAGE4_FULL_FROZEN_EVALUATION_COMPLETE":
        raise ValueError("Stage 4 frozen evaluation is not complete")
    rows = []
    for candidate_id, candidate in data["candidates"].items():
        rows.append({
            "legacy_id": candidate_id,
            "cluster": candidate["cluster"],
            "four_context_status": "complete_unresolved",
            "functional_prediction": "unresolved",
            "experimental_pam_validation": "not_performed",
            "stage4_claim_boundary": data["claim_boundary"],
        })
    return rows


def export_evidence(summary_path: Path, output_path: Path) -> pd.DataFrame:
    data = json.loads(summary_path.read_text(encoding="utf-8-sig"))
    rows = _prospective_rows(data) if "results" in data else _legacy_rows(data)
    result = pd.DataFrame(rows).sort_values("legacy_id").reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Export frozen four-context evidence")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "results/module4_function.csv")
    args = parser.parse_args()
    result = export_evidence(args.summary, args.output)
    print(result[["legacy_id", "cluster", "four_context_status"]].to_string(index=False))


if __name__ == "__main__":
    main()
