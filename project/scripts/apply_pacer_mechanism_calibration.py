"""Apply three-anchor, assay-specific alpha-beta/tauB calibration to candidates.

The calibration is deliberately unavailable for pKB: current cross-publication
evidence shows that affinity does not transfer reliably. Alpha-beta and tauB
are calibrated only inside the lead medicinal-chemistry series and only when
the same probe/readout operational model passes QC for all three anchors.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


P = Path(__file__).resolve().parents[1]
RELATIVE = P / "results" / "pacer_assay_closed_loop_v01" / "candidate_pacer_fs_predictions.csv"
ROUND1 = P / "results" / "pacer_assay_closed_loop_v01" / "round1_results_template.csv"
SELECTION = P / "results" / "pacer_assay_closed_loop_v01" / "selection_audit.json"
OUT = P / "results" / "pacer_assay_closed_loop_v01"
ENDPOINTS = ("log_alpha_beta", "log_tauB")


def qc_true(value) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "pass", "passed", "y"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROUND1)
    args = parser.parse_args()
    candidates = pd.read_csv(RELATIVE)
    selection = json.loads(SELECTION.read_text(encoding="utf-8"))
    lead = selection["lead_reference_series"]
    results = pd.read_csv(args.results, dtype=str).fillna("") if args.results.exists() else pd.DataFrame()
    anchors = results[results.assay_role == "calibration_anchor"].copy() if len(results) else pd.DataFrame()

    candidates["predicted_log_alpha_beta"] = np.nan
    candidates["predicted_log_tauB"] = np.nan
    candidates["predicted_functional_pKB"] = np.nan
    candidates["mechanism_calibration_status"] = "awaiting_three_operational_model_anchors"
    endpoint_audit = {}
    common_assay = False
    if len(anchors) == 3:
        probes = set(anchors.assay_probe.str.strip()) - {""}
        readouts = set(anchors.assay_readout.str.strip()) - {""}
        common_assay = len(probes) == 1 and len(readouts) == 1
    if len(anchors) == 3 and common_assay:
        anchor_ids = set(selection["anchor_ids"])
        valid_ids = set(anchors.candidate_id)
        for endpoint in ENDPOINTS:
            valid = []
            for _, row in anchors.iterrows():
                try:
                    value = float(row[endpoint])
                except (TypeError, ValueError):
                    continue
                if row.candidate_id in anchor_ids and qc_true(row.qc_pass) and qc_true(row.operational_model_qc):
                    relative = float(candidates.loc[candidates.candidate_id == row.candidate_id,
                                                    "pacer_relative_sar"].iloc[0])
                    valid.append((row.candidate_id, value - relative))
            endpoint_audit[endpoint] = {"valid_anchor_ids": [item[0] for item in valid]}
            if len(valid) == 3 and valid_ids == anchor_ids:
                offsets = np.asarray([item[1] for item in valid], float)
                offset = float(offsets.mean())
                residual_sd = float(offsets.std(ddof=1)) if len(offsets) > 1 else np.nan
                mask = candidates.reference_series == lead
                candidates.loc[mask, f"predicted_{endpoint}"] = candidates.loc[mask, "pacer_relative_sar"] + offset
                candidates.loc[mask, "mechanism_calibration_status"] = "three_anchor_same_assay_calibrated_hypothesis"
                candidates.loc[~mask, "mechanism_calibration_status"] = "different_chemotype_not_transferable"
                endpoint_audit[endpoint].update({"offset": offset, "anchor_residual_sd": residual_sd})

    output = OUT / "candidate_mechanism_predictions.csv"
    candidates.to_csv(output, index=False)
    audit = {
        "lead_reference_series": lead,
        "same_probe_and_readout": common_assay,
        "endpoint_audit": endpoint_audit,
        "pKB_policy": "not predicted; cross-publication and Monash evidence show non-transferability",
        "status": "calibrated" if candidates.predicted_log_alpha_beta.notna().any() else "awaiting_valid_operational_model_anchors",
        "claim_boundary": (
            "Alpha-beta/tauB values are assay-specific hypotheses sharing the frozen relative-SAR rank. "
            "They do not confirm PAM identity and cannot transfer across chemotypes, probes, or readouts."
        ),
    }
    (OUT / "candidate_mechanism_calibration_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
