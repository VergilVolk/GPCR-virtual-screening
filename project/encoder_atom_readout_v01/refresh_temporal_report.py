#!/usr/bin/env python
"""Refresh derived temporal/report artifacts without running encoder inference."""

from __future__ import annotations

import json
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[2]
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))

from project.encoder_atom_readout_v01.experiment_d import (
    EXPERIMENT_C_ROOT,
    file_record,
    render_report,
    sha256,
    temporal_diagnostics,
    write_json,
)


def main() -> None:
    output = WORKSPACE / "project/results/encoder_atom_readout_D_v01/run_001"
    extraction = json.loads((output / "extraction_audit.json").read_text(encoding="utf-8"))
    c_extraction = json.loads((EXPERIMENT_C_ROOT / "extraction_provenance.json").read_text(encoding="utf-8"))
    c_lookup = {row["relative_path"]: row for row in c_extraction["rows"]}
    rows = []
    for item in extraction["rows"]:
        rows.append(
            {
                "relative_path": item["relative_path"],
                "replica": item["replica"],
                "split": item["split"],
                "condition": item["condition"],
                "window": item["window"],
                "cache": Path(item["cache_path"]),
                "c_cache": Path(c_lookup[item["relative_path"]]["cache_path"]),
            }
        )
    numerical = json.loads((output / "numerical_health.json").read_text(encoding="utf-8"))
    temporal = temporal_diagnostics(rows, numerical, output)
    provenance = json.loads((output / "provenance.json").read_text(encoding="utf-8"))
    disk = json.loads((output / "disk_preflight.json").read_text(encoding="utf-8"))
    smoke = json.loads((output / "smoke_test.json").read_text(encoding="utf-8"))
    backbone = json.loads((output / "backbone_torsion_probe.json").read_text(encoding="utf-8"))["metrics"]
    chi = json.loads((output / "chi1_probe.json").read_text(encoding="utf-8"))["metrics"]
    comparison = json.loads((output / "representation_comparison.json").read_text(encoding="utf-8"))
    replay = json.loads((output / "replay_audit.json").read_text(encoding="utf-8"))
    integrity_path = output / "integrity_completion_audit.json"
    integrity = json.loads(integrity_path.read_text(encoding="utf-8"))
    integrity["experiment_d_implementation"] = file_record(
        WORKSPACE / "project/encoder_atom_readout_v01/experiment_d.py"
    )
    integrity["derived_temporal_artifacts_refreshed_without_encoder_inference"] = True
    write_json(integrity_path, integrity)
    report_path = output / "REPORT.md"
    report_path.write_text(
        render_report(
            provenance, disk, smoke, extraction, numerical, temporal,
            backbone, chi, comparison, replay, integrity,
        ),
        encoding="utf-8",
    )
    write_json(
        output / "EXPERIMENT_D_COMPLETE.json",
        {
            "status": integrity["status"],
            "report": file_record(report_path),
            "integrity_audit_sha256": sha256(integrity_path),
            "followup_started": False,
        },
    )


if __name__ == "__main__":
    main()
