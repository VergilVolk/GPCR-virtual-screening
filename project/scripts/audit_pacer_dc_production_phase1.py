#!/usr/bin/env python
"""Frozen 5-ns production QC gate for the 18 PACER-DC system/replica jobs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "config" / "pacer_dc_gpu_phase1_manifest.csv"
DEFAULT_RUNS = ROOT / "results" / "pacer_dc_production_v01"
DEFAULT_OUT = ROOT / "results" / "pacer_dc_production_phase1_audit_v01"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--runs", type=Path, default=DEFAULT_RUNS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    jobs = pd.read_csv(args.manifest)
    rows = []
    for job in jobs.itertuples(index=False):
        run = args.runs / job.system / f"replica_{int(job.replica):02d}"
        ledger_path = run / "progress.json"
        state_path = run / "state.csv"
        row = {
            "array_id": int(job.array_id),
            "system": job.system,
            "replica": int(job.replica),
            "target_ns": float(job.target_ns),
            "status": "missing",
            "completed_ns": 0.0,
            "n_report_frames": 0,
            "median_temperature_K": np.nan,
            "temperature_fraction_270_330K": np.nan,
            "median_density_g_ml": np.nan,
            "volume_drift_percent": np.nan,
            "finite": False,
            "qc_pass": False,
        }
        if ledger_path.exists():
            ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
            row["status"] = ledger.get("status", "unknown")
            row["completed_ns"] = float(ledger.get("completed_ns", 0.0))
        if state_path.exists():
            state = pd.read_csv(state_path)
            numeric = state.select_dtypes(include=[np.number]).to_numpy()
            row["n_report_frames"] = len(state)
            row["finite"] = bool(np.isfinite(numeric).all())
            temp = state["Temperature (K)"].to_numpy(float)
            density = state["Density (g/mL)"].to_numpy(float)
            volume = state["Box Volume (nm^3)"].to_numpy(float)
            row["median_temperature_K"] = float(np.median(temp))
            row["temperature_fraction_270_330K"] = float(np.mean((temp >= 270) & (temp <= 330)))
            row["median_density_g_ml"] = float(np.median(density))
            row["volume_drift_percent"] = float((volume[-1] / volume[0] - 1.0) * 100.0)
        row["qc_pass"] = bool(
            row["status"] == "complete"
            and row["completed_ns"] >= float(job.target_ns) - 1e-9
            and row["n_report_frames"] >= max(10, int(float(job.target_ns) * 1000 / 10) - 1)
            and row["finite"]
            and 285 <= row["median_temperature_K"] <= 315
            and row["temperature_fraction_270_330K"] >= 0.95
            and 0.90 <= row["median_density_g_ml"] <= 1.10
            and abs(row["volume_drift_percent"]) <= 5.0
        )
        rows.append(row)

    frame = pd.DataFrame(rows).sort_values("array_id")
    args.out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out / "job_qc.csv", index=False)
    complete = int((frame.status == "complete").sum())
    passed = int(frame.qc_pass.sum())
    paired_ready = []
    for replica in (1, 2, 3):
        subset = frame[frame.replica == replica]
        paired_ready.append({
            "replica": replica,
            "systems_passed": int(subset.qc_pass.sum()),
            "paired_context_ready": bool(len(subset) == 6 and subset.qc_pass.all()),
        })
    audit = {
        "expected_jobs": len(frame),
        "complete_jobs": complete,
        "qc_pass_jobs": passed,
        "all_jobs_pass": passed == len(frame),
        "paired_replica_readiness": paired_ready,
        "extension_to_20ns_allowed": passed == len(frame),
        "claim_boundary": (
            "This gate checks numerical/thermodynamic trajectory integrity only. Passing 5 ns does not "
            "establish conformational convergence, functional PAM activity, or model performance."
        ),
    }
    (args.out / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if complete == len(frame) and passed != len(frame):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
