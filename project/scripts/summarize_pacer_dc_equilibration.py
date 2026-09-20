#!/usr/bin/env python
"""Summarize matched PACER-DC short-equilibration audits without overstating them."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "pacer_dc_short_equilibration_v01"


def main() -> None:
    rows = []
    for path in sorted(SOURCE.glob("*/replica_*/audit.json")):
        audit = json.loads(path.read_text(encoding="utf-8"))
        final = audit["stages"][-1]
        # The OpenMM reporter uses the engine's constraint-aware degrees of freedom.
        # Prefer its last recorded temperature, including for legacy pilot audits.
        stage_name = str(final["stage"]).lower()
        state_csv = path.parent / f"{stage_name}_state.csv"
        with state_csv.open(newline="", encoding="utf-8") as handle:
            records = list(csv.DictReader(handle))
        temperature = float(records[-1]["Temperature (K)"])
        numerical = bool(audit.get("numerically_stable", audit.get("passed", False)))
        thermal = bool(temperature is not None and 270.0 <= float(temperature) <= 330.0)
        rows.append({
            "system": audit["system"],
            "replica": audit["replica"],
            "total_ps": float(audit["nvt_ps"]) + float(audit["npt_ps"]),
            "final_temperature_K": float(temperature),
            "volume_change_percent": float(audit["volume_change_percent"]),
            "numerically_stable": numerical,
            "thermalized_270_330K": thermal,
            "equilibration_gate": numerical and thermal,
        })
    if not rows:
        raise SystemExit("No equilibration audits found.")
    with (SOURCE / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "n_runs": len(rows),
        "n_numerically_stable": sum(x["numerically_stable"] for x in rows),
        "n_thermalized": sum(x["thermalized_270_330K"] for x in rows),
        "n_equilibration_gate": sum(x["equilibration_gate"] for x in rows),
        "rows": rows,
        "claim_boundary": (
            "Numerical stability and target-temperature attainment are necessary but not sufficient for "
            "equilibrium, convergence, replica independence, or PAM functional inference."
        ),
    }
    (SOURCE / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    lines = [
        "# PACER-DC short-equilibration audit",
        "",
        f"Numerically stable: **{summary['n_numerically_stable']}/{summary['n_runs']}**  ",
        f"Temperature gate (270–330 K): **{summary['n_thermalized']}/{summary['n_runs']}**  ",
        f"Combined gate: **{summary['n_equilibration_gate']}/{summary['n_runs']}**",
        "",
        "| System | ps | Final T (K) | Volume change | Numerical | Thermal |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['system']} | {row['total_ps']:.2f} | {row['final_temperature_K']:.1f} | "
            f"{row['volume_change_percent']:.2f}% | {row['numerically_stable']} | "
            f"{row['thermalized_270_330K']} |"
        )
    lines.extend(["", "## Claim boundary", "", summary["claim_boundary"], ""])
    (SOURCE / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in (
        "n_runs", "n_numerically_stable", "n_thermalized", "n_equilibration_gate"
    )}))


if __name__ == "__main__":
    main()
