#!/usr/bin/env python
"""Freeze the six-system PACER-DC numerical-stability gate."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "pacer_dc_membrane_smoke_v01"
SYSTEMS = (
    "apo",
    "probe_only",
    "LY2119620__candidate_no_probe",
    "LY2119620__candidate_probe",
    "compound110__candidate_no_probe",
    "compound110__candidate_probe",
)


def main() -> None:
    rows = []
    for name in SYSTEMS:
        run_dir = SOURCE / name
        audit = json.loads((run_dir / "audit.json").read_text(encoding="utf-8"))
        pre = json.loads((run_dir / "preflight.json").read_text(encoding="utf-8"))
        max_force = float(pre["max_force_kj_mol_nm"])
        constraint = float(pre["max_constraint_error_nm"])
        atom = str(pre["max_force_atom"])
        passed = bool(
            audit["finite"]
            and max_force < 5000.0
            and constraint < 1e-5
            and any(x in atom for x in (":HOH", ":WAT", ":SOL"))
        )
        rows.append({
            "system": name,
            "finite": bool(audit["finite"]),
            "simulated_ps": float(audit["simulated_ps"]),
            "max_force_kj_mol_nm": max_force,
            "max_force_atom": atom,
            "max_constraint_error_nm": constraint,
            "passed": passed,
        })

    overall = all(row["passed"] for row in rows)
    summary = {
        "gate": "PACER-DC six-system numerical stability",
        "criteria": {
            "finite_energy_and_coordinates": True,
            "max_force_kj_mol_nm_lt": 5000.0,
            "max_constraint_error_nm_lt": 1e-5,
            "largest_residual_force_on_solvent": True,
        },
        "n_passed": sum(row["passed"] for row in rows),
        "n_total": len(rows),
        "passed": overall,
        "systems": rows,
        "claim_boundary": (
            "This gate establishes numerical stability for 0.2 ps only. It does not establish "
            "equilibrium, converged receptor dynamics, PAM activity, or predictive performance."
        ),
    }
    (SOURCE / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (SOURCE / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        "# PACER-DC six-system smoke-test gate",
        "",
        f"**Result: {'PASS' if overall else 'FAIL'} ({summary['n_passed']}/{summary['n_total']})**",
        "",
        "| System | Finite | Max force | Max-force atom | Constraint error (nm) | Gate |",
        "|---|---:|---:|---|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['system']} | {row['finite']} | {row['max_force_kj_mol_nm']:.1f} | "
            f"{row['max_force_atom']} | {row['max_constraint_error_nm']:.2e} | "
            f"{'PASS' if row['passed'] else 'FAIL'} |"
        )
    lines.extend([
        "",
        "## Claim boundary",
        "",
        summary["claim_boundary"],
        "",
    ])
    (SOURCE / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"passed": overall, "n_passed": summary["n_passed"], "n_total": len(rows)}))
    if not overall:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
