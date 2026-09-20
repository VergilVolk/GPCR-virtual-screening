#!/usr/bin/env python
"""Audit whether this host can prepare and run the frozen PACER-DC MD ledger."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil

from openmm import Platform


ROOT = Path(__file__).resolve().parents[1]


def available(module: str) -> bool:
    try:
        return importlib.util.find_spec(module) is not None
    except ModuleNotFoundError:
        return False


def main() -> None:
    commands = {name: shutil.which(name) for name in
                ["gmx", "gmx_mpi", "pmemd.cuda", "pmemd", "sander", "tleap"]}
    modules = {name: available(name) for name in
               ["openmm", "pdbfixer", "openmmforcefields", "openff.toolkit",
                "parmed", "MDAnalysis", "rdkit"]}
    platforms = [Platform.getPlatform(i).getName()
                 for i in range(Platform.getNumPlatforms())]
    parameterization_audit_path = (
        ROOT / "results" / "pacer_dc_ligand_parameterization_audit.json"
    )
    parameterization_smoke_passed = False
    if parameterization_audit_path.exists():
        try:
            parameterization_smoke_passed = bool(
                json.loads(parameterization_audit_path.read_text(encoding="utf-8"))
                .get("all_passed", False)
            )
        except (json.JSONDecodeError, OSError):
            parameterization_smoke_passed = False
    ligand_parameterization_ready = bool(
        (modules["openff.toolkit"] or commands["tleap"])
        and parameterization_smoke_passed
    )
    audit = {
        "commands": commands,
        "python_modules": modules,
        "openmm_platforms": platforms,
        "trajectory_analysis_ready": bool(modules["MDAnalysis"] and modules["openmm"]),
        "protein_system_prep_partial": bool(modules["pdbfixer"] and modules["openmm"]),
        "ligand_parameterization_ready": ligand_parameterization_ready,
        "parameterization_smoke_passed": parameterization_smoke_passed,
        "parameterization_audit": str(parameterization_audit_path),
        "production_ready": bool(ligand_parameterization_ready and modules["openmm"]),
        "blocking_reason": (
            "none" if ligand_parameterization_ready else
            "Reference-ligand parameterization smoke test has not passed; imports alone are insufficient."
        ),
        "recommended_environment": "project/environment_pacer_dc_md.yml",
    }
    out = ROOT / "results" / "pacer_dc_md_environment_audit.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
