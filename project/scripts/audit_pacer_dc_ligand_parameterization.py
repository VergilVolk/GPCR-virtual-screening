#!/usr/bin/env python
"""Audit deterministic OpenFF parameterization of PACER-DC reference ligands.

This is a force-field readiness test, not an MD result.  It records exact
chemical identity, stereochemistry warnings and whether an OpenMM System can
be constructed for each ligand.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import traceback

from openff.toolkit import Molecule
from openmm import app
from openmmforcefields.generators import SMIRNOFFTemplateGenerator


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ligand-dir", type=Path,
        default=ROOT / "data" / "pdb" / "m4_ligands",
    )
    parser.add_argument(
        "--force-field", default="openff_unconstrained-2.2.1.offxml",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "results" / "pacer_dc_ligand_parameterization_audit.json",
    )
    args = parser.parse_args()

    records: list[dict[str, object]] = []
    for path in sorted(args.ligand_dir.glob("*_ideal.sdf")):
        record: dict[str, object] = {
            "file": str(path), "sha256": sha256(path), "success": False,
        }
        try:
            molecule = Molecule.from_file(path, allow_undefined_stereo=True)
            # Use OpenMM's residue-template route.  Calling
            # ForceField.create_openmm_system() through Interchange ties this
            # audit to a second versioned API and previously produced a false
            # environment-ready result.
            template = SMIRNOFFTemplateGenerator(
                molecules=[molecule], forcefield=args.force_field
            )
            openmm_force_field = app.ForceField()
            openmm_force_field.registerTemplateGenerator(template.generator)
            system = openmm_force_field.createSystem(
                molecule.to_topology().to_openmm(),
                nonbondedMethod=app.NoCutoff,
                constraints=None,
            )
            record.update({
                "success": True,
                "name": molecule.name,
                "n_atoms": molecule.n_atoms,
                "n_particles": system.getNumParticles(),
                "formal_charge": int(molecule.total_charge.m),
                "mapped_smiles": molecule.to_smiles(mapped=True),
                "canonical_isomeric_smiles": molecule.to_smiles(
                    isomeric=True, explicit_hydrogens=False
                ),
            })
        except Exception as exc:  # audit must retain every failure
            record["error"] = f"{type(exc).__name__}: {exc}"
            record["traceback"] = traceback.format_exc()
        records.append(record)

    audit = {
        "method": "PACER-DC ligand parameterization readiness",
        "force_field": args.force_field,
        "ligand_count": len(records),
        "success_count": sum(bool(x["success"]) for x in records),
        "all_passed": bool(records) and all(bool(x["success"]) for x in records),
        "records": records,
        "claim_boundary": "Parameterization QC only; not evidence of binding or PAM function.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
