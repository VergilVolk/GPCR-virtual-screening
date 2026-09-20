#!/usr/bin/env python
"""Prepare the common 7TRS protein assembly without inventing missing loops."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from openmm import app
from openmm.app import PDBFile
from pdbfixer import PDBFixer


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=ROOT / "data" / "pdb" / "7TRS.pdb")
    parser.add_argument(
        "--outdir", type=Path,
        default=ROOT / "results" / "pacer_dc_common_protein_v01",
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    fixer = PDBFixer(filename=str(args.input))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    detected_missing_residue_blocks = {
        f"{chain_index}:{residue_index}": list(names)
        for (chain_index, residue_index), names in fixer.missingResidues.items()
    }
    # Large unresolved GPCR loops and termini must not be hallucinated by an
    # automated preparation step.  Their treatment requires a separately
    # justified model; the pilot keeps the experimental construct boundaries.
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    nonstandard = [
        {"residue": str(residue), "replacement": replacement}
        for residue, replacement in fixer.nonstandardResidues
    ]
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    missing_atoms_before = sum(len(v) for v in fixer.missingAtoms.values())
    missing_terminals_before = sum(len(v) for v in fixer.missingTerminals.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    output = args.outdir / "7TRS_common_protein_pH74.pdb"
    with output.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    force_field = app.ForceField("amber14/protein.ff14SB.xml")
    system = force_field.createSystem(
        fixer.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds,
    )
    audit = {
        "input": str(args.input), "output": str(output),
        "chains": sum(1 for _ in fixer.topology.chains()),
        "residues": sum(1 for _ in fixer.topology.residues()),
        "atoms": sum(1 for _ in fixer.topology.atoms()),
        "system_particles": system.getNumParticles(),
        "detected_missing_residue_blocks_not_built": detected_missing_residue_blocks,
        "nonstandard_replacements": nonstandard,
        "missing_atoms_added": missing_atoms_before,
        "missing_terminal_atoms_added": missing_terminals_before,
        "forcefield_template_test": "passed",
        "claim_boundary": "Protein-coordinate preparation QC only; not a production membrane system.",
    }
    (args.outdir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
