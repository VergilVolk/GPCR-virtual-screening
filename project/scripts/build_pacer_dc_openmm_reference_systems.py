#!/usr/bin/env python
"""Build and minimize unsolvated PACER-DC reference systems for topology QC.

The resulting systems are not production trajectories.  They verify chemical
graph/coordinate mapping and common force-field construction before the much
more expensive membrane-solvation stage.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import traceback

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from openff.toolkit import Molecule
from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, XmlSerializer
from openmm import app, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator


ROOT = Path(__file__).resolve().parents[1]
LIGAND_DIR = ROOT / "data" / "pdb" / "m4_ligands"
COMPLEX_DIR = ROOT / "results" / "pacer_dc_reference_complexes_v01"
PROTEIN = ROOT / "results" / "pacer_dc_common_protein_v01" / "7TRS_common_protein_pH74.pdb"


def pdb_ligand_block(path: Path, residue_name: str | None = None, chain: str | None = None) -> str:
    selected = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("HETATM"):
            continue
        if residue_name is not None and line[17:20].strip() != residue_name:
            continue
        if chain is not None and line[21:22].strip() != chain:
            continue
        selected.append(line)
    if not selected:
        raise ValueError(f"No ligand atoms selected from {path}")
    return "\n".join(selected + ["END"]) + "\n"


def posed_from_pdb(path: Path, template_sdf: Path, *, residue_name=None, chain=None) -> Chem.Mol:
    observed = Chem.MolFromPDBBlock(
        pdb_ligand_block(path, residue_name=residue_name, chain=chain),
        sanitize=False, removeHs=False, proximityBonding=True,
    )
    if observed is None:
        raise ValueError(f"RDKit could not parse posed ligand from {path}")
    template = Chem.RemoveHs(Chem.SDMolSupplier(str(template_sdf), removeHs=False)[0])
    observed = Chem.RemoveHs(observed)
    assigned = AllChem.AssignBondOrdersFromTemplate(template, observed)
    assigned = Chem.AddHs(assigned, addCoords=True)
    Chem.SanitizeMol(assigned)
    return assigned


def posed_from_sdf(path: Path) -> Chem.Mol:
    molecule = Chem.SDMolSupplier(str(path), removeHs=False)[0]
    if molecule is None:
        raise ValueError(f"RDKit could not parse {path}")
    Chem.SanitizeMol(molecule)
    return molecule


def openff_molecule(rdkit_molecule: Chem.Mol, name: str) -> Molecule:
    # PDB-derived heavy atoms retain residue metadata while newly added
    # hydrogens do not.  If retained, OpenFF may split one ligand into an
    # incomplete named residue plus orphan atoms, preventing template matching.
    rdkit_molecule = Chem.Mol(rdkit_molecule)
    for atom in rdkit_molecule.GetAtoms():
        atom.SetMonomerInfo(None)
    molecule = Molecule.from_rdkit(
        rdkit_molecule, allow_undefined_stereo=True, hydrogens_are_explicit=True,
    )
    molecule.name = name
    return molecule


def add_positional_restraints(system, positions, protein_atom_count: int, k=1000.0):
    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", k * unit.kilojoule_per_mole / unit.nanometer**2)
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    xyz = positions.value_in_unit(unit.nanometer)
    for index in range(protein_atom_count):
        force.addParticle(index, [float(x) for x in xyz[index]])
    system.addForce(force)


def build_one(name: str, ligands: list[Molecule], output_dir: Path) -> dict[str, object]:
    pdb = app.PDBFile(str(PROTEIN))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    for ligand in ligands:
        modeller.add(ligand.to_topology().to_openmm(), ligand.conformers[0].to_openmm())

    force_field = app.ForceField("amber14/protein.ff14SB.xml")
    for ligand in ligands:
        template = SMIRNOFFTemplateGenerator(
            molecules=[ligand], forcefield="openff_unconstrained-2.2.1.offxml",
        )
        force_field.registerTemplateGenerator(template.generator)
    system = force_field.createSystem(
        modeller.topology, nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds,
    )
    protein_atom_count = sum(1 for _ in pdb.topology.atoms())
    add_positional_restraints(system, modeller.positions, protein_atom_count)

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
    )
    platform = Platform.getPlatformByName("CPU")
    simulation = app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    initial = simulation.context.getState(getEnergy=True)
    initial_energy = initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    simulation.minimizeEnergy(maxIterations=250)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_energy = final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    system_dir = output_dir / name
    system_dir.mkdir(parents=True, exist_ok=True)
    with (system_dir / "minimized.pdb").open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(modeller.topology, final_state.getPositions(), handle)
    (system_dir / "system.xml").write_text(XmlSerializer.serialize(system), encoding="utf-8")
    (system_dir / "state.xml").write_text(
        XmlSerializer.serialize(final_state), encoding="utf-8"
    )
    return {
        "system": name,
        "ligands": [ligand.name for ligand in ligands],
        "atoms": sum(1 for _ in modeller.topology.atoms()),
        "particles": system.getNumParticles(),
        "initial_energy_kj_mol": initial_energy,
        "minimized_energy_kj_mol": final_energy,
        "energy_decreased": final_energy < initial_energy,
        "output_dir": str(system_dir),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outdir", type=Path,
        default=ROOT / "results" / "pacer_dc_openmm_reference_qc_v01",
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    ach = openff_molecule(posed_from_pdb(
        ROOT / "data" / "pdb" / "7TRS.pdb", LIGAND_DIR / "ACH_ideal.sdf",
        residue_name="ACH",
    ), "ACH")
    ly = openff_molecule(posed_from_pdb(
        COMPLEX_DIR / "LY2119620__candidate_probe.pdb", LIGAND_DIR / "2CU_ideal.sdf",
        chain="L",
    ), "LY2119620")
    compound110 = openff_molecule(posed_from_sdf(
        COMPLEX_DIR / "compound110_7TRS_redocked.sdf"
    ), "compound110")

    definitions = {
        "probe_only": [ach],
        "apo": [],
        "LY2119620__candidate_probe": [ach, ly],
        "LY2119620__candidate_no_probe": [ly],
        "compound110__candidate_probe": [ach, compound110],
        "compound110__candidate_no_probe": [compound110],
    }
    records = []
    for name, ligands in definitions.items():
        try:
            records.append({"success": True, **build_one(name, ligands, args.outdir)})
        except Exception as exc:
            records.append({
                "system": name, "success": False,
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            })
    audit = {
        "method": "PACER-DC unsolvated restrained topology/minimization QC",
        "system_count": len(records),
        "success_count": sum(bool(x["success"]) for x in records),
        "all_passed": all(bool(x["success"]) for x in records),
        "records": records,
        "claim_boundary": "Topology and minimization QC only; not production MD or biological evidence.",
    }
    (args.outdir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
