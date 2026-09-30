#!/usr/bin/env python
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from openff.toolkit import Molecule
from openff.units import unit as off_unit
from openmm import (
    CustomExternalForce,
    LangevinMiddleIntegrator,
    MonteCarloMembraneBarostat,
    Platform,
    XmlSerializer,
    app,
    unit,
)
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "pdb" / "7TRS.pdb"
POSE = ROOT / "results" / "pacer_dc_reference_complexes_v01" / "CM00734_7TRS_redocked.sdf"
OUT = ROOT / "results" / "pacer_dc_membrane_reference_v01"

SRCROOT = Path(
    "/mnt/c/projects/GPCR-virtual-screening-fkg-v02/"
    "project/results/pacer_dc_membrane_reference_v01"
)

DEFS = {
    "CM00734__candidate_no_probe": {
        "source": SRCROOT / "LY2119620__candidate_no_probe" / "minimized.pdb",
        "drop_chain": "J",
        "probe": False,
    },
    "CM00734__candidate_probe": {
        "source": SRCROOT / "LY2119620__candidate_probe" / "minimized.pdb",
        "drop_chain": "K",
        "probe": True,
    },
}

CENTER = np.array([110.2121366, 107.8369051, 68.1181557], dtype=float)
ACH_SMILES = "CC(=O)OCC[N+](C)(C)C"
PROTEIN_CHAINS = set("ABCDE")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def ca_xyz(path: Path, chain: str) -> np.ndarray:
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        if (
            line.startswith("ATOM")
            and len(line) >= 54
            and line[21] == chain
            and line[12:16].strip() == "CA"
            and line[16] in (" ", "A")
        ):
            rows.append(
                [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ]
            )
    return np.asarray(rows, dtype=float)


def kabsch_row(P: np.ndarray, Q: np.ndarray):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt
    t = qc - pc @ R
    fit = P @ R + t
    rmsd = float(np.sqrt(np.mean(np.sum((fit - Q) ** 2, axis=1))))
    return R, t, rmsd


def load_cm() -> Molecule:
    rd = Chem.SDMolSupplier(str(POSE), removeHs=False)[0]
    if rd is None:
        raise RuntimeError(f"Could not read {POSE}")
    Chem.SanitizeMol(rd)

    n_h = sum(1 for atom in rd.GetAtoms() if atom.GetAtomicNum() == 1)
    if n_h == 0:
        raise RuntimeError(
            "CM00734 pose SDF has no explicit hydrogens; "
            "refusing silent hydrogen reconstruction."
        )

    rd = Chem.Mol(rd)
    for atom in rd.GetAtoms():
        atom.SetMonomerInfo(None)

    mol = Molecule.from_rdkit(
        rd,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True,
    )
    mol.name = "CM00734"

    if len(mol.conformers) != 1:
        raise RuntimeError(
            f"Expected one CM00734 conformer, got {len(mol.conformers)}"
        )
    return mol


def forcefield(ligands):
    ff = app.ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/lipid17.xml",
        "amber14/tip3p.xml",
    )
    for ligand in ligands:
        generator = SMIRNOFFTemplateGenerator(
            molecules=[ligand],
            forcefield="openff_unconstrained-2.2.1.offxml",
        )
        ff.registerTemplateGenerator(generator.generator)
    return ff


def remove_cm_solvent_clashes(
    modeller,
    cm_atom_count: int,
    cutoff_nm: float = 0.25,
) -> int:
    atoms = list(modeller.topology.atoms())
    xyz = np.asarray(
        modeller.positions.value_in_unit(unit.nanometer),
        dtype=float,
    )
    cm_start = len(atoms) - cm_atom_count

    tree = cKDTree(xyz[:cm_start])
    clashing_indices = set()

    for point in xyz[cm_start:]:
        clashing_indices.update(tree.query_ball_point(point, cutoff_nm))

    removable = {
        atoms[i].residue
        for i in clashing_indices
        if atoms[i].residue.name in {"HOH", "NA", "CL"}
    }

    if removable:
        modeller.delete(list(removable))

    return len(removable)


def add_protein_restraints(
    system,
    topology,
    positions,
    k: float = 1000.0,
) -> int:
    restraint = CustomExternalForce(
        "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
    )
    restraint.addGlobalParameter(
        "k",
        k * unit.kilojoule_per_mole / unit.nanometer**2,
    )

    for name in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(name)

    xyz = positions.value_in_unit(unit.nanometer)
    count = 0

    for atom in topology.atoms():
        if atom.residue.chain.id not in PROTEIN_CHAINS:
            continue
        if atom.element is not None and atom.element.symbol == "H":
            continue
        restraint.addParticle(
            atom.index,
            [float(v) for v in xyz[atom.index]],
        )
        count += 1

    system.addForce(restraint)
    return count


def chain_atom_count(topology, chain_id: str) -> int:
    return sum(
        1
        for atom in topology.atoms()
        if atom.residue.chain.id == chain_id
    )


def chain_bond_count(topology, chain_id: str) -> int:
    ids = {
        atom.index
        for atom in topology.atoms()
        if atom.residue.chain.id == chain_id
    }
    return sum(
        1
        for a, b in topology.bonds()
        if a.index in ids and b.index in ids
    )


def build_one(
    name: str,
    spec: dict,
    cm_template: Molecule,
    raw_ca: np.ndarray,
) -> dict:
    source = spec["source"]
    pdb = app.PDBFile(str(source))

    mem_ca = ca_xyz(source, "E")
    if raw_ca.shape != (270, 3) or mem_ca.shape != (270, 3):
        raise RuntimeError(
            f"{name}: receptor CA count mismatch "
            f"{raw_ca.shape} vs {mem_ca.shape}"
        )

    R, t, rmsd = kabsch_row(raw_ca, mem_ca)
    mapped_center = CENTER @ R + t

    observed_drop_atoms = chain_atom_count(
        pdb.topology,
        spec["drop_chain"],
    )
    if observed_drop_atoms != 53:
        raise RuntimeError(
            f"{name}: expected 53 LY atoms in chain "
            f"{spec['drop_chain']}, got {observed_drop_atoms}"
        )

    retained_ach = None
    if spec["probe"]:
        ach_atoms = chain_atom_count(pdb.topology, "J")
        ach_bonds = chain_bond_count(pdb.topology, "J")
        if ach_atoms != 26 or ach_bonds == 0:
            raise RuntimeError(
                f"{name}: retained ACh chain J unexpected: "
                f"atoms={ach_atoms}, bonds={ach_bonds}"
            )

        retained_ach = Molecule.from_smiles(ACH_SMILES)
        retained_ach.name = "ACH"
        if retained_ach.n_atoms != 26:
            raise RuntimeError(
                f"ACh template atom count unexpected: "
                f"{retained_ach.n_atoms}"
            )

    modeller = app.Modeller(pdb.topology, pdb.positions)

    drop_residues = [
        residue
        for residue in modeller.topology.residues()
        if residue.chain.id == spec["drop_chain"]
    ]
    if not drop_residues:
        raise RuntimeError(
            f"{name}: source LY chain {spec['drop_chain']} missing"
        )

    modeller.delete(drop_residues)

    cm = copy.deepcopy(cm_template)
    xyz = cm.conformers[0].m_as(off_unit.angstrom)
    xyz = xyz @ R + t
    cm.conformers[0] = xyz * off_unit.angstrom

    modeller.add(
        cm.to_topology().to_openmm(),
        cm.conformers[0].to_openmm(),
    )

    removed_solvent = remove_cm_solvent_clashes(
        modeller,
        cm.n_atoms,
    )

    ligands = [cm]
    if retained_ach is not None:
        ligands.insert(0, retained_ach)

    ff = forcefield(ligands)

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=5e-4,
    )

    restrained = add_protein_restraints(
        system,
        modeller.topology,
        modeller.positions,
        k=1000.0,
    )

    system.addForce(
        MonteCarloMembraneBarostat(
            1.0 * unit.bar,
            0.0 * unit.bar * unit.nanometer,
            300 * unit.kelvin,
            MonteCarloMembraneBarostat.XYIsotropic,
            MonteCarloMembraneBarostat.ZFree,
            25,
        )
    )

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    simulation = app.Simulation(
        modeller.topology,
        system,
        integrator,
        Platform.getPlatformByName("CPU"),
    )

    simulation.context.setPositions(modeller.positions)

    initial = simulation.context.getState(getEnergy=True)
    initial_energy = initial.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )

    simulation.minimizeEnergy(maxIterations=100)

    final = simulation.context.getState(
        getEnergy=True,
        getPositions=True,
    )
    final_energy = final.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )

    output = OUT / name
    output.mkdir(parents=True, exist_ok=True)

    with (output / "input_with_cm.pdb").open("w") as handle:
        app.PDBFile.writeFile(
            modeller.topology,
            modeller.positions,
            handle,
        )

    with (output / "minimized.pdb").open("w") as handle:
        app.PDBFile.writeFile(
            modeller.topology,
            final.getPositions(),
            handle,
        )

    (output / "system.xml").write_text(
        XmlSerializer.serialize(system)
    )
    (output / "state.xml").write_text(
        XmlSerializer.serialize(final)
    )

    box = final.getPeriodicBoxVectors(
        asNumpy=True
    ).value_in_unit(unit.nanometer)

    return {
        "system": name,
        "source_geometry": str(source),
        "source_geometry_sha256": sha256(source),
        "source_removed_chain": spec["drop_chain"],
        "probe_retained": bool(spec["probe"]),
        "raw_to_membrane_ca_rmsd_A": rmsd,
        "mapped_frozen_center_A": mapped_center.tolist(),
        "cm_atoms": cm.n_atoms,
        "cm_explicit_hydrogens": sum(
            1 for atom in cm.atoms if atom.atomic_number == 1
        ),
        "removed_clashing_solvent_or_ions": removed_solvent,
        "protein_heavy_atoms_restrained": restrained,
        "particles": system.getNumParticles(),
        "initial_energy_kj_mol": initial_energy,
        "minimized_energy_kj_mol": final_energy,
        "box_nm": np.asarray(box).tolist(),
        "files": {
            "minimized_pdb": str(output / "minimized.pdb"),
            "system_xml": str(output / "system.xml"),
            "state_xml": str(output / "state.xml"),
            "minimized_pdb_sha256": sha256(
                output / "minimized.pdb"
            ),
            "system_xml_sha256": sha256(
                output / "system.xml"
            ),
            "state_xml_sha256": sha256(
                output / "state.xml"
            ),
        },
    }


def main():
    for path in (RAW, POSE):
        if not path.exists():
            raise SystemExit(
                f"Missing required input: {path}"
            )

    for spec in DEFS.values():
        if not spec["source"].exists():
            raise SystemExit(
                f"Missing historical membrane geometry: "
                f"{spec['source']}"
            )

    raw_ca = ca_xyz(RAW, "R")
    cm = load_cm()

    OUT.mkdir(parents=True, exist_ok=True)

    records = []

    for name, spec in DEFS.items():
        print(f"BUILDING {name}", flush=True)
        records.append(
            build_one(
                name,
                spec,
                cm,
                raw_ca,
            )
        )

    audit = {
        "schema": "pacer_dc.cm00734.stage_b_membrane.v1",
        "status": (
            "CM00734_STAGE_B_MEMBRANE_REFERENCE_COMPLETE"
        ),
        "cm_pose_sdf": str(POSE),
        "cm_pose_sdf_sha256": sha256(POSE),
        "force_fields": [
            "amber14/protein.ff14SB.xml",
            "amber14/lipid17.xml",
            "amber14/tip3p.xml",
            "openff_unconstrained-2.2.1.offxml",
        ],
        "nonbonded": {
            "method": "PME",
            "cutoff_nm": 1.0,
            "ewald_error_tolerance": 5e-4,
        },
        "constraints": "HBonds",
        "rigid_water": True,
        "protein_restraint_k_kj_mol_nm2": 1000.0,
        "membrane_barostat": {
            "pressure_bar": 1.0,
            "surface_tension_bar_nm": 0.0,
            "temperature_K": 300,
            "frequency": 25,
            "mode": "XYIsotropic/ZFree",
        },
        "minimization_max_iterations": 100,
        "systems": records,
        "claim_boundary": (
            "Historical membrane geometry reuse plus new "
            "CM00734 SMIRNOFF parameterization and restrained "
            "minimization QC only; not production dynamics or "
            "efficacy evidence."
        ),
    }

    audit_path = (
        OUT / "CM00734_STAGE_B_MEMBRANE_AUDIT_v01.json"
    )
    audit_path.write_text(
        json.dumps(audit, indent=2)
    )

    print(
        json.dumps(
            {
                "status": audit["status"],
                "audit": str(audit_path),
                "systems": [
                    {
                        "system": row["system"],
                        "ca_rmsd_A": (
                            row["raw_to_membrane_ca_rmsd_A"]
                        ),
                        "removed_solvent": (
                            row[
                                "removed_clashing_solvent_or_ions"
                            ]
                        ),
                        "initial_energy_kj_mol": (
                            row["initial_energy_kj_mol"]
                        ),
                        "minimized_energy_kj_mol": (
                            row["minimized_energy_kj_mol"]
                        ),
                        "particles": row["particles"],
                    }
                    for row in records
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
