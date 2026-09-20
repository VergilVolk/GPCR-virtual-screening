#!/usr/bin/env python
"""Build a shared OPM-oriented POPC/water box for PACER-DC reference systems."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from openmm import CustomExternalForce, LangevinMiddleIntegrator, MonteCarloMembraneBarostat
from openmm import Platform, XmlSerializer, Vec3, app, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.units import unit as off_unit

import build_pacer_dc_openmm_reference_systems as qc


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "pacer_dc_membrane_reference_v01"
OPM = ROOT / "data" / "pdb" / "7TRS_OPM.pdb"
RAW = ROOT / "data" / "pdb" / "7TRS.pdb"
MINIMIZED_PROTEIN = ROOT / "results" / "pacer_dc_openmm_reference_qc_v01" / "apo" / "minimized.pdb"


def ca_map(path: Path):
    result = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ATOM") and line[21:22] == "R" and line[12:16].strip() == "CA":
            key = (line[21:22], int(line[22:26]))
            result[key] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return result


def ca_order(path: Path, chain: str):
    rows = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ATOM") and line[21:22] == chain and line[12:16].strip() == "CA":
            rows.append((int(line[22:26]), np.array([
                float(line[30:38]), float(line[38:46]), float(line[46:54])
            ])))
    return rows


def opm_transform(source_path: Path = MINIMIZED_PROTEIN):
    source_order = ca_order(source_path, "E")
    common_order = ca_order(qc.PROTEIN, "R")
    target = ca_map(OPM)
    if len(source_order) != len(common_order):
        raise RuntimeError(f"receptor CA mismatch: {len(source_order)} vs {len(common_order)}")
    pairs = [(src_xyz, target[("R", common_resid)])
             for (_, src_xyz), (common_resid, _) in zip(source_order, common_order)
             if ("R", common_resid) in target]
    p, q = np.vstack([x[0] for x in pairs]), np.vstack([x[1] for x in pairs])
    pc, qc_ = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc_))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = u @ vt
    fitted = (p - pc) @ rotation + qc_
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - q) ** 2, axis=1))))
    return rotation, pc, qc_, len(pairs), rmsd


def transform_positions(positions, rotation, source_center, target_center):
    xyz = np.array([[v.x, v.y, v.z] for v in positions.value_in_unit(unit.angstrom)])
    xyz = (xyz - source_center) @ rotation + target_center
    return unit.Quantity([Vec3(*row) for row in xyz], unit.angstrom)


def transform_ligand(molecule, rotation, source_center, target_center):
    xyz = molecule.conformers[0].m_as(off_unit.angstrom)
    transformed = (xyz - source_center) @ rotation + target_center
    molecule.conformers[0] = transformed * off_unit.angstrom
    return molecule


def align_ligand_to_base(molecule, source_complex: Path):
    source = ca_order(source_complex, "E")
    target = ca_order(MINIMIZED_PROTEIN, "E")
    if len(source) != len(target):
        raise RuntimeError("source/base receptor CA mismatch")
    p = np.vstack([x[1] for x in source])
    q = np.vstack([x[1] for x in target])
    pc, qc_ = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc_))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = u @ vt
    return transform_ligand(molecule, rotation, pc, qc_)


def replace_conformer_from_chain(molecule, pdb_path: Path, chain: str = "F"):
    xyz = []
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("HETATM") and line[21:22] == chain:
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if len(xyz) != molecule.n_atoms:
        raise RuntimeError(f"ligand atom count mismatch: {len(xyz)} vs {molecule.n_atoms}")
    molecule.conformers[0] = np.asarray(xyz) * off_unit.angstrom
    return molecule


def forcefield(ligands=()):
    ff = app.ForceField("amber14/protein.ff14SB.xml", "amber14/lipid17.xml", "amber14/tip3p.xml")
    for ligand in ligands:
        generator = SMIRNOFFTemplateGenerator(
            molecules=[ligand], forcefield="openff_unconstrained-2.2.1.offxml"
        )
        ff.registerTemplateGenerator(generator.generator)
    return ff


def restrain_protein(system, topology, positions, count, k=1000.0):
    restraint = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", k * unit.kilojoule_per_mole / unit.nanometer**2)
    for name in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(name)
    xyz = positions.value_in_unit(unit.nanometer)
    atoms = list(topology.atoms())
    for i in range(count):
        if atoms[i].element is not None and atoms[i].element.symbol != "H":
            restraint.addParticle(i, [float(v) for v in xyz[i]])
    system.addForce(restraint)


def remove_ligand_solvent_clashes(modeller, ligand_atom_count: int, cutoff_nm: float = 0.25):
    if ligand_atom_count == 0:
        return 0
    atoms = list(modeller.topology.atoms())
    xyz = np.asarray(modeller.positions.value_in_unit(unit.nanometer), dtype=float)
    ligand_start = len(atoms) - ligand_atom_count
    tree = cKDTree(xyz[:ligand_start])
    clashing_indices = set()
    for point in xyz[ligand_start:]:
        clashing_indices.update(tree.query_ball_point(point, cutoff_nm))
    removable = {atoms[i].residue for i in clashing_indices
                 if atoms[i].residue.name in {"HOH", "NA", "CL"}}
    if removable:
        modeller.delete(list(removable))
    return len(removable)


def hydrogen_geometry_qc(modeller):
    atoms = list(modeller.topology.atoms())
    xyz = np.asarray(modeller.positions.value_in_unit(unit.nanometer), dtype=float)
    bonded = {frozenset((a.index, b.index)) for a, b in modeller.topology.bonds()}
    worst = (float("inf"), None)
    for atom in atoms:
        if atom.element is None or atom.element.symbol != "H":
            continue
        for other in atom.residue.atoms():
            if other.element is None or other.element.symbol == "H":
                continue
            if frozenset((atom.index, other.index)) in bonded:
                continue
            distance = float(np.linalg.norm(xyz[atom.index] - xyz[other.index]))
            if distance < worst[0]:
                worst = (distance, (atom, other))
    distance, pair = worst
    label = None if pair is None else (
        f"{pair[0].residue.chain.id}:{pair[0].residue.name}{pair[0].residue.id}:"
        f"{pair[0].name}-{pair[1].name}"
    )
    return {"min_nonbonded_H_heavy_nm": distance, "closest_pair": label,
            "passed": bool(distance >= 0.06)}


def rebuild_protein_hydrogens(pdb, oriented, ff):
    attempts = []
    for seed in range(1701, 1721):
        random.seed(seed)
        np.random.seed(seed)
        modeller = app.Modeller(pdb.topology, oriented)
        modeller.topology.setPeriodicBoxVectors(None)
        modeller.delete([atom for atom in modeller.topology.atoms()
                         if atom.element is not None and atom.element.symbol == "H"])
        modeller.addHydrogens(ff, pH=7.4)
        result = {"seed": seed, **hydrogen_geometry_qc(modeller)}
        attempts.append(result)
        if result["passed"]:
            return modeller, result, attempts
    raise RuntimeError(f"no hydrogen placement passed QC: {attempts}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--systems", nargs="*", default=["probe_only", "apo"])
    parser.add_argument("--outdir", type=Path, default=OUT)
    parser.add_argument("--hydrogen-qc-only", action="store_true")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rotation, source_center, target_center, n_ca, alignment_rmsd = opm_transform()
    pdb = app.PDBFile(str(MINIMIZED_PROTEIN))
    oriented = transform_positions(pdb.positions, rotation, source_center, target_center)
    base, hydrogen_qc, hydrogen_attempts = rebuild_protein_hydrogens(pdb, oriented, forcefield())
    if args.hydrogen_qc_only:
        print(json.dumps({"selected": hydrogen_qc, "attempts": hydrogen_attempts}, indent=2))
        return
    protein_atoms = sum(1 for _ in base.topology.atoms())
    base.addMembrane(
        forcefield(), lipidType="POPC", minimumPadding=1.0 * unit.nanometer,
        ionicStrength=0.15 * unit.molar, neutralize=True,
    )
    base_pdb = args.outdir / "shared_OPM_POPC_water.pdb"
    with base_pdb.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(base.topology, base.positions, handle)

    qc_root = ROOT / "results" / "pacer_dc_openmm_reference_qc_v01"

    # Ligands are resolved lazily so a missing pose fails only the systems that
    # need it.  Loading them up front aborted the whole run -- including systems
    # whose inputs were present -- and hid which file was actually absent.
    def load_ach():
        source = qc_root / "probe_only" / "minimized.pdb"
        posed = qc.openff_molecule(qc.posed_from_pdb(
            ROOT / "data" / "pdb" / "7TRS.pdb", qc.LIGAND_DIR / "ACH_ideal.sdf",
            residue_name="ACH"
        ), "ACH")
        posed = align_ligand_to_base(replace_conformer_from_chain(posed, source), source)
        return transform_ligand(posed, rotation, source_center, target_center)

    def load_ly2119620():
        source = qc_root / "LY2119620__candidate_no_probe" / "minimized.pdb"
        posed = qc.openff_molecule(qc.posed_from_pdb(
            qc.COMPLEX_DIR / "LY2119620__candidate_probe.pdb", qc.LIGAND_DIR / "2CU_ideal.sdf",
            chain="L"
        ), "LY2119620")
        posed = align_ligand_to_base(replace_conformer_from_chain(posed, source), source)
        return transform_ligand(posed, rotation, source_center, target_center)

    def load_compound110():
        # Produced by the 7TRS pocket re-docking step, not by
        # build_pacer_dc_reference_complexes.py.
        source = qc_root / "compound110__candidate_no_probe" / "minimized.pdb"
        posed = qc.openff_molecule(qc.posed_from_sdf(
            qc.COMPLEX_DIR / "compound110_7TRS_redocked.sdf"
        ), "compound110")
        posed = align_ligand_to_base(replace_conformer_from_chain(posed, source), source)
        return transform_ligand(posed, rotation, source_center, target_center)

    definitions = {
        "probe_only": [load_ach], "apo": [],
        "LY2119620__candidate_probe": [load_ach, load_ly2119620],
        "LY2119620__candidate_no_probe": [load_ly2119620],
        "compound110__candidate_probe": [load_ach, load_compound110],
        "compound110__candidate_no_probe": [load_compound110],
    }
    unknown = set(args.systems) - set(definitions)
    if unknown:
        raise ValueError(f"unknown systems: {sorted(unknown)}")

    records = []
    for name in args.systems:
        ligands = [loader() for loader in definitions[name]]
        modeller = app.Modeller(base.topology, base.positions)
        for ligand in ligands:
            modeller.add(ligand.to_topology().to_openmm(), ligand.conformers[0].to_openmm())
        removed_solvent = remove_ligand_solvent_clashes(
            modeller, sum(ligand.n_atoms for ligand in ligands)
        )
        ff = forcefield(ligands)
        system = ff.createSystem(
            modeller.topology, nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds,
            rigidWater=True, ewaldErrorTolerance=5e-4,
        )
        restrain_protein(system, modeller.topology, modeller.positions, protein_atoms)
        system.addForce(MonteCarloMembraneBarostat(
            1.0 * unit.bar, 0.0 * unit.bar * unit.nanometer,
            300 * unit.kelvin, MonteCarloMembraneBarostat.XYIsotropic,
            MonteCarloMembraneBarostat.ZFree, 25,
        ))
        integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        simulation = app.Simulation(modeller.topology, system, integrator, Platform.getPlatformByName("CPU"))
        simulation.context.setPositions(modeller.positions)
        initial = simulation.context.getState(getEnergy=True)
        simulation.minimizeEnergy(maxIterations=100)
        final = simulation.context.getState(getEnergy=True, getPositions=True)
        system_dir = args.outdir / name
        system_dir.mkdir(exist_ok=True)
        with (system_dir / "minimized.pdb").open("w", encoding="utf-8") as handle:
            app.PDBFile.writeFile(modeller.topology, final.getPositions(), handle)
        (system_dir / "system.xml").write_text(XmlSerializer.serialize(system), encoding="utf-8")
        (system_dir / "state.xml").write_text(XmlSerializer.serialize(final), encoding="utf-8")
        records.append({
            "system": name, "ligands": [x.name for x in ligands],
            "atoms": system.getNumParticles(),
            "removed_clashing_solvent_residues": removed_solvent,
            "initial_energy_kj_mol": initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
            "minimized_energy_kj_mol": final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
        })
    audit = {
        "opm_source": str(OPM), "opm_sha256": "2f9060a0281bd5ba787fecfc24d244b695d9a86919a634c5cf12cef1c2a1aba",
        "receptor_ca_pairs": n_ca, "orientation_alignment_rmsd_A": alignment_rmsd,
        "hydrogen_geometry_qc": hydrogen_qc, "hydrogen_attempts": hydrogen_attempts,
        "shared_base_atoms": sum(1 for _ in base.topology.atoms()), "protein_atoms": protein_atoms,
        "systems": records,
        "claim_boundary": "OPM-oriented membrane setup and restrained minimization QC only; not production MD.",
    }
    (args.outdir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
