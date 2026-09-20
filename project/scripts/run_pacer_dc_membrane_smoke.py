#!/usr/bin/env python
"""Very short membrane-MD stability smoke test; not a production trajectory."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from openmm import LangevinMiddleIntegrator, MonteCarloMembraneBarostat, Platform, XmlSerializer, app, unit


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "pacer_dc_membrane_reference_v01"
OUT = ROOT / "results" / "pacer_dc_membrane_smoke_v01"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("system")
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--seed", type=int, default=1701)
    args = parser.parse_args()
    source = SOURCE / args.system
    out = OUT / args.system
    out.mkdir(parents=True, exist_ok=True)

    source_state_path = source / "state.xml"
    source_sha256 = hashlib.sha256(source_state_path.read_bytes()).hexdigest()
    refined_state_path = out / "refined_state.xml"
    refined_hash_path = out / "refined_source_sha256.txt"
    use_cached_refinement = (
        refined_state_path.exists() and refined_hash_path.exists()
        and refined_hash_path.read_text(encoding="utf-8").strip() == source_sha256
    )

    pdb = app.PDBFile(str(source / "minimized.pdb"))
    system = XmlSerializer.deserialize((source / "system.xml").read_text(encoding="utf-8"))
    state_path = refined_state_path if use_cached_refinement else source_state_path
    state = XmlSerializer.deserialize(state_path.read_text(encoding="utf-8"))
    for force in system.getForces():
        if isinstance(force, MonteCarloMembraneBarostat):
            force.setFrequency(0)
    integrator = LangevinMiddleIntegrator(50 * unit.kelvin, 1 / unit.picosecond, 0.0005 * unit.picoseconds)
    integrator.setRandomNumberSeed(args.seed)
    simulation = app.Simulation(pdb.topology, system, integrator, Platform.getPlatformByName("CPU"))
    simulation.context.setState(state)
    raw = simulation.context.getState(getForces=True)
    raw_forces = raw.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    pre_refinement_max_force = float(np.max(np.linalg.norm(raw_forces, axis=1)))
    performed_deep_minimization = bool(not use_cached_refinement and pre_refinement_max_force > 5000.0)
    if performed_deep_minimization:
        simulation.minimizeEnergy(tolerance=10 * unit.kilojoule_per_mole / unit.nanometer, maxIterations=1000)
    pre = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
    if not use_cached_refinement:
        refined_state_path.write_text(XmlSerializer.serialize(pre), encoding="utf-8")
        refined_hash_path.write_text(source_sha256 + "\n", encoding="utf-8")
        with (out / "refined.pdb").open("w", encoding="utf-8") as handle:
            app.PDBFile.writeFile(pdb.topology, pre.getPositions(), handle)
    forces = pre.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    positions = pre.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    max_force_index = int(np.argmax(np.linalg.norm(forces, axis=1)))
    atoms = list(pdb.topology.atoms())
    atom = atoms[max_force_index]
    max_constraint_error = 0.0
    max_constraint = None
    for i in range(system.getNumConstraints()):
        a, b, target = system.getConstraintParameters(i)
        observed = float(np.linalg.norm(positions[int(a)] - positions[int(b)]))
        target_nm = target.value_in_unit(unit.nanometer)
        error = abs(observed - target_nm)
        if error > max_constraint_error:
            max_constraint_error = error
            max_constraint = [int(a), int(b), target_nm, observed]
    preflight = {
        "potential_kj_mol": pre.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
        "max_force_kj_mol_nm": float(np.linalg.norm(forces[max_force_index])),
        "max_force_atom_index": max_force_index,
        "max_force_atom": f"{atom.residue.chain.id}:{atom.residue.name}{atom.residue.id}:{atom.name}",
        "max_constraint_error_nm": max_constraint_error,
        "max_constraint": max_constraint,
        "source_state_sha256": source_sha256,
        "used_cached_refinement": use_cached_refinement,
        "pre_refinement_max_force_kj_mol_nm": pre_refinement_max_force,
        "performed_deep_minimization": performed_deep_minimization,
    }
    (out / "preflight.json").write_text(json.dumps(preflight, indent=2), encoding="utf-8")
    print(json.dumps({"preflight": preflight}, indent=2), flush=True)
    if args.steps == 0:
        return
    simulation.context.setVelocitiesToTemperature(50 * unit.kelvin, args.seed)
    initial = simulation.context.getState(getEnergy=True)
    simulation.reporters.append(app.StateDataReporter(
        str(out / "state.csv"), max(1, args.steps // 10), step=True,
        potentialEnergy=True, kineticEnergy=True, temperature=True,
        density=True, speed=True, separator=",",
    ))
    simulation.reporters.append(app.DCDReporter(str(out / "trajectory.dcd"), max(1, args.steps // 10)))
    remaining = args.steps
    for temperature in (50, 100, 200, 300):
        block = remaining if temperature == 300 else args.steps // 4
        remaining -= block if temperature != 300 else 0
        integrator.setTemperature(temperature * unit.kelvin)
        simulation.step(block)
    final = simulation.context.getState(getEnergy=True, getPositions=True)
    potential = final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    kinetic = final.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    coords = final.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    finite = bool(math.isfinite(potential) and math.isfinite(kinetic) and bool((coords == coords).all()))
    with (out / "final.pdb").open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(pdb.topology, final.getPositions(), handle)
    (out / "final_state.xml").write_text(XmlSerializer.serialize(final), encoding="utf-8")
    audit = {
        "system": args.system, "steps": args.steps, "timestep_fs": 0.5,
        "simulated_ps": args.steps * 0.0005,
        "barostat_enabled": False, "heating_K": [50, 100, 200, 300],
        "initial_potential_kj_mol": initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
        "final_potential_kj_mol": potential, "final_kinetic_kj_mol": kinetic,
        "finite": finite,
        "claim_boundary": "Numerical stability smoke test only; no equilibrium or biological inference.",
    }
    (out / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    if not finite:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
