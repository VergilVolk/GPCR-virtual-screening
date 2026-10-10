#!/usr/bin/env python
"""Matched-seed NVT/NPT equilibration pilot for one PACER-DC context.

This is a protocol/QC run, not a production trajectory and not functional evidence.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from openmm import LangevinMiddleIntegrator, MonteCarloMembraneBarostat, Platform, XmlSerializer, app, unit


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "results" / "pacer_dc_membrane_reference_v01"
SMOKE = ROOT / "results" / "pacer_dc_membrane_smoke_v01"
OUT = ROOT / "results" / "pacer_dc_short_equilibration_v01"


def set_barostat_frequency(system, frequency: int) -> bool:
    found = False
    for force in system.getForces():
        if isinstance(force, MonteCarloMembraneBarostat):
            force.setFrequency(frequency)
            found = True
    return found


def effective_dof(system) -> int:
    moving = sum(
        1 for index in range(system.getNumParticles())
        if system.getParticleMass(index).value_in_unit(unit.dalton) > 0
    )
    return max(1, 3 * moving - system.getNumConstraints() - 3)


def assert_finite(state, system, stage: str) -> dict:
    potential = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
    finite = bool(
        math.isfinite(potential)
        and math.isfinite(kinetic)
        and np.isfinite(positions).all()
        and np.isfinite(box).all()
    )
    if not finite:
        raise RuntimeError(f"Non-finite state after {stage}")
    temperature = (
        2 * state.getKineticEnergy() / (effective_dof(system) * unit.MOLAR_GAS_CONSTANT_R)
    ).value_in_unit(unit.kelvin)
    return {
        "stage": stage,
        "potential_kj_mol": potential,
        "kinetic_kj_mol": kinetic,
        "temperature_K": temperature,
        "box_nm": np.asarray(box).tolist(),
        "finite": finite,
    }


def make_simulation(
    pdb, system, state, seed: int, out: Path, label: str, report_interval: int, timestep_ps: float
):
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, timestep_ps * unit.picoseconds
    )
    integrator.setRandomNumberSeed(seed)
    simulation = app.Simulation(pdb.topology, system, integrator, Platform.getPlatformByName("CPU"))
    simulation.context.setState(state)
    simulation.reporters.append(app.StateDataReporter(
        str(out / f"{label}_state.csv"), report_interval, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, temperature=True, density=True,
        volume=True, speed=True, separator=",",
    ))
    simulation.reporters.append(app.DCDReporter(str(out / f"{label}.dcd"), report_interval))
    return simulation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("system")
    parser.add_argument("--replica", type=int, default=1)
    parser.add_argument("--seed", type=int, default=27101)
    parser.add_argument("--nvt-ps", type=float, default=1.0)
    parser.add_argument("--npt-ps", type=float, default=1.0)
    parser.add_argument("--report-ps", type=float, default=0.1)
    parser.add_argument("--timestep-fs", type=float, default=0.5)
    args = parser.parse_args()
    if args.nvt_ps < 0 or args.npt_ps < 0 or args.nvt_ps + args.npt_ps <= 0:
        raise SystemExit("Equilibration durations must be non-negative with a positive total.")

    source = REFERENCE / args.system
    refined = SMOKE / args.system / "refined_state.xml"
    if not refined.exists():
        raise SystemExit(f"Missing smoke-refined state: {refined}")
    out = OUT / args.system / f"replica_{args.replica:02d}"
    out.mkdir(parents=True, exist_ok=True)

    pdb = app.PDBFile(str(source / "minimized.pdb"))
    initial_state = XmlSerializer.deserialize(refined.read_text(encoding="utf-8"))
    timestep_ps = args.timestep_fs / 1000.0
    if not (0 < args.timestep_fs <= 2.0):
        raise SystemExit("timestep-fs must be in (0, 2].")
    report_steps = max(1, round(args.report_ps / timestep_ps))
    stages = []

    nvt_system = XmlSerializer.deserialize((source / "system.xml").read_text(encoding="utf-8"))
    set_barostat_frequency(nvt_system, 0)
    nvt = make_simulation(
        pdb, nvt_system, initial_state, args.seed, out, "nvt", report_steps, timestep_ps
    )
    nvt.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)
    nvt_steps = round(args.nvt_ps / timestep_ps)
    if nvt_steps:
        nvt.step(nvt_steps)
    nvt_state = nvt.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
    stages.append(assert_finite(nvt_state, nvt_system, "NVT"))
    (out / "nvt_final_state.xml").write_text(XmlSerializer.serialize(nvt_state), encoding="utf-8")
    with (out / "nvt_final.pdb").open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(pdb.topology, nvt_state.getPositions(), handle)

    final_state = nvt_state
    has_barostat = False
    npt_steps = round(args.npt_ps / 0.0005)
    if npt_steps:
        npt_system = XmlSerializer.deserialize((source / "system.xml").read_text(encoding="utf-8"))
        has_barostat = set_barostat_frequency(npt_system, 25)
        if not has_barostat:
            raise RuntimeError("No MonteCarloMembraneBarostat found; refusing mislabeled NPT stage.")
        npt_seed = args.seed + 1_000_003
        npt = make_simulation(
            pdb, npt_system, nvt_state, npt_seed, out, "npt", report_steps, timestep_ps
        )
        npt_steps = round(args.npt_ps / timestep_ps)
        npt.step(npt_steps)
        final_state = npt.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
        stages.append(assert_finite(final_state, npt_system, "NPT"))
        (out / "npt_final_state.xml").write_text(XmlSerializer.serialize(final_state), encoding="utf-8")
        with (out / "npt_final.pdb").open("w", encoding="utf-8") as handle:
            app.PDBFile.writeFile(pdb.topology, final_state.getPositions(), handle)

    box0 = np.asarray(initial_state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    box1 = np.asarray(final_state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    volume0 = float(abs(np.linalg.det(box0)))
    volume1 = float(abs(np.linalg.det(box1)))
    volume_change = (volume1 / volume0 - 1.0) * 100.0
    audit = {
        "system": args.system,
        "replica": args.replica,
        "matched_seed": args.seed,
        "timestep_fs": args.timestep_fs,
        "nvt_ps": args.nvt_ps,
        "npt_ps": args.npt_ps,
        "membrane_barostat_present": has_barostat,
        "volume_initial_nm3": volume0,
        "volume_final_nm3": volume1,
        "volume_change_percent": volume_change,
        "stages": stages,
        "numerically_stable": bool(all(x["finite"] for x in stages) and abs(volume_change) < 20.0),
        "thermalized": bool(270.0 <= stages[-1]["temperature_K"] <= 330.0),
        "passed": bool(
            all(x["finite"] for x in stages)
            and abs(volume_change) < 20.0
            and 270.0 <= stages[-1]["temperature_K"] <= 330.0
        ),
        "claim_boundary": (
            "Short equilibration protocol validation only; not converged dynamics, independent sampling, "
            "or evidence of PAM function."
        ),
    }
    (out / "audit.json").write_text(json.dumps(audit, indent=2, default=str), encoding="utf-8")
    print(json.dumps({
        "system": args.system,
        "replica": args.replica,
        "passed": audit["passed"],
        "volume_change_percent": volume_change,
        "final_temperature_K": stages[-1]["temperature_K"],
    }, indent=2))
    if not audit["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
