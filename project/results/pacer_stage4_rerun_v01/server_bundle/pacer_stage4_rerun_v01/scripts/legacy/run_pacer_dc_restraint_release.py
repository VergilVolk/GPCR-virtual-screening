#!/usr/bin/env python
"""Staged protein-heavy-atom restraint release for PACER-DC membrane systems."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from openmm import LangevinMiddleIntegrator, MonteCarloMembraneBarostat, Platform, XmlSerializer, app, unit


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "results" / "pacer_dc_membrane_reference_v01"
EQUIL = ROOT / "results" / "pacer_dc_short_equilibration_v01"
OUT = ROOT / "results" / "pacer_dc_restraint_release_v01"


def set_barostat(system, frequency: int) -> None:
    bars = [f for f in system.getForces() if isinstance(f, MonteCarloMembraneBarostat)]
    if len(bars) != 1:
        raise RuntimeError(f"Expected one membrane barostat, found {len(bars)}")
    bars[0].setFrequency(frequency)


def finite_state(state) -> bool:
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    return bool(math.isfinite(pe) and math.isfinite(ke) and np.isfinite(xyz).all())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("system")
    parser.add_argument("--replica", type=int, default=1)
    parser.add_argument("--seed", type=int, default=27101)
    parser.add_argument("--stage-ps", type=float, default=0.5)
    parser.add_argument("--timestep-fs", type=float, default=1.0)
    parser.add_argument("--schedule", default="500,100,10,0")
    args = parser.parse_args()
    schedule = [float(x) for x in args.schedule.split(",")]
    if not schedule or schedule[-1] != 0 or any(a <= b for a, b in zip(schedule, schedule[1:])):
        raise SystemExit("Schedule must be strictly decreasing and end at zero.")

    source = REFERENCE / args.system
    prior = EQUIL / args.system / f"replica_{args.replica:02d}"
    prior_audit = json.loads((prior / "audit.json").read_text(encoding="utf-8"))
    if not prior_audit.get("passed", False):
        raise RuntimeError("Refusing restraint release: short-equilibration gate has not passed.")
    state = XmlSerializer.deserialize((prior / "npt_final_state.xml").read_text(encoding="utf-8"))
    pdb = app.PDBFile(str(source / "minimized.pdb"))
    out = OUT / args.system / f"replica_{args.replica:02d}"
    out.mkdir(parents=True, exist_ok=True)
    timestep_ps = args.timestep_fs / 1000.0
    steps = round(args.stage_ps / timestep_ps)
    report_steps = max(1, steps // 10)
    volume_initial = float(abs(np.linalg.det(
        np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    )))
    stages = []

    for stage_index, k in enumerate(schedule, start=1):
        system = XmlSerializer.deserialize((source / "system.xml").read_text(encoding="utf-8"))
        set_barostat(system, 25)
        integrator = LangevinMiddleIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, timestep_ps * unit.picoseconds
        )
        integrator.setRandomNumberSeed(args.seed + stage_index * 1_000_003)
        simulation = app.Simulation(pdb.topology, system, integrator, Platform.getPlatformByName("CPU"))
        simulation.context.setState(state)
        try:
            simulation.context.setParameter("k", k)
        except Exception as exc:
            raise RuntimeError("Protein-restraint global parameter 'k' was not found.") from exc
        label = f"stage_{stage_index:02d}_k_{k:g}"
        simulation.reporters.append(app.StateDataReporter(
            str(out / f"{label}.csv"), report_steps, step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, temperature=True, density=True,
            volume=True, speed=True, separator=",",
        ))
        simulation.reporters.append(app.DCDReporter(str(out / f"{label}.dcd"), report_steps))
        simulation.step(steps)
        state = simulation.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
        finite = finite_state(state)
        box = np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
        volume = float(abs(np.linalg.det(box)))
        stage = {
            "stage": stage_index,
            "restraint_k_kj_mol_nm2": k,
            "simulated_ps": args.stage_ps,
            "potential_kj_mol": state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
            "volume_nm3": volume,
            "volume_change_from_start_percent": (volume / volume_initial - 1.0) * 100.0,
            "finite": finite,
        }
        stages.append(stage)
        (out / f"{label}_state.xml").write_text(XmlSerializer.serialize(state), encoding="utf-8")
        if not finite or abs(stage["volume_change_from_start_percent"]) >= 20.0:
            raise RuntimeError(f"Restraint-release gate failed at {label}: {stage}")

    final_csv = out / f"stage_{len(schedule):02d}_k_0.csv"
    import csv
    with final_csv.open(newline="", encoding="utf-8") as handle:
        records = list(csv.DictReader(handle))
    final_temperature = float(records[-1]["Temperature (K)"])
    audit = {
        "system": args.system,
        "replica": args.replica,
        "matched_seed": args.seed,
        "schedule_k_kj_mol_nm2": schedule,
        "stage_ps": args.stage_ps,
        "total_ps": args.stage_ps * len(schedule),
        "timestep_fs": args.timestep_fs,
        "final_temperature_K": final_temperature,
        "stages": stages,
        "passed": bool(270.0 <= final_temperature <= 330.0 and all(x["finite"] for x in stages)),
        "production_start_eligible": bool(
            270.0 <= final_temperature <= 330.0 and all(x["finite"] for x in stages)
        ),
        "claim_boundary": (
            "Passing creates an unrestrained production starting state. It does not demonstrate "
            "equilibrium convergence or any PAM-related biological effect."
        ),
    }
    (out / "production_start_state.xml").write_text(XmlSerializer.serialize(state), encoding="utf-8")
    (out / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps({
        "system": args.system,
        "passed": audit["passed"],
        "final_temperature_K": final_temperature,
        "total_ps": audit["total_ps"],
    }, indent=2))
    if not audit["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
