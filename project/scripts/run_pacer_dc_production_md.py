#!/usr/bin/env python
"""Checkpointed unrestrained PACER-DC production MD for CPU/CUDA/OpenCL.

Each invocation is one independent system/replica.  The script starts from the
audited zero-restraint state, explicitly sets restraint k=0, and writes an
atomic progress ledger.  Short runs are sampling pilots, not PAM evidence.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

from openmm import LangevinMiddleIntegrator, MonteCarloMembraneBarostat, Platform, XmlSerializer, app, unit


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "results" / "pacer_dc_membrane_reference_v01"
START = ROOT / "results" / "pacer_dc_restraint_release_v01"
OUT = ROOT / "results" / "pacer_dc_production_v01"
VALID_SYSTEMS = (
    "apo",
    "probe_only",
    "LY2119620__candidate_no_probe",
    "LY2119620__candidate_probe",
    "compound110__candidate_no_probe",
    "compound110__candidate_probe",
)
SEEDS = {1: 27101, 2: 38201, 3: 49301}


def atomic_json(path: Path, payload: dict) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(temp, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("system", choices=VALID_SYSTEMS)
    parser.add_argument("--replica", type=int, choices=(1, 2, 3), required=True)
    parser.add_argument("--ns", type=float, default=5.0)
    parser.add_argument("--platform", choices=("CUDA", "OpenCL", "CPU"), default="CUDA")
    parser.add_argument("--device-index", default="0")
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--trajectory-ps", type=float, default=10.0)
    parser.add_argument("--checkpoint-ps", type=float, default=100.0)
    parser.add_argument("--output-root", type=Path, default=OUT)
    args = parser.parse_args()
    if args.ns <= 0 or args.timestep_fs <= 0:
        raise SystemExit("ns and timestep-fs must be positive")

    source = REFERENCE / args.system
    start_dir = START / args.system / "replica_01"
    start_audit = json.loads((start_dir / "audit.json").read_text(encoding="utf-8"))
    if not start_audit.get("production_start_eligible", False):
        raise RuntimeError("Audited zero-restraint production start is absent.")
    out = args.output_root / args.system / f"replica_{args.replica:02d}"
    out.mkdir(parents=True, exist_ok=True)
    ledger_path = out / "progress.json"
    checkpoint_path = out / "checkpoint.chk"

    pdb = app.PDBFile(str(source / "minimized.pdb"))
    system = XmlSerializer.deserialize((source / "system.xml").read_text(encoding="utf-8"))
    bars = [force for force in system.getForces() if isinstance(force, MonteCarloMembraneBarostat)]
    if len(bars) != 1:
        raise RuntimeError(f"Expected one membrane barostat, found {len(bars)}")
    bars[0].setFrequency(25)
    timestep_ps = args.timestep_fs / 1000.0
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, timestep_ps * unit.picoseconds
    )
    seed = SEEDS[args.replica]
    integrator.setRandomNumberSeed(seed)
    platform = Platform.getPlatformByName(args.platform)
    properties = {}
    if args.platform in ("CUDA", "OpenCL"):
        properties["DeviceIndex"] = args.device_index
        properties["Precision"] = "mixed"
    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
    try:
        simulation.context.setParameter("k", 0.0)
    except Exception as exc:
        raise RuntimeError("Protein restraint parameter k is missing; refusing production.") from exc

    resumed = checkpoint_path.exists() or (out / "latest_state.xml").exists()
    resume_mode = "new"
    if checkpoint_path.exists():
        try:
            simulation.loadCheckpoint(str(checkpoint_path))
            resume_mode = "checkpoint"
        except Exception:
            # Binary checkpoints can be hardware/platform specific.  The XML state
            # and progress ledger provide a portable fallback across GPU nodes.
            portable = XmlSerializer.deserialize((out / "latest_state.xml").read_text(encoding="utf-8"))
            simulation.context.setState(portable)
            prior = json.loads(ledger_path.read_text(encoding="utf-8"))
            simulation.currentStep = round(float(prior["completed_ns"]) * 1000.0 / timestep_ps)
            resume_mode = "portable_xml"
    elif (out / "latest_state.xml").exists():
        portable = XmlSerializer.deserialize((out / "latest_state.xml").read_text(encoding="utf-8"))
        simulation.context.setState(portable)
        prior = json.loads(ledger_path.read_text(encoding="utf-8"))
        simulation.currentStep = round(float(prior["completed_ns"]) * 1000.0 / timestep_ps)
        resume_mode = "portable_xml"
    else:
        state = XmlSerializer.deserialize(
            (start_dir / "production_start_state.xml").read_text(encoding="utf-8")
        )
        simulation.context.setState(state)
        simulation.currentStep = 0
        simulation.context.setVelocitiesToTemperature(300 * unit.kelvin, seed)

    target_steps = round(args.ns / timestep_ps * 1000.0)
    current_steps = int(simulation.currentStep)
    if current_steps >= target_steps:
        print(json.dumps({"status": "already_complete", "steps": current_steps}))
        return
    trajectory_steps = max(1, round(args.trajectory_ps / timestep_ps))
    checkpoint_steps = max(1, round(args.checkpoint_ps / timestep_ps))
    append = resumed and (out / "trajectory.dcd").exists()
    simulation.reporters.append(app.DCDReporter(
        str(out / "trajectory.dcd"), trajectory_steps, append=append, enforcePeriodicBox=False
    ))
    simulation.reporters.append(app.StateDataReporter(
        str(out / "state.csv"), trajectory_steps, append=append, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, temperature=True, density=True,
        volume=True, speed=True, remainingTime=True, totalSteps=target_steps, separator=",",
    ))
    started = time.time()
    while simulation.currentStep < target_steps:
        block = min(checkpoint_steps, target_steps - simulation.currentStep)
        simulation.step(block)
        simulation.saveCheckpoint(str(checkpoint_path))
        state = simulation.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
        (out / "latest_state.xml").write_text(XmlSerializer.serialize(state), encoding="utf-8")
        atomic_json(ledger_path, {
            "status": "running",
            "system": args.system,
            "replica": args.replica,
            "seed": seed,
            "platform": args.platform,
            "device_index": args.device_index,
            "target_ns": args.ns,
            "completed_ns": simulation.currentStep * timestep_ps / 1000.0,
            "timestep_fs": args.timestep_fs,
            "restraint_k_kj_mol_nm2": 0.0,
            "resumed": resumed,
            "resume_mode": resume_mode,
            "wall_seconds_this_invocation": time.time() - started,
            "claim_boundary": "Production sampling progress only; no convergence or PAM inference.",
        })
    final = json.loads(ledger_path.read_text(encoding="utf-8"))
    final["status"] = "complete"
    final["completed_ns"] = args.ns
    final["claim_boundary"] = (
        "Completed trajectory is eligible for QC. Biological inference still requires convergence, "
        "three independent replicas, matched contexts, and frozen comparisons."
    )
    atomic_json(ledger_path, final)
    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    main()
