#!/usr/bin/env python
"""CM00734 Stage-B checkpointed unrestrained production MD."""
from __future__ import annotations
import argparse, json, os, time
from pathlib import Path
from openmm import LangevinMiddleIntegrator, MonteCarloMembraneBarostat, Platform, XmlSerializer, app, unit

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "results" / "pacer_dc_membrane_reference_v01"
START = ROOT / "results" / "pacer_dc_restraint_release_v01"
OUT = ROOT / "results" / "pacer_dc_cm00734_stage_b_20ns_v01" / "production"
VALID_SYSTEMS = ("CM00734__candidate_no_probe","CM00734__candidate_probe")
SEEDS = {1:27101,2:38201,3:49301}

def atomic_json(path,payload):
    tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(payload,indent=2),encoding="utf-8")
    os.replace(tmp,path)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("system",choices=VALID_SYSTEMS)
    ap.add_argument("--replica",type=int,choices=(1,2,3),required=True)
    ap.add_argument("--ns",type=float,default=20.0)
    ap.add_argument("--platform",choices=("CUDA","OpenCL","CPU"),default="CUDA")
    ap.add_argument("--device-index",default="0")
    ap.add_argument("--timestep-fs",type=float,default=2.0)
    ap.add_argument("--trajectory-ps",type=float,default=50.0)
    ap.add_argument("--checkpoint-ps",type=float,default=100.0)
    ap.add_argument("--output-root",type=Path,default=OUT)
    a=ap.parse_args()
    if a.ns<=0 or a.timestep_fs<=0: raise SystemExit("ns and timestep-fs must be positive")

    source=REFERENCE/a.system
    start_dir=START/a.system/f"replica_{a.replica:02d}"
    start_audit_path=start_dir/"audit.json"
    start_state_path=start_dir/"production_start_state.xml"
    for p in (source/"minimized.pdb",source/"system.xml",start_audit_path,start_state_path):
        if not p.exists(): raise SystemExit(f"Missing required input: {p}")

    audit=json.loads(start_audit_path.read_text(encoding="utf-8"))
    if not audit.get("production_start_eligible",False):
        raise RuntimeError("Audited zero-restraint production start is absent")
    seed=SEEDS[a.replica]
    if int(audit.get("matched_seed",-1))!=seed:
        raise RuntimeError(f"Matched-seed mismatch: audit={audit.get('matched_seed')} expected={seed}")

    out=a.output_root/a.system/f"replica_{a.replica:02d}"
    out.mkdir(parents=True,exist_ok=True)
    ledger=out/"progress.json"; chk=out/"checkpoint.chk"; latest=out/"latest_state.xml"

    pdb=app.PDBFile(str(source/"minimized.pdb"))
    system=XmlSerializer.deserialize((source/"system.xml").read_text(encoding="utf-8"))
    bars=[f for f in system.getForces() if isinstance(f,MonteCarloMembraneBarostat)]
    if len(bars)!=1: raise RuntimeError(f"Expected one membrane barostat, found {len(bars)}")
    bars[0].setFrequency(25)

    timestep_ps=a.timestep_fs/1000.0
    integ=LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,timestep_ps*unit.picoseconds)
    integ.setRandomNumberSeed(seed)
    platform=Platform.getPlatformByName(a.platform)
    props={}
    if a.platform in ("CUDA","OpenCL"):
        props["DeviceIndex"]=a.device_index; props["Precision"]="mixed"
    sim=app.Simulation(pdb.topology,system,integ,platform,props)
    try: sim.context.setParameter("k",0.0)
    except Exception as e: raise RuntimeError("Protein restraint parameter k is missing") from e

    resumed=chk.exists() or latest.exists(); mode="new"
    if chk.exists():
        try:
            sim.loadCheckpoint(str(chk)); mode="checkpoint"
        except Exception:
            if not latest.exists() or not ledger.exists(): raise
            sim.context.setState(XmlSerializer.deserialize(latest.read_text(encoding="utf-8")))
            prior=json.loads(ledger.read_text(encoding="utf-8"))
            sim.currentStep=round(float(prior["completed_ns"])*1000.0/timestep_ps); mode="portable_xml"
    elif latest.exists():
        if not ledger.exists(): raise RuntimeError("latest_state.xml exists but progress.json is missing")
        sim.context.setState(XmlSerializer.deserialize(latest.read_text(encoding="utf-8")))
        prior=json.loads(ledger.read_text(encoding="utf-8"))
        sim.currentStep=round(float(prior["completed_ns"])*1000.0/timestep_ps); mode="portable_xml"
    else:
        sim.context.setState(XmlSerializer.deserialize(start_state_path.read_text(encoding="utf-8")))
        sim.currentStep=0
        sim.context.setVelocitiesToTemperature(300*unit.kelvin,seed)

    target=round(a.ns/timestep_ps*1000.0)
    if sim.currentStep>=target:
        print(json.dumps({"status":"already_complete","system":a.system,"replica":a.replica,"steps":int(sim.currentStep)},indent=2)); return

    traj_steps=max(1,round(a.trajectory_ps/timestep_ps))
    chk_steps=max(1,round(a.checkpoint_ps/timestep_ps))
    append=resumed and (out/"trajectory.dcd").exists()
    sim.reporters.append(app.DCDReporter(str(out/"trajectory.dcd"),traj_steps,append=append,enforcePeriodicBox=False))
    sim.reporters.append(app.StateDataReporter(str(out/"state.csv"),traj_steps,append=append,step=True,time=True,potentialEnergy=True,kineticEnergy=True,temperature=True,density=True,volume=True,speed=True,remainingTime=True,totalSteps=target,separator=","))

    started=time.time()
    while sim.currentStep<target:
        block=min(chk_steps,target-sim.currentStep)
        sim.step(block)
        sim.saveCheckpoint(str(chk))
        st=sim.context.getState(getEnergy=True,getPositions=True,getVelocities=True)
        latest.write_text(XmlSerializer.serialize(st),encoding="utf-8")
        atomic_json(ledger,{
            "status":"running","system":a.system,"replica":a.replica,"seed":seed,
            "platform":a.platform,"device_index":a.device_index,"target_ns":a.ns,
            "completed_ns":sim.currentStep*timestep_ps/1000.0,
            "timestep_fs":a.timestep_fs,"trajectory_ps":a.trajectory_ps,"checkpoint_ps":a.checkpoint_ps,
            "restraint_k_kj_mol_nm2":0.0,"resumed":resumed,"resume_mode":mode,
            "wall_seconds_this_invocation":time.time()-started,
            "claim_boundary":"Production sampling progress only; no convergence or efficacy inference."
        })

    final=json.loads(ledger.read_text(encoding="utf-8"))
    final["status"]="complete"; final["completed_ns"]=a.ns
    final["claim_boundary"]="Completed trajectory is eligible for QC and frozen FKG analysis; interpretation requires all three matched replicas and both contexts."
    atomic_json(ledger,final)
    print(json.dumps(final,indent=2))

if __name__=="__main__":
    main()
