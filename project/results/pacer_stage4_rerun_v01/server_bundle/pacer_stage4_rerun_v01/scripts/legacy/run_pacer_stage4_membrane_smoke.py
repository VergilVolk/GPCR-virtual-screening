#!/usr/bin/env python
"""CUDA-only 500-step PACER-DC membrane smoke; never equilibration/production."""
import argparse
import json
import math
import numpy as np
import pacer_stage4_prospective as s4

def finite_metrics(state,system,positions,forces):
    from openmm import unit
    potential=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    xyz=positions.value_in_unit(unit.nanometer)
    f=forces.value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
    s4.require(math.isfinite(potential) and np.isfinite(xyz).all() and np.isfinite(f).all(),
               'Nonfinite potential/coordinates/forces')
    errors=[]
    for i in range(system.getNumConstraints()):
        a,b,d=system.getConstraintParameters(i)
        errors.append(abs(float(np.linalg.norm(xyz[int(a)]-xyz[int(b)]))-d.value_in_unit(unit.nanometer)))
    s4.require(np.isfinite(errors).all(),'Nonfinite constraint error')
    return dict(potential_kj_mol=float(potential),max_force_kj_mol_nm=float(np.linalg.norm(f,axis=1).max()),
                max_constraint_error_nm=max(errors,default=0),finite=True)

def smoke(name,device_index):
    from openmm import LangevinMiddleIntegrator,MonteCarloMembraneBarostat,Platform,XmlSerializer,app,unit
    source=s4.OUT/'systems'/name; output=s4.OUT/'smoke'/name
    build=json.loads((source/'build_audit.json').read_text())
    s4.require(build['status']=='COMPLETE' and build['finite'],'Build is not eligible')
    for filename,digest in build['outputs'].items():
        s4.require(s4.sha(source/filename)==digest,'Built asset hash mismatch: '+filename)
    cid=build['candidate']; expected=s4.systems_for(next(r for r in s4.frozen_rows() if r['candidate_id']==cid))[name]
    s4.require(build['ligands']==expected,'Context ligand mismatch')
    s4.require((build['candidate_atoms']>0)==(cid in expected) and (build['probe_atoms']>0)==('ACH' in expected),
               'Context presence/absence mismatch')
    output.mkdir(parents=True,exist_ok=True)
    pdb=app.PDBFile(str(source/'minimized.pdb'))
    system=XmlSerializer.deserialize((source/'system.xml').read_text())
    s4.require(system.getNumParticles()==pdb.topology.getNumAtoms(),'Topology particle count mismatch')
    bars=[f for f in system.getForces() if isinstance(f,MonteCarloMembraneBarostat)]
    s4.require(len(bars)==1,'Expected one membrane barostat'); bars[0].setFrequency(0)
    integrator=LangevinMiddleIntegrator(50*unit.kelvin,1/unit.picosecond,0.0005*unit.picoseconds)
    integrator.setRandomNumberSeed(1701)
    simulation=app.Simulation(pdb.topology,system,integrator,Platform.getPlatformByName('CUDA'),
                              {'DeviceIndex':device_index,'Precision':'mixed'})
    digest=s4.sha(source/'state.xml')+s4.sha(source/'system.xml')+s4.sha(source/'minimized.pdb')
    cache=output/'refined_source_sha256.txt'
    cached=(output/'refined_state.xml').exists() and cache.exists() and cache.read_text().strip()==digest
    simulation.context.setState(XmlSerializer.deserialize((output/'refined_state.xml' if cached else source/'state.xml').read_text()))
    raw=simulation.context.getState(getEnergy=True,getForces=True,getPositions=True)
    initial=finite_metrics(raw,system,raw.getPositions(asNumpy=True),raw.getForces(asNumpy=True))
    deep=not cached and initial['max_force_kj_mol_nm']>5000
    if deep: simulation.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=1000)
    pre=simulation.context.getState(getEnergy=True,getForces=True,getPositions=True)
    metrics=finite_metrics(pre,system,pre.getPositions(asNumpy=True),pre.getForces(asNumpy=True))
    (output/'refined_state.xml').write_text(XmlSerializer.serialize(pre)); cache.write_text(digest+'\n')
    with (output/'refined.pdb').open('w') as f: app.PDBFile.writeFile(pdb.topology,pre.getPositions(),f)
    preflight=dict(metrics,system=name,platform='CUDA',device_index=device_index,precision='mixed',
        source_state_sha256=s4.sha(source/'state.xml'),build_audit_sha256=s4.sha(source/'build_audit.json'),
        used_cached_refinement=cached,performed_deep_minimization=deep,
        pre_refinement_max_force_kj_mol_nm=initial['max_force_kj_mol_nm'])
    s4.dump(output/'preflight.json',preflight)
    simulation.context.setVelocitiesToTemperature(50*unit.kelvin,1701)
    simulation.reporters.append(app.StateDataReporter(str(output/'state.csv'),50,step=True,time=True,
        potentialEnergy=True,kineticEnergy=True,temperature=True,density=True,speed=True,separator=','))
    simulation.reporters.append(app.DCDReporter(str(output/'trajectory.dcd'),50))
    for temperature in [50,100,200,300]:
        integrator.setTemperature(temperature*unit.kelvin); simulation.step(125)
        stage=simulation.context.getState(getEnergy=True,getForces=True,getPositions=True)
        finite_metrics(stage,system,stage.getPositions(asNumpy=True),stage.getForces(asNumpy=True))
    checkpoint=output/'smoke_checkpoint.chk'
    simulation.saveCheckpoint(str(checkpoint))
    final=simulation.context.getState(getEnergy=True,getForces=True,getPositions=True,getVelocities=True)
    final_metrics=finite_metrics(final,system,final.getPositions(asNumpy=True),final.getForces(asNumpy=True))
    s4.require(math.isfinite(final.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)),'Nonfinite kinetic energy')
    # Checkpoint roundtrip without extending the frozen 500-step smoke length.
    before=final.getPositions(asNumpy=True).value_in_unit(unit.nanometer).copy()
    simulation.loadCheckpoint(str(checkpoint))
    after=simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    s4.require(simulation.currentStep==500 and np.array_equal(before,after),'Checkpoint roundtrip mismatch')
    s4.require((output/'trajectory.dcd').stat().st_size>0 and (output/'state.csv').stat().st_size>0,'Reporter output absent')
    (output/'final_state.xml').write_text(XmlSerializer.serialize(final))
    with (output/'final.pdb').open('w') as f: app.PDBFile.writeFile(pdb.topology,final.getPositions(),f)
    audit=dict(final_metrics,system=name,steps=500,timestep_fs=0.5,simulated_ps=0.25,
        barostat_enabled=False,heating_K=[50,100,200,300],checkpoint_roundtrip=True,
        platform='CUDA',device_index=device_index,production_started=False,
        build_audit_sha256=s4.sha(source/'build_audit.json'),
        outputs={f:s4.sha(output/f) for f in ['preflight.json','refined_state.xml','refined.pdb','state.csv',
                 'trajectory.dcd','final_state.xml','smoke_checkpoint.chk']})
    s4.dump(output/'audit.json',audit); print(json.dumps(audit),flush=True)
    if (output/'failure.json').exists(): (output/'failure.json').unlink()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--candidate',choices=list(s4.EXPECTED)); group.add_argument('--all',action='store_true')
    p.add_argument('--device-index',default='0')
    args=p.parse_args(); rows=s4.frozen_rows()
    chosen=rows if args.all else [r for r in rows if r['candidate_id']==args.candidate]
    for row in chosen:
        for name in s4.systems_for(row):
            try: smoke(name,args.device_index)
            except Exception as exc:
                s4.dump(s4.OUT/'smoke'/name/'failure.json',dict(status='FAILED',error=str(exc),production_started=False))
                raise
    s4.update_manifest(rows)

if __name__=='__main__': main()
