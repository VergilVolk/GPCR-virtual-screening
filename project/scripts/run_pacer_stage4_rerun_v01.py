"""One authenticated CUDA job. Immutable committed segments make resume unambiguous."""
import argparse
import fcntl
import json
import os
import uuid
from pathlib import Path
from build_pacer_stage4_rerun_v01 import ROOT,sha,table,verify_inputs


def dump(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_suffix(p.suffix+'.tmp'); temp.write_text(json.dumps(obj,indent=2)+'\n'); os.replace(temp,p)


def fingerprint(job):
    verify_inputs()
    d=ROOT/'systems'/job['system_id']
    gate=json.loads((ROOT/'logs/server_preflight_pacer_stage4_rerun_v01.json').read_text())
    assert gate['status']=='PASS' and gate['config_sha256']==sha(ROOT/'config/pacer_stage4_rerun_v01.json')
    assert gate['builds'][job['system_id']]==sha(d/'build_audit.json')
    assert gate['smokes'][job['system_id']]==sha(ROOT/'smoke'/job['system_id']/'audit.json')
    audit=json.loads((d/'build_audit.json').read_text())
    for file,digest in audit['outputs'].items(): assert sha(d/file)==digest
    return dict(job=job,config_sha256=gate['config_sha256'],system_sha256=sha(d/'system.xml'),
                topology_sha256=sha(d/'minimized.pdb'),smoke_sha256=gate['smokes'][job['system_id']])


def simulation(job,device,state_path,dt,barostat,seed,k,temperature=300):
    from openmm import LangevinMiddleIntegrator,MonteCarloMembraneBarostat,Platform,XmlSerializer,app,unit
    d=ROOT/'systems'/job['system_id']; pdb=app.PDBFile(str(d/'minimized.pdb'))
    system=XmlSerializer.deserialize((d/'system.xml').read_text())
    bars=[f for f in system.getForces() if isinstance(f,MonteCarloMembraneBarostat)]
    assert len(bars)==1; bars[0].setFrequency(barostat); bars[0].setRandomNumberSeed(seed)
    integ=LangevinMiddleIntegrator(temperature*unit.kelvin,1/unit.picosecond,dt*unit.picoseconds)
    integ.setRandomNumberSeed(seed)
    sim=app.Simulation(pdb.topology,system,integ,Platform.getPlatformByName('CUDA'),{'DeviceIndex':str(device),'Precision':'mixed'})
    sim.context.setState(XmlSerializer.deserialize(state_path.read_text())); sim.context.setParameter('k',k)
    return sim,system


def prepare(job,device):
    from openmm import XmlSerializer,unit,app
    import numpy as np
    import sys
    sys.path.insert(0,str(ROOT/'scripts/legacy'))
    from run_pacer_dc_short_equilibration import assert_finite
    identity=fingerprint(job); seed=int(job['seed'])
    base=ROOT/'equilibration'/job['system_id']/f"replica_{int(job['replica']):02d}"
    base.mkdir(parents=True,exist_ok=True)
    lock=(base/'task.lock').open('a'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    initial=ROOT/'smoke'/job['system_id']/'refined_state.xml'
    identity['initial_state_sha256']=sha(initial)
    previous=initial; previous_hash=sha(initial); stages=[]
    definitions=[('nvt',1,.0005,0,seed,1000),('npt',1,.0005,25,seed+1000003,1000)]
    definitions += [(f'release_{i:02d}',.5,.001,25,seed+i*1000003,k) for i,k in enumerate([500,100,10,0],1)]
    for label,ps,dt,bfreq,stage_seed,k in definitions:
        d=base/label; receipt=d/'receipt.json'
        stage_identity=dict(identity,stage=label,previous_state_sha256=previous_hash)
        if receipt.exists():
            r=json.loads(receipt.read_text()); assert r['identity']==stage_identity and r['passed']
            for file,digest in r['outputs'].items(): assert sha(d/file)==digest
        else:
            if d.exists(): d.rename(base/(label+'.interrupted.'+uuid.uuid4().hex))
            d.mkdir()
            sim,system=simulation(job,device,previous,dt,bfreq,stage_seed,k)
            if label=='nvt': sim.context.setVelocitiesToTemperature(300*unit.kelvin,seed)
            prior=sim.context.getState(getPositions=True)
            volume0=abs(np.linalg.det(prior.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)))
            steps=round(ps/dt); interval=200 if label in ['nvt','npt'] else max(1,steps//10)
            sim.reporters.append(app.StateDataReporter(str(d/'state.csv'),interval,step=True,time=True,potentialEnergy=True,kineticEnergy=True,temperature=True,volume=True))
            sim.reporters.append(app.DCDReporter(str(d/'trajectory.dcd'),interval))
            sim.step(steps); state=sim.context.getState(getEnergy=True,getPositions=True,getVelocities=True)
            metrics=assert_finite(state,system,label)
            volume1=abs(np.linalg.det(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)))
            change=(volume1/volume0-1)*100
            # Release follows the inherited total-volume reference, not a fresh reference at each stage.
            release_start=base/'release_volume_reference.json'
            if label=='release_01': dump(release_start,dict(volume_nm3=float(volume0),identity=identity))
            if label.startswith('release'):
                vr=json.loads(release_start.read_text()); assert vr['identity']==identity
                change=(volume1/vr['volume_nm3']-1)*100
            assert abs(change)<20, 'Equilibration/release volume gate failed'
            if label in ['npt','release_04']: assert 270<=metrics['temperature_K']<=330, 'Temperature gate failed'
            (d/'state.xml').write_text(XmlSerializer.serialize(state))
            r=dict(identity=stage_identity,passed=True,metrics=metrics,volume_change_percent=float(change),
                   outputs={f:sha(d/f) for f in ['state.xml','state.csv','trajectory.dcd']})
            dump(receipt,r); del sim,system
        stages.append(sha(receipt)); previous=d/'state.xml'; previous_hash=sha(previous)
    dump(base/'ready_pacer_stage4_rerun_v01.json',dict(identity=identity,passed=True,stage_receipts=stages,
        start_state=str(previous.relative_to(ROOT)),start_state_sha256=previous_hash))


def production(job,device,authorized):
    assert authorized, 'Manual --authorize-production required'
    from openmm import unit,app
    identity=fingerprint(job)
    prep=ROOT/'equilibration'/job['system_id']/f"replica_{int(job['replica']):02d}"
    ready=json.loads((prep/'ready_pacer_stage4_rerun_v01.json').read_text()); assert ready['passed']
    # Validate preparation identity including the smoke-refined source and every stage output.
    expected=dict(identity,initial_state_sha256=sha(ROOT/'smoke'/job['system_id']/'refined_state.xml'))
    assert ready['identity']==expected
    for index,label in enumerate(['nvt','npt','release_01','release_02','release_03','release_04']):
        rpath=prep/label/'receipt.json'; assert sha(rpath)==ready['stage_receipts'][index]
        r=json.loads(rpath.read_text()); assert r['passed']
        for file,digest in r['outputs'].items(): assert sha(rpath.parent/file)==digest
    start=ROOT/ready['start_state']; assert sha(start)==ready['start_state_sha256']
    identity.update(start_state_sha256=sha(start),preparation_receipt_sha256=sha(prep/'ready_pacer_stage4_rerun_v01.json'))
    out=ROOT/job['output_path']; out.mkdir(parents=True,exist_ok=True)
    lock=(out/'task.lock').open('a'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    sim,system=simulation(job,device,start,.002,25,int(job['seed']),0)
    committed=0; last_checkpoint=None
    for end in range(50000,5000001,50000):
        d=out/'segments'/f'{end:08d}'; receipt=d/'receipt.json'
        if receipt.exists():
            r=json.loads(receipt.read_text()); assert r['identity']==identity and r['start_step']==committed and r['end_step']==end
            for file,digest in r['outputs'].items(): assert sha(d/file)==digest
            committed=end; last_checkpoint=d/'checkpoint.chk'
            continue
        # A committed later segment with a missing predecessor is not recoverable by guessing.
        assert not any((out/'segments'/f'{later:08d}'/'receipt.json').exists() for later in range(end+50000,5000001,50000))
        if last_checkpoint:
            sim.loadCheckpoint(str(last_checkpoint)); assert sim.currentStep==committed
        elif committed==0:
            sim.currentStep=0; sim.context.setTime(0*unit.picoseconds)
            sim.context.setVelocitiesToTemperature(300*unit.kelvin,int(job['seed']))
        if d.exists(): d.rename(d.with_name(d.name+'.interrupted.'+uuid.uuid4().hex))
        d.mkdir(parents=True)
        sim.reporters.append(app.DCDReporter(str(d/'trajectory.dcd'),5000,enforcePeriodicBox=False))
        sim.reporters.append(app.StateDataReporter(str(d/'state.csv'),5000,step=True,time=True,potentialEnergy=True,kineticEnergy=True,temperature=True,volume=True))
        sim.step(end-committed)
        state=sim.context.getState(getEnergy=True,getPositions=True,getVelocities=True)
        import sys
        sys.path.insert(0,str(ROOT/'scripts/legacy'))
        from run_pacer_dc_short_equilibration import assert_finite
        metrics=assert_finite(state,system,'production')
        sim.saveCheckpoint(str(d/'checkpoint.chk'))
        # Close trajectory handles before hashing; committed segments never receive appended frames.
        sim.reporters.clear()
        dump(receipt,dict(identity=identity,start_step=committed,end_step=end,metrics=metrics,
             outputs={f:sha(d/f) for f in ['trajectory.dcd','state.csv','checkpoint.chk']}))
        committed=end; last_checkpoint=d/'checkpoint.chk'
        dump(out/'progress.json',dict(identity=identity,status='COMPLETE' if end==5000000 else 'RUNNING',
                                     completed_steps=end,completed_ns=end*.002/1000,last_segment=str(d.relative_to(ROOT))))
    assert committed==5000000
    dump(out/'progress.json',dict(identity=identity,status='COMPLETE',completed_steps=committed,completed_ns=10))


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--job',required=True); p.add_argument('--device',choices=['0','1'],required=True)
    p.add_argument('--phase',choices=['prepare','production'],required=True); p.add_argument('--authorize-production',action='store_true'); a=p.parse_args()
    assert not os.environ.get('CUDA_VISIBLE_DEVICES'), 'Unset CUDA_VISIBLE_DEVICES'
    rows=[r for r in table('STAGE4_RERUN_JOB_MATRIX.csv') if r['job_id']==a.job]; assert len(rows)==1
    try:
        prepare(rows[0],a.device) if a.phase=='prepare' else production(rows[0],a.device,a.authorize_production)
    except Exception as exc:
        dump(ROOT/'logs'/f'{a.job}__{a.phase}__failure.json',dict(status='FAILED',error=str(exc))); raise
