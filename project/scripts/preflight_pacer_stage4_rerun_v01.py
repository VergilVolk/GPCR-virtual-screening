"""Actual two-device CUDA integration and authenticated 12-system smoke gates."""
import argparse
import json
import os
import subprocess
from pathlib import Path
from build_pacer_stage4_rerun_v01 import ROOT,sha,setup,table,verify_inputs


def environment():
    assert not os.environ.get('CUDA_VISIBLE_DEVICES'), 'Unset CUDA_VISIBLE_DEVICES; physical DeviceIndex 0/1 is required'
    verify_inputs()
    import openmm
    from openmm import System,Vec3,Platform,LangevinMiddleIntegrator,Context,unit,CustomExternalForce
    import numpy as np
    inventory=subprocess.check_output(['nvidia-smi','--query-gpu=index,uuid,name,driver_version','--format=csv,noheader'],text=True)
    lines=[x.strip() for x in inventory.strip().splitlines()]
    assert len(lines)>=2 and all(any(line.startswith(str(i)+',') for line in lines) for i in [0,1])
    # Imports alone do not verify the frozen force-field/charge stack.
    from openff.toolkit import Molecule
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    import pdbfixer
    from openff.interchange import Interchange
    platform=Platform.getPlatformByName('CUDA'); results=[]
    for device in ['0','1']:
        system=System(); system.addParticle(12*unit.dalton)
        force=CustomExternalForce('0.5*(x*x+y*y+z*z)'); force.addParticle(0,[]); system.addForce(force)
        integrator=LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,0.002*unit.picoseconds)
        context=Context(system,integrator,platform,{'DeviceIndex':device,'Precision':'mixed'})
        context.setPositions([Vec3(0.1,0.2,0.3)]*unit.nanometer)
        context.setVelocitiesToTemperature(300*unit.kelvin,27101); integrator.step(10)
        state=context.getState(getEnergy=True,getPositions=True,getForces=True)
        assert np.isfinite(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        assert np.isfinite(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)).all()
        results.append(dict(device=device,finite=True,actual_device=platform.getPropertyValue(context,'DeviceIndex')))
        del context,integrator
    receipt=dict(status='PASS',openmm_version=openmm.__version__,nvidia_smi=lines,devices=results,
        config_sha256=sha(ROOT/'config/pacer_stage4_rerun_v01.json'))
    (ROOT/'logs/environment_pacer_stage4_rerun_v01.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def systems():
    environment()
    build=json.loads((ROOT/'logs/build_pacer_stage4_rerun_v01.json').read_text()); assert build['status']=='PASS'
    assert build['config_sha256']==sha(ROOT/'config/pacer_stage4_rerun_v01.json')
    setup()
    import run_pacer_stage4_membrane_smoke as inherited
    receipts={}
    for i,row in enumerate(table('STAGE4_RERUN_SYSTEM_MANIFEST.csv')):
        sid=row['system_id']; assert sha(ROOT/row['system_path']/'build_audit.json')==build['builds'][sid]
        receipt=ROOT/'smoke'/sid/'audit.json'
        if receipt.exists():
            audit=json.loads(receipt.read_text()); assert audit['build_audit_sha256']==build['builds'][sid]
            assert audit['steps']==500 and audit['finite'] and audit['checkpoint_roundtrip']
            for name,digest in audit['outputs'].items(): assert sha(receipt.parent/name)==digest
        else: inherited.smoke(sid,str(i%2))
        receipts[sid]=sha(receipt)
    (ROOT/'logs/server_preflight_pacer_stage4_rerun_v01.json').write_text(json.dumps(dict(status='PASS',
        config_sha256=sha(ROOT/'config/pacer_stage4_rerun_v01.json'),builds=build['builds'],smokes=receipts),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(); g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--environment-only',action='store_true'); g.add_argument('--systems',action='store_true'); a=p.parse_args()
    systems() if a.systems else environment()
