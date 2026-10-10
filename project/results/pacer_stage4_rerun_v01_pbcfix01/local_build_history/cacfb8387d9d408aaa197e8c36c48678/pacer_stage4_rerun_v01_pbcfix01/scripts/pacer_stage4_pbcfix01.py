"""Versioned, read-only-source recovery of the Stage 4 PBC restraint defect."""
import argparse
import copy
import csv
import hashlib
import itertools
import json
import math
import os
import shutil
from pathlib import Path

VERSION = 'pacer_stage4_rerun_v01_pbcfix01'
OLD = 'pacer_stage4_rerun_v01'
BAD = '0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)'
GOOD = '0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2'
ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation: failed runs and successful evidence are never overwritten.
    with path.open('x', encoding='utf-8') as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write('\n')


def table(name):
    with (ROOT / 'manifests' / name).open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def verify_manifest(root, filename):
    lines = (root / filename).read_text().splitlines()
    require(bool(lines), 'Empty integrity manifest')
    for line in lines:
        digest, name = line.split('  ', 1)
        p = (root / name).resolve()
        require(root.resolve() in p.parents and p.is_file(), 'Unsafe/missing member: ' + name)
        require(sha(p) == digest, 'Integrity mismatch: ' + str(p))


def source():
    verify_manifest(ROOT, 'SHA256SUMS_PBCFIX01.txt')
    info = json.loads((ROOT / 'source.json').read_text())
    old = Path(info['source_root']).resolve()
    require(old != ROOT.resolve() and old not in ROOT.resolve().parents, 'Output must be a sibling of source')
    require(sha(old / 'SHA256SUMS.txt') == info['legacy_manifest_sha256'], 'Legacy manifest changed')
    verify_manifest(old, 'SHA256SUMS.txt')
    return old


def snapshot_source():
    old = source()
    systems = {}
    for row in table('STAGE4_RERUN_SYSTEM_MANIFEST.csv'):
        d = old / row['source_system_path']
        audit = json.loads((d/'build_audit.json').read_text())
        require(audit['status']=='COMPLETE', 'Incomplete source build')
        for name, digest in audit['outputs'].items():
            require(sha(d/name)==digest, 'Legacy built asset mismatch: '+str(d/name))
        systems[row['source_system_id']] = {p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()}
    all_files = {p.relative_to(old).as_posix():sha(p) for p in sorted(old.rglob('*')) if p.is_file()}
    dump(ROOT/'evidence/source_snapshot.json', dict(source_root=str(old), systems=systems, all_source_files=all_files,
         legacy_manifest_entries=len((old/'SHA256SUMS.txt').read_text().splitlines()),
         note='Hashes retain failed refinement evidence; none of its states are used as starting coordinates'))
    for folder in ['logs','smoke']:
        for p in (old/folder).rglob('*'):
            if p.is_file() and p.suffix in ['.json','.log','.txt']:
                dest = ROOT/'evidence/legacy'/p.relative_to(old)
                dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(p,dest)
    write_sums(ROOT/'evidence')


def fix(system):
    from openmm import CustomExternalForce, XmlSerializer
    indices = [i for i, f in enumerate(system.getForces())
               if isinstance(f, CustomExternalForce) and f.getEnergyFunction().replace(' ', '') == BAD]
    require(len(indices) == 1, 'Expected exactly one known defective restraint')
    index = indices[0]
    force = system.getForce(index)
    require([force.getPerParticleParameterName(i) for i in range(force.getNumPerParticleParameters())]
            == ['x0', 'y0', 'z0'], 'Unexpected restraint parameters')
    require(force.getNumParticles() > 0, 'Empty restraint')
    before = XmlSerializer.serialize(system)
    force.setEnergyFunction(GOOD)
    require(force.usesPeriodicBoundaryConditions(), 'Restraint must use PBC')
    # Prove all parameters, bonds, constraints and other forces remain byte-identical.
    restored = XmlSerializer.deserialize(XmlSerializer.serialize(system))
    restored.getForce(index).setEnergyFunction(BAD)
    require(XmlSerializer.serialize(restored) == before, 'Unexpected scientific parameter mutation')
    return index


def periodic_squared(delta, box):
    """Independent nearest-lattice calculation for OpenMM reduced triclinic boxes.

    Enumerate all candidates within a mathematically bounded radius, rather than
    fractional component rounding (which is wrong for general triclinic boxes).
    """
    import numpy as np
    delta, box = np.asarray(delta, float), np.asarray(box, float)
    frac = np.linalg.solve(box.T, delta)
    center = np.rint(frac).astype(int)
    trial = delta - center @ box
    radius = np.linalg.norm(trial) / np.linalg.svd(box, compute_uv=False).min() + 1e-10
    ranges = [range(math.ceil(x-radius), math.floor(x+radius)+1) for x in frac]
    return min(float(np.dot(delta-np.asarray(n)@box, delta-np.asarray(n)@box))
               for n in itertools.product(*ranges))


def energy_comparison(system, index, positions, platform='Reference', device='0'):
    from openmm import System, XmlSerializer, Context, VerletIntegrator, Platform, unit
    import numpy as np
    force = system.getForce(index)
    xyz = np.asarray(positions.value_in_unit(unit.nanometer))
    box = np.asarray([v.value_in_unit(unit.nanometer) for v in system.getDefaultPeriodicBoxVectors()])
    k = next(force.getGlobalParameterDefaultValue(i) for i in range(force.getNumGlobalParameters())
             if force.getGlobalParameterName(i) == 'k')
    analytical = {'before': 0., 'after': 0.}
    for i in range(force.getNumParticles()):
        atom, ref = force.getParticleParameters(i)
        delta = xyz[atom] - np.asarray(ref)
        analytical['before'] += 0.5*k*float(np.dot(delta, delta))
        analytical['after'] += 0.5*k*periodic_squared(delta, box)
    measured = {}
    for label, expression in [('before', BAD), ('after', GOOD)]:
        isolated = System()
        for i in range(system.getNumParticles()):
            isolated.addParticle(system.getParticleMass(i))
        isolated.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        f = XmlSerializer.deserialize(XmlSerializer.serialize(force))
        f.setEnergyFunction(expression)
        isolated.addForce(f)
        integ = VerletIntegrator(0.001)
        props = {'DeviceIndex': device, 'Precision': 'mixed'} if platform == 'CUDA' else {}
        ctx = Context(isolated, integ, Platform.getPlatformByName(platform), props)
        ctx.setPositions(positions)
        measured[label] = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        require(math.isclose(measured[label], analytical[label], rel_tol=2e-5, abs_tol=0.05),
                'Independent restraint energy mismatch: ' + label)
        del ctx, integ
    return dict(openmm_kj_mol=measured, independent_kj_mol=analytical, platform=platform,
                method='bounded nearest-lattice enumeration; original reference coordinates unchanged')


def metrics(state):
    import numpy as np
    from openmm import unit
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
    require(math.isfinite(pe) and math.isfinite(ke) and np.isfinite(xyz).all() and np.isfinite(forces).all(),
            'Nonfinite energy/position/force')
    return dict(potential_kj_mol=float(pe), kinetic_kj_mol=float(ke),
                max_force_kj_mol_nm=float(np.linalg.norm(forces, axis=1).max()), finite=True)


def identity(sid):
    source()
    d = ROOT / 'systems' / sid
    audit = json.loads((d / 'build_audit.json').read_text())
    require(audit['protocol_id'] == VERSION and audit['status'] == 'COMPLETE', 'Wrong build')
    for name, digest in audit['outputs'].items():
        require(sha(d / name) == digest, 'Build changed: ' + name)
    return dict(protocol_id=VERSION, system_id=sid, build_audit_sha256=sha(d/'build_audit.json'),
                config_sha256=sha(ROOT/'config'/f'{VERSION}.json'),
                package_manifest_sha256=sha(ROOT/'SHA256SUMS_PBCFIX01.txt'),
                system_sha256=sha(d/'system.xml'), topology_sha256=sha(d/'minimized.pdb'),
                state_sha256=sha(d/'state.xml'), source_snapshot_sha256=sha(ROOT/'evidence/source_snapshot.json'))


def recover(row, device):
    from openmm import XmlSerializer, app, unit, Platform, LangevinMiddleIntegrator, MonteCarloMembraneBarostat, NonbondedForce
    old = source()
    sid = row['system_id']
    src = old / row['source_system_path']
    snapshot = json.loads((ROOT/'evidence/source_snapshot.json').read_text())
    for name, digest in snapshot['systems'][row['source_system_id']].items():
        require(sha(src/name) == digest, 'Source changed after snapshot: '+name)
    out = ROOT/'systems'/sid
    out.mkdir(parents=True, exist_ok=False)
    # Never deserialize old state.xml/refined_state.xml. Only system parameters and raw coordinates.
    pdb = app.PDBFile(str(src/'input.pdb'))
    system = XmlSerializer.deserialize((src/'system.xml').read_text())
    require(pdb.topology.getNumAtoms() == system.getNumParticles(), 'Particle count mismatch')
    index = fix(system)
    comparison = energy_comparison(system, index, pdb.positions, 'CUDA', device)
    dump(out/'restraint_regression.json', comparison)
    shutil.copyfile(src/'input.pdb', out/'input.pdb')
    shutil.copyfile(src/'parameters.json', out/'parameters.json')
    # Keep original barostat parameters in the serialized scientific system; disable only in minimization context.
    (out/'system.xml').write_text(XmlSerializer.serialize(system))
    bars = [f for f in system.getForces() if isinstance(f, MonteCarloMembraneBarostat)]
    require(len(bars) == 1, 'Expected one membrane barostat')
    bars[0].setFrequency(0)
    integ = LangevinMiddleIntegrator(50*unit.kelvin, 1/unit.picosecond, .0005*unit.picoseconds)
    sim = app.Simulation(pdb.topology, system, integ, Platform.getPlatformByName('CUDA'),
                         {'DeviceIndex': device, 'Precision': 'mixed'})
    sim.context.setPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    sim.context.setPositions(pdb.positions)
    initial = metrics(sim.context.getState(getEnergy=True, getPositions=True, getForces=True))
    sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer, maxIterations=1000)
    state = sim.context.getState(getEnergy=True, getPositions=True, getForces=True, getVelocities=True)
    final = metrics(state)
    (out/'state.xml').write_text(XmlSerializer.serialize(state))
    with (out/'minimized.pdb').open('x') as f:
        app.PDBFile.writeFile(pdb.topology, state.getPositions(), f)
    nb = [f for f in system.getForces() if isinstance(f, NonbondedForce)]
    require(len(nb) == 1, 'Expected one nonbonded force')
    charge = sum(nb[0].getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                 for i in range(system.getNumParticles()))
    require(min(abs(charge-x) for x in [0, 1, 2, 3]) < 1e-4, 'Unexpected charge; scientific review required')
    outputs = {n: sha(out/n) for n in ['input.pdb','parameters.json','system.xml','state.xml','minimized.pdb','restraint_regression.json']}
    dump(out/'build_audit.json', dict(protocol_id=VERSION, status='COMPLETE', source_system_id=row['source_system_id'],
         source_hashes=snapshot['systems'][row['source_system_id']], outputs=outputs, initial=initial, minimized=final,
         net_charge_e=charge, ions={n: sum(r.name.upper()==n for r in pdb.topology.residues()) for n in ['NA','CL']},
         ion_policy='Unchanged: baseline neutralized before ligand insertion; production scientific review required',
         restraint_expression=GOOD, only_force_expression_changed=True, reused_old_state=False,
         platform='CUDA', device=device, precision='mixed', production_started=False))
    write_sums(out)


def write_sums(directory):
    paths = sorted(p for p in directory.rglob('*') if p.is_file() and p.name != 'SHA256SUMS_OUTPUT.txt')
    with (directory/'SHA256SUMS_OUTPUT.txt').open('x') as f:
        f.writelines(f'{sha(p)}  {p.relative_to(directory).as_posix()}\n' for p in paths)


def smoke(sid, device, single=False):
    from openmm import XmlSerializer, app, unit, Platform, LangevinMiddleIntegrator, MonteCarloMembraneBarostat
    import numpy as np
    ident = identity(sid)
    d = ROOT/'systems'/sid
    out = ROOT/('single_apo' if single else 'smoke')/sid
    out.mkdir(parents=True, exist_ok=False)
    pdb = app.PDBFile(str(d/'minimized.pdb'))
    system = XmlSerializer.deserialize((d/'system.xml').read_text())
    for f in system.getForces():
        if isinstance(f, MonteCarloMembraneBarostat):
            f.setFrequency(0)
    integ = LangevinMiddleIntegrator(50*unit.kelvin, 1/unit.picosecond, .0005*unit.picoseconds)
    integ.setRandomNumberSeed(1701)
    sim = app.Simulation(pdb.topology, system, integ, Platform.getPlatformByName('CUDA'),
                         {'DeviceIndex':device, 'Precision':'mixed'})
    sim.context.setState(XmlSerializer.deserialize((d/'state.xml').read_text()))
    stages = [dict(label='minimized', **metrics(sim.context.getState(getEnergy=True,getPositions=True,getForces=True)))]
    sim.context.setVelocitiesToTemperature(50*unit.kelvin, 1701)
    for temp in [50,100,200,300]:
        integ.setTemperature(temp*unit.kelvin)
        sim.step(125)
        stages.append(dict(temperature_K=temp, **metrics(sim.context.getState(getEnergy=True,getPositions=True,getForces=True))))
    state = sim.context.getState(getEnergy=True,getPositions=True,getForces=True,getVelocities=True)
    (out/'final_state.xml').write_text(XmlSerializer.serialize(state))
    sim.saveCheckpoint(str(out/'smoke_checkpoint.chk'))
    sim.loadCheckpoint(str(out/'smoke_checkpoint.chk'))
    after = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    require(sim.currentStep == 500 and np.array_equal(after, state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)),
            'Checkpoint roundtrip failed')
    dump(out/'audit.json', dict(protocol_id=VERSION, identity=ident, steps=500, stages=stages, finite=True,
         checkpoint_roundtrip=True, platform='CUDA', device=device, production_started=False,
         outputs={n:sha(out/n) for n in ['final_state.xml','smoke_checkpoint.chk']}))
    write_sums(out)


def check_smoke(sid, single=False):
    path = ROOT/('single_apo' if single else 'smoke')/sid/'audit.json'
    audit = json.loads(path.read_text())
    require(audit['identity'] == identity(sid) and audit['steps']==500 and audit['finite']
            and audit['checkpoint_roundtrip'] and audit['platform']=='CUDA', 'Invalid smoke')
    for name, digest in audit['outputs'].items():
        require(sha(path.parent/name)==digest, 'Smoke output changed')
    return sha(path)


def storage_budget():
    """Measured binary checkpoints + 10 full-system frames per 100 ps segment."""
    from openmm import XmlSerializer
    rows = table('STAGE4_RERUN_SYSTEM_MANIFEST.csv')
    detail = []
    for row in rows:
        sid = row['system_id']
        check_smoke(sid)
        n = XmlSerializer.deserialize((ROOT/'systems'/sid/'system.xml').read_text()).getNumParticles()
        checkpoint = (ROOT/'smoke'/sid/'smoke_checkpoint.chk').stat().st_size
        # 12 bytes/atom/frame DCD, generous headers/CSV/receipt; 3 replicas * 100 segments.
        per = checkpoint + 10*n*12 + 65536
        detail.append(dict(system_id=sid, particles=n, measured_checkpoint_bytes=checkpoint,
                           segments=300, estimated_bytes=300*per))
    raw = sum(x['estimated_bytes'] for x in detail)
    required = math.ceil(raw*1.5) + 10*1024**3
    usage = shutil.disk_usage(ROOT)
    return dict(protocol_id=VERSION, trajectories=36, total_ns=360, checkpoint_ps=100, checkpoints=3600,
                trajectory_ps=10, systems=detail, raw_bytes=raw, safety_factor=1.5, reserve_bytes=10*1024**3,
                required_free_bytes=required, data_path=str(ROOT), disk_free_bytes=usage.free,
                disk_total_bytes=usage.total, sufficient=usage.free>=required,
                limitations='Binary checkpoint size measured on CUDA; 50% margin plus 10 GiB reserve; excludes unrelated campaigns')


def preflight():
    rows = table('STAGE4_RERUN_SYSTEM_MANIFEST.csv')
    require(len(rows)==12 and len({r['system_id'] for r in rows})==12, 'Expected 12 distinct systems')
    first = next(r for r in rows if r['context']=='apo')['system_id']
    single = check_smoke(first, True)
    smokes = {r['system_id']:check_smoke(r['system_id']) for r in rows}
    budget = storage_budget()
    dump(ROOT/'logs/storage_budget.json', budget)
    require(budget['sufficient'], 'Insufficient data disk capacity')
    dump(ROOT/'logs'/f'server_preflight_{VERSION}.json', dict(status='PASS', protocol_id=VERSION,
         config_sha256=sha(ROOT/'config'/f'{VERSION}.json'),
         builds={r['system_id']:sha(ROOT/'systems'/r['system_id']/'build_audit.json') for r in rows},
         smokes=smokes, single_apo_sha256=single, storage_budget_sha256=sha(ROOT/'logs/storage_budget.json'),
         scientific_review='PENDING', user_authorization='PENDING', production_started=False))


def authorized_gate(job):
    """Called by both prepare and production; authorization is a separate manual artifact."""
    ident = identity(job['system_id'])
    gatepath = ROOT/'logs'/f'server_preflight_{VERSION}.json'
    gate = json.loads(gatepath.read_text())
    require(gate['status']=='PASS' and gate['builds'][job['system_id']]==ident['build_audit_sha256'], 'Preflight mismatch')
    require(gate['config_sha256']==ident['config_sha256'], 'Preflight config mismatch')
    for row in table('STAGE4_RERUN_SYSTEM_MANIFEST.csv'):
        require(check_smoke(row['system_id'])==gate['smokes'][row['system_id']], 'Smoke gate changed')
    first = next(r for r in table('STAGE4_RERUN_SYSTEM_MANIFEST.csv') if r['context']=='apo')['system_id']
    require(check_smoke(first, True)==gate['single_apo_sha256'], 'Single Apo gate changed')
    require(sha(ROOT/'logs/storage_budget.json')==gate['storage_budget_sha256'], 'Storage receipt changed')
    require(storage_budget()['sufficient'], 'Current disk capacity insufficient')
    approval = json.loads((ROOT/'authorization.json').read_text())
    require(approval['protocol_id']==VERSION and approval['preflight_sha256']==sha(gatepath)
            and approval['scientific_review_approved'] is True and approval['user_explicit_authorization'] is True
            and bool(approval['reviewer']) and bool(approval['user_authorization_reference']), 'Manual review/authorization missing')
    return dict(ident, preflight_sha256=sha(gatepath), authorization_sha256=sha(ROOT/'authorization.json'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['snapshot','recover-single-apo','smoke-single-apo','recover-all','smoke-all','preflight'])
    p.add_argument('--device', choices=['0','1'], default='0')
    a = p.parse_args()
    require(not os.environ.get('CUDA_VISIBLE_DEVICES'), 'Unset CUDA_VISIBLE_DEVICES')
    source()
    rows = table('STAGE4_RERUN_SYSTEM_MANIFEST.csv')
    first = next(r for r in rows if r['context']=='apo')
    try:
        if a.action=='snapshot': snapshot_source()
        elif a.action=='recover-single-apo': recover(first, a.device)
        elif a.action=='smoke-single-apo':
            smoke(first['system_id'], a.device, True)
            print('PBC_FIX_SINGLE_APO_SMOKE_PASS', flush=True)
        elif a.action=='recover-all':
            check_smoke(first['system_id'], True)
            for row in rows:
                if row['system_id']==first['system_id']:
                    identity(row['system_id'])
                else: recover(row, a.device)
        elif a.action=='smoke-all':
            check_smoke(first['system_id'], True)
            for i,row in enumerate(rows): smoke(row['system_id'], str(i%2))
        else: preflight()
    except Exception as exc:
        import uuid, traceback
        dump(ROOT/'logs'/f'{a.action}_failure_{uuid.uuid4().hex}.json',
             dict(status='FAILED', error=str(exc), traceback=traceback.format_exc(), production_started=False))
        raise


if __name__=='__main__': main()
