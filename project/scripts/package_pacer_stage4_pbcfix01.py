"""Create a minimal additive recovery package; never modify the v01 bundle."""
import ast
import csv
import hashlib
import io
import json
import shutil
import subprocess
import tarfile
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
OLD = 'pacer_stage4_rerun_v01'
NEW = OLD+'_pbcfix01'
LEGACY = PROJECT/'results'/OLD/'server_bundle'/OLD
OUT = PROJECT/'results'/NEW
BUNDLE = OUT/'server_bundle'/NEW


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def replace_checked(text, before, after):
    if before not in text:
        raise RuntimeError('Upstream adapter anchor missing: '+before)
    return text.replace(before,after)


def write(rel, text):
    p=BUNDLE/rel
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8',newline='\n') as f: f.write(text)


def main():
    branch=subprocess.check_output(['git','branch','--show-current'],text=True).strip()
    if branch!='results/pacer-stage4-rerun-preparation-v01': raise RuntimeError('Wrong branch')
    if BUNDLE.exists():
        import uuid
        if OUT.resolve() not in BUNDLE.resolve().parents: raise RuntimeError('Unsafe output path')
        history=OUT/'local_build_history'/uuid.uuid4().hex
        history.mkdir(parents=True)
        BUNDLE.rename(history/NEW)
        for name in [f'{NEW}_repair_bundle.tar.gz',f'{NEW}_repair_bundle.tar.gz.sha256','LOCAL_TEST_RESULTS.txt']:
            p=OUT/name
            if p.exists(): p.rename(history/name)
    BUNDLE.mkdir(parents=True,exist_ok=False)
    for folder in ['logs','manifests','config','scripts']:
        (BUNDLE/folder).mkdir()
    write('source.json',json.dumps(dict(source_root='/data/pacer_stage4_rerun_v01_deployment/'+OLD,
         legacy_manifest_sha256=sha(LEGACY/'SHA256SUMS.txt'), source_branch=branch,
         source_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()),indent=2)+'\n')
    config=json.loads((LEGACY/'config'/f'{OLD}.json').read_text())
    config.update(protocol_id=NEW,restraint_expression='0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2',
                  initial_state_policy='Re-minimized raw input.pdb; never old refined_state.xml',
                  scientific_review_required=True,user_authorization_required_before_prepare=True)
    write(f'config/{NEW}.json',json.dumps(config,indent=2)+'\n')
    for filename in ['STAGE4_RERUN_SYSTEM_MANIFEST.csv','STAGE4_RERUN_JOB_MATRIX.csv']:
        with (LEGACY/'manifests'/filename).open(newline='') as f: rows=list(csv.DictReader(f))
        old_fields=list(rows[0])
        source_fields=['source_system_id','source_system_path'] if 'SYSTEM_MANIFEST' in filename else ['source_job_id']
        for row in rows:
            originals=dict(row)
            for k,v in row.items(): row[k]=v.replace(OLD,NEW)
            if 'SYSTEM_MANIFEST' in filename:
                row.update(source_system_id=originals['system_id'],source_system_path=originals['system_path'])
            else: row['source_job_id']=originals['job_id']
        stream=io.StringIO(newline=''); w=csv.DictWriter(stream,fieldnames=old_fields+source_fields)
        w.writeheader(); w.writerows(rows); write('manifests/'+filename,stream.getvalue())
    shutil.copyfile(PROJECT/'scripts/pacer_stage4_pbcfix01.py', BUNDLE/'scripts/pacer_stage4_pbcfix01.py')
    shutil.copyfile(PROJECT/'tests/test_pacer_stage4_pbcfix01.py', BUNDLE/'scripts/test_pacer_stage4_pbcfix01.py')
    # Reuse established dynamics and segment handling via a separately materialized adapter.
    runner=(PROJECT/'scripts'/f'run_{OLD}.py').read_text().replace(OLD,NEW)
    runner=replace_checked(runner,
        f'from build_{NEW} import ROOT,sha,table,verify_inputs',
        'from pacer_stage4_pbcfix01 import ROOT,sha,table,source as verify_inputs,authorized_gate')
    runner=replace_checked(runner,'    verify_inputs()\n    d=ROOT/',
        '    authorization_identity=authorized_gate(job)\n    verify_inputs()\n    d=ROOT/')
    runner=replace_checked(runner,'    return dict(job=job,config_sha256=',
        '    return dict(authorization_identity=authorization_identity,job=job,config_sha256=')
    runner=replace_checked(runner,"ROOT/'smoke'/job['system_id']/'refined_state.xml'",
        "ROOT/'systems'/job['system_id']/'state.xml'")
    # Both phases require the explicit flag; authorization file additionally binds review to preflight.
    runner=replace_checked(runner,"    assert not os.environ.get('CUDA_VISIBLE_DEVICES')",
        "    assert a.authorize_production, 'Explicit authorization required for prepare AND production'\n    assert not os.environ.get('CUDA_VISIBLE_DEVICES')")
    runner=replace_checked(runner,"sys.path.insert(0,str(ROOT/'scripts/legacy'))", "sys.path.insert(0,str(ROOT/'scripts'))")
    write(f'scripts/run_{NEW}.py',runner)
    # Only this finite-metrics dependency is needed; copy intact and authenticate it.
    shutil.copyfile(PROJECT/'scripts/run_pacer_dc_short_equilibration.py',BUNDLE/'scripts/run_pacer_dc_short_equilibration.py')
    scheduler=(PROJECT/'scripts'/f'schedule_{OLD}.py').read_text().replace(OLD,NEW)
    scheduler=replace_checked(scheduler,f'from build_{NEW} import ROOT,table',
                              'from pacer_stage4_pbcfix01 import ROOT,table,authorized_gate')
    scheduler=replace_checked(scheduler,"assert a.phase!='production' or a.authorize_production",
                              'assert a.authorize_production')
    scheduler=replace_checked(scheduler,
        f"    subprocess.run([sys.executable,str(ROOT/'scripts/preflight_{NEW}.py'),'--environment-only'],check=True)",
        "    for job in table('STAGE4_RERUN_JOB_MATRIX.csv'): authorized_gate(job)")
    write(f'scripts/schedule_{NEW}.py',scheduler)
    guide=(PROJECT/'scripts/PBCFIX01_DEPLOYMENT_GUIDE.md').read_text(encoding='utf-8')
    write('DEPLOYMENT_GUIDE.md',guide)
    write('evidence/server_experiment_report.json',json.dumps(dict(
        evidence_type='User supplied server experiment; not independently rerun locally',
        before_restraint_kj_mol=269484748,after_restraint_kj_mol=322.347,
        minimized_potential_kj_mol=-1857111.364,max_force_kj_mol_nm=2877.257,
        steps=500,heating_K=[50,100,200,300],result='PBC_FIX_SINGLE_APO_SMOKE_PASS'),indent=2)+'\n')
    for p in BUNDLE.rglob('*.py'): ast.parse(p.read_text(encoding='utf-8'))
    members=sorted(p for p in BUNDLE.rglob('*') if p.is_file())
    write('SHA256SUMS_PBCFIX01.txt',''.join(f'{sha(p)}  {p.relative_to(BUNDLE).as_posix()}\n' for p in members))
    archive=OUT/f'{NEW}_repair_bundle.tar.gz'
    with tarfile.open(archive,'w:gz') as tar: tar.add(BUNDLE,arcname=NEW)
    (OUT/(archive.name+'.sha256')).write_text(f'{sha(archive)}  {archive.name}\n')
    shutil.copyfile(BUNDLE/'DEPLOYMENT_GUIDE.md',OUT/'DEPLOYMENT_GUIDE.md')
    print(json.dumps(dict(archive=str(archive),sha256=sha(archive),files=len(members)+1,bytes=archive.stat().st_size),indent=2))


if __name__=='__main__': main()
