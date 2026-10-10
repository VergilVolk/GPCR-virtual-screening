"""Read-only local deployment QA; write a QC receipt in the new rerun log directory."""
import ast
import csv
import hashlib
import json
import subprocess
import sys
import tarfile
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
OUT=REPO/'project/results/pacer_stage4_rerun_v01'
ROOT=OUT/'server_bundle/pacer_stage4_rerun_v01'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    if '--attempt-local-build' in sys.argv:
        result=subprocess.run([sys.executable,str(ROOT/'scripts/build_pacer_stage4_rerun_v01.py'),'--all'],capture_output=True,text=True)
        error=result.stdout+result.stderr
        (OUT/'logs/local_build_attempt_pacer_stage4_rerun_v01.log').write_text(error,encoding='utf-8')
        attempt=dict(status='BLOCKED_LOCAL_DEPENDENCIES' if result.returncode else 'PASS',exit_code=result.returncode,
                     command='python server_bundle/pacer_stage4_rerun_v01/scripts/build_pacer_stage4_rerun_v01.py --all',
                     production_started=False)
        (OUT/'logs/local_build_attempt_pacer_stage4_rerun_v01.json').write_text(json.dumps(attempt,indent=2)+'\n')
        assert result.returncode==0 or "ModuleNotFoundError: No module named 'openff'" in error, error
    def csvrows(name):
        with (OUT/name).open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))
    jobs=csvrows('STAGE4_RERUN_JOB_MATRIX.csv'); systems=csvrows('STAGE4_RERUN_SYSTEM_MANIFEST.csv')
    assert len(systems)==12 and len(jobs)==36
    assert len({r['system_id'] for r in systems})==12
    assert len({r['job_id'] for r in jobs})==len({r['output_path'] for r in jobs})==36
    for candidate in ['PACERGEN02432','PACERGEN00123','PACERGEN00462']:
        family=[j for j in jobs if j['candidate_id']==candidate]
        assert len(family)==12
        for context in ['apo','probe_only','candidate_no_probe','candidate_probe']:
            triple=[j for j in family if j['context']==context]
            assert {(int(j['replica']),int(j['seed'])) for j in triple}=={(1,27101),(2,38201),(3,49301)}
    for row in systems:
        for key in ['receptor','candidate_sdf','probe_sdf']:
            if row[key]: assert (ROOT/row[key]).is_file() and not Path(row[key]).is_absolute()
        assert row['has_candidate']==str(row['context'].startswith('candidate'))
        assert row['has_ACH']==str(row['context'] in ['probe_only','candidate_probe'])
    config=json.loads((ROOT/'config/pacer_stage4_rerun_v01.json').read_text())
    assert config['production_steps']==5000000 and config['trajectory_ps']==10 and config['timestep_fs']==2
    sources=list((ROOT/'scripts').rglob('*.py'))+list((REPO/'project/scripts').glob('*pacer_stage4_rerun_v01.py'))
    for p in sources:
        text=p.read_text(encoding='utf-8'); ast.parse(text,filename=str(p))
        if ROOT in p.parents:
            assert 'C:/' not in text and 'C:\\' not in text, 'Windows execution path in package: '+str(p)
    archive=OUT/'pacer_stage4_rerun_v01_server_bundle.tar.gz'
    checks={line.split('  ',1)[1]:line.split('  ',1)[0] for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines()}
    with tarfile.open(archive) as t:
        for name,digest in checks.items():
            assert sha(ROOT/name)==digest
            assert hashlib.sha256(t.extractfile(ROOT.name+'/'+name).read()).hexdigest()==digest
        members=t.getmembers()
        assert not any(m.name.endswith(('.dcd','.chk','.pyc','.npz','.pt')) for m in members)
        assert all(not m.issym() and not m.islnk() and '..' not in Path(m.name).parts for m in members)
    before=json.loads((OUT/'logs/protected_before.json').read_text())
    for name,meta in before.items():
        p=REPO/name; assert p.is_file() and [p.stat().st_size,p.stat().st_mtime_ns]==meta
    assert sha(REPO/'project/results/pacer_xr_rerun_v01/pacer_xr_rerun_full200_ranking.csv')=='6fff7e34d94466ff288b52281b1c3a64fef4ea6cd9da84ca45bc703429a04559'
    receipt=dict(status='PASS',systems=12,jobs=36,total_ns=360,python_sources_syntax_checked=len(sources),
                 package_hashes_verified=len(checks),package_sha256=sha(archive),
                 protected_files_checked=len(before),production_started=False,
                 limitations=['Linux CUDA and MD stack execution unavailable locally; server steps pending'])
    (OUT/'logs/local_QC_pacer_stage4_rerun_v01.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__': main()
