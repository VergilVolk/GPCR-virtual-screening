"""Local Reference tests and non-executing package/adapter validation."""
import ast
import hashlib
import io
import json
import sys
import unittest
from pathlib import Path

PROJECT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT/'.pbcfix01-test-deps'))
sys.path.insert(0,str(PROJECT/'tests'))
import test_pacer_stage4_pbcfix01


def main():
    stream=io.StringIO()
    suite=unittest.defaultTestLoader.loadTestsFromModule(test_pacer_stage4_pbcfix01)
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    output=stream.getvalue()
    print(output)
    if not result.wasSuccessful(): raise SystemExit(1)
    import package_pacer_stage4_pbcfix01 as pack
    pack.main()
    bundle=pack.BUNDLE
    for line in (bundle/'SHA256SUMS_PBCFIX01.txt').read_text().splitlines():
        digest,name=line.split('  ',1)
        if pack.sha(bundle/name)!=digest: raise RuntimeError('Package hash mismatch')
    runner=(bundle/'scripts'/f'run_{pack.NEW}.py').read_text()
    scheduler=(bundle/'scripts'/f'schedule_{pack.NEW}.py').read_text()
    for p in bundle.rglob('*.py'): ast.parse(p.read_text(encoding='utf-8'))
    if 'refined_state.xml' in runner or 'scripts/legacy' in runner or 'from build_' in runner:
        raise RuntimeError('Old runtime reference remains')
    if 'authorized_gate(job)' not in runner or 'authorized_gate(job)' not in scheduler:
        raise RuntimeError('Missing authorization binding')
    import csv
    with (bundle/'manifests/STAGE4_RERUN_SYSTEM_MANIFEST.csv').open(newline='') as f: systems=list(csv.DictReader(f))
    with (bundle/'manifests/STAGE4_RERUN_JOB_MATRIX.csv').open(newline='') as f: jobs=list(csv.DictReader(f))
    if len(systems)!=12 or len(jobs)!=36: raise RuntimeError('Wrong campaign size')
    for job in jobs:
        if not all(pack.NEW in job[k] for k in ['job_id','system_id','output_path','checkpoint_path']):
            raise RuntimeError('Old job/path binding')
    output+='\nPackage: syntax, SHA256 members, 12 systems, 36 versioned jobs, runtime state and authorization anchors PASS.\n'
    output+='Local platform: OpenMM 8.1.1 Reference; isolated workspace dependency. CUDA server validation PENDING.\n'
    (pack.OUT/'LOCAL_TEST_RESULTS.txt').write_text(output,encoding='utf-8')


if __name__=='__main__': main()
