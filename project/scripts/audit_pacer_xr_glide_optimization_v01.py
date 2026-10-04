"""Read-only baseline inventory and preservation snapshot for Glide optimization."""
import json,os
from collections import Counter
from pathlib import Path
from pacer_xr_production_v01 import *
ROOT=OUTPUT/'glide_optimization_v01'

def main():
    ROOT.mkdir(exist_ok=True); baseline=ROOT/'baseline.json'
    require(not baseline.exists(),'Baseline already exists; do not overwrite')
    roster=load_roster(); assets=load_asset_lock(); jobs=OUTPUT/'glide/jobs'
    hashes={str(p):sha(p) for p in jobs.rglob('*') if p.is_file()}
    states=[]; prep=[]
    for row in roster.itertuples():
        p=jobs/('prepare_'+row.candidate_id)/'status.json'
        record=json.loads(p.read_text()) if p.exists() else {}
        prep.append(dict(candidate_id=row.candidate_id,status=record.get('status','NOT_STARTED'),artifact_hash_valid=JobStore.valid(None,record,record.get('fingerprint'))))
        for receptor in assets['receptors']:
            p=jobs/(row.candidate_id+'_'+receptor)/'status.json'; record=json.loads(p.read_text()) if p.exists() else {}
            if record.get('status')=='PASS' and JobStore.valid(None,record,record.get('fingerprint')): state='COMPLETE'
            elif record.get('status')=='FAILED': state='FAILED_RETRYABLE'
            elif p.parent.exists(): state='INCOMPLETE'
            else: state='NOT_STARTED'
            states.append(dict(candidate_id=row.candidate_id,receptor_id=receptor,state=state,error=record.get('error'),status_path=str(p)))
    frozen={str(P/a['path']):sha(P/a['path']) for a in assets['files']}
    config=P/'config/pacer_xr_references_v01.json'; frozen[str(config)]=sha(config)
    for a in json.loads(config.read_text())['channels'].values(): frozen[str(P/a['path'])]=sha(P/a['path'])
    scripts={str(p):sha(p) for p in (P/'scripts').glob('*pacer_xr*v01.py')}
    write_json(baseline,dict(source_sha256=SOURCE_SHA,logical_cpus=os.cpu_count(),prepared_statuses=prep,evaluations=states,counts=dict(Counter(x['state'] for x in states)),production_job_file_hashes=hashes,frozen_asset_hashes=frozen,initial_script_hashes=scripts))
    print(json.dumps(dict(preparations_complete=sum(x['artifact_hash_valid'] for x in prep),evaluations=dict(Counter(x['state'] for x in states)),files_preserved=len(hashes),logical_cpus=os.cpu_count()),indent=2))

if __name__=='__main__': main()
