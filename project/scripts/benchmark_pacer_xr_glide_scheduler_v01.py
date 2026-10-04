"""Bounded ranks1–10 benchmarking; never invokes production or Vina."""
import argparse,json,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from pacer_xr_production_v01 import *
from pacer_xr_glide_scheduler_v01 import *

def comparison(record,expected):
    checks=[]
    for cid,result in record['result']['evaluations'].items():
        key=cid+'_'+result['receptor_id']
        if key in expected:
            delta=abs(result['glide_gscore']-expected[key]['glide_gscore'])
            require(delta<=1e-12,'Numerical equivalence FAIL: '+key+' difference '+str(delta))
            checks.append(dict(candidate_id=cid,receptor_id=result['receptor_id'],expected=expected[key]['glide_gscore'],observed=result['glide_gscore'],difference=delta))
    return checks

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--njobs',type=int,choices=[4,8,16]); ap.add_argument('--dual',action='store_true'); ap.add_argument('--serial16',action='store_true'); ap.add_argument('--equivalence',action='store_true'); args=ap.parse_args()
    assets=load_asset_lock(); roster=load_roster(); science_contract()
    prepared={r.candidate_id:valid_preparation(r,OUTPUT/'glide',assets,science_contract()) for r in roster.itertuples()}
    baseline=json.loads((OPT/'baseline.json').read_text()); expected={}
    for item in baseline['evaluations']:
        if item['state']=='COMPLETE': expected[item['candidate_id']+'_'+item['receptor_id']]=json.loads(Path(item['status_path']).read_text())['result']
    # Acceptance is a directly comparable verified source for ranks1–3 on all grids.
    for path in (ACCEPTANCE/'glide/jobs').glob('PACERGEN*/status.json'):
        record=json.loads(path.read_text())
        if record.get('status')=='PASS' and JobStore.valid(None,record,record['fingerprint']):
            key=record['job_key']; result=record['result']
            if key not in expected: expected[key]=result
    # All ten scores, rather than only historical ranks1–3, must agree across parallelism settings.
    for path in (OPT/'experiment_single4/batches').glob('*/status.json'):
        reference=json.loads(path.read_text())
        if reference.get('status')=='PASS' and batch_valid(reference,assets,science_contract(),prepared):
            for cid,result in reference['result']['evaluations'].items():
                expected[cid+'_'+result['receptor_id']]=result
    parents=[prepared[x] for x in roster.head(10).candidate_id]
    name='dual8' if args.dual else 'serial16_two_receptors' if args.serial16 else 'equivalence' if args.equivalence else 'single'+str(args.njobs)
    benchmark_dir=OPT/'benchmarks'; benchmark_dir.mkdir(exist_ok=True)
    # Configurations need distinct immutable records even though their scientific fingerprint is identical.
    scoped=OPT/('experiment_'+name)/'batches'
    start=time.perf_counter(); records=[]
    if args.dual:
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(campaign,scoped,rec,parents,8,assets,1) for rec in ['7TRS','cluster_00']]
            records=[f.result() for f in futures]
    elif args.serial16:
        records=[campaign(scoped,rec,parents,16,assets,1) for rec in ['7TRS','cluster_00']]
    elif args.equivalence:
        # One candidate across all remaining grids; rank2 additionally tests failed jobs and prior successes.
        with ThreadPoolExecutor(max_workers=2) as pool:
            subset=[prepared[x] for x in roster.head(2).candidate_id]
            futures=[pool.submit(campaign,scoped,receptor,subset,8,assets,1) for receptor in list(assets['receptors'])[1:]]
            records=[f.result() for f in futures]
    else:
        require(args.njobs is not None,'Specify --njobs, --dual or --equivalence')
        records=[campaign(scoped,'7TRS',parents,args.njobs,assets,1)]
    elapsed=time.perf_counter()-start; checks=[]
    for record in records:
        checks.extend(comparison(record,expected))
        # Retain the successful, equivalent benchmark record in the completion catalog, without touching legacy jobs.
        key=record['inputs']['receptor_id']+'_'+name
        write_json(benchmark_dir/key/'status.json',record)
    row=dict(configuration=name,campaign_concurrency=2 if args.dual or args.equivalence else 1,njobs_per_campaign=8 if args.dual or args.equivalence else 16 if args.serial16 else args.njobs,candidates_per_campaign=2 if args.equivalence else 10,receptor_campaigns=len(records),wall_seconds=elapsed,successful_parent_evaluations=sum(len(r['result']['evaluations']) for r in records),failures=0,cpu_mean=float(np.mean([r['result']['metrics']['system_cpu_percent_mean'] for r in records])),checks=checks,records=[str(Path(r['work']).parent/'status.json') for r in records])
    write_json(OPT/('benchmark_'+name+'.json'),row)
    print(json.dumps(row,indent=2))

if __name__=='__main__': main()
