"""Execution-only receptor campaigns; immutable legacy caches and bounded retries."""
from __future__ import annotations
import ast,json,os,subprocess,time,uuid
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import pandas as pd
from pacer_xr_production_v01 import *

OPT=OUTPUT/'glide_optimization_v01'
BRIDGE=P/'scripts/pacer_xr_glide_batch_bridge_v01.py'
SCIENCE=P/'config/pacer_xr_glide_science_v01.json'
SCHEDULER=P/'config/pacer_xr_glide_scheduler_v01.json'

def scheduler_contract():
    config=json.loads(SCHEDULER.read_text())
    require(1<=config['campaign_concurrency']<=2 and 1<=config['njobs_per_campaign']<=16 and config['campaign_concurrency']*config['njobs_per_campaign']<=16,'Glide scheduler exceeds tested resource budget')
    require(config['batch_bridge_sha256']==sha(BRIDGE),'Batch identity/parser bridge changed since equivalence testing')
    return config

def function_hash(path,names):
    tree=ast.parse(Path(path).read_text())
    functions={node.name:ast.dump(node,include_attributes=False) for node in tree.body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef)) and node.name in names}
    require(set(functions)==set(names),'Scientific function missing')
    return fingerprint(functions)

def science_contract():
    lock=json.loads(SCIENCE.read_text())
    require(function_hash(P/'scripts/run_pacer_xr_glide_rerun_v01.py',['prepare_candidate'])==lock['preparation_function_ast_sha256'],'Preparation chemistry code changed')
    require(sha(P/'scripts/pacer_xr_production_v01.py')==lock['shared_contract_sha256'],'Scoring/reduction contract changed')
    require(sha(P/'scripts/pacer_xr_schrodinger_v01.py')==lock['original_parser_sha256'],'Scientific parser changed')
    require(sha(P/'config/pacer_xr_assets_v01.json')==lock['assets_lock_sha256'],'Frozen assets changed')
    return lock

def freeze_science():
    baseline=json.loads((OPT/'baseline.json').read_text())
    data=dict(preparation_function_ast_sha256=function_hash(P/'scripts/run_pacer_xr_glide_rerun_v01.py',['prepare_candidate']),shared_contract_sha256=sha(P/'scripts/pacer_xr_production_v01.py'),original_parser_sha256=sha(P/'scripts/pacer_xr_schrodinger_v01.py'),assets_lock_sha256=sha(P/'config/pacer_xr_assets_v01.json'),legacy_runner_sha256=baseline['initial_script_hashes'][str(P/'scripts/run_pacer_xr_glide_rerun_v01.py')],version='Glide SP production v01; scheduling excluded from scientific identity')
    if SCIENCE.exists(): require(json.loads(SCIENCE.read_text())==data,'Science lock differs')
    else: write_json(SCIENCE,data)
    return data

def bridge(work,mode,meta):
    write_json(work/'metadata.json',meta)
    run_command(work,mode,[SCH/'run.exe',BRIDGE,mode,work/'metadata.json'],True)
    return json.loads((work/(mode+'_result.json' if mode!='parse' else 'batch_poses.json')).read_text())

def receptor_meta(name,assets):
    rec=assets['receptors'][name]
    return dict(receptor_id=name,grid_sha256=sha(P/rec['grid']),center=rec['center'],inner=rec['inner'],outer=rec['outer'])

def valid_preparation(row,root,assets,science):
    status=root/'jobs'/('prepare_'+row.candidate_id)/'status.json'
    require(status.exists(),'Missing prepared manifest: '+row.candidate_id)
    record=json.loads(status.read_text()); inputs=record['inputs']; result=record['result']
    require(JobStore.valid(None,record,record['fingerprint']) and fingerprint(inputs)==record['fingerprint'],'Prepared artifact integrity failed: '+row.candidate_id)
    require(inputs['candidate_id']==row.candidate_id and inputs['canonical_smiles']==row.canonical_smiles and inputs['rank']==int(row.rerun_final_rank) and inputs['roster_sha256']==SOURCE_SHA,'Prepared parent identity differs')
    require(inputs['protocol']==assets['protocol'] and inputs['ligprep_executable_sha256']==assets['schrodinger_executables']['ligprep.exe'],'Preparation contract differs')
    require(inputs['codes']['pacer_xr_production_v01.py']==science['shared_contract_sha256'] and inputs['codes']['pacer_xr_schrodinger_v01.py']==science['original_parser_sha256'],'Preparation science provenance differs')
    # New preparation attempts use the same frozen function; legacy full-run caches retain original runner hash.
    require(inputs['codes']['run_pacer_xr_glide_rerun_v01.py'] in [science['legacy_runner_sha256'],sha(P/'scripts/run_pacer_xr_glide_rerun_v01.py')],'Unapproved preparation runner provenance')
    require(result['candidate_id']==row.candidate_id and result['canonical_smiles']==row.canonical_smiles and result['rerun_final_rank']==int(row.rerun_final_rank),'Prepared manifest mapping differs')
    require(bool(result['variants']) and len({v['prepared_variant_id'] for v in result['variants']})==len(result['variants']),'Duplicate or empty prepared variants')
    require(all(v['candidate_id']==row.candidate_id and v['canonical_smiles']==row.canonical_smiles and v['rerun_final_rank']==int(row.rerun_final_rank) and v['prepared_title']==v['prepared_variant_id'] and v['prepared_variant_id'].startswith(v['molscrub_state_id']+'__v') for v in result['variants']),'Invalid prepared state/variant mapping')
    require(Path(result['ligand']).resolve().as_posix() in {Path(x).resolve().as_posix() for x in record['artifacts']},'Prepared ligand lacks hash provenance')
    return result

def legacy_compatible(record,row,receptor,prepared,assets,science):
    require(JobStore.valid(None,record,record.get('fingerprint')) and fingerprint(record['inputs'])==record['fingerprint'],'Docking artifacts invalid')
    inputs=record['inputs']; result=record['result']; rec=receptor_meta(receptor,assets)
    require(inputs['candidate_id']==row.candidate_id and inputs['canonical_smiles']==row.canonical_smiles and inputs['receptor_id']==receptor,'Docking parent/receptor mismatch')
    require(inputs['grid_sha256']==rec['grid_sha256'] and inputs['center']==rec['center'] and inputs['ligand_sha256']==sha(prepared['ligand']),'Docking input hashes differ')
    require(inputs['protocol']==assets['protocol'] and inputs['glide_executable_sha256']==assets['schrodinger_executables']['glide.exe'],'Docking method differs')
    require(inputs['codes']['run_pacer_xr_glide_rerun_v01.py']==science['legacy_runner_sha256'] and inputs['codes']['pacer_xr_production_v01.py']==science['shared_contract_sha256'] and inputs['codes']['pacer_xr_schrodinger_v01.py']==science['original_parser_sha256'],'Legacy science provenance differs')
    require(result['command']['exit_code']==0 and result['status']=='PASS' and math.isfinite(result['glide_gscore']),'Docking exit/score invalid')
    pv=next(Path(x) for x in record['artifacts'] if Path(x).name=='dock_pv.maegz')
    return dict(candidate_id=row.candidate_id,poseviewer=str(pv),variants=prepared['variants'],receptor=rec),result

def audit_index(root,roster,assets,allow_prepare=False):
    science=science_contract(); prepared={}; invalid=[]
    for row in roster.itertuples():
        try: prepared[row.candidate_id]=valid_preparation(row,root,assets,science)
        except (ValueError,FileNotFoundError,KeyError,json.JSONDecodeError) as exc: invalid.append((row,str(exc)))
    if invalid and allow_prepare:
        from run_pacer_xr_glide_rerun_v01 import prepare_candidate
        store=JobStore(root/'jobs')
        for row,_ in invalid:
            prepare_candidate(row,store,assets,code_hashes())
            prepared[row.candidate_id]=valid_preparation(row,root,assets,science)
    else: require(not invalid,'Invalid preparation caches: '+str([(x.candidate_id,e) for x,e in invalid]))
    print(f'Glide preparation cached: {len(roster)-len(invalid)}/{len(roster)}; new preparations: {len(invalid) if allow_prepare else 0}',flush=True)
    audit_work=OPT/('cache_audit_'+uuid.uuid4().hex[:10]); audit_work.mkdir(parents=True)
    index={}; checks=[]; old_results={}
    for row in roster.itertuples():
        for receptor in assets['receptors']:
            key=row.candidate_id+'_'+receptor; status=root/'jobs'/key/'status.json'
            record=json.loads(status.read_text()) if status.exists() else {}
            state='NOT_STARTED'; error=None
            if record.get('status')=='PASS':
                try:
                    check,result=legacy_compatible(record,row,receptor,prepared[row.candidate_id],assets,science)
                    checks.append(check); old_results[key]=result; state='COMPLETE'
                except (ValueError,KeyError,StopIteration) as exc: state='INCOMPLETE'; error=str(exc)
            elif record.get('status')=='FAILED': state='FAILED_RETRYABLE'; error=record.get('error')
            elif status.parent.exists(): state='INCOMPLETE'
            index[key]=dict(candidate_id=row.candidate_id,receptor_id=receptor,state=state,error=error,status_path=str(status))
    audit=bridge(audit_work,'audit',dict(prepared=list(prepared.values()),completed=checks))
    for key,poses in audit['completed_poses'].items():
        require(abs(reduce_glide(poses)['glide_gscore']-old_results[key]['glide_gscore'])<=1e-12,'Cached score differs from actual poseviewer')
    # A batch may cover only the pending parents for its receptor; all covered hashes must still validate.
    roots=[root/'batches',OPT/'benchmarks']
    batch_records=[]
    for batchroot in roots:
        for path in sorted(batchroot.glob('*/status.json')) if batchroot.exists() else []:
            record=json.loads(path.read_text())
            if record.get('status')!='PASS': continue
            if not batch_valid(record,assets,science,prepared): continue
            batch_records.append(record)
            for cid,result in record['result']['evaluations'].items():
                key=cid+'_'+result['receptor_id']
                if key not in index: continue
                if key in old_results:
                    require(abs(result['glide_gscore']-old_results[key]['glide_gscore'])<=1e-12,'Batch/legacy score differs')
                    continue
                old_results[key]=result
                index[key].update(state='COMPLETE',batch_status_path=str(path),status_path=str(path))
    write_json(OPT/'completion_index.json',dict(root=str(root),counts=dict(Counter(x['state'] for x in index.values())),preparations_cached=len(prepared),new_preparations=len(invalid) if allow_prepare else 0,evaluations=index,scientific_signature=fingerprint(science),audit_output=str(audit_work/'audit_result.json')))
    return prepared,index,old_results

def batch_valid(record,assets,science,prepared):
    try:
        if not JobStore.valid(None,record,record.get('fingerprint')) or fingerprint(record['inputs'])!=record['fingerprint']: return False
        inputs=record['inputs']; rec=receptor_meta(inputs['receptor_id'],assets)
        if inputs['science']!=fingerprint(science) or inputs['grid_sha256']!=rec['grid_sha256']: return False
        if any(cid not in prepared or sha(prepared[cid]['ligand'])!=digest for cid,digest in inputs['ligands'].items()): return False
        if set(record['result']['evaluations'])!=set(inputs['ligands']): return False
        return all(x['status']=='PASS' and math.isfinite(x['glide_gscore']) and x['command']['exit_code']==0 and x['candidate_id']==cid and x['receptor_id']==inputs['receptor_id'] and x['canonical_smiles']==prepared[cid]['canonical_smiles'] for cid,x in record['result']['evaluations'].items())
    except (KeyError,ValueError,FileNotFoundError): return False

def host_identity(work):
    """Unique Job Control host name per campaign work directory.

    Execution metadata only.  Two campaigns that both register a host named
    ``localhost`` with different ``tmpdir`` values collide in the shared Job
    Control database (``%LOCALAPPDATA%/Schrodinger/.jobdb2``); the later
    campaign's subjobs are then dispatched against the other campaign's temp
    directory and die in ``[transfer] [local_copy]``.  ``host: localhost`` still
    targets this machine and no scientific setting changes.
    """
    return 'pacer_host_'+fingerprint(('host',str(Path(work).resolve())))[:12]

def run_batch_command(work,njobs):
    env=os.environ.copy(); env.update(PYTHONDONTWRITEBYTECODE='1',SCHRODINGER=str(SCH),SCHRODINGER_TMPDIR=str(work/'tmp'),SCHRODINGER_HOSTS=str(work/'schrodinger.hosts'),OMP_NUM_THREADS='1')
    (work/'tmp').mkdir(exist_ok=True)
    hostname=host_identity(work)
    (work/'schrodinger.hosts').write_text('name: '+hostname+'\nhost: localhost\nprocessors: '+str(njobs)+'\ntmpdir: '+str(work/'tmp')+'\nschrodinger: '+str(SCH)+'\n')
    jobname='pacer_'+work.parent.name+'_'+work.name+'_'+uuid.uuid4().hex[:10]
    argv=[str(SCH/'glide.exe'),'campaign.in','-WAIT','-LOCAL','-HOST',f'{hostname}:{njobs}','-NJOBS',str(njobs),'-JOBNAME',jobname,'-max_retries','0','-no_cleanup','-noforce']
    import psutil
    start=time.perf_counter(); samples=[]; peak=0; proc=None; interrupted=False
    try:
        with (work/'glide.console.log').open('w',encoding='utf-8') as log:
            proc=subprocess.Popen(argv,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT)
            psutil.cpu_percent(interval=None)
            while proc.poll() is None:
                samples.append(psutil.cpu_percent(interval=1))
                try: peak=max(peak,len(psutil.Process(proc.pid).children(recursive=True)))
                except psutil.Error: pass
                if time.perf_counter()-start>7200: proc.terminate(); raise TimeoutError('Batch exceeded 7200 seconds')
            code=proc.returncode
    except BaseException:
        interrupted=True
        if proc is not None and proc.poll() is None: proc.terminate()
        raise
    finally:
        record=dict(argv=argv,cwd=str(work),exit_code=None if proc is None else proc.poll(),elapsed_seconds=time.perf_counter()-start,system_cpu_percent_mean=float(np.mean(samples)) if samples else None,system_cpu_percent_peak=max(samples) if samples else None,peak_descendant_processes=peak,njobs=njobs,interrupted=interrupted,log=str(work/'glide.console.log'))
        write_json(work/'glide.command.json',record)
    require(code==0,'Batch failed: '+str(work/'glide.console.log'))
    return record,jobname

def remaining_campaign_limit(batchroot,legacy_root,receptor,parents,inputs):
    """Three total production attempts per parent/receptor, including interrupted legacy attempts."""
    this_key=receptor+'_'+fingerprint(inputs)[:14]; limits=[]
    for parent in parents:
        cid=parent['candidate_id']; used=len(list((legacy_root/'jobs'/(cid+'_'+receptor)).glob('attempt_*')))
        for status in batchroot.glob('*/status.json') if batchroot.exists() else []:
            old=json.loads(status.read_text())
            if status.parent.name!=this_key and old.get('inputs',{}).get('receptor_id')==receptor and cid in old.get('inputs',{}).get('ligands',{}): used+=int(old.get('attempt',0))
        limits.append(3-used)
    require(min(limits)>0,'Retry budget exhausted for parent/receptor; manual diagnosis required')
    return min(limits)

def campaign(root,receptor,parents,njobs,assets,maximum_attempts=3,legacy_root=None):
    require(1<=njobs<=16 and parents,'Invalid campaign size or parallelism')
    science=science_contract(); rec=receptor_meta(receptor,assets)
    inputs=dict(science=fingerprint(science),source_sha256=SOURCE_SHA,receptor_id=receptor,grid_sha256=rec['grid_sha256'],ligands={p['candidate_id']:sha(p['ligand']) for p in parents})
    key=receptor+'_'+fingerprint(inputs)[:14]; jobroot=Path(root)/key; jobroot.mkdir(parents=True,exist_ok=True); status=jobroot/'status.json'
    if status.exists():
        old=json.loads(status.read_text())
        if batch_valid(old,assets,science,{p['candidate_id']:p for p in parents}): return old
    if legacy_root is not None: maximum_attempts=min(maximum_attempts,remaining_campaign_limit(Path(root),legacy_root,receptor,parents,inputs))
    attempted=max([int(p.name[-4:]) for p in jobroot.glob('attempt_*')]+[0])
    require(attempted<maximum_attempts,'Campaign retry budget exhausted: '+key)
    for number in range(attempted+1,maximum_attempts+1):
        work=jobroot/f'attempt_{number:04d}'; work.mkdir()
        record=dict(status='RUNNING',fingerprint=fingerprint(inputs),inputs=inputs,attempt=number,work=str(work),execution=dict(njobs=max(1,njobs//(2**(number-attempted-1))),maximum_attempts=maximum_attempts),started_utc=time.time())
        write_json(status,record)
        phase='pack'
        try:
            pack=bridge(work,'pack',dict(prepared=parents))
            (work/'campaign.in').write_text('GRIDFILE '+(P/assets['receptors'][receptor]['grid']).as_posix()+'\nLIGANDFILE '+(work/'library.maegz').as_posix()+'\nPRECISION SP\nPOSES_PER_LIG 5\nPOSTDOCK_NPOSE 5\nFORCEFIELD OPLS_2005\n')
            phase='dock'; cmd,jobname=run_batch_command(work,record['execution']['njobs'])
            phase='parse'
            pv=work/(jobname+'_pv.maegz'); require(pv.is_file(),'Successful campaign missing poseviewer: '+str(pv))
            poses=bridge(work,'parse',dict(prepared=parents,poseviewer=str(pv),receptor=rec))
            evaluations={}
            for parent in parents:
                cid=parent['candidate_id']; best=reduce_glide(poses[cid])
                evaluations[cid]=dict(candidate_id=cid,canonical_smiles=parent['canonical_smiles'],receptor_id=receptor,cluster_id=None if receptor=='7TRS' else int(receptor[-2:]),glide_gscore=best['glide_gscore'],selected_variant_id=best['prepared_variant_id'],selected_entry_index=best['entry_index'],poses=poses[cid],selected_pose=best,n_poses=len(poses[cid]),command=cmd,status='PASS',poseviewer=str(pv))
            write_json(work/'evaluations.json',evaluations)
            artifacts=[pv,work/'library.maegz',work/'batch_poses.json',work/'evaluations.json',work/'glide.command.json',work/'campaign.in',work/'metadata.json']
            record.update(status='PASS',completed_utc=time.time(),result=dict(evaluations=evaluations,metrics=cmd),artifacts={str(p):sha(p) for p in artifacts})
            write_json(status,record); return record
        except BaseException as exc:
            interrupted=isinstance(exc,(KeyboardInterrupt,SystemExit))
            record.update(status='INCOMPLETE' if interrupted else 'FAILED',error=str(exc),completed_utc=time.time())
            write_json(status,record)
            if interrupted or number==maximum_attempts: raise
            log=(work/'glide.console.log').read_text(errors='replace') if (work/'glide.console.log').exists() else ''
            transient=any(x in log.lower() for x in ['control-c','interrupt','license','job control','jobcontrol','connection','lock','temporar','resource']) or isinstance(exc,TimeoutError)
            require(phase=='dock' and transient,'Scientific/parser failure is not automatically retried: '+str(exc))
            print('Retry execution-only '+key+' after '+str(exc),flush=True)
    raise RuntimeError('Retry budget exhausted')

def collect_complete(root,roster,assets):
    prepared,index,results=audit_index(root,roster,assets,False)
    require(len(results)==11*len(roster) and all(x['state']=='COMPLETE' for x in index.values()),'Glide coverage incomplete; fusion prohibited')
    return prepared,index,results

def run_optimized(acceptance=False,workers=None):
    roster=select_scope(acceptance); assets=load_asset_lock(); verify_environment(assets,'glide'); science_contract()
    config=scheduler_contract(); root=output_root(acceptance)/'glide'
    prepared,index,results=audit_index(root,roster,assets,True)
    print('Glide docking cached: '+str(len(results))+'/'+str(len(roster)*11),flush=True)
    plans=[(rec,[prepared[r.candidate_id] for r in roster.itertuples() if index[r.candidate_id+'_'+rec]['state']!='COMPLETE']) for rec in assets['receptors']]
    plans=[x for x in plans if x[1]]; errors=[]
    with ThreadPoolExecutor(max_workers=config['campaign_concurrency']) as pool:
        futures={pool.submit(campaign,root/'batches',rec,parents,config['njobs_per_campaign'],assets,3,root):rec for rec,parents in plans}
        for future in as_completed(futures):
            try: record=future.result(); print('Glide receptor campaign PASS '+futures[future]+' '+str(len(record['result']['evaluations']))+' parents',flush=True)
            except Exception as exc: errors.append(dict(receptor=futures[future],error=str(exc))); print('Glide campaign FAILED '+futures[future]+': '+str(exc),flush=True)
    write_json(root/'execution_summary.json',dict(scope='acceptance_rank1_3' if acceptance else 'full200',execution_model='receptor_major',scheduler=config,errors=errors,source_sha256=SOURCE_SHA))
    require(not errors,'Glide campaign failed; resume from intact completed batches')
    prepared,index,results=collect_complete(root,roster,assets)
    emit_channels(root,roster,results,config)
    return pd.read_csv(root/'glide_channels.csv')

def emit_channels(root,roster,results,config):
    channels=[]; clusters=[]; poses=[]
    for row in roster.itertuples():
        pdb=results[row.candidate_id+'_7TRS']; ensemble=[results[row.candidate_id+f'_cluster_{i:02d}'] for i in range(10)]; be=ensemble_score(ensemble,'glide_gscore')
        channels.append(dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,rerun_final_rank=int(row.rerun_final_rank),Glide_PDB_raw=pdb['glide_gscore'],Glide_BEmin_raw=be['BEmin'],Glide_BEavg_raw=be['BEavg'],Glide_BEmin_cluster=be['winning_cluster'],n_clusters=10,status='PASS'))
        for c in ensemble: clusters.append(dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,cluster_id=c['cluster_id'],glide_gscore=c['glide_gscore'],selected_variant_id=c['selected_variant_id'],selected_entry_index=c['selected_entry_index'],pmf=load_pmf()[c['cluster_id']],BE_i=be['BE_i'][c['cluster_id']]))
        for c in [pdb]+ensemble: poses.extend(c['poses'])
    frame=pd.DataFrame(channels).sort_values('rerun_final_rank'); validate_scores(frame,roster,CHANNELS[:3])
    write_csv(root/'glide_all_poses.csv',pd.DataFrame(poses).sort_values(['candidate_id','receptor_id','pose_rank']))
    write_csv(root/'glide_cluster_scores.csv',pd.DataFrame(clusters)); write_csv(root/'glide_channels.csv',frame)
    write_json(root/'channels_provenance.json',dict(source_sha256=SOURCE_SHA,codes=code_hashes(),channel_sha256=sha(root/'glide_channels.csv'),complete=True,n_candidates=len(roster),evaluations=len(results),assets_lock_sha256=sha(P/'config/pacer_xr_assets_v01.json'),execution_model='receptor_major',scientific_signature=fingerprint(science_contract()),scheduler=config))
