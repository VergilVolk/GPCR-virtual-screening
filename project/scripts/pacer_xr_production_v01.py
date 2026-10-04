"""Shared frozen contracts; no docking or result writes occur at import time."""
from __future__ import annotations
import hashlib, json, math, os, subprocess, threading, time
from pathlib import Path
from typing import Callable
import numpy as np
import pandas as pd
from rdkit import Chem
from apply_pacer_xr_candidates import empirical_rank

P = Path(__file__).resolve().parents[1]
REPO = P.parent
SOURCE = P / 'results/pacer_rerun_v01/module2_top200_output.csv'
SOURCE_REF = 'origin/codex/drugclip-pacer-handoff'
SOURCE_SHA = '41c5a69e6c54825d0cbe9149a2f1fe9dbe5ea4ae2de137c161450607c3fa28eb'
OUTPUT = P / 'results/pacer_xr_rerun_v01'
ACCEPTANCE = OUTPUT / 'acceptance_rank1_3'
SCH = Path(r'C:\Program Files\Schrodinger2023-1')
CHANNELS = ('Glide_PDB','Glide_BEmin','Glide_BEavg','Vina_PDB','Vina_BEmin','Vina_BEavg')
ANCHORS = {'PACER0010':('PACERGEN01755',17),'PACER0073':('PACERGEN00094',60),'PACER0027':('PACERGEN00123',114)}

def require(condition, message):
    if not condition: raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def fingerprint(data): return hashlib.sha256(json.dumps(data,sort_keys=True,allow_nan=False).encode()).hexdigest()

def write_json(path,data):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.writing')
    temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8'); os.replace(temp,path)

def write_csv(path, frame):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.writing'); frame.to_csv(temp,index=False); os.replace(temp,path)

def canonical(smiles):
    require(isinstance(smiles,str) and bool(smiles.strip()),'Missing SMILES')
    mol=Chem.MolFromSmiles(smiles)
    require(mol is not None,'Invalid SMILES')
    return Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)

def validate_roster(frame,expected=200):
    require(set(['gen_id','smiles','score','rankpct','final_rank']).issubset(frame),'Rerun schema missing fields')
    require(len(frame)==expected,'Rerun row count mismatch')
    out=frame.copy()
    require(out.gen_id.notna().all() and out.gen_id.astype(str).str.fullmatch(r'PACERGEN\d{5}').all(),'Invalid candidate IDs')
    require(not out.gen_id.duplicated().any(),'Duplicate candidate IDs')
    out['candidate_id']=out.gen_id
    out['canonical_smiles']=out.smiles.map(canonical)
    require(not out.canonical_smiles.duplicated().any(),'Duplicate canonical structures')
    ranks=pd.to_numeric(out.final_rank,errors='raise')
    require(sorted(ranks.tolist())==list(range(1,expected+1)),'Rerun ranks must be exactly 1..200')
    out['rerun_final_rank']=ranks.astype(int)
    require(np.isfinite(out[['score','rankpct']].to_numpy(float)).all(),'Nonfinite Module-2 metadata')
    return out.sort_values('rerun_final_rank',kind='stable').reset_index(drop=True)

def load_roster():
    if not SOURCE.exists():
        data=subprocess.check_output(['git','-c','safe.directory='+REPO.as_posix(),'show',SOURCE_REF+':project/results/pacer_rerun_v01/module2_top200_output.csv'],cwd=REPO)
        require(hashlib.sha256(data).hexdigest()==SOURCE_SHA,'Rerun source SHA256 mismatch')
        SOURCE.parent.mkdir(parents=True,exist_ok=True); SOURCE.write_bytes(data)
    require(sha(SOURCE)==SOURCE_SHA,'Rerun source SHA256 mismatch; FAIL CLOSED')
    return validate_roster(pd.read_csv(SOURCE))

def verify_anchors(roster):
    historical=pd.read_csv(P/'results/pacer_candidates_v01/predock_portfolio.csv')
    result=[]
    for legacy,(candidate,rank) in ANCHORS.items():
        old=historical[historical.candidate_id==legacy]; new=roster[roster.candidate_id==candidate]
        require(len(old)==len(new)==1,'Legacy anchor absent or ambiguous')
        require(canonical(old.iloc[0].canonical_smiles)==new.iloc[0].canonical_smiles,'Legacy structural mapping mismatch')
        require(int(new.iloc[0].rerun_final_rank)==rank,'Legacy rerun rank mismatch')
        result.append(dict(legacy_id=legacy,candidate_id=candidate,rerun_final_rank=rank,canonical_smiles=new.iloc[0].canonical_smiles))
    return result

def load_pmf():
    from analyze_m4_gamd_ensemble import load_pmf as existing_load, PMF_FILE
    pmf=existing_load(PMF_FILE)
    audit=json.loads((P/'results/m4_gamd_ensemble/ensemble_preparation_audit.json').read_text())
    require(set(pmf)==set(range(10)),'PMF must contain exactly ten clusters')
    for c in audit['clusters']:
        require(math.isfinite(pmf[c['cluster']]) and abs(pmf[c['cluster']]-c['pmf_kcal_mol'])<1e-12,'PMF differs from frozen audit')
    return pmf

def load_asset_lock():
    path=P/'config/pacer_xr_assets_v01.json'; lock=json.loads(path.read_text())
    for asset in lock['files']:
        require(sha(P/asset['path'])==asset['sha256'],'Frozen asset hash mismatch: '+asset['path'])
    load_pmf()
    return lock

def load_references(lock_path=None):
    lock=json.loads(Path(lock_path or P/'config/pacer_xr_references_v01.json').read_text())
    require(set(lock['channels'])==set(CHANNELS),'Six frozen reference distributions required')
    arrays={}
    for channel in CHANNELS:
        rec=lock['channels'][channel]; path=P/rec['path']
        require(path.is_file() and sha(path)==rec['sha256'],'REFERENCE_RESTORE_BLOCKED: '+channel)
        values=pd.read_csv(path,usecols=[rec['column']])[rec['column']].to_numpy(float)
        values=values[np.isfinite(values)]
        require(len(values)>0 and len(values)==rec['finite_rows'],'Empty or changed reference array: '+channel)
        arrays[channel]=values
    return arrays,lock

def verify_environment(assets,engine):
    import importlib.metadata as metadata
    for package,version in assets['preparation_versions'].items():
        require(metadata.version(package)==version,'Preparation version mismatch: '+package)
    if engine=='glide':
        for name,digest in assets['schrodinger_executables'].items():
            require(sha(SCH/name)==digest,'Schrodinger executable hash changed: '+name)
    else:
        for name,digest in assets['vina_executables'].items():
            require(sha(P/'tools'/name)==digest,'Vina executable hash changed: '+name)

def validate_scores(frame,roster,channels=CHANNELS):
    required=['candidate_id','canonical_smiles']+[x+'_raw' for x in channels]
    require(set(required).issubset(frame),'Missing raw channel columns')
    require(len(frame)==len(roster),'Incomplete candidate coverage')
    require(frame.candidate_id.notna().all() and not frame.candidate_id.duplicated().any(),'Duplicate or missing candidate IDs')
    require(set(frame.candidate_id)==set(roster.candidate_id),'Candidate universe differs from frozen roster')
    expected=roster.set_index('candidate_id').canonical_smiles
    require(all(canonical(r.canonical_smiles)==expected[r.candidate_id] for r in frame.itertuples()),'Canonical structure mismatch')
    require(np.isfinite(frame[[x+'_raw' for x in channels]].to_numpy(float)).all(),'Missing or nonfinite raw scores')
    return frame

def assemble(roster,glide,vina):
    validate_scores(glide,roster,CHANNELS[:3]); validate_scores(vina,roster,CHANNELS[3:])
    result=roster.copy()
    for scores,channels in [(glide,CHANNELS[:3]),(vina,CHANNELS[3:])]:
        result=result.merge(scores[['candidate_id']+[x+'_raw' for x in channels]],on='candidate_id',validate='one_to_one',how='left')
    validate_scores(result,roster)
    return result

def ensemble_score(rows,score_key):
    require(len(rows)==10 and {int(r['cluster_id']) for r in rows}==set(range(10)),'Incomplete or duplicated cluster coverage; no imputation')
    pmf=load_pmf(); values=[]
    for r in sorted(rows,key=lambda x:x['cluster_id']):
        value=float(r[score_key]); require(math.isfinite(value),'Nonfinite cluster score')
        values.append(value+pmf[int(r['cluster_id'])])
    return dict(BEmin=min(values),BEavg=float(np.mean(values)),winning_cluster=int(np.argmin(values)),BE_i=values)

def reduce_glide(poses):
    require(bool(poses),'No returned Glide poses')
    require(all(math.isfinite(float(p['glide_gscore'])) for p in poses),'Nonfinite Glide score')
    return min(poses,key=lambda p:(float(p['glide_gscore']),p['prepared_variant_id'],int(p['entry_index'])))

def parse_vina(path):
    values=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('REMARK VINA RESULT:'): values.append(float(line.split()[3]))
    require(bool(values) and all(math.isfinite(x) for x in values),'Missing or nonfinite Vina affinity')
    # Existing protocol reads rank-1 affinity, not an auxiliary energy component.
    return values[0]

def reduce_vina(states):
    require(bool(states) and all(math.isfinite(float(s['vina_affinity'])) for s in states),'Invalid Vina state coverage')
    return min(states,key=lambda s:(float(s['vina_affinity']),s['state_id']))

def fuse(scores,roster,references):
    validate_scores(scores,roster)
    require(set(references)==set(CHANNELS),'All six references required')
    out=scores.set_index('candidate_id').loc[roster.candidate_id].reset_index()
    for channel in CHANNELS:
        ref=np.asarray(references[channel],float)
        require(len(ref)>0 and np.isfinite(ref).all(),'Invalid reference array')
        out[channel+'_rank']=empirical_rank(out[channel+'_raw'],ref)
    out['PACER_XR']=out[[x+'_rank' for x in CHANNELS]].mean(axis=1)
    threshold=float(np.quantile(references['Glide_BEmin'],.01,method='higher'))
    out['cascade_head_top1pct']=out.Glide_BEmin_raw<=threshold
    out['cascade_stage']=np.where(out.cascade_head_top1pct,1,2)
    out['cascade_stage_score']=np.where(out.cascade_head_top1pct,out.Glide_BEmin_rank,out.PACER_XR)
    out=out.sort_values(['cascade_stage','cascade_stage_score'],ascending=[True,False],kind='stable').reset_index(drop=True)
    out['candidate_cascade_rank']=np.arange(1,len(out)+1)
    return out,threshold

def natural_top20(ranked):
    require(len(ranked)>=20 and sorted(ranked.candidate_cascade_rank)==list(range(1,len(ranked)+1)),'Invalid complete ranking')
    require(not ranked.candidate_id.duplicated().any() and not ranked.canonical_smiles.map(canonical).duplicated().any(),'Duplicate ranking identity')
    return ranked.sort_values('candidate_cascade_rank',kind='stable').head(20).copy()

def shortlist20(ranked,anchors=ANCHORS):
    natural=natural_top20(ranked); selected=set(natural.candidate_id); anchor_ids={x[0] for x in anchors.values()}
    require(anchor_ids.issubset(set(ranked.candidate_id)),'Required legacy anchor missing from full ranking')
    cascade=ranked.set_index('candidate_id').candidate_cascade_rank.to_dict()
    for legacy,(candidate,_) in anchors.items():
        if candidate not in selected:
            selected.add(candidate)
            removable=selected-anchor_ids
            require(bool(removable),'No non-anchor available for removal')
            selected.remove(max(removable,key=lambda cid:cascade[cid]))
    out=ranked[ranked.candidate_id.isin(selected)].sort_values('candidate_cascade_rank',kind='stable').copy()
    require(len(out)==20 and out.canonical_smiles.map(canonical).nunique()==20,'Shortlist must contain exactly 20 unique canonical structures')
    reverse={v[0]:k for k,v in anchors.items()}; original=set(natural.candidate_id)
    out.insert(0,'shortlist_position',np.arange(1,21))
    out['legacy_id']=out.candidate_id.map(reverse).fillna('')
    out['anchor_required']=out.candidate_id.isin(anchor_ids)
    out['selection_reason']=['legacy_anchor_forced' if cid not in original else 'natural_pacer_xr_top20_legacy_anchor' if cid in anchor_ids else 'natural_pacer_xr_top20' for cid in out.candidate_id]
    return out

def code_hashes():
    paths=list((P/'scripts').glob('*pacer_xr*v01.py'))+[P/'scripts/dock_m4_gamd_ensemble.py',P/'scripts/apply_pacer_xr_candidates.py']
    return {p.name:sha(p) for p in paths}

class JobStore:
    """Each attempt is immutable; resume requires matching inputs and all output hashes."""
    def __init__(self,root):
        self.root=Path(root); self.root.mkdir(parents=True,exist_ok=True); self.lock=threading.Lock()
        self.counts={'resumed':0,'executed':0,'failed':0}

    def valid(self,record,token):
        if record.get('status')!='PASS' or record.get('fingerprint')!=token: return False
        for p,digest in record.get('artifacts',{}).items():
            if not Path(p).is_file() or sha(p)!=digest: return False
        return bool(record.get('artifacts')) and 'result' in record

    def execute(self,key,inputs,operation:Callable):
        require(all(c.isalnum() or c in '_-' for c in key),'Unsafe job key')
        token=fingerprint(inputs); job=self.root/key; job.mkdir(exist_ok=True); status=job/'status.json'
        old=json.loads(status.read_text()) if status.exists() else {}
        if self.valid(old,token):
            with self.lock: self.counts['resumed']+=1
            return old['result']
        attempt_number=max([int(p.name.split('_')[-1]) for p in job.glob('attempt_*')]+[0])+1
        work=job/f'attempt_{attempt_number:04d}'; work.mkdir()
        record=dict(job_key=key,fingerprint=token,inputs=inputs,status='RUNNING',attempt=attempt_number,started_utc=time.time(),work=str(work))
        def publish():
            with self.lock:
                write_json(status,record)
                with (self.root/'ledger.jsonl').open('a',encoding='utf-8') as f: f.write(json.dumps(record,allow_nan=False)+'\n')
        publish()
        try:
            result,artifacts=operation(work)
            require(bool(artifacts) and all(Path(p).is_file() for p in artifacts),'Job missing output artifacts')
            record.update(status='PASS',result=result,artifacts={str(Path(p).resolve()):sha(p) for p in artifacts},completed_utc=time.time())
            with self.lock: self.counts['executed']+=1
        except Exception as e:
            record.update(status='FAILED',error=str(e),completed_utc=time.time())
            with self.lock: self.counts['failed']+=1
            publish(); raise
        publish(); return result

def run_command(work,label,argv,schrodinger=False):
    work=Path(work); env=os.environ.copy(); env['PYTHONDONTWRITEBYTECODE']='1'
    if schrodinger:
        tmp=work/'tmp'; tmp.mkdir(exist_ok=True)
        hosts=work/'schrodinger.hosts'
        hosts.write_text('name: localhost\nhost: localhost\nprocessors: 1\ntmpdir: '+str(tmp)+'\nschrodinger: '+str(SCH)+'\n')
        env.update(SCHRODINGER=str(SCH),SCHRODINGER_HOSTS=str(hosts),SCHRODINGER_TMPDIR=str(tmp))
    log=work/(label+'.console.log'); start=time.time()
    with log.open('w',encoding='utf-8') as f:
        code=subprocess.run([str(x) for x in argv],cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800).returncode
    record=dict(argv=[str(x) for x in argv],cwd=str(work),exit_code=code,elapsed_seconds=time.time()-start,log=str(log))
    write_json(work/(label+'.command.json'),record)
    require(code==0,'Command failed: '+str(log))
    return record

def select_scope(acceptance):
    roster=load_roster(); verify_anchors(roster)
    return roster.head(3).copy() if acceptance else roster

def output_root(acceptance): return ACCEPTANCE if acceptance else OUTPUT

def runner_cli(description):
    import argparse
    ap=argparse.ArgumentParser(description=description)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--run',action='store_true',help='Full roster; execute only with explicit user authorization')
    mode.add_argument('--acceptance',action='store_true',help='Real docking for rerun ranks 1, 2, 3 only')
    ap.add_argument('--workers',type=int,default=3)
    return ap
