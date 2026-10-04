"""Dedicated rerun Vina runner; legacy results are never reused without provenance."""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from pathlib import Path
import pandas as pd
from dock_m4_gamd_ensemble import enumerate_ph7_states,prepare_ligand
from pacer_xr_production_v01 import *

def prepare(row,ensemble,store,assets,codes):
    cid=row.candidate_id; protocol=assets['protocol']['vina_ensemble' if ensemble else 'vina_pdb']
    inputs=dict(candidate_id=cid,canonical_smiles=row.canonical_smiles,ensemble=ensemble,protocol=protocol,versions=assets['preparation_versions'],codes=codes,source_sha256=SOURCE_SHA)
    def operation(work):
        smiles=enumerate_ph7_states(row.canonical_smiles,max_states=16) if ensemble else [row.canonical_smiles]
        require(0<len(smiles)<=16,'No accepted Vina ligand state')
        states=[]
        for sid,smi in enumerate(smiles):
            lp=work/f'state_{sid:02d}.pdbqt'; prepare_ligand(smi,lp,seed=42+sid)
            states.append(dict(candidate_id=cid,canonical_smiles=row.canonical_smiles,state_id=cid+f'__s{sid:02d}',state_smiles=smi,embedding_seed=42+sid,ligand=str(lp),ligand_sha256=sha(lp)))
        write_json(work/'state_manifest.json',states)
        return states,[work/'state_manifest.json']+[Path(s['ligand']) for s in states]
    return store.execute('prepare_'+cid+('_ensemble' if ensemble else '_pdb'),inputs,operation)

def dock(row,receptor,state,store,assets,codes):
    ensemble=receptor!='7TRS'; rec=assets['receptors'][receptor]
    protocol=assets['protocol']['vina_ensemble' if ensemble else 'vina_pdb']; exe=P/'tools'/('vina_1.2.7.exe' if ensemble else 'vina.exe')
    inputs=dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,state=state,receptor_id=receptor,receptor_sha256=sha(P/rec['vina_receptor']),center=rec['vina_center'],size=rec['vina_size'],protocol=protocol,executable_sha256=sha(exe),codes=codes,source_sha256=SOURCE_SHA)
    def operation(work):
        argv=[exe,'--receptor',P/rec['vina_receptor'],'--ligand',state['ligand']]
        for axis,c,size in zip('xyz',rec['vina_center'],rec['vina_size']): argv+=['--center_'+axis,str(c),'--size_'+axis,str(size)]
        argv+=['--exhaustiveness',str(protocol['exhaustiveness']),'--num_modes',str(protocol['num_modes']),'--cpu','1','--seed','42','--out',work/'pose.pdbqt']
        if ensemble: argv+=['--energy_range','3.0']
        command=run_command(work,'vina',argv)
        affinity=parse_vina(work/'pose.pdbqt')
        require('AutoDock Vina v1.2.7' in (work/'vina.console.log').read_text(errors='replace'),'Unexpected Vina backend version')
        result=dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,receptor_id=receptor,cluster_id=None if not ensemble else int(receptor[-2:]),state_id=state['state_id'],state_smiles=state['state_smiles'],vina_affinity=affinity,pose=str(work/'pose.pdbqt'),command=command,status='PASS')
        return result,[work/'pose.pdbqt',work/'vina.command.json']
    return store.execute(row.candidate_id+'_'+receptor+'_'+state['state_id'].split('__')[-1],inputs,operation)

def run_vina(acceptance=False,workers=8):
    roster=select_scope(acceptance); assets=load_asset_lock(); verify_environment(assets,'vina'); codes=code_hashes()
    root=output_root(acceptance)/'vina'; store=JobStore(root/'jobs'); states={}
    for row in roster.itertuples():
        for ensemble in [False,True]: states[(row.candidate_id,ensemble)]=prepare(row,ensemble,store,assets,codes)
        print('Vina preparation PASS '+row.candidate_id,flush=True)
    write_json(root/'ligand_state_manifest.json',[dict(candidate_id=cid,ensemble=ensemble,states=values) for (cid,ensemble),values in states.items()])
    results=[]; errors=[]; task_count=0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures={}
        for row in roster.itertuples():
            for receptor in assets['receptors']:
                for state in states[(row.candidate_id,receptor!='7TRS')]:
                    f=pool.submit(dock,row,receptor,state,store,assets,codes); futures[f]=(row.candidate_id,receptor,state['state_id']); task_count+=1
        for future in as_completed(futures):
            cid,rec,sid=futures[future]
            try: results.append(future.result()); print('Vina PASS '+cid+' '+rec+' '+sid,flush=True)
            except Exception as exc: errors.append(dict(candidate_id=cid,receptor_id=rec,state_id=sid,error=str(exc))); print('Vina FAILED '+cid+' '+rec+' '+sid+': '+str(exc),flush=True)
    write_json(root/'execution_summary.json',dict(scope='acceptance_rank1_3' if acceptance else 'full200',n_candidates=len(roster),required_state_jobs=task_count,completed_state_jobs=len(results),errors=errors,resume_counts=store.counts,source_sha256=SOURCE_SHA,historical_reuse=False,historical_reuse_reason='Protocol and preparation fingerprints not independently established; rerun all entries.'))
    require(not errors and len(results)==task_count,'Vina incomplete; no channel output')
    channels=[]; cluster_rows=[]
    for row in roster.itertuples():
        reduced=[]
        for receptor in assets['receptors']:
            part=[r for r in results if r['candidate_id']==row.candidate_id and r['receptor_id']==receptor]
            require(len(part)==len(states[(row.candidate_id,receptor!='7TRS')]),'Missing Vina state evaluation')
            reduced.append(reduce_vina(part))
        pdb=[x for x in reduced if x['receptor_id']=='7TRS'][0]; ensemble=[x for x in reduced if x['cluster_id'] is not None]; be=ensemble_score(ensemble,'vina_affinity')
        channels.append(dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,rerun_final_rank=int(row.rerun_final_rank),Vina_PDB_raw=pdb['vina_affinity'],Vina_BEmin_raw=be['BEmin'],Vina_BEavg_raw=be['BEavg'],Vina_BEmin_cluster=be['winning_cluster'],n_clusters=10,status='PASS'))
        for r in sorted(ensemble,key=lambda x:x['cluster_id']): cluster_rows.append({k:r[k] for k in ['candidate_id','canonical_smiles','cluster_id','state_id','vina_affinity']}|dict(pmf=load_pmf()[r['cluster_id']],BE_i=be['BE_i'][r['cluster_id']]))
    frame=pd.DataFrame(channels).sort_values('rerun_final_rank'); validate_scores(frame,roster,CHANNELS[3:])
    write_csv(root/'vina_statewise_scores.csv',pd.DataFrame(results).drop(columns=['command']).sort_values(['candidate_id','receptor_id','state_id']))
    write_csv(root/'vina_cluster_scores.csv',pd.DataFrame(cluster_rows)); write_csv(root/'vina_channels.csv',frame)
    write_json(root/'channels_provenance.json',dict(source_sha256=SOURCE_SHA,codes=codes,channel_sha256=sha(root/'vina_channels.csv'),complete=True,n_candidates=len(roster),evaluations=len(roster)*11,state_jobs=len(results),assets_lock_sha256=sha(P/'config/pacer_xr_assets_v01.json')))
    return frame

def main():
    args=runner_cli(__doc__).parse_args(); require(args.workers>=1,'workers must be positive'); run_vina(args.acceptance,args.workers)

if __name__=='__main__': main()
