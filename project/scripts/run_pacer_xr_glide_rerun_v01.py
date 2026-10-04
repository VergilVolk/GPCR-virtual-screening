"""Hash-aware, resumable SP production runner for the exact rerun roster."""
from __future__ import annotations
import json,math
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from dock_m4_gamd_ensemble import enumerate_ph7_states
from pacer_xr_production_v01 import *

BRIDGE=P/'scripts/pacer_xr_schrodinger_v01.py'

def prepare_candidate(row,store,assets,codes):
    cid=row.candidate_id
    inputs=dict(candidate_id=cid,canonical_smiles=row.canonical_smiles,rank=int(row.rerun_final_rank),roster_sha256=SOURCE_SHA,protocol=assets['protocol'],codes=codes,ligprep_executable_sha256=assets['schrodinger_executables']['ligprep.exe'])
    if row.rerun_final_rank==1: inputs['validated_rank1_ligand_sha256']=sha(P/assets['rank1_ligand'])
    def operation(work):
        if row.rerun_final_rank==1:
            smoke=P/'results/pacer_xr_rerun_glide_v01/smoke_7trs'
            old=json.loads((smoke/'ligand/prepared_state_properties.json').read_text())
            states=json.loads((smoke/'ligand/state_manifest.json').read_text())
            variants=[dict(candidate_id=cid,canonical_smiles=row.canonical_smiles,rerun_final_rank=1,molscrub_state_id=v['original_title'],prepared_variant_id=v['title'],prepared_title=v['title'],prepared_smiles=v['smiles']) for v in old]
            result=dict(candidate_id=cid,canonical_smiles=row.canonical_smiles,rerun_final_rank=1,ligand=str(P/assets['rank1_ligand']),states=states,variants=variants,provenance='Byte-exact docking-ready file from successful pH7 Molscrub/LigPrep smoke; no ligand regeneration')
            write_json(work/'variants.json',variants); write_json(work/'ligand_manifest.json',result)
            return result,[P/assets['rank1_ligand'],work/'variants.json',work/'ligand_manifest.json']
        states=enumerate_ph7_states(row.canonical_smiles,max_states=16)
        require(0<len(states)<=16,'Invalid Molscrub state coverage')
        records=[]
        with Chem.SDWriter(str(work/'ph7_3d.sdf')) as writer:
            for sid,smi in enumerate(states):
                state=cid+f'__s{sid:02d}'; mol=Chem.AddHs(Chem.MolFromSmiles(smi)); params=AllChem.ETKDGv3(); params.randomSeed=42+sid
                require(AllChem.EmbedMolecule(mol,params)==0,'RDKit embedding failed: '+state)
                status=AllChem.MMFFOptimizeMolecule(mol,maxIters=300)
                for key,value in dict(_Name=state,candidate_id=cid,parent_smiles=row.canonical_smiles,ligand_state_id=state,state_smiles=smi).items(): mol.SetProp(key,value)
                writer.write(mol)
                records.append(dict(molscrub_state_id=state,state_smiles=smi,embedding_seed=42+sid,mmff_status=status))
        meta=dict(candidate_id=cid,canonical_smiles=row.canonical_smiles,rerun_final_rank=int(row.rerun_final_rank),states={s['molscrub_state_id']:s for s in records})
        write_json(work/'bridge_metadata.json',meta)
        run_command(work,'ligprep',[SCH/'ligprep.exe','-isd','ph7_3d.sdf','-omae','prepared.maegz','-i','0','-nt','-ns','-bff','14','-WAIT','-LOCAL','-HOST','localhost'],True)
        run_command(work,'stamp',[SCH/'run.exe',BRIDGE,'stamp',work/'bridge_metadata.json'],True)
        variants=json.loads((work/'variants.json').read_text())
        result=dict(candidate_id=cid,canonical_smiles=row.canonical_smiles,rerun_final_rank=int(row.rerun_final_rank),ligand=str(work/'docking_ready.maegz'),states=records,variants=variants,provenance='Molscrub 0.1.1 pH7 cap16; RDKit ETKDGv3/MMFF300; LigPrep -i0 -nt -ns -bff14; all returned variants')
        write_json(work/'ligand_manifest.json',result)
        return result,[work/'docking_ready.maegz',work/'ph7_3d.sdf',work/'variants.json',work/'ligand_manifest.json']
    return store.execute('prepare_'+cid,inputs,operation)

def dock_one(row,receptor,prepared,store,assets,codes):
    grid=assets['receptors'][receptor]; gridpath=P/grid['grid']; ligand=Path(prepared['ligand'])
    inputs=dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,receptor_id=receptor,grid_sha256=sha(gridpath),ligand_sha256=sha(ligand),center=grid['center'],protocol=assets['protocol'],codes=codes,glide_executable_sha256=assets['schrodinger_executables']['glide.exe'])
    def operation(work):
        (work/'dock.in').write_text('GRIDFILE '+gridpath.as_posix()+'\nLIGANDFILE '+ligand.as_posix()+'\nPRECISION SP\nPOSES_PER_LIG 5\nPOSTDOCK_NPOSE 5\nFORCEFIELD OPLS_2005\n')
        command=run_command(work,'glide',[SCH/'glide.exe','dock.in','-WAIT','-LOCAL','-HOST','localhost','-NJOBS','1'],True)
        pv=work/'dock_pv.maegz'
        require(pv.exists(),'Glide completed without poseviewer output')
        require('Glide version 98128' in (work/'dock.log').read_text(errors='replace'),'Glide backend version mismatch')
        meta=dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,receptor_id=receptor,grid_sha256=inputs['grid_sha256'],center=grid['center'],variants=prepared['variants'],poseviewer=str(pv))
        write_json(work/'bridge_metadata.json',meta)
        run_command(work,'parse',[SCH/'run.exe',BRIDGE,'parse',work/'bridge_metadata.json'],True)
        poses=json.loads((work/'poses.json').read_text()); best=reduce_glide(poses)
        result=dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,receptor_id=receptor,cluster_id=None if receptor=='7TRS' else int(receptor[-2:]),glide_gscore=best['glide_gscore'],selected_variant_id=best['prepared_variant_id'],selected_entry_index=best['entry_index'],n_poses=len(poses),n_returned_variants=len({p['prepared_variant_id'] for p in poses}),accepted_variants=len(prepared['variants']),selected_pose=best,poses=poses,command=command,status='PASS')
        return result,[pv,work/'poses.json',work/'dock.in',work/'bridge_metadata.json']
    return store.execute(row.candidate_id+'_'+receptor,inputs,operation)

def run_glide_legacy(acceptance=False,workers=3):
    roster=select_scope(acceptance); assets=load_asset_lock(); verify_environment(assets,'glide'); codes=code_hashes()
    root=output_root(acceptance)/'glide'; store=JobStore(root/'jobs')
    prepared={}
    # Enumeration is sequential: preserve the validated Molscrub state order.
    for row in roster.itertuples():
        prepared[row.candidate_id]=prepare_candidate(row,store,assets,codes)
        print('Glide preparation PASS '+row.candidate_id,flush=True)
    write_json(root/'ligand_state_manifest.json',dict(source_sha256=SOURCE_SHA,scope='acceptance_rank1_3' if acceptance else 'full200',candidates=list(prepared.values())))
    results=[]; errors=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(dock_one,row,receptor,prepared[row.candidate_id],store,assets,codes):(row.candidate_id,receptor) for row in roster.itertuples() for receptor in assets['receptors']}
        for future in as_completed(futures):
            cid,receptor=futures[future]
            try: results.append(future.result()); print('Glide PASS '+cid+' '+receptor,flush=True)
            except Exception as exc: errors.append(dict(candidate_id=cid,receptor_id=receptor,error=str(exc))); print('Glide FAILED '+cid+' '+receptor+': '+str(exc),flush=True)
    write_json(root/'execution_summary.json',dict(scope='acceptance_rank1_3' if acceptance else 'full200',n_candidates=len(roster),required_evaluations=len(roster)*11,completed_evaluations=len(results),errors=errors,resume_counts=store.counts,source_sha256=SOURCE_SHA))
    require(not errors and len(results)==11*len(roster),'Glide incomplete; inspect execution_summary.json; no channel output')
    channels=[]; cluster_rows=[]; poses=[]
    for row in roster.itertuples():
        part=[x for x in results if x['candidate_id']==row.candidate_id]
        pdb=[x for x in part if x['receptor_id']=='7TRS']; ensemble=[x for x in part if x['cluster_id'] is not None]
        require(len(pdb)==1,'Incomplete Glide PDB coverage'); be=ensemble_score(ensemble,'glide_gscore')
        channels.append(dict(candidate_id=row.candidate_id,canonical_smiles=row.canonical_smiles,rerun_final_rank=int(row.rerun_final_rank),Glide_PDB_raw=pdb[0]['glide_gscore'],Glide_BEmin_raw=be['BEmin'],Glide_BEavg_raw=be['BEavg'],Glide_BEmin_cluster=be['winning_cluster'],n_clusters=10,status='PASS'))
        for result in sorted(ensemble,key=lambda x:x['cluster_id']):
            cluster_rows.append({k:result[k] for k in ['candidate_id','canonical_smiles','cluster_id','glide_gscore','selected_variant_id','selected_entry_index']}|dict(pmf=load_pmf()[result['cluster_id']],BE_i=be['BE_i'][result['cluster_id']]))
        for result in part: poses.extend(result['poses'])
    frame=pd.DataFrame(channels).sort_values('rerun_final_rank'); validate_scores(frame,roster,CHANNELS[:3])
    write_csv(root/'glide_all_poses.csv',pd.DataFrame(poses).sort_values(['candidate_id','receptor_id','pose_rank']))
    write_csv(root/'glide_cluster_scores.csv',pd.DataFrame(cluster_rows)); write_csv(root/'glide_channels.csv',frame)
    write_json(root/'channels_provenance.json',dict(source_sha256=SOURCE_SHA,codes=codes,channel_sha256=sha(root/'glide_channels.csv'),complete=True,n_candidates=len(roster),evaluations=len(results),assets_lock_sha256=sha(P/'config/pacer_xr_assets_v01.json')))
    return frame

def run_glide(acceptance=False,workers=3):
    # Scheduling changes do not invalidate frozen preparation or completed scientific scores.
    from pacer_xr_glide_scheduler_v01 import run_optimized
    return run_optimized(acceptance,workers)

def main():
    args=runner_cli(__doc__).parse_args(); require(args.workers>=1,'workers must be positive')
    run_glide(args.acceptance,args.workers)

if __name__=='__main__': main()
