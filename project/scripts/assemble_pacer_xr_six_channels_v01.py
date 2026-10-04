"""Fail-closed assembly of complete protocol-matched rerun Glide/Vina channels."""
import json
import pandas as pd
from pacer_xr_production_v01 import *

def verify_channel_provenance(root,engine,roster):
    path=root/engine/(engine+'_channels.csv'); prov=json.loads((root/engine/'channels_provenance.json').read_text())
    require(prov['complete'] and prov['n_candidates']==len(roster) and prov['evaluations']==len(roster)*11,'Incomplete engine provenance')
    require(prov['source_sha256']==SOURCE_SHA and prov['channel_sha256']==sha(path),'Raw channel provenance mismatch')
    if engine=='glide' and prov.get('execution_model')=='receptor_major':
        from pacer_xr_glide_scheduler_v01 import science_contract,scheduler_contract,collect_complete
        scheduler_contract()
        require(prov.get('scientific_signature')==fingerprint(science_contract()),'Glide scientific contract changed')
        require(prov['assets_lock_sha256']==sha(P/'config/pacer_xr_assets_v01.json'),'Asset contract changed')
        _,_,results=collect_complete(root/engine,roster,load_asset_lock())
        require(len(results)==len(roster)*11,'Glide completion index missing evaluations')
        frame=pd.read_csv(path); validate_scores(frame,roster,CHANNELS[:3])
        for row in frame.itertuples():
            ensemble=[results[row.candidate_id+f'_cluster_{i:02d}'] for i in range(10)]
            be=ensemble_score(ensemble,'glide_gscore')
            require(abs(row.Glide_PDB_raw-results[row.candidate_id+'_7TRS']['glide_gscore'])<=1e-12 and abs(row.Glide_BEmin_raw-be['BEmin'])<=1e-12 and abs(row.Glide_BEavg_raw-be['BEavg'])<=1e-12,'Glide channel scores differ from validated completion records')
        return frame
    current=code_hashes()
    relevant={'pacer_xr_production_v01.py','run_pacer_xr_'+engine+'_rerun_v01.py','dock_m4_gamd_ensemble.py'}
    if engine=='glide':
        relevant.add('pacer_xr_schrodinger_v01.py')
        # Preserve scientifically valid legacy channels across execution-layer updates.
        from pacer_xr_glide_scheduler_v01 import science_contract
        science=science_contract()
        require(prov['codes'].get('run_pacer_xr_glide_rerun_v01.py') in [science['legacy_runner_sha256'],current['run_pacer_xr_glide_rerun_v01.py']],'Unapproved legacy Glide runner')
        relevant.remove('run_pacer_xr_glide_rerun_v01.py')
    require(all(prov['codes'].get(name)==current[name] for name in relevant),'Runner/parser code changed; rerun or resume the engine stage')
    require(prov['assets_lock_sha256']==sha(P/'config/pacer_xr_assets_v01.json'),'Asset contract changed')
    statuses=list((root/engine/'jobs').glob('*/status.json'))
    require(bool(statuses),'No completed docking job ledger')
    receptor_names=['7TRS']+[f'cluster_{i:02d}' for i in range(10)]
    expected_jobs=set()
    if engine=='glide':
        for cid in roster.candidate_id:
            expected_jobs.add('prepare_'+cid)
            expected_jobs.update(cid+'_'+rec for rec in receptor_names)
    else:
        manifest=json.loads((root/engine/'ligand_state_manifest.json').read_text())
        require(len(manifest)==2*len(roster),'Incomplete Vina preparation manifest')
        state_map={(x['candidate_id'],x['ensemble']):x['states'] for x in manifest}
        require(len(state_map)==2*len(roster),'Duplicate Vina preparation identity')
        for cid in roster.candidate_id:
            for ensemble in [False,True]:
                expected_jobs.add('prepare_'+cid+('_ensemble' if ensemble else '_pdb'))
                require(bool(state_map.get((cid,ensemble))),'Missing Vina prepared states')
            for rec in receptor_names:
                expected_jobs.update(cid+'_'+rec+'_'+s['state_id'].split('__')[-1] for s in state_map[(cid,rec!='7TRS')])
    require({path.parent.name for path in statuses}==expected_jobs,'Missing, unexpected, or stale docking job identities')
    for status in statuses:
        record=json.loads(status.read_text())
        require(JobStore.valid(None,record,record.get('fingerprint')),'Invalid resumable job artifact: '+str(status))
    frame=pd.read_csv(path); validate_scores(frame,roster,CHANNELS[:3] if engine=='glide' else CHANNELS[3:])
    return frame

def assemble_production(acceptance=False):
    roster=select_scope(acceptance); root=output_root(acceptance); load_asset_lock()
    glide=verify_channel_provenance(root,'glide',roster); vina=verify_channel_provenance(root,'vina',roster)
    frame=assemble(roster,glide,vina)
    target=root/'candidate_six_channel_raw_scores.csv'; write_csv(target,frame)
    write_json(root/'six_channel_provenance.json',dict(source_sha256=SOURCE_SHA,n_candidates=len(roster),complete=True,scope='acceptance_rank1_3' if acceptance else 'full200',raw_sha256=sha(target),glide_sha256=sha(root/'glide/glide_channels.csv'),vina_sha256=sha(root/'vina/vina_channels.csv'),codes=code_hashes()))
    return frame

def main():
    args=runner_cli(__doc__).parse_args(); assemble_production(args.acceptance)

if __name__=='__main__': main()
