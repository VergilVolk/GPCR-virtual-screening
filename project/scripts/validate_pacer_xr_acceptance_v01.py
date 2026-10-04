"""Validate only rank1–3 acceptance; optional two-job Vina repeatability check."""
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd
from dock_m4_gamd_ensemble import prepare_ligand
from pacer_xr_production_v01 import *
from assemble_pacer_xr_six_channels_v01 import verify_channel_provenance

def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--repeat-vina',action='store_true'); args=ap.parse_args()
    roster=select_scope(True); assets=load_asset_lock(); refs,ref_lock=load_references()
    glide=verify_channel_provenance(ACCEPTANCE,'glide',roster); vina=verify_channel_provenance(ACCEPTANCE,'vina',roster)
    raw=pd.read_csv(ACCEPTANCE/'candidate_six_channel_raw_scores.csv'); validate_scores(raw,roster)
    rank1=glide[glide.rerun_final_rank==1].iloc[0]
    expected=dict(Glide_PDB_raw=-6.29971543794357,Glide_BEmin_raw=-7.18155879467582,Glide_BEavg_raw=-5.025427724052952)
    reproduction={name:dict(expected=value,observed=float(rank1[name]),absolute_difference=abs(float(rank1[name])-value)) for name,value in expected.items()}
    require(all(x['absolute_difference']<=1e-6 for x in reproduction.values()),'Rank1 Glide smoke reproduction failed')
    repeats=[]; preparation_checks=[]
    if args.repeat_vina:
        from dock_candidate_portfolio import ligand as legacy_pdb_prepare
        work=ACCEPTANCE/'determinism_rechecks'; work.mkdir(exist_ok=True)
        manifest=json.loads((ACCEPTANCE/'vina/ligand_state_manifest.json').read_text())
        for entry in manifest:
            for state in entry['states']:
                target=work/(state['state_id']+('_ensemble' if entry['ensemble'] else '_pdb')+'.pdbqt')
                if entry['ensemble']: prepare_ligand(state['state_smiles'],target,state['embedding_seed'])
                else: legacy_pdb_prepare(state['state_smiles'],target,state['embedding_seed'])
                require(sha(target)==state['ligand_sha256'],'Vina ligand preparation not byte deterministic')
                preparation_checks.append(dict(candidate_id=entry['candidate_id'],state_id=state['state_id'],ensemble=entry['ensemble'],sha256=sha(target),PASS=True))
        cid=roster.iloc[0].candidate_id
        for rec in ['7TRS','cluster_00']:
            status=json.loads((ACCEPTANCE/'vina/jobs'/(cid+'_'+rec+'_s00')/'status.json').read_text()); result=status['result']
            argv=result['command']['argv'].copy(); target=work/(rec+'_repeat.pdbqt'); argv[argv.index('--out')+1]=str(target)
            run_command(work,rec+'_repeat',argv)
            affinity=parse_vina(target); require(affinity==result['vina_affinity'],'Vina repeat affinity differs')
            repeats.append(dict(candidate_id=cid,receptor_id=rec,expected_from_first_real_run=result['vina_affinity'],repeat_affinity=affinity,PASS=True))
    # Rerun contract tests with stdout retained under acceptance, not production ranking.
    with (ACCEPTANCE/'contract_tests.log').open('w',encoding='utf-8') as log:
        result=subprocess.run([sys.executable,'-B','-m','unittest','discover','-s',str(P/'tests'),'-p','test_pacer_xr_production_v01.py','-v'],cwd=REPO,stdout=log,stderr=subprocess.STDOUT)
    require(result.returncode==0,'Contract test suite failed')
    statuses={engine:json.loads((ACCEPTANCE/engine/'execution_summary.json').read_text()) for engine in ['glide','vina']}
    for engine in ['glide','vina']:
        records=[json.loads(p.read_text()) for p in (ACCEPTANCE/engine/'jobs').glob('*/status.json') if p.parent.name.startswith('PACERGEN')]
        require({r['result']['candidate_id'] for r in records}==set(roster.candidate_id),'Acceptance exceeds ranks1–3')
    for name in ['pacer_xr_rerun_full200_ranking.csv','pacer_xr_natural_top20.csv','pacer_xr_final_shortlist20.csv']:
        require(not (OUTPUT/name).exists() and not (ACCEPTANCE/name).exists(),'Unauthorized scientific ranking output exists')
    # Frozen locks include source/prepared assets and both historical smoke manifests.
    for name in ['smoke_7trs/smoke_manifest.json','smoke_ensemble/smoke_ensemble_manifest.json']:
        root=P/'results/pacer_xr_rerun_glide_v01'/Path(name).parent; old=json.loads((root/'smoke_manifest.json' if 'smoke_7trs' in name else root/'smoke_ensemble_manifest.json').read_text())
        for relative,digest in old['artifact_sha256'].items(): require(sha(root/relative)==digest,'Historical smoke artifact changed: '+relative)
    report=dict(status='PASS',scope='rerun_ranks_1_2_3_only',source_sha256=SOURCE_SHA,glide_counts=dict(PDB=3,ensemble=30,BEmin=3,BEavg=3),vina_counts=dict(PDB=3,ensemble=30,BEmin=3,BEavg=3),rank1_glide_reproduction=reproduction,vina_channels=vina.to_dict('records'),glide_channels=glide.to_dict('records'),vina_preparation_determinism=preparation_checks,vina_repeat_jobs=repeats,reference_commit=ref_lock['reference_commit'],reference_channels=ref_lock['channels'],contract_tests='PASS; see contract_tests.log',historical_smoke_hashes_unchanged=True,full_campaign_started=False,full_ranking_produced=False,final_top20_produced=False,execution_summaries=statuses)
    write_json(ACCEPTANCE/'acceptance_report.json',report)
    print(json.dumps(dict(status='PASS',rank1_glide_reproduction=reproduction,vina_channels=report['vina_channels'],repeat_jobs=repeats),indent=2))

if __name__=='__main__': main()
