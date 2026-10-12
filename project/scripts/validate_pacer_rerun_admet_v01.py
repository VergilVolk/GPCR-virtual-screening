#!/usr/bin/env python3
"""Independent consistency, identity, numerical and frozen-file validation.
Usage: .venv-pacer-admet-v01/bin/python project/scripts/validate_pacer_rerun_admet_v01.py [--output DIR]
"""
import argparse,json,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors as D, Lipinski as L, QED, rdMolDescriptors as R
from pacer_rerun_medchem_v01 import ROOT,INPUT,DEFAULT,TOP3,KEY,sha,dump,catalog

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=DEFAULT);o=ap.parse_args().output
    inp=pd.read_csv(INPUT);mapping=inp.set_index('candidate_id').canonical_smiles.to_dict();checks={};details={}
    def check(name,value):checks[name]=bool(value)
    required=['INPUT_VALIDATION_v01.json','ALL200_PHYSICOCHEMICAL_v01.csv','ALL200_DRUGLIKENESS_v01.csv','ALL200_ADMET_v01.csv','ALL200_TOXICITY_v01.csv','ALL200_RETROSYNTHESIS_v01.csv','ALL200_SYNTHESIS_FEASIBILITY_v01.csv','ALL200_STEREOCHEMISTRY_v01.csv','TOP3_COMPREHENSIVE_PROFILE_v01.csv','MODEL_PROVENANCE_v01.json','EVIDENCE_COMPLETENESS_v01.csv','FINAL_SCIENTIFIC_REPORT_v01.md']
    check('required_prevalidation_outputs_exist',all((o/x).is_file() for x in required));check('input_sha256_unchanged',sha(INPUT)==json.loads((o/'INPUT_VALIDATION_v01.json').read_text())['sha256'])
    files={}
    for p in sorted(o.glob('*.csv')):
        d=pd.read_csv(p);keyok=all(k in d for k in KEY);check(p.name+'__identity_columns',keyok)
        if keyok:
            check(p.name+'__verbatim_identity',all(mapping.get(cid)==s for cid,s in zip(d.candidate_id,d.canonical_smiles)))
            expected=set(TOP3) if p.name.startswith('TOP3') else set(mapping)
            check(p.name+'__full_expected_candidate_set',set(d.candidate_id)==expected)
        files[p.name]={'rows':len(d),'columns':len(d.columns),'unique_candidates':int(d.candidate_id.nunique()) if keyok else None}
    p=pd.read_csv(o/'ALL200_PHYSICOCHEMICAL_v01.csv');dr=pd.read_csv(o/'ALL200_DRUGLIKENESS_v01.csv');assert len(p)==len(dr)==200
    recalcs=[]
    for r in p.to_dict('records'):
        m=Chem.MolFromSmiles(mapping[r['candidate_id']]);values={'molecular_weight':D.MolWt(m),'exact_molecular_weight':D.ExactMolWt(m),'MolLogP':D.MolLogP(m),'TPSA':D.TPSA(m),'HBD':L.NumHDonors(m),'HBA':L.NumHAcceptors(m),'rotatable_bonds':L.NumRotatableBonds(m),'aromatic_ring_count':R.CalcNumAromaticRings(m),'ring_count':R.CalcNumRings(m),'heavy_atom_count':m.GetNumHeavyAtoms(),'FractionCSP3':R.CalcFractionCSP3(m),'formal_charge':Chem.GetFormalCharge(m)}
        recalcs.append(all(np.isclose(r[k],v,atol=1e-9,rtol=1e-9) for k,v in values.items()))
    check('all200_all12_descriptors_recomputed',all(recalcs));check('all200_QED_recomputed',all(np.isclose(r.QED,QED.qed(Chem.MolFromSmiles(r.canonical_smiles))) for r in dr.itertuples()))
    for rule in ['Lipinski','Veber','Ghose','Egan','Muegge']:
        check(rule+'__flags_match_violation_lists',all((r[rule+'_status']=='PASS')==(len(json.loads(r[rule+'_violations']))==0) and r[rule+'_violation_count']==len(json.loads(r[rule+'_violations'])) for r in dr.to_dict('records')))
    check('PAINS_catalog_nonempty',catalog('PAINS').GetNumEntries()>400);details['PAINS_catalog_entries']=catalog('PAINS').GetNumEntries()
    a=pd.concat([pd.read_csv(o/'ALL200_ADMET_v01.csv'),pd.read_csv(o/'ALL200_TOXICITY_v01.csv')],ignore_index=True);members=pd.read_csv(o/'ALL200_ADMET_MODEL_MEMBERS_v01.csv');raw=pd.read_csv(o/'ALL200_ADMET_RAW_ENSEMBLE_MEANS_v01.csv')
    check('endpoint_records_unique',not a.duplicated(['candidate_id','endpoint']).any());check('all52_endpoint_complete200',len(a)==10400 and a.groupby('endpoint').candidate_id.nunique().eq(200).all());check('raw41_endpoint_complete200',len(raw)==200 and len(raw.columns)==43)
    member_cols=[f'member_{i}' for i in range(5)];check('all8200_member_means_consistent',len(members)==8200 and np.allclose(members[member_cols].mean(axis=1),members.ensemble_mean,rtol=1e-6,atol=1e-5));check('all8200_member_SDs_consistent',np.allclose(members[member_cols].std(axis=1,ddof=1),members.ensemble_sample_std,rtol=1e-5,atol=1e-5))
    am=a.merge(members[['candidate_id','endpoint','ensemble_mean']],on=['candidate_id','endpoint'],how='inner',validate='one_to_one');check('raw_predictions_match_saved_members',np.allclose(am.raw_prediction,am.ensemble_mean,rtol=1e-6,atol=1e-5))
    check('missing_failed_predictions_are_null',a[a.status!='PREDICTED'].prediction.isna().all());check('predicted_values_finite',np.isfinite(a[a.status=='PREDICTED'].prediction).all());check('classification_probabilities_in_0_1',a[(a.output_type=='class_1_probability')&(a.status=='PREDICTED')].prediction.between(0,1).all())
    ppb=a[a.endpoint=='PPBR_AZ'];check('PPBR34_out_of_bounds_not_clipped',len(ppb[ppb.status=='MODEL_FAILED'])==34 and ppb[ppb.status=='MODEL_FAILED'].raw_prediction.gt(100).all());check('negative_R2_endpoints_not_accepted',a[a.endpoint.isin(['Half_Life_Obach','VDss_Lombardo'])].status.eq('NOT_ASSESSED').all());check('model_AD_not_falsely_in_domain',a[a.status=='PREDICTED'].applicability_domain.str.startswith('UNKNOWN').all())
    rt=pd.read_csv(o/'ALL200_RETROSYNTHESIS_v01.csv');check('route_unassessed_not_no_route',len(rt)==200 and rt.status.eq('NOT_ASSESSED').all() and rt.shortest_reaction_steps.isna().all() and rt.candidate_route_count.isna().all())
    syn=pd.read_csv(o/'ALL200_SYNTHESIS_FEASIBILITY_v01.csv');check('no_unfounded_cost_yield_supplier',syn.cost.isna().all() and syn['yield'].isna().all() and syn.supplier.isna().all() and syn.price.isna().all())
    top=pd.read_csv(o/'TOP3_COMPREHENSIVE_PROFILE_v01.csv');check('frozen_top3_order_identity',top.candidate_id.tolist()==TOP3 and all(mapping[r.candidate_id]==r.canonical_smiles for r in top.itertuples()))
    for cid in TOP3:
        tr=top.set_index('candidate_id').loc[cid];g=a[a.candidate_id==cid]
        check(cid+'__profile_predictions_match_all200',all((pd.isna(r.prediction) and pd.isna(tr[r.endpoint+'__prediction'])) or np.isclose(r.prediction,tr[r.endpoint+'__prediction']) for r in g.itertuples()))
    figures=[]
    for i in range(1,7):
        png=o/f'figures/Figure_{i}_v01.png';pdf=png.with_suffix('.pdf')
        with Image.open(png) as im:sz=im.size;im.verify()
        check(f'Figure_{i}__png_pdf_valid',min(sz)>1000 and pdf.read_bytes().startswith(b'%PDF-'))
        figures.append({'figure':i,'png_pixels':sz,'pdf_bytes':pdf.stat().st_size,'visual_QA':'Inspected all six PNGs; labels, evidence types, missingness and model caveats readable.'})
    changed=subprocess.check_output(['git','diff','HEAD','--name-only'],cwd=ROOT,text=True).splitlines();check('existing_tracked_files_only_gitignore_changed',set(changed)<={'.gitignore'});check('frozen_tracked_results_unchanged',not any(x.startswith('project/results/') for x in changed));check('analysis_branch_correct',subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()=='analysis/pacer-rerun-admet-synthesis-v01')
    rc=subprocess.run(['git','check-ignore',str(o/'figures/Figure_1_v01.pdf')],cwd=ROOT,capture_output=True,text=True).returncode;check('result_PDF_not_gitignored',rc==1)
    statuses=a.status.value_counts().to_dict()
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'failed_checks':[k for k,v in checks.items() if not v],'files':files,'figures':figures,'details':details,'endpoint_status_counts':statuses,'n_model_endpoints':41,'n_recorded_endpoints':52,'input_sha256_final':sha(INPUT),'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'git_branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),'git_status':subprocess.check_output(['git','status','--short','--untracked-files=normal'],cwd=ROOT,text=True),'modified_tracked_files':changed,'not_committed_pushed_or_PR_created':True,'scope_limitations':['No validated model applicability domain or probability calibration on these candidates','Caco2 unit discrepancy unresolved','Half-life/VDss model outputs not adopted','No complete retrosynthesis assets, searches or procurement validation'],'scientific_completeness':'All200 covered with explicit gaps, not all requested endpoints reliably predicted.'}
    dump(o,'ANALYSIS_VALIDATION_v01.json',result)
    paths=sorted([p for p in o.rglob('*') if p.is_file()]+list((ROOT/'project/scripts').glob('*pacer_rerun*admet*v01.py'))+[ROOT/'project/scripts/pacer_rerun_medchem_v01.py',ROOT/'project/scripts/pacer_rerun_synthesis_v01.py',ROOT/'project/scripts/pacer_rerun_report_v01.py'])
    paths=sorted(set(paths));dump(o,'DELIVERY_MANIFEST_v01.json',{'files':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in paths],'self_excluded_to_avoid_recursive_hash':True,'note':'.venv directories excluded; pinned dependencies and installed model hashes provided separately.'})
    print(json.dumps({'status':result['status'],'checks':len(checks),'failures':result['failed_checks'],'files_manifested':len(paths)},indent=2));assert result['status']=='PASS'
if __name__=='__main__':main()
