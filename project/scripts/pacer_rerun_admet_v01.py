#!/usr/bin/env python3
"""Phase C: local ADMET-AI 2.0.1 inference, member outputs, explicit missing endpoints.
Usage: .venv-pacer-admet-v01/bin/python project/scripts/pacer_rerun_admet_v01.py [--output DIR]
"""
import argparse, json, math, importlib.metadata as md
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from lightning import pytorch as pl
from rdkit import Chem
from admet_ai import ADMETModel
from admet_ai.constants import DEFAULT_MODELS_DIR, DEFAULT_ADMET_PATH
from pacer_rerun_medchem_v01 import INPUT, DEFAULT, KEY, sha, csv, dump, identity, encode

class AuditedADMET(ADMETModel):
    def _make_ensemble_predictions(self,data_loader):
        means={};self.member_predictions={};self.stds={}
        with torch.inference_mode():
            for tasks,models in zip(self.task_lists,self.model_lists):
                trainer=pl.Trainer(logger=False,enable_checkpointing=False,enable_progress_bar=False,accelerator='cpu',devices=1)
                a=torch.stack([torch.cat(trainer.predict(model=m,dataloaders=data_loader),dim=0) for m in models]).numpy()
                for i,t in enumerate(tasks):
                    self.member_predictions[t]=a[:,:,i].T; means[t]=a[:,:,i].mean(axis=0);self.stds[t]=a[:,:,i].std(axis=0,ddof=1)
        return means

MISSING={
    'Pgp_substrate':'ADMET-AI Pgp_Broccatelli is inhibition only; no verified local substrate model.',
    'CYP_metabolism_site':'No validated atom-level site model/assets; inhibition and substrate probabilities do not locate metabolism sites.',
    'Human_microsomal_half_life':'Clearance_Microsome_AZ is in-vitro intrinsic clearance, not microsomal half-life; no conversion performed.',
    'Fraction_unbound':'No independent fraction-unbound model; PPBR retained as predicted binding %, no additional derived endpoint.',
    'Systemic_clearance':'Hepatocyte/microsome intrinsic clearance is not in-vivo systemic clearance; no physiological scaling model.',
    'Oral_bioavailability_fraction':'Bioavailability_Ma is binary classification; no calibrated continuous F model.',
    'CNS_brain_exposure':'BBB classification is not measured brain exposure or Kp,uu; no exposure model.',
    'CNS_MPO_full':'No verified pKa/logD(pH7.4) inputs; full six-component CNS MPO not calculated.',
    'Hepatotoxicity_general':'DILI model is a specific dataset endpoint; no separate generic hepatotoxicity model.',
    'Reactive_metabolite_model':'No validated bioactivation model/assets; separate motif structural alerts are available.',
    'Human_microsomal_stability_class':'No validated microsomal stability classification model; retain direct microsomal clearance prediction separately.'}
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=DEFAULT);o=ap.parse_args().output
    df=pd.read_csv(INPUT);meta=pd.read_csv(DEFAULT_ADMET_PATH).set_index('id');torch.set_num_threads(4);torch.manual_seed(20261012);np.random.seed(20261012)
    model=AuditedADMET(include_physchem=False,drugbank_path=None,num_workers=0)
    pred=model.predict(df.canonical_smiles.tolist());assert list(pred.index)==df.canonical_smiles.tolist() and len(pred)==200
    csv(o,'ALL200_ADMET_RAW_ENSEMBLE_MEANS_v01.csv',pd.concat([df[KEY].reset_index(drop=True),pred.reset_index(drop=True)],axis=1))
    members=[];admet=[];tox=[]
    for ix,r in enumerate(df.to_dict('records')):
        k=identity(r)
        for ep in pred.columns:
            x=meta.loc[ep];raw=float(pred.iloc[ix][ep]);std=float(model.stds[ep][ix]);classification=x.task_type=='classification'
            vals=[float(v) for v in model.member_predictions[ep][ix]]
            members.append({**k,'endpoint':ep,**{f'member_{j}':v for j,v in enumerate(vals)},'ensemble_mean':raw,'ensemble_sample_std':std})
            finite=math.isfinite(raw);physical=finite and ((0<=raw<=1) if classification else float(x.minimum)<=raw<=float(x.maximum))
            # A negative external R2 in shipped metadata does not support an actionable regression.
            reject=not classification and pd.notna(x['R^2']) and x['R^2']<0
            reason='Negative R2 in shipped model metadata; raw output archived, no reliable endpoint assessment.' if reject else ('Output violates endpoint physical bounds; retained raw, not clipped.' if not physical else '')
            unit=x.units
            unote=''
            if ep=='Caco2_Wang': unote='Package declares log(10^-6 cm/s), TDC leaderboard declares cm/s. Raw model-scale output preserved; scale discrepancy unresolved, no conversion or permeability threshold used.'
            rec={**k,'endpoint':ep,'endpoint_name':x['name'],'category':x.category,'model_name':'ADMET-AI Chemprop ensemble','model_version':md.version('admet-ai'),'source':x.url,'output_type':'class_1_probability' if classification else 'regression','prediction':raw if physical and not reject else None,'raw_prediction':raw,'unit':'probability' if classification else unit,'unit_note':unote,'applicability_domain':'UNKNOWN: checkpoints do not include endpoint training structures or validated per-molecule AD estimator','uncertainty':'Five-member sample SD; not a calibrated confidence interval; training/test error does not establish error on PACERGEN scaffolds','ensemble_sample_std':std,'ensemble_size':len(vals),'status':'NOT_ASSESSED' if reject else ('PREDICTED' if physical else 'MODEL_FAILED'),'status_reason':reason,'evidence_grade':'D_UNRELIABLE_MODEL' if reject else ('C_EXPLORATORY_MODEL_UNVERIFIED_AD' if physical else 'D_INVALID_OUTPUT'),'classification_threshold':None,'threshold_source':'No threshold applied; class-1 probability retained without binary safety label' if classification else None,'reported_AUROC':x.AUROC if pd.notna(x.AUROC) else None,'reported_AUPRC':x.AUPRC if pd.notna(x.AUPRC) else None,'reported_R2':x['R^2'] if pd.notna(x['R^2']) else None,'reported_MAE':x.MAE if pd.notna(x.MAE) else None,'performance_provenance':'Shipped admet.csv; version-specific external replication not established; not PACERGEN validation','training_size_reported':x['size'],'species_reported':x.species,'assay_conditions':'Dataset-specific; no uniform pH, crystal form, dose or assay conditions implied','limitations':'No experimental validation on these 200; unspecified stereo and protonation retained; probabilities may be miscalibrated.'}
            (tox if x.category=='Toxicity' else admet).append(rec)
        for ep,reason in MISSING.items():
            rec={**k,'endpoint':ep,'model_name':None,'model_version':None,'source':None,'output_type':None,'prediction':None,'raw_prediction':None,'unit':None,'applicability_domain':'NOT_ASSESSED','uncertainty':'No reliable estimate','status':'NOT_ASSESSED','status_reason':reason,'evidence_grade':'D_MISSING'}
            (tox if ep in ['Hepatotoxicity_general','Reactive_metabolite_model'] else admet).append(rec)
    csv(o,'ALL200_ADMET_MODEL_MEMBERS_v01.csv',members);csv(o,'ALL200_ADMET_v01.csv',admet);csv(o,'ALL200_TOXICITY_v01.csv',tox)
    # Deliberately report only known CNS MPO components, never a counterfeit six-point score.
    phy=pd.read_csv(o/'ALL200_PHYSICOCHEMICAL_v01.csv');cn=[];reac=[]
    dec=lambda v,a,b:float(np.clip((b-v)/(b-a),0,1))
    reactive={'thiophene_bioactivation_motif':'[s;r5]1[c;r5][c;r5][c;r5][c;r5]1','furan_bioactivation_motif':'[o;r5]1[c;r5][c;r5][c;r5][c;r5]1','nitroaromatic_motif':'[a][N+](=O)[O-]','phenol_motif':'[c][OX2H]','epoxide_reactive_group':'[O;r3]1[C;r3][C;r3]1','aldehyde_reactive_group':'[CX3H1](=O)[#6]','michael_acceptor_group':'[C;!a]=[C;!a]-[C,S](=O)'}
    for r in phy.to_dict('records'):
        k=identity(r);p=r['TPSA'];scorepsa=dec(p,90,120) if p>90 else (float(np.clip((p-20)/20,0,1)) if p<40 else 1.)
        scores={'CNS_MPO_MW_component':dec(r['molecular_weight'],360,500),'CNS_MPO_MolLogP_component_proxy':dec(r['MolLogP'],3,5),'CNS_MPO_TPSA_component':scorepsa,'CNS_MPO_HBD_component':dec(r['HBD'],.5,3.5)}
        cn.append({**k,**scores,'partial_component_sum_0_to_4':sum(scores.values()),'logD_pH7_4':None,'most_basic_pKa':None,'full_CNS_MPO':None,'status':'PARTIAL_HEURISTIC_ONLY','interpretation':'Wager 2010 functions with RDKit descriptors, MolLogP proxy for ClogP; not BBB or brain-exposure evidence.'})
        m=Chem.MolFromSmiles(r['canonical_smiles']);hits=[]
        for name,q in reactive.items():
            matches=m.GetSubstructMatches(Chem.MolFromSmarts(q))
            if matches:hits.append({'rule_id':name,'smarts':q,'atom_indices_0based':[list(t) for t in matches],'matched_fragments':[Chem.MolFragmentToSmiles(m,list(t)) for t in matches]})
        reac.append({**k,'endpoint':'reactive_or_bioactivation_motif_screen','prediction':bool(hits),'status':'STRUCTURAL_ALERT' if hits else 'SCREENED_NO_MATCH','model_name':'Explicit motif screen v01 (non-exhaustive, not trained)','model_version':'v01','source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5871347/','output_type':'structural_motif_matches','unit':None,'applicability_domain':'Structure matching only; no formation probability or metabolic pathway assignment','uncertainty':'Substitution, fusion, alternative metabolism and detoxification affect relevance; no experimental formation established','matches':encode(hits),'evidence_grade':'B_STRUCTURAL_ALERT','interpretation':'These flags nominate confirmatory studies, not demonstrated reactive-metabolite formation or toxicity.'})
    csv(o,'ALL200_CNS_COMPONENTS_v01.csv',cn);csv(o,'ALL200_REACTIVITY_ALERTS_v01.csv',reac)
    assets=[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(DEFAULT_MODELS_DIR.rglob('*.pt'))]
    meta_clean=json.loads(meta.reset_index().to_json(orient='records'))
    dump(o,'MODEL_PROVENANCE_v01.json',{'inference_package':'admet-ai','version':md.version('admet-ai'),'chemprop':md.version('chemprop'),'torch':md.version('torch'),'source':'https://github.com/swansonk14/admet_ai','license':'MIT (package); underlying dataset licenses vary, see TDC','assets':assets,'metadata_sha256':sha(DEFAULT_ADMET_PATH),'endpoint_metadata':meta_clean,'domain_policy':'UNKNOWN for all molecules; no false in-domain declarations or arbitrary similarity thresholds','uncertainty_method':'Sample standard deviation ddof=1 over 5 trained members; no calibrated coverage claim','classification_policy':'Raw P(class1), no binary threshold or combined development probability','regression_policy':'Raw units, no clipping, physical bound violations MODEL_FAILED, negative metadata R2 NOT_ASSESSED and raw output archived','v2_warning':'Published 2024 paper describes v1; v2 retrained models differ. Package metadata metrics not independently established as v2 validation.','Caco2_unit_warning':'Package unit and TDC leaderboard conflict. Raw output only, no physical conversion.','RDKit_version':md.version('rdkit'),'SA_score_source':'https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score','SA_assets':[{'path':str(p),'sha256':sha(p)} for p in sorted(Path(__import__('rdkit').__file__).parent.joinpath('Contrib/SA_Score').glob('*')) if p.is_file() and p.suffix in ['.py','.gz']],'rules_source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5335600/','rules_convention':'All PASS means zero violations. Lipinski <=1 additionally reported; Ghose atoms includes hydrogens; Muegge MolLogP proxy; Egan common rectangular approximation. QED default weighted mean.','CNS_source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC3368654/','reactivity_source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5871347/','missing_endpoints':MISSING})
    dump(o,'MODEL_ENDPOINT_METADATA_v01.json',meta_clean)
    dump(o,'PHASE_C_COMPLETE_v01.json',{'status':'COMPLETE_WITH_EXPLICIT_GAPS','n_molecules':200,'n_model_endpoints':len(pred.columns),'n_raw_predictions':pred.size,'input_sha256':sha(INPUT),'all_model_members_saved':True})
    print('Phase C complete:',pred.shape)
if __name__=='__main__':main()
