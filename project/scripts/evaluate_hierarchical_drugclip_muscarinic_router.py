#!/usr/bin/env python3
"""DrugCLIP family router plus muscarinic ECFP3 subtype adapter."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from torch.nn import functional as F

from evaluate_muscarinic_ecfp_pocket_baseline import fps
from select_molecular_adapter_on_non_m4 import Adapter,train
from evaluate_gpcr_drugclip_target_loso import TARGETS,metrics,scaffold_bootstrap_delta


def scaffold(s):
    m=Chem.MolFromSmiles(s);return MurckoScaffold.MurckoScaffoldSmiles(mol=m) or Chem.MolToSmiles(m,canonical=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--muscarinic-representations',type=Path,required=True);p.add_argument('--muscarinic-triplets',type=Path,required=True);p.add_argument('--calibration-representations',type=Path,required=True);p.add_argument('--calibration-pairs',type=Path,required=True);p.add_argument('--retrieval-representations',type=Path,required=True);p.add_argument('--retrieval-benchmark',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seeds',default='20260924,20260925,20260926');p.add_argument('--epochs',type=int,default=100);p.add_argument('--lr',type=float,default=1e-3);p.add_argument('--weight-decay',type=float,default=1e-2);p.add_argument('--margin',type=float,default=.1);p.add_argument('--bootstrap',type=int,default=1000);a=p.parse_args();torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    ma=np.load(a.muscarinic_representations,allow_pickle=False);ca=np.load(a.calibration_representations,allow_pickle=False);ra=np.load(a.retrieval_representations,allow_pickle=False);table0=pd.read_csv(a.retrieval_benchmark);ridx={v:i for i,v in enumerate(map(str,ra['molecule_ids']))};table=table0[table0.canonical_smiles.astype(str).isin(ridx)].copy().reset_index(drop=True);order=np.asarray([ridx[v] for v in table.canonical_smiles.astype(str)]);official=ra['scores'].T[order].astype(np.float32);labels=np.asarray([TARGETS.index(v) for v in table.target])
    calibration=pd.read_csv(a.calibration_pairs);calibration=calibration[((calibration.positive_subtype=='M2')&(calibration.negative_subtype=='M4'))|((calibration.positive_subtype=='M4')&(calibration.negative_subtype=='M2'))].reset_index(drop=True);retrieval_scaffolds=set(table.murcko_scaffold.astype(str));calibration_scaffolds={scaffold(v) for v in calibration.canonical_smiles};mids=list(map(str,ma['molecule_ids']));mi={v:i for i,v in enumerate(mids)};mpids=[str(v).split('_')[0] for v in ma['pocket_ids']];mpi={v:i for i,v in enumerate(mpids)};trip0=pd.read_csv(a.muscarinic_triplets);trip=trip0[trip0.canonical_smiles.astype(str).isin(mi)].copy();trip['murcko_scaffold']=[scaffold(v) for v in trip.canonical_smiles];trip=trip[~trip.murcko_scaffold.isin(retrieval_scaffolds|calibration_scaffolds)].reset_index(drop=True)
    mx=fps(mids,radius=3);mpocket=F.normalize(torch.as_tensor(ma['pocket_embeddings'].astype(np.float32)),dim=-1);rx=fps(table.canonical_smiles.astype(str).tolist(),radius=3);rpocket=F.normalize(torch.as_tensor(ra['pocket_embeddings'].astype(np.float32)),dim=-1);rnames=[str(v).split('_')[0] for v in ra['pocket_ids']];rpi={v:i for i,v in enumerate(rnames)}
    cids=list(map(str,ca['molecule_ids']));cmi={v:i for i,v in enumerate(cids)};cx=fps(cids,radius=3);cpocket=F.normalize(torch.as_tensor(ca['pocket_embeddings'].astype(np.float32)),dim=-1);cnames=[str(v).split('_')[0] for v in ca['pocket_ids']];cpi={v:i for i,v in enumerate(cnames)};local=[];cal_scores=[]
    for seed in map(int,a.seeds.split(',')):
        model=train(mx,mpocket,trip,mi,mpi,seed,a)
        with torch.inference_mode():local.append((model(rx)@rpocket.T).numpy());cal_scores.append((model(cx)@cpocket.T).numpy())
    subtype=np.mean(local,axis=0);cal_matrix=np.mean(cal_scores,axis=0);musc=[rpi['M2R'],rpi['M4R']];other=[rpi['B2AR'],rpi['CCR2']];route=np.max(official[:,musc],axis=1)>np.max(official[:,other],axis=1);hier=official.copy();hier[route][:,musc]=subtype[route][:,musc]
    cal_delta=np.asarray([cal_matrix[cmi[str(r.canonical_smiles)],cpi['M4']]-cal_matrix[cmi[str(r.canonical_smiles)],cpi['M2']] for r in calibration.itertuples()]);cal_y=(calibration.positive_subtype.to_numpy()=='M4');candidates=np.r_[cal_delta.min()-1e-6,(np.sort(np.unique(cal_delta))[:-1]+np.sort(np.unique(cal_delta))[1:])/2,cal_delta.max()+1e-6];best=None
    for threshold in candidates:
        pred=cal_delta>threshold;bal=.5*(pred[cal_y].mean()+(~pred[~cal_y]).mean());key=(float(bal),-abs(float(threshold)))
        if best is None or key>best[0]:best=(key,float(threshold))
    threshold=best[1];local_delta=subtype[:,rpi['M4R']]-subtype[:,rpi['M2R']]
    cofficial=ca['scores'].T.astype(np.float32);coff_delta=np.asarray([cofficial[cmi[str(r.canonical_smiles)],cpi['M4']]-cofficial[cmi[str(r.canonical_smiles)],cpi['M2']] for r in calibration.itertuples()]);Xcal=np.c_[coff_delta,cal_delta];groups=np.asarray([scaffold(v) for v in calibration.canonical_smiles]);oof=np.full(len(calibration),np.nan);splitter=GroupKFold(n_splits=min(5,len(np.unique(groups))))
    for tr_idx,te_idx in splitter.split(Xcal,cal_y,groups):
        fold_model=make_pipeline(StandardScaler(),LogisticRegression(C=1.0,class_weight='balanced',max_iter=2000));fold_model.fit(Xcal[tr_idx],cal_y[tr_idx]);oof[te_idx]=fold_model.decision_function(Xcal[te_idx])
    oof_pred=oof>0;oof_bal=.5*(oof_pred[cal_y].mean()+(~oof_pred[~cal_y]).mean());gate=make_pipeline(StandardScaler(),LogisticRegression(C=1.0,class_weight='balanced',max_iter=2000));gate.fit(Xcal,cal_y);retrieval_official_delta=official[:,rpi['M4R']]-official[:,rpi['M2R']];gate_decision=gate.decision_function(np.c_[retrieval_official_delta,local_delta])
    # Chained indexing above would write a copy; assign the selected block explicitly.
    rr=np.flatnonzero(route);hier[np.ix_(rr,musc)]=subtype[np.ix_(rr,musc)];hier[rr,rpi['M4R']]-=threshold;hier[np.ix_(rr,other)]=-1e6
    oo=np.flatnonzero(~route);hier[np.ix_(oo,musc)]=-1e6
    hier_gate=official.copy();hier_gate[np.ix_(rr,other)]=-1e6;hier_gate[rr,rpi['M2R']]=0.;hier_gate[rr,rpi['M4R']]=gate_decision[rr];hier_gate[np.ix_(oo,musc)]=-1e6
    true_musc=np.isin(labels,musc);family_acc=float(np.mean(route==true_musc))
    report={'method':'Frozen DrugCLIP broad family router plus scaffold-purged muscarinic ECFP3-to-pocket subtype adapter','training_audit':{'muscarinic_triplets_before_purge':int(len(trip0)),'usable_after_embedding_and_both_scaffold_purges':int(len(trip)),'retrieval_scaffold_overlap':0,'calibration_scaffold_overlap':0},'calibration':{'m2_m4_pairs':int(len(calibration)),'threshold_m4_minus_m2':threshold,'threshold_balanced_accuracy':best[0][0],'two_feature_scaffold_oof_balanced_accuracy':float(oof_bal),'two_features':['official_drugclip_M4_minus_M2','ecfp3_adapter_M4_minus_M2']},'family_routing':{'accuracy':family_acc,'muscarinic_recall':float(route[true_musc].mean()),'non_muscarinic_specificity':float((~route[~true_musc]).mean())},'metrics':{'official_drugclip':metrics(labels,official),'hierarchical_threshold_adapter':metrics(labels,hier),'hierarchical_two_feature_gate':metrics(labels,hier_gate)},'threshold_minus_official_macro_recall1_scaffold_bootstrap_95ci':scaffold_bootstrap_delta(table,labels,hier,official,a.bootstrap,20260926),'two_feature_gate_minus_official_macro_recall1_scaffold_bootstrap_95ci':scaffold_bootstrap_delta(table,labels,hier_gate,official,a.bootstrap,20260927),'claim_boundary':'Four-target retrospective GPCR retrieval with external muscarinic calibration and scaffold purge; not binding, potency, PAM efficacy, or broad GPCR SOTA.'}
    out=table[['canonical_molecule_id','target','murcko_scaffold']].copy();out['family_route_muscarinic']=route
    for name,val in [('official',official),('hierarchical',hier),('hierarchical_gate',hier_gate),('local_subtype',subtype)]:
        for j,t in enumerate(TARGETS):out[f'{name}_{t}']=val[:,j]
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');out.to_csv(a.output.with_suffix('.csv'),index=False);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
