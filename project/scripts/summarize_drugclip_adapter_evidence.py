#!/usr/bin/env python3
"""Paired uncertainty audit for DrugCLIP molecular-adapter experiments."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

KEYS=["held_subtype","canonical_smiles","positive_subtype","negative_subtype"]
M4KEY=["canonical_smiles","positive_subtype","negative_subtype"]

def collapse(path,col,keys=KEYS):
    return pd.read_csv(path).groupby(keys,as_index=False)[col].mean()

def cluster_ci(frame,left,right,group="canonical_smiles",draws=5000,seed=20260926):
    rng=np.random.default_rng(seed); groups=frame[group].unique(); by={g:frame.index[frame[group]==g].to_numpy() for g in groups}; vals=[]
    for _ in range(draws):
        sampled=rng.choice(groups,len(groups),replace=True); idx=np.concatenate([by[g] for g in sampled]); vals.append(float((frame.loc[idx,left]>0).mean()-(frame.loc[idx,right]>0).mean()))
    return list(map(float,np.quantile(vals,[.025,.5,.975])))

def fold_metrics(frame,cols):
    out={}
    for held,g in frame.groupby('held_subtype'):
        out[held]={c:float((g[c]>0).mean()) for c in cols}
        out[held]['ecfp3_minus_drugclip_ci']=cluster_ci(g,'ecfp3','drugclip')
        out[held]['ecfp3_minus_ecfp2_ci']=cluster_ci(g,'ecfp3','ecfp2')
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--drugclip-loso',type=Path,required=True);p.add_argument('--ecfp2-loso',type=Path,required=True);p.add_argument('--ecfp3-loso',type=Path,required=True);p.add_argument('--concat-loso',type=Path,required=True);p.add_argument('--drugclip-m4',type=Path,required=True);p.add_argument('--ecfp2-m4',type=Path,required=True);p.add_argument('--selected-m4',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    d=collapse(a.drugclip_loso,'targeted_delta').rename(columns={'targeted_delta':'drugclip'});e2=collapse(a.ecfp2_loso,'targeted_delta').rename(columns={'targeted_delta':'ecfp2'});e3=collapse(a.ecfp3_loso,'targeted_delta').rename(columns={'targeted_delta':'ecfp3'});co=collapse(a.concat_loso,'concat_delta').rename(columns={'concat_delta':'concat'})
    loso=d.merge(e2,on=KEYS).merge(e3,on=KEYS).merge(co,on=KEYS);cols=['drugclip','ecfp2','ecfp3','concat'];folds=fold_metrics(loso,cols)
    macro={c:float(np.mean([folds[h][c] for h in folds])) for c in cols}
    dm=pd.read_csv(a.drugclip_m4)[M4KEY+['gpcr_triplet_delta']].rename(columns={'gpcr_triplet_delta':'drugclip'});e2m=pd.read_csv(a.ecfp2_m4)[M4KEY+['ecfp_pocket_delta']].rename(columns={'ecfp_pocket_delta':'ecfp2'});e3m=pd.read_csv(a.selected_m4)[M4KEY+['selected_adapter_delta']].rename(columns={'selected_adapter_delta':'ecfp3'});m4=dm.merge(e2m,on=M4KEY).merge(e3m,on=M4KEY)
    m4_metrics={c:float((m4[c]>0).mean()) for c in ['drugclip','ecfp2','ecfp3']};m4_metrics['ecfp3_minus_drugclip_ci']=cluster_ci(m4,'ecfp3','drugclip');m4_metrics['ecfp3_minus_ecfp2_ci']=cluster_ci(m4,'ecfp3','ecfp2')
    report={'loso_protocol':'Molecule-disjoint and held-pocket-disjoint target LOSO; ECFP3 was selected later on a separate non-M4 calibration panel, so non-M4 LOSO comparisons are developmental.','loso_fold_metrics':folds,'loso_macro_accuracy':macro,'frozen_m4_protocol':'ECFP3 selected without M4 labels, then evaluated once on 102 molecules / 153 pairs.','frozen_m4':m4_metrics,'claim_boundary':'Subtype activity ranking only; no PAM, affinity, efficacy, or prospective experimental claim.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
