# -*- coding: utf-8 -*-
"""Few-shot medicinal-chemistry-series adaptation benchmark for M4 PAM potency.

For every held-out source component, train a zero-shot RF on all other series.
Then reveal k potency anchors from the held-out series and correct query
predictions using similarity-weighted support residuals.

This simulates the intended workflow: known PAM anchors or a tiny functional
assay batch are available before ranking generated analogues in that series.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error

RDLogger.DisableLog("rdApp.*")
PROJECT=Path(__file__).resolve().parents[1]
DATA=PROJECT/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv"
OUT=PROJECT/"results"/"pacer_fewshot_v01"
SHOTS=[0,1,3,5]
REPEATS=30
SEED=42

DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,
      Descriptors.NumHAcceptors,Descriptors.TPSA,Descriptors.NumRotatableBonds,
      Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]

def features(smiles):
    fps=[]; desc=[]; bitvectors=[]
    for s in smiles:
        m=Chem.MolFromSmiles(s); fp=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048)
        a=np.zeros(2048,np.float32); DataStructs.ConvertToNumpyArray(fp,a)
        fps.append(a); bitvectors.append(fp); desc.append([f(m) for f in DESC])
    return np.asarray(fps),np.asarray(desc,np.float32),bitvectors

def prep(tr_fp,te_fp,tr_d,te_d):
    bits=np.argsort(tr_fp.var(0))[-512:]; mean=tr_d.mean(0); std=tr_d.std(0); std[std<1e-8]=1
    return np.hstack([tr_fp[:,bits],(tr_d-mean)/std]),np.hstack([te_fp[:,bits],(te_d-mean)/std])

def adapt(base, support_idx, query_idx, y, bitvectors, bandwidth=0.25, shrinkage=2.0):
    residual=y[support_idx]-base[support_idx]; output=[]
    for q in query_idx:
        sims=np.asarray(DataStructs.BulkTanimotoSimilarity(bitvectors[q],[bitvectors[i] for i in support_idx]))
        weights=np.exp((sims-sims.max())/bandwidth)
        correction=float(np.sum(weights*residual)/(np.sum(weights)+shrinkage))
        output.append(base[q]+correction)
    return np.asarray(output)

def main():
    OUT.mkdir(parents=True,exist_ok=True); d=pd.read_csv(DATA).reset_index(drop=True)
    fp,desc,bits=features(d.canonical_smiles); y=d.pEC50.to_numpy(float); groups=d.source_component.to_numpy(str)
    groups_eval=[g for g,n in d.source_component.value_counts().items() if n>=8]
    rows=[]; rng=np.random.RandomState(SEED)
    for group in groups_eval:
        target=np.where(groups==group)[0]; train=np.where(groups!=group)[0]
        train_x,target_x=prep(fp[train],fp[target],desc[train],desc[target])
        model=RandomForestRegressor(n_estimators=600,min_samples_leaf=2,random_state=SEED,n_jobs=-1)
        model.fit(train_x,y[train]); base_all=np.full(len(d),np.nan); base_all[target]=model.predict(target_x)
        for shot in SHOTS:
            reps=1 if shot==0 else REPEATS
            if shot>=len(target)-2: continue
            for repeat in range(reps):
                support=np.asarray([],dtype=int) if shot==0 else rng.choice(target,shot,replace=False)
                query=np.asarray([i for i in target if i not in set(support)],dtype=int)
                pred=base_all[query] if shot==0 else adapt(base_all,support,query,y,bits)
                rows.append({"group":group,"n_group":len(target),"shot":shot,"repeat":repeat,
                             "n_query":len(query),"Spearman":float(spearmanr(y[query],pred).statistic),
                             "MAE":float(mean_absolute_error(y[query],pred)),
                             "support_ids":"|".join(d.loc[support,'canonical_molecule_id'])})
        print(group,"done",flush=True)
    result=pd.DataFrame(rows); result.to_csv(OUT/"fewshot_episodes.csv",index=False)
    # Equal weight per series: average repeats within series, then average series.
    per_series=result.groupby(["shot","group"])[["Spearman","MAE"]].mean().reset_index()
    summary=per_series.groupby("shot").agg(
        macro_Spearman=("Spearman","mean"),median_Spearman=("Spearman","median"),
        worst_Spearman=("Spearman","min"),positive_series=("Spearman",lambda x:int((x>0).sum())),
        n_series=("Spearman","size"),macro_MAE=("MAE","mean")
    ).reset_index()
    summary.to_csv(OUT/"fewshot_summary.csv",index=False)
    payload={str(int(r.shot)):{k:(int(v) if k in ['positive_series','n_series'] else float(v))
                              for k,v in r.drop('shot').items()} for _,r in summary.iterrows()}
    with open(OUT/"fewshot_metrics.json","w",encoding="utf-8") as f:json.dump(payload,f,ensure_ascii=False,indent=2)
    print(json.dumps(payload,ensure_ascii=False,indent=2))

if __name__=="__main__":main()
