# -*- coding: utf-8 -*-
"""LightGBM strong baseline on the frozen 12-series M4 potency LOSO task."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

RDLogger.DisableLog("rdApp.*")
PROJECT=Path(__file__).resolve().parents[1]
DATA=PROJECT/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv"
OUT=PROJECT/"results"/"pacer_lightgbm_loso_v01"
DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,
      Descriptors.TPSA,Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,
      Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]

def features(smiles):
    fp=[]; ds=[]
    for s in smiles:
        m=Chem.MolFromSmiles(s); b=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048)
        a=np.zeros(2048,np.float32); DataStructs.ConvertToNumpyArray(b,a); fp.append(a)
        ds.append([f(m) for f in DESC])
    return np.asarray(fp),np.asarray(ds,np.float32)

def main():
    OUT.mkdir(parents=True,exist_ok=True); d=pd.read_csv(DATA); fp,ds=features(d.canonical_smiles)
    y=d.pEC50.to_numpy(float); groups=d.source_component.astype(str).to_numpy(); pred=np.full(len(d),np.nan); folds=[]
    for g,n in d.source_component.value_counts().items():
        if n<8: continue
        te=groups==g; tr=~te; selected=np.argsort(fp[tr].var(0))[-1024:]
        mean,std=ds[tr].mean(0),ds[tr].std(0); std[std<1e-8]=1
        xtr=np.hstack([fp[tr][:,selected],(ds[tr]-mean)/std]); xte=np.hstack([fp[te][:,selected],(ds[te]-mean)/std])
        model=LGBMRegressor(n_estimators=500,learning_rate=.025,num_leaves=15,max_depth=5,
                            min_child_samples=12,subsample=.8,colsample_bytree=.7,
                            reg_alpha=.1,reg_lambda=1.0,random_state=42,n_jobs=4,verbosity=-1)
        model.fit(xtr,y[tr]); pred[te]=model.predict(xte)
        rho=float(spearmanr(y[te],pred[te]).statistic); mae=float(mean_absolute_error(y[te],pred[te]))
        folds.append({"group":g,"n":int(te.sum()),"Spearman":rho,"MAE":mae}); print(g,rho,mae,flush=True)
    r=np.asarray([x["Spearman"] for x in folds]); m=np.asarray([x["MAE"] for x in folds])
    report={"model":"LightGBM fixed strong baseline","folds":folds,"aggregate":{
        "macro_Spearman":float(np.mean(r)),"median_Spearman":float(np.median(r)),
        "worst_Spearman":float(np.min(r)),"positive_rho_series":int(np.sum(r>0)),
        "n_series":len(folds),"macro_MAE":float(np.mean(m))}}
    o=d[["canonical_molecule_id","canonical_smiles","pEC50","source_component"]].copy(); o["LightGBM"]=pred
    o.to_csv(OUT/"loso_predictions.csv",index=False)
    (OUT/"loso_metrics.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(report["aggregate"],ensure_ascii=False,indent=2))

if __name__=="__main__":main()
