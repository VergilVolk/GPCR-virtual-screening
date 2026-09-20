# -*- coding: utf-8 -*-
"""PACER-StateMoE v0.1: nested, group-balanced structural residual experts.

The ligand-only LightGBM is the null expert. Low-dimensional receptor-state
blocks may only correct its residual when selected by inner leave-one-source
validation. The outer medicinal-chemistry series remains untouched.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem,Descriptors
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1]; DATA=P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv"
FILES={"Q":P/"results"/"pacer_structure_loso_v01"/"docking_features.csv",
       "P":P/"results"/"pacer_structure_7trp_v01"/"docking_features.csv",
       "S":P/"results"/"pacer_structure_7trs_v01"/"docking_features.csv"}
OUT=P/"results"/"pacer_state_residual_moe_v01"
DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,
Descriptors.TPSA,Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]

def molecular(smiles):
    fp=[];ds=[]
    for s in smiles:
        m=Chem.MolFromSmiles(s);b=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048);a=np.zeros(2048,np.float32)
        DataStructs.ConvertToNumpyArray(b,a);fp.append(a);ds.append([f(m) for f in DESC])
    return np.asarray(fp),np.asarray(ds,np.float32)

def lgbm():return LGBMRegressor(n_estimators=500,learning_rate=.025,num_leaves=15,max_depth=5,min_child_samples=12,
subsample=.8,colsample_bytree=.7,reg_alpha=.1,reg_lambda=1,random_state=42,n_jobs=4,verbosity=-1)

def x2(fp,ds,tr,te):
    sel=np.argsort(fp[tr].var(0))[-1024:];mean,std=ds[tr].mean(0),ds[tr].std(0);std[std<1e-8]=1
    return np.hstack([fp[tr][:,sel],(ds[tr]-mean)/std]),np.hstack([fp[te][:,sel],(ds[te]-mean)/std])

def group_weights(g):
    counts=pd.Series(g).value_counts();return np.asarray([1/counts[x] for x in g],float)*len(g)/len(counts)

def main():
    OUT.mkdir(parents=True,exist_ok=True);d=pd.read_csv(DATA)
    block_cols={}
    for tag,path in FILES.items():
        f=pd.read_csv(path);f=f[(f.seed==42)&(f.error.fillna("")=="")].drop_duplicates("canonical_molecule_id")
        cols=["vina_affinity","crystal_ifp_jaccard","crystal_centroid_distance_A","unique_pocket_residue_contacts",
              "pocket_atom_pair_contacts","tyr439_min_distance_A"]
        names=[f"{tag}_{c}" for c in cols];block_cols[tag]=names
        d=d.merge(f[["canonical_molecule_id"]+cols].rename(columns=dict(zip(cols,names))),on="canonical_molecule_id")
    configs={"S_IFP":["S_crystal_ifp_jaccard","S_crystal_centroid_distance_A","S_unique_pocket_residue_contacts"],
             "All_IFP":sum([[f"{t}_crystal_ifp_jaccard",f"{t}_crystal_centroid_distance_A",f"{t}_unique_pocket_residue_contacts"] for t in "QPS"],[]),
             "All_State_Summary":sum(block_cols.values(),[])}
    fp,ds=molecular(d.canonical_smiles);y=d.pEC50.to_numpy(float);groups=d.source_component.astype(str).to_numpy()
    outer=[g for g,n in d.source_component.value_counts().items() if n>=8];base_all=np.full(len(d),np.nan);moe_all=np.full(len(d),np.nan);folds=[]
    for outer_g in outer:
        outer_te=groups==outer_g;outer_tr=~outer_te
        a,b=x2(fp,ds,outer_tr,outer_te);base=lgbm();base.fit(a,y[outer_tr]);base_test=base.predict(b);base_all[outer_te]=base_test
        train_idx=np.where(outer_tr)[0];inner_groups=[g for g,n in pd.Series(groups[outer_tr]).value_counts().items() if n>=5]
        oof=np.full(len(train_idx),np.nan)
        for ig in inner_groups:
            ite=groups[train_idx]==ig;itr=~ite;xa,xb=x2(fp[train_idx],ds[train_idx],itr,ite)
            m=lgbm();m.fit(xa,y[train_idx][itr]);oof[ite]=m.predict(xb)
        valid=np.isfinite(oof);idx=train_idx[valid];resid=y[idx]-oof[valid];igroups=groups[idx]
        best=None
        for cname,cols in configs.items():
            xs=d.loc[idx,cols].to_numpy(float)
            for alpha in (0.1,1.0,10.0,100.0):
                for beta in (0.25,0.5,1.0):
                    rhos=[]
                    for vg in sorted(set(igroups)):
                        va=igroups==vg;tr=~va
                        if va.sum()<5:continue
                        r=make_pipeline(StandardScaler(),Ridge(alpha=alpha));r.fit(xs[tr],resid[tr],ridge__sample_weight=group_weights(igroups[tr]))
                        pred=oof[valid][va]+beta*r.predict(xs[va]);rhos.append(spearmanr(y[idx][va],pred).statistic)
                    score=float(np.nanmean(rhos))
                    if best is None or score>best[0]:best=(score,cname,alpha,beta)
        _,cname,alpha,beta=best;cols=configs[cname];r=make_pipeline(StandardScaler(),Ridge(alpha=alpha))
        r.fit(d.loc[idx,cols].to_numpy(float),resid,ridge__sample_weight=group_weights(igroups))
        moe=base_test+beta*r.predict(d.loc[outer_te,cols].to_numpy(float));moe_all[outer_te]=moe
        folds.append({"group":outer_g,"n":int(outer_te.sum()),"selected":{"config":cname,"alpha":alpha,"beta":beta,"inner_macro":best[0]},
          "Base":{"Spearman":float(spearmanr(y[outer_te],base_test).statistic),"MAE":float(mean_absolute_error(y[outer_te],base_test))},
          "StateResidualMoE":{"Spearman":float(spearmanr(y[outer_te],moe).statistic),"MAE":float(mean_absolute_error(y[outer_te],moe))}})
        print(folds[-1],flush=True)
    agg={}
    for name in ("Base","StateResidualMoE"):
        r=np.array([f[name]["Spearman"] for f in folds]);ma=np.array([f[name]["MAE"] for f in folds])
        agg[name]={"macro_Spearman":float(r.mean()),"median_Spearman":float(np.median(r)),"worst_Spearman":float(r.min()),
                   "positive_rho_series":int((r>0).sum()),"macro_MAE":float(ma.mean())}
    delta=np.array([f["StateResidualMoE"]["Spearman"]-f["Base"]["Spearman"] for f in folds]);rng=np.random.default_rng(42)
    boot=np.array([rng.choice(delta,len(delta),replace=True).mean() for _ in range(100000)])
    report={"method":"nested group-balanced structural residual experts","folds":folds,"aggregate":agg,
            "paired_delta":{"mean":float(delta.mean()),"median":float(np.median(delta)),"ci95":np.quantile(boot,[.025,.975]).tolist(),"p_gt_0":float((boot>0).mean())}}
    (OUT/"metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    out=d[["canonical_molecule_id","source_component","pEC50"]].copy();out["Base"]=base_all;out["StateResidualMoE"]=moe_all;out.to_csv(OUT/"predictions.csv",index=False)
    print(json.dumps({"aggregate":agg,"paired_delta":report["paired_delta"]},indent=2))

if __name__=="__main__":main()
