# -*- coding: utf-8 -*-
"""Retrospective temporal benchmark: tune <=2021, evaluate 2022+ sources."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem,Descriptors
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error

RDLogger.DisableLog("rdApp.*");P=Path(__file__).resolve().parents[1]
DATA=P/"data"/"benchmarks"/"m4_pam_v1"/"temporal"/"potency_molecules_temporal.csv"
SFEAT=P/"results"/"pacer_structure_7trs_v01"/"docking_features.csv";OUT=P/"results"/"pacer_temporal_v01"
D=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,Descriptors.TPSA,
Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]

def features(smiles):
 fp=[];ds=[]
 for s in smiles:
  m=Chem.MolFromSmiles(s);b=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048);a=np.zeros(2048,np.float32);DataStructs.ConvertToNumpyArray(b,a)
  fp.append(a);ds.append([f(m) for f in D])
 return np.asarray(fp),np.asarray(ds,np.float32)
def prep(fp,ds,tr,te,nbit):
 sel=np.argsort(fp[tr].var(0))[-nbit:];mu,sd=ds[tr].mean(0),ds[tr].std(0);sd[sd<1e-8]=1
 return np.hstack([fp[tr][:,sel],(ds[tr]-mu)/sd]),np.hstack([fp[te][:,sel],(ds[te]-mu)/sd])
def score(y,p,g):
 groups=[]
 for x in sorted(set(g)):
  z=g==x
  if z.sum()>=8:groups.append({"group":x,"n":int(z.sum()),"Spearman":float(spearmanr(y[z],p[z]).statistic),"MAE":float(mean_absolute_error(y[z],p[z]))})
 return {"pooled_Spearman":float(spearmanr(y,p).statistic),"MAE":float(mean_absolute_error(y,p)),"groups":groups,
         "macro_group_Spearman":float(np.mean([x["Spearman"] for x in groups])) if groups else None}
def main():
 OUT.mkdir(parents=True,exist_ok=True);d=pd.read_csv(DATA);s=pd.read_csv(SFEAT);s=s[(s.seed==42)&(s.error.fillna("")=="")].drop_duplicates("canonical_molecule_id")
 scols=["crystal_ifp_jaccard","crystal_centroid_distance_A","unique_pocket_residue_contacts"]
 d=d.merge(s[["canonical_molecule_id"]+scols],on="canonical_molecule_id");fp,ds=features(d.canonical_smiles);y=d.pEC50.to_numpy(float);g=d.source_component.astype(str).to_numpy()
 early=d.temporal_split_enriched.eq("train_le_2018").to_numpy();val=d.temporal_split_enriched.eq("validation_2019_2021").to_numpy();test=d.temporal_split_enriched.eq("test_ge_2022").to_numpy()
 configs=[]
 for nbit in (512,1024):
  a,b=prep(fp,ds,early,val,nbit)
  for leaves,depth,child in [(7,3,12),(15,5,12),(31,6,8)]:
   m=LGBMRegressor(n_estimators=500,learning_rate=.025,num_leaves=leaves,max_depth=depth,min_child_samples=child,subsample=.8,colsample_bytree=.7,reg_alpha=.1,reg_lambda=1,random_state=42,n_jobs=4,verbosity=-1)
   m.fit(a,y[early]);p=m.predict(b);r=score(y[val],p,g[val]);configs.append({"family":"2D","params":[nbit,leaves,depth,child],"validation":r})
   m.fit(np.hstack([a,d.loc[early,scols]]),y[early]);p=m.predict(np.hstack([b,d.loc[val,scols]]));r=score(y[val],p,g[val]);configs.append({"family":"2D_plus_ACh_IFP","params":[nbit,leaves,depth,child],"validation":r})
 # Structure-only low-capacity control.
 rf=RandomForestRegressor(n_estimators=600,min_samples_leaf=3,random_state=42,n_jobs=-1);rf.fit(d.loc[early,scols],y[early]);p=rf.predict(d.loc[val,scols])
 configs.append({"family":"ACh_IFP_only","params":[],"validation":score(y[val],p,g[val])})
 best=max(configs,key=lambda x:x["validation"]["macro_group_Spearman"]);train=early|val
 if best["family"]=="ACh_IFP_only":
  model=RandomForestRegressor(n_estimators=600,min_samples_leaf=3,random_state=42,n_jobs=-1);model.fit(d.loc[train,scols],y[train]);pred=model.predict(d.loc[test,scols])
 else:
  nbit,leaves,depth,child=best["params"];a,b=prep(fp,ds,train,test,nbit);model=LGBMRegressor(n_estimators=500,learning_rate=.025,num_leaves=leaves,max_depth=depth,min_child_samples=child,subsample=.8,colsample_bytree=.7,reg_alpha=.1,reg_lambda=1,random_state=42,n_jobs=4,verbosity=-1)
  if best["family"]=="2D_plus_ACh_IFP":a=np.hstack([a,d.loc[train,scols]]);b=np.hstack([b,d.loc[test,scols]])
  model.fit(a,y[train]);pred=model.predict(b)
 report={"split_counts":{"early":int(early.sum()),"validation":int(val.sum()),"test":int(test.sum())},"all_validation_configs":configs,"selected":best,"test":score(y[test],pred,g[test]),
 "caution":"Retrospective pseudo-temporal test. The 2022+ series appeared in earlier LOSO analyses, so this is not a never-seen prospective validation."}
 (OUT/"metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8");o=d.loc[test,["canonical_molecule_id","source_component","publication_year_enriched","pEC50"]].copy();o["prediction"]=pred;o.to_csv(OUT/"test_predictions.csv",index=False)
 print(json.dumps({"selected":best,"test":report["test"]},indent=2))
if __name__=="__main__":main()
