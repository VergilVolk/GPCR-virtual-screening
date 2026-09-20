# -*- coding: utf-8 -*-
"""Frozen ligand-only LightGBM v2 baseline (descriptor scaling ablation)."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem,Descriptors
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error
RDLogger.DisableLog("rdApp.*");P=Path(__file__).resolve().parents[1];OUT=P/"results"/"pacer_lightgbm_v2_loso";OUT.mkdir(parents=True,exist_ok=True)
D=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,Descriptors.TPSA,Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]
def feat(smiles):
 z=[]
 for s in smiles:
  m=Chem.MolFromSmiles(s);f=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048);a=np.zeros(2048,np.float32);DataStructs.ConvertToNumpyArray(f,a);z.append(np.r_[a,[q(m) for q in D]])
 return np.asarray(z,np.float32)
def main():
 d=pd.read_csv(P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv");X=feat(d.canonical_smiles);y=d.pEC50.to_numpy(float);g=d.source_component.astype(str).to_numpy();pred=np.full(len(d),np.nan);folds=[]
 for group,n in d.source_component.value_counts().items():
  if n<8:continue
  tr=g!=group;te=~tr;sel=np.argsort(X[tr,:2048].var(0))[-1024:];xx=np.c_[X[:,sel],X[:,2048:]];m=LGBMRegressor(n_estimators=350,learning_rate=.03,num_leaves=15,max_depth=5,min_child_samples=12,reg_alpha=.1,reg_lambda=1,random_state=42,verbosity=-1,n_jobs=4).fit(xx[tr],y[tr]);pred[te]=m.predict(xx[te]);folds.append({"group":group,"n":int(te.sum()),"Spearman":float(spearmanr(y[te],pred[te]).statistic),"MAE":float(mean_absolute_error(y[te],pred[te]))})
 r=np.asarray([x["Spearman"] for x in folds]);report={"model":"LightGBM v2 fixed ligand-only","folds":folds,"aggregate":{"macro_Spearman":float(r.mean()),"worst_Spearman":float(r.min()),"positive_series":int((r>0).sum()),"macro_MAE":float(np.mean([x["MAE"] for x in folds]))}}
 pd.DataFrame({"canonical_molecule_id":d.canonical_molecule_id,"source_component":g,"pEC50":y,"LightGBM_v2":pred}).to_csv(OUT/"predictions.csv",index=False);(OUT/"metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8");print(json.dumps(report["aggregate"],indent=2))
if __name__=="__main__":main()
