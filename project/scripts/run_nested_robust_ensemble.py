# -*- coding: utf-8 -*-
"""Nested source-LOSO robust rank ensemble for PACER-M4.

Inner source folds select a convex combination of ligand LightGBM, a
source-balanced LightGBM, RF, and ACh-state IFP.  The outer source is never
used for expert training or weight selection.
"""
from __future__ import annotations
import itertools,json
from pathlib import Path
import numpy as np,pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem,Descriptors
from scipy.stats import rankdata,spearmanr
from sklearn.ensemble import RandomForestRegressor

RDLogger.DisableLog("rdApp.*");P=Path(__file__).resolve().parents[1];OUT=P/"results"/"pacer_nested_robust_ensemble_v01";OUT.mkdir(parents=True,exist_ok=True)
DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,Descriptors.TPSA,Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]
def ligx(smiles):
 z=[]
 for s in smiles:
  m=Chem.MolFromSmiles(s);f=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048);a=np.zeros(2048,np.float32);DataStructs.ConvertToNumpyArray(f,a);z.append(np.r_[a,[q(m) for q in DESC]])
 return np.asarray(z,np.float32)
def rho(y,p):return float(spearmanr(y,p).statistic) if len(y)>2 and np.std(y)>0 and np.std(p)>0 else 0.
def weights(groups):
 u,c=np.unique(groups,return_counts=True);mp=dict(zip(u,c));return np.asarray([1/mp[x] for x in groups])*len(groups)/len(u)
def fit_predict(X,S,y,g,tr,te):
 selected=np.argsort(X[tr,:2048].var(0))[-1024:];xx=np.c_[X[:,selected],X[:,2048:]]
 kw=dict(n_estimators=350,learning_rate=.03,num_leaves=15,max_depth=5,min_child_samples=12,reg_alpha=.1,reg_lambda=1,random_state=42,verbosity=-1,n_jobs=4)
 a=LGBMRegressor(**kw).fit(xx[tr],y[tr]).predict(xx[te]);b=LGBMRegressor(**kw).fit(xx[tr],y[tr],sample_weight=weights(g[tr])).predict(xx[te])
 c=RandomForestRegressor(n_estimators=400,min_samples_leaf=2,max_features=.35,random_state=42,n_jobs=-1).fit(X[tr],y[tr]).predict(X[te])
 q=RandomForestRegressor(n_estimators=400,min_samples_leaf=3,max_features=.7,random_state=42,n_jobs=-1).fit(S[tr],y[tr],sample_weight=weights(g[tr])).predict(S[te])
 return np.c_[a,b,c,q]
def norm_rows(pred,groups):
 out=np.zeros_like(pred,float)
 for g in np.unique(groups):
  ix=np.where(groups==g)[0]
  for j in range(pred.shape[1]):out[ix,j]=rankdata(pred[ix,j])/len(ix)
 return out
def grid():
 return [np.asarray(x)/4 for x in itertools.product(range(5),repeat=4) if sum(x)==4]
def main():
 d=pd.read_csv(P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv");st=pd.read_csv(P/"results"/"pacer_structure_7trs_v01"/"docking_features.csv")
 cols=[x for x in st.columns if x in {"crystal_ifp_jaccard","unique_pocket_residue_contacts","pocket_atom_pair_contacts","ligand_atoms_contacting_pocket"} or x.startswith("res")]
 st=st.set_index("canonical_molecule_id").loc[d.canonical_molecule_id];S=st[cols].fillna(st[cols].median()).to_numpy(float);X=ligx(d.canonical_smiles);y=d.pEC50.to_numpy(float);g=d.source_component.astype(str).to_numpy();pred=np.full(len(d),np.nan);base=np.full(len(d),np.nan);folds=[]
 eligible=[x for x,n in d.source_component.value_counts().items() if n>=8]
 for outer in eligible:
  otr=np.where(g!=outer)[0];ote=np.where(g==outer)[0];inner_pred=[];inner_y=[];inner_g=[]
  for inner in [x for x in eligible if x!=outer]:
   ite=np.where(g==inner)[0];itr=np.where((g!=outer)&(g!=inner))[0]
   z=fit_predict(X,S,y,g,itr,ite);inner_pred.append(z);inner_y.extend(y[ite]);inner_g.extend(g[ite])
  ip=norm_rows(np.vstack(inner_pred),np.asarray(inner_g));iy=np.asarray(inner_y);ig=np.asarray(inner_g);best=None
  for w in grid():
   p=ip@w;rs=[rho(iy[ig==x],p[ig==x]) for x in np.unique(ig)];score=np.mean(rs)+.25*np.min(rs)
   if best is None or score>best[0]:best=(score,w,np.mean(rs),np.min(rs))
  op=fit_predict(X,S,y,g,otr,ote);op=norm_rows(op,g[ote]);pred[ote]=op@best[1];base[ote]=op[:,0]
  folds.append({"group":outer,"n":len(ote),"Spearman":rho(y[ote],pred[ote]),"base_Spearman":rho(y[ote],base[ote]),"weights":best[1].tolist(),"inner_macro":best[2],"inner_worst":best[3]});print(outer,folds[-1],flush=True)
 rs=np.asarray([x["Spearman"] for x in folds]);rb=np.asarray([x["base_Spearman"] for x in folds]);delta=rs-rb;rng=np.random.default_rng(42);boot=np.asarray([np.mean(rng.choice(delta,len(delta),replace=True)) for _ in range(20000)])
 report={"protocol":"nested source LOSO; inner objective macro Spearman + 0.25*worst; convex step 0.25; experts LGBM/LGBM-balanced/RF/ACh-IFP","folds":folds,"aggregate":{"macro_Spearman":float(rs.mean()),"worst_Spearman":float(rs.min()),"positive_series":int((rs>0).sum()),"base_macro_in_run":float(rb.mean()),"delta_macro":float(delta.mean()),"delta_95CI":[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))],"P_delta_gt_0":float((boot>0).mean())}}
 pd.DataFrame({"canonical_molecule_id":d.canonical_molecule_id,"source_component":g,"pEC50":y,"robust_ensemble":pred,"base_LGBM_rank":base}).to_csv(OUT/"predictions.csv",index=False);(OUT/"metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8");print(json.dumps(report["aggregate"],indent=2))
if __name__=="__main__":main()
