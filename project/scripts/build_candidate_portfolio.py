# -*- coding: utf-8 -*-
"""Build a leakage-aware, pre-docking candidate portfolio without a weighted mega-score."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem,Descriptors,QED,rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog,FilterCatalogParams
from sklearn.ensemble import RandomForestClassifier

RDLogger.DisableLog("rdApp.*");P=Path(__file__).resolve().parents[1]
POT=P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv";CLS=P/"data"/"benchmarks"/"m4_pam_v1"/"pam_vs_inactive.csv"
GEN=P/"results"/"generated";OUT=P/"results"/"pacer_candidates_v01";OUT.mkdir(parents=True,exist_ok=True)
DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,Descriptors.TPSA,
Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]

def fp(m):return AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048)
def arr(f):
 a=np.zeros(2048,np.float32);DataStructs.ConvertToNumpyArray(f,a);return a
def xmat(smiles):
 out=[]
 for s in smiles:
  m=Chem.MolFromSmiles(s);out.append(np.r_[arr(fp(m)),[f(m) for f in DESC]])
 return np.asarray(out,np.float32)
def main():
 pot=pd.read_csv(POT);known_m=[Chem.MolFromSmiles(s) for s in pot.canonical_smiles]
 known_smiles={Chem.MolToSmiles(m,isomericSmiles=True) for m in known_m};known_fp=[fp(m) for m in known_m];y=pot.pEC50.to_numpy(float)
 candidates={}
 for file,source in [(GEN/"generated_pam_analogs.csv","fragment"),(GEN/"lstm_generated_pam_analogs.csv","lstm"),(GEN/"gpt_generated_pam_analogs.csv","gpt"),(GEN/"diffusion_generated_pam_analogs.csv","diffusion")]:
  if file.exists():
   for s in pd.read_csv(file).smiles:
    m=Chem.MolFromSmiles(s)
    if m is None:continue
    can=Chem.MolToSmiles(m,isomericSmiles=True)
    candidates.setdefault(can,set()).add(source)
 rows=[];params=FilterCatalogParams();params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS);catalog=FilterCatalog(params)
 for s,sources in candidates.items():
  m=Chem.MolFromSmiles(s)
  if m is None:continue
  can=Chem.MolToSmiles(m,isomericSmiles=True)
  if can in known_smiles:continue
  f=fp(m);sims=np.asarray(DataStructs.BulkTanimotoSimilarity(f,known_fp));order=np.argsort(sims)[::-1];top=order[:5];w=np.maximum(sims[top],1e-4)**3
  local=float(np.average(y[top],weights=w));local_sd=float(np.sqrt(np.average((y[top]-local)**2,weights=w)))
  maxsim=float(sims[top[0]]);nneigh=int((sims>=.45).sum())
  tier="local" if .55<=maxsim<=.85 and nneigh>=3 else ("exploratory" if .35<=maxsim<.55 else "reject_domain")
  mw=Descriptors.MolWt(m);logp=Descriptors.MolLogP(m);tpsa=Descriptors.TPSA(m);qed=QED.qed(m);pains=int(catalog.HasMatch(m))
  druglike=int(250<=mw<=550 and -1<=logp<=5.5 and tpsa<=140 and Descriptors.NumRotatableBonds(m)<=10 and pains==0)
  rows.append({"canonical_smiles":can,"generator":"+".join(sorted(sources)),"domain_tier":tier,"max_tanimoto":maxsim,"n_neighbors_045":nneigh,
   "local_knn_pEC50":local,"local_neighbor_sd":local_sd,"MW":mw,"logP":logp,"TPSA":tpsa,"QED":qed,"PAINS":pains,"druglike":druglike,
   "nearest_known_id":pot.iloc[top[0]].canonical_molecule_id})
 d=pd.DataFrame(rows).drop_duplicates("canonical_smiles")
 # Strict-inactive probability is a risk flag only; its series-held-out calibration is weak.
 cls=pd.read_csv(CLS);xc=xmat(cls.canonical_smiles);model=RandomForestClassifier(n_estimators=600,min_samples_leaf=2,class_weight="balanced",random_state=42,n_jobs=-1);model.fit(xc,cls.target)
 d["strict_inactive_risk_ref"]=1-model.predict_proba(xmat(d.canonical_smiles))[:,1]
 # Full-data potency model is reference-only; local kNN is the promoted in-domain estimate.
 reg=LGBMRegressor(n_estimators=500,learning_rate=.025,num_leaves=15,max_depth=5,min_child_samples=12,reg_alpha=.1,reg_lambda=1,random_state=42,verbosity=-1);reg.fit(xmat(pot.canonical_smiles),y);d["global_LGBM_pEC50_ref"]=reg.predict(xmat(d.canonical_smiles))
 eligible=d[(d.domain_tier.isin(["local","exploratory"]))&(d.druglike==1)&(d.strict_inactive_risk_ref<.45)].copy()
 # Pre-docking portfolio: local quality/diversity + a bounded exploratory allocation.
 local=eligible[eligible.domain_tier=="local"].sort_values(["local_knn_pEC50","local_neighbor_sd","QED"],ascending=[False,True,False]).head(160)
 explore=eligible[eligible.domain_tier=="exploratory"].sort_values(["QED","strict_inactive_risk_ref"],ascending=[False,True]).head(40)
 portfolio=pd.concat([local,explore]).drop_duplicates("canonical_smiles").reset_index(drop=True);portfolio.insert(0,"candidate_id",[f"PACER{i:04d}" for i in range(1,len(portfolio)+1)])
 d.to_csv(OUT/"all_generated_audit.csv",index=False);portfolio.to_csv(OUT/"predock_portfolio.csv",index=False)
 report={"n_raw_unique":len(d),"domain_counts":d.domain_tier.value_counts().to_dict(),"n_eligible":len(eligible),"n_predock":len(portfolio),
 "generator_counts":portfolio.generator.value_counts().to_dict(),"tier_counts":portfolio.domain_tier.value_counts().to_dict(),
 "policy":"No weighted total score. Exact known molecules excluded; local potency used only inside similarity domain; classifier probability is a risk flag."}
 (OUT/"audit.json").write_text(json.dumps(report,indent=2),encoding="utf-8");print(json.dumps(report,indent=2))
if __name__=="__main__":main()
