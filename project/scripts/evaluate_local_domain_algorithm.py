# -*- coding: utf-8 -*-
"""Predeclared risk/coverage curves for local-SAR prediction under source LOSO.

No threshold is selected on held-out labels.  Four fixed similarity thresholds
are reported so that abstention cannot be tuned after seeing test performance.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1];DATA=P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv"
OUT=P/"results"/"pacer_local_domain_v01";THRESHOLDS=(.35,.45,.55,.65)

def fp(s):return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s),2,nBits=2048)
def safe_rho(y,p):return float(spearmanr(y,p).statistic) if len(y)>=3 and np.std(y)>0 and np.std(p)>0 else None

def main():
 OUT.mkdir(parents=True,exist_ok=True);d=pd.read_csv(DATA);fps=[fp(s) for s in d.canonical_smiles];y=d.pEC50.to_numpy(float);groups=d.source_component.astype(str).to_numpy()
 rows=[]
 for g,n in d.source_component.value_counts().items():
  if n<8:continue
  te=np.where(groups==g)[0];tr=np.where(groups!=g)[0]
  for i in te:
   sims=np.asarray(DataStructs.BulkTanimotoSimilarity(fps[i],[fps[j] for j in tr]));o=np.argsort(sims)[::-1][:5];w=np.maximum(sims[o],1e-4)**3;vals=y[tr[o]]
   pred=float(np.average(vals,weights=w));sd=float(np.sqrt(np.average((vals-pred)**2,weights=w)))
   rows.append({"canonical_molecule_id":d.iloc[i].canonical_molecule_id,"group":g,"true":y[i],"pred":pred,"max_similarity":float(sims[o[0]]),"neighbor_sd":sd,"abs_error":abs(y[i]-pred)})
 r=pd.DataFrame(rows);r.to_csv(OUT/"loso_predictions.csv",index=False);curves=[]
 for t in THRESHOLDS:
  q=r[r.max_similarity>=t];folds=[]
  for g,x in q.groupby("group"):
   folds.append({"group":g,"n":len(x),"coverage":len(x)/int((r.group==g).sum()),"Spearman":safe_rho(x.true,x.pred),"MAE":float(mean_absolute_error(x.true,x.pred))})
  valid=[x for x in folds if x["n"]>=5 and x["Spearman"] is not None]
  curves.append({"threshold":t,"n":len(q),"coverage":len(q)/len(r),"pooled_Spearman":safe_rho(q.true,q.pred),"pooled_MAE":float(mean_absolute_error(q.true,q.pred)) if len(q) else None,"macro_Spearman_nge5":float(np.mean([x["Spearman"] for x in valid])) if valid else None,"worst_Spearman_nge5":float(np.min([x["Spearman"] for x in valid])) if valid else None,"n_evaluable_series":len(valid),"folds":folds})
 # Uncertainty audit: neighbor spread should rise with actual error if useful.
 unc={"spearman_neighbor_sd_vs_abs_error":safe_rho(r.neighbor_sd,r.abs_error),"error_by_neighbor_sd_quartile":[]}
 for name,x in r.groupby(pd.qcut(r.neighbor_sd,4,duplicates="drop")):
  unc["error_by_neighbor_sd_quartile"].append({"bin":str(name),"n":len(x),"MAE":float(x.abs_error.mean())})
 report={"protocol":"source-component LOSO; top-5 ECFP4 similarity^3 local mean; fixed thresholds; no held-out threshold selection","n":len(r),"curves":curves,"uncertainty":unc}
 (OUT/"metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8");print(json.dumps({"curves":[{k:v for k,v in x.items() if k!="folds"} for x in curves],"uncertainty":unc},indent=2))
if __name__=="__main__":main()
