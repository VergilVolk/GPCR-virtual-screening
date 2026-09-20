# -*- coding: utf-8 -*-
"""Scaffold-held-out M4R early-enrichment baselines and LambdaRank model."""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
from lightgbm import LGBMClassifier,LGBMRanker
from rdkit import Chem,RDLogger
from rdkit.Chem import Descriptors,rdFingerprintGenerator
from rdkit.ML.Scoring import Scoring
from sklearn.metrics import roc_auc_score,average_precision_score
RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1];B=P/"data"/"benchmarks"/"m4r_enrichment_v1";OUT=P/"results"/"pacer_m4r_enrichment_v01";CACHE=B/"feature_cache"
DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,Descriptors.TPSA,Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]
def features(d):
 xp=CACHE/"ecfp4_2048.npz";dp=CACHE/"descriptors.npy";ip=CACHE/"molecule_ids.txt"
 if xp.exists() and dp.exists() and ip.exists() and ip.read_text().splitlines()==d.molecule_id.tolist():return sparse.load_npz(xp),np.load(dp)
 gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048);rows=[];cols=[];ds=[];t=time.time()
 for i,s in enumerate(d.canonical_smiles):
  m=Chem.MolFromSmiles(s);bits=gen.GetFingerprint(m).GetOnBits();rows.extend([i]*len(bits));cols.extend(bits);ds.append([q(m) for q in DESC])
  if (i+1)%20000==0:print(f"features {i+1}/{len(d)} elapsed={time.time()-t:.0f}s",flush=True)
 X=sparse.csr_matrix((np.ones(len(rows),np.float32),(rows,cols)),shape=(len(d),2048));D=np.asarray(ds,np.float32);sparse.save_npz(xp,X);np.save(dp,D);ip.write_text("\n".join(d.molecule_id));return X,D
def metrics(y,score):
 order=np.argsort(score)[::-1];ys=np.asarray(y,int)[order];n=len(ys);p=ys.sum();prev=p/n
 out={"n":n,"actives":int(p),"ROC_AUC":float(roc_auc_score(y,score)),"PR_AUC":float(average_precision_score(y,score))}
 ranked=[(float(score[i]),int(y[i])) for i in order]
 for alpha in (20,80):out[f"BEDROC{alpha}"]=float(Scoring.CalcBEDROC(ranked,1,alpha))
 for frac in (.005,.01):
  k=max(1,int(np.ceil(frac*n)));hits=int(ys[:k].sum());out[f"EF{frac*100:g}%"]=float((hits/k)/prev);out[f"hits_top_{k}"]=hits
 for k in (100,500):out[f"hits_top{k}"]=int(ys[:min(k,n)].sum())
 return out
def aggregate(folds):
 keys=[k for k in folds[0] if k not in {"fold","n","actives"} and not k.startswith("hits_top_")];return {"macro_"+k:float(np.mean([x[k] for x in folds])) for k in keys}|{"worst_EF0.5%":float(np.min([x["EF0.5%"] for x in folds]))}
def main():
 global OUT,CACHE
 ap=argparse.ArgumentParser();ap.add_argument("--dataset",default="molecules.csv");ap.add_argument("--tag",default="full");a=ap.parse_args();OUT=P/"results"/("pacer_m4r_enrichment_v01" if a.tag=="full" else f"pacer_m4r_enrichment_{a.tag}_v01");OUT.mkdir(parents=True,exist_ok=True);stem=Path(a.dataset).stem;CACHE=B/("feature_cache" if stem=="molecules" else f"feature_cache_{stem}");CACHE.mkdir(exist_ok=True)
 d=pd.read_csv(B/a.dataset);X,D=features(d);XD=sparse.hstack([X,sparse.csr_matrix(D)],format="csr");y=d.is_active.to_numpy(int);fold=d.scaffold_fold.to_numpy(int);reports={};pred=[]
 methods=("PDB_Vina","Ensemble_BEmin","Ensemble_BEavg","DescriptorClassifier","ECFP_Classifier","ECFP_LambdaRank")
 for method in methods:reports[method]={"folds":[]}
 for f in sorted(np.unique(fold)):
  tr=fold!=f;te=~tr;yt=y[te]
  static={"PDB_Vina":-d.loc[te,"pdb_vina"].to_numpy(float),"Ensemble_BEmin":-d.loc[te,"ensemble_BEmin"].fillna(d.ensemble_BEmin.median()).to_numpy(float),"Ensemble_BEavg":-d.loc[te,"ensemble_BEavg"].fillna(d.ensemble_BEavg.median()).to_numpy(float)}
  for name,sc in static.items():reports[name]["folds"].append({"fold":int(f),**metrics(yt,sc)})
  clfkw=dict(n_estimators=350,learning_rate=.04,num_leaves=31,max_depth=7,min_child_samples=50,subsample=.8,colsample_bytree=.6,reg_alpha=.1,reg_lambda=1,class_weight="balanced",random_state=42,n_jobs=6,verbosity=-1)
  dc=LGBMClassifier(**clfkw).fit(D[tr],y[tr]);sc=dc.predict_proba(D[te])[:,1];reports["DescriptorClassifier"]["folds"].append({"fold":int(f),**metrics(yt,sc)})
  ec=LGBMClassifier(**clfkw).fit(XD[tr],y[tr]);sc2=ec.predict_proba(XD[te])[:,1];reports["ECFP_Classifier"]["folds"].append({"fold":int(f),**metrics(yt,sc2)})
  # LightGBM limits a query to 10k rows. Split each natural training screen
  # into four deterministic ID-hash shards without changing outer folds.
  tri=np.where(tr)[0];shard=np.asarray([int(x[-6:])%4 for x in d.molecule_id]);query=fold*4+shard;o=np.argsort(query[tri],kind="stable");tri=tri[o];groups=pd.Series(query[tri]).value_counts(sort=False).sort_index().to_list()
  rank=LGBMRanker(objective="lambdarank",metric="ndcg",label_gain=[0,1],n_estimators=350,learning_rate=.04,num_leaves=31,max_depth=7,min_child_samples=50,subsample=.8,colsample_bytree=.6,reg_alpha=.1,reg_lambda=1,random_state=42,n_jobs=6,verbosity=-1).fit(XD[tri],y[tri],group=groups);sc3=rank.predict(XD[te]);reports["ECFP_LambdaRank"]["folds"].append({"fold":int(f),**metrics(yt,sc3)})
  for i,a,b,c in zip(np.where(te)[0],sc,sc2,sc3):pred.append({"molecule_id":d.iloc[i].molecule_id,"fold":int(f),"is_active":int(y[i]),"DescriptorClassifier":a,"ECFP_Classifier":b,"ECFP_LambdaRank":c})
  print("fold",f,{m:round(reports[m]["folds"][-1]["EF0.5%"],2) for m in methods},flush=True)
 for name in methods:reports[name]["aggregate"]=aggregate(reports[name]["folds"])
 (OUT/"metrics.json").write_text(json.dumps(reports,indent=2),encoding="utf-8");pd.DataFrame(pred).to_csv(OUT/"oof_predictions.csv",index=False);print(json.dumps({k:v["aggregate"] for k,v in reports.items()},indent=2))
if __name__=="__main__":main()
