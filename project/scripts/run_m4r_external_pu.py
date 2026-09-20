# -*- coding: utf-8 -*-
"""External strict-inactive test of M4R positive-unlabeled recognition.

All external molecule scaffolds are removed from M4R positives and the
unlabeled corpus before training.  The 463 PAM / 66 experimental-inactive set
is never used for fitting or tuning.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
from lightgbm import LGBMClassifier
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import Descriptors,rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics import roc_auc_score,average_precision_score
RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1];M=P/"data"/"benchmarks"/"m4r_enrichment_v1";E=P/"data"/"benchmarks"/"m4_pam_v1"/"pam_vs_inactive.csv";CORP=P/"data"/"generated"/"corpus.smi";OUT=P/"results"/"pacer_m4r_external_pu_v01";OUT.mkdir(parents=True,exist_ok=True)
DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,Descriptors.TPSA,Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount];GEN=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
def scaf(m):
 x=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m),isomericSmiles=True);return x
def featurize(smiles):
 rows=[];cols=[];ds=[];fps=[];valid=[]
 for i,s in enumerate(smiles):
  m=Chem.MolFromSmiles(s)
  if m is None:continue
  f=GEN.GetFingerprint(m);bits=f.GetOnBits();rows.extend([len(valid)]*len(bits));cols.extend(bits);ds.append([q(m) for q in DESC]);fps.append(f);valid.append(s)
 X=sparse.csr_matrix((np.ones(len(rows),np.float32),(rows,cols)),shape=(len(valid),2048));return X,np.asarray(ds,np.float32),fps,valid
def evals(frame,scores):
 out={"pooled_ROC_AUC":float(roc_auc_score(frame.target,scores)),"pooled_PR_AUC":float(average_precision_score(frame.target,scores)),"groups":[]}
 for g,x in frame.assign(score=scores).groupby("source_component"):
  if x.target.nunique()==2 and min(x.target.value_counts())>=5:out["groups"].append({"group":g,"n":len(x),"ROC_AUC":float(roc_auc_score(x.target,x.score)),"PR_AUC":float(average_precision_score(x.target,x.score))})
 out["macro_group_AUC"]=float(np.mean([x["ROC_AUC"] for x in out["groups"]]));out["worst_group_AUC"]=float(np.min([x["ROC_AUC"] for x in out["groups"]]));return out
def model(seed=42):return LGBMClassifier(n_estimators=350,learning_rate=.04,num_leaves=31,max_depth=7,min_child_samples=30,subsample=.8,colsample_bytree=.6,reg_alpha=.1,reg_lambda=1,class_weight="balanced",random_state=seed,n_jobs=6,verbosity=-1)
def main():
 ext=pd.read_csv(E);ext_sc=set(ext.murcko_scaffold.fillna(""));ext_can=set(ext.canonical_smiles);m4=pd.read_csv(M/"molecules.csv");pos=m4[(m4.is_active==1)&(~m4.murcko_scaffold.isin(ext_sc))].copy();dec=m4[(m4.is_active==0)&(~m4.murcko_scaffold.isin(ext_sc))].copy()
 # Build a scaffold-excluded unlabeled drug-like background.
 u=[];seen=set()
 for s in CORP.read_text(encoding="utf-8").splitlines():
  mol=Chem.MolFromSmiles(s)
  if mol is None:continue
  can=Chem.MolToSmiles(mol,isomericSmiles=True)
  if can in seen or can in ext_can or scaf(mol) in ext_sc:continue
  seen.add(can);u.append(can)
 print(f"train positives={len(pos)} unlabeled={len(u)} DUD-E decoys={len(dec)} external={len(ext)}",flush=True)
 Xp,Dp,Fp,_=featurize(pos.canonical_smiles);Xu,Du,_,_=featurize(u);Xe,De,Fe,_=featurize(ext.canonical_smiles);Xd,Dd,_,_=featurize(dec.canonical_smiles)
 rng=np.random.default_rng(42);preds=[];predsD=[]
 for b in range(10):
  ix=rng.choice(len(u),min(5*len(pos),len(u)),replace=False);X=sparse.vstack([Xp,Xu[ix]]);y=np.r_[np.ones(len(pos),int),np.zeros(len(ix),int)];preds.append(model(100+b).fit(sparse.hstack([X,sparse.csr_matrix(np.vstack([Dp,Du[ix]]))]),y).predict_proba(sparse.hstack([Xe,sparse.csr_matrix(De)]))[:,1]);predsD.append(model(200+b).fit(np.vstack([Dp,Du[ix]]),y).predict_proba(De)[:,1])
 pu=np.mean(preds,axis=0);pud=np.mean(predsD,axis=0)
 ix=rng.choice(len(dec),min(5*len(pos),len(dec)),replace=False);Xs=sparse.vstack([Xp,Xd[ix]]);ys=np.r_[np.ones(len(pos),int),np.zeros(len(ix),int)];sup=model(42).fit(sparse.hstack([Xs,sparse.csr_matrix(np.vstack([Dp,Dd[ix]]))]),ys).predict_proba(sparse.hstack([Xe,sparse.csr_matrix(De)]))[:,1]
 # Label-free positive-neighborhood baseline.
 sim=[]
 for f in Fe:
  s=np.sort(DataStructs.BulkTanimotoSimilarity(f,Fp))[-5:];sim.append(float(np.mean(s)))
 methods={"PositiveSimilarity":np.asarray(sim),"Descriptor_PU":pud,"ECFP_PU":pu,"DUD_E_Supervised":sup};report={k:evals(ext,v) for k,v in methods.items()}
 audit={"train_positives_after_external_scaffold_exclusion":len(pos),"unlabeled_after_exclusion":len(u),"training_DUD_E_decoys_after_exclusion":len(dec),"external_n":len(ext),"external_PAM":int(ext.target.sum()),"external_inactive":int((ext.target==0).sum()),"rule":"All exact external molecules and all external Murcko scaffolds excluded from every training pool."}
 pd.DataFrame({"canonical_molecule_id":ext.canonical_molecule_id,"source_component":ext.source_component,"target":ext.target,**methods}).to_csv(OUT/"external_predictions.csv",index=False);(OUT/"metrics.json").write_text(json.dumps({"audit":audit,"methods":report},indent=2),encoding="utf-8");print(json.dumps({k:{"macro":v["macro_group_AUC"],"worst":v["worst_group_AUC"],"pooled":v["pooled_ROC_AUC"]} for k,v in report.items()},indent=2))
if __name__=="__main__":main()
