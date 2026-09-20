# -*- coding: utf-8 -*-
"""Construct fold-pure adversarial decoys that resist property and ECFP shortcuts."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
from lightgbm import LGBMClassifier
P=Path(__file__).resolve().parents[1];B=P/"data"/"benchmarks"/"m4r_enrichment_v1";K=10
def main():
 d=pd.read_csv(B/"molecules.csv");X=sparse.load_npz(B/"feature_cache"/"ecfp4_2048.npz");D=np.load(B/"feature_cache"/"descriptors.npy");y=d.is_active.to_numpy(int);fold=d.scaffold_fold.to_numpy(int);rng=np.random.default_rng(42);chosen=[];assign=[]
 for f in sorted(np.unique(fold)):
  tr=fold!=f;ai=np.where((fold==f)&(y==1))[0];di=np.where((fold==f)&(y==0))[0]
  prop=LGBMClassifier(n_estimators=300,learning_rate=.04,num_leaves=31,max_depth=6,min_child_samples=50,class_weight="balanced",random_state=42,n_jobs=6,verbosity=-1).fit(D[tr],y[tr]).predict_proba(D[di])[:,1]
  # Retain the 15k most active-like decoys by a model that never saw this fold.
  pool=di[np.argsort(prop)[::-1][:min(15000,len(di))]];inter=(X[ai]@X[pool].T).toarray().astype(np.float32);na=np.asarray(X[ai].sum(1)).ravel()[:,None];nd=np.asarray(X[pool].sum(1)).ravel()[None,:];sim=inter/(na+nd-inter+1e-8)
  used=set();order=rng.permutation(len(ai))
  for q in order:
   picks=[]
   for z in np.argsort(sim[q])[::-1]:
    idx=int(pool[z])
    if idx not in used:picks.append((idx,float(sim[q,z])));used.add(idx)
    if len(picks)==K:break
   for idx,s in picks:chosen.append(idx);assign.append({"active_molecule_id":d.iloc[ai[q]].molecule_id,"decoy_molecule_id":d.iloc[idx].molecule_id,"scaffold_fold":int(f),"ecfp4_tanimoto":s})
  print(f"fold={f} actives={len(ai)} adversarial_decoys={len(used)} mean_pair_sim={np.mean([x['ecfp4_tanimoto'] for x in assign if x['scaffold_fold']==f]):.3f}",flush=True)
 keep=np.sort(np.r_[np.where(y==1)[0],np.asarray(chosen,int)]);h=d.iloc[keep].copy();h.to_csv(B/"adversarial_molecules.csv",index=False);pd.DataFrame(assign).to_csv(B/"adversarial_assignments.csv",index=False)
 names=["MW","logP","HBD","HBA","TPSA","RotB","AromRings","FracCSP3","Heavy"];DD=D[keep];yy=h.is_active.to_numpy(int);smd={n:float((DD[yy==1,i].mean()-DD[yy==0,i].mean())/np.sqrt((DD[yy==1,i].var()+DD[yy==0,i].var())/2)) for i,n in enumerate(names)}
 audit={"actives":int(h.is_active.sum()),"adversarial_decoys":int((h.is_active==0).sum()),"decoys_per_active":K,"max_abs_SMD":float(max(abs(x) for x in smd.values())),"standardized_mean_differences":smd,"mean_assigned_ECFP4_similarity":float(pd.DataFrame(assign).ecfp4_tanimoto.mean()),"folds":h.groupby("scaffold_fold").is_active.agg(n="size",actives="sum").reset_index().to_dict("records"),"rules":"For each outer fold, a descriptor classifier trained only on other folds identifies active-like decoys; each active is then assigned ten unique nearest ECFP4 decoys from that fold. All scaffolds remain fold-pure."};(B/"adversarial_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8");print(json.dumps(audit,indent=2))
if __name__=="__main__":main()
