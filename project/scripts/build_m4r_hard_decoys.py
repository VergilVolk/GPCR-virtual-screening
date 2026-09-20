# -*- coding: utf-8 -*-
"""Create fold-pure, descriptor-matched hard decoys for M4R enrichment."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.neighbors import NearestNeighbors
P=Path(__file__).resolve().parents[1];B=P/"data"/"benchmarks"/"m4r_enrichment_v1";K=10
def main():
 d=pd.read_csv(B/"molecules.csv");D=np.load(B/"feature_cache"/"descriptors.npy");names=["MW","logP","HBD","HBA","TPSA","RotB","AromRings","FracCSP3","Heavy"]
 chosen=[];rng=np.random.default_rng(42)
 for f in sorted(d.scaffold_fold.unique()):
  ai=np.where((d.scaffold_fold==f)&(d.is_active==1))[0];di=np.where((d.scaffold_fold==f)&(d.is_active==0))[0];mean=D[di].mean(0);sd=D[di].std(0)+1e-8;importance=np.asarray([1,1,1.5,1,1,1.5,2,3,1],float);A=(D[ai]-mean)/sd*importance;Z=(D[di]-mean)/sd*importance
  nn=NearestNeighbors(n_neighbors=min(400,len(di)),algorithm="auto",n_jobs=-1).fit(Z);_,ix=nn.kneighbors(A);used=set();order=rng.permutation(len(ai))
  for q in order:
   pick=[]
   for z in ix[q]:
    idx=int(di[z])
    if idx not in used:pick.append(idx);used.add(idx)
    if len(pick)==K:break
   if len(pick)<K:
    rest=[int(x) for x in rng.permutation(di) if int(x) not in used][:K-len(pick)];pick.extend(rest);used.update(rest)
   chosen.extend(pick)
  print(f"fold={f} actives={len(ai)} hard_decoys={len(used)}",flush=True)
 keep=np.sort(np.r_[np.where(d.is_active==1)[0],np.asarray(chosen,int)]);h=d.iloc[keep].copy();h.to_csv(B/"hard_molecules.csv",index=False)
 # Standardized mean differences after matching; values near zero indicate property balance.
 smd={}
 DD=D[keep];yy=h.is_active.to_numpy(int)
 for i,n in enumerate(names):smd[n]=float((DD[yy==1,i].mean()-DD[yy==0,i].mean())/np.sqrt((DD[yy==1,i].var()+DD[yy==0,i].var())/2))
 audit={"actives":int(h.is_active.sum()),"hard_decoys":int((h.is_active==0).sum()),"decoys_per_active":K,"folds":h.groupby("scaffold_fold").is_active.agg(n="size",actives="sum").reset_index().to_dict("records"),"standardized_mean_differences":smd,"max_abs_SMD":float(max(abs(x) for x in smd.values())),"rules":f"Each active receives {K} unique nearest decoys in weighted nine-descriptor space from the same scaffold fold; decoy scaffolds therefore remain fold-pure."}
 (B/"hard_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8");print(json.dumps(audit,indent=2))
if __name__=="__main__":main()
