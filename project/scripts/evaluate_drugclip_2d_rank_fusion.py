#!/usr/bin/env python3
"""Exploratory parameter-free rank fusion of 2D and M4-pocket OOF scores."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--adapter",type=Path,required=True); ap.add_argument("--baseline",type=Path,required=True)
    ap.add_argument("--benchmark",type=Path,required=True); ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--column",default="pocket_rank"); ap.add_argument("--weight",type=float,default=0.5)
    ap.add_argument("--bootstrap",type=int,default=5000); ap.add_argument("--seed",type=int,default=20260924)
    a=ap.parse_args()
    adapter=pd.read_csv(a.adapter).groupby(["canonical_molecule_id","target"],as_index=False)[a.column].mean()
    base=pd.read_csv(a.baseline)[["canonical_molecule_id","target","prediction"]].rename(columns={"target":"base_target"})
    groups=pd.read_csv(a.benchmark)[["canonical_molecule_id","source_component"]]
    d=adapter.merge(base,on="canonical_molecule_id",validate="one_to_one").merge(groups,on="canonical_molecule_id",validate="one_to_one")
    if not np.array_equal(d.target,d.base_target): raise ValueError("target mismatch")
    y=d.target.to_numpy(int); n=len(d)
    dr=rankdata(d[a.column].to_numpy())/n; br=rankdata(d.prediction.to_numpy())/n
    fusion=a.weight*dr+(1-a.weight)*br
    base_auc=float(roc_auc_score(y,br)); fusion_auc=float(roc_auc_score(y,fusion))
    rng=np.random.default_rng(a.seed); molecule=[]; cluster=[]; sources=d.source_component.unique()
    for _ in range(a.bootstrap):
        idx=rng.integers(0,n,n)
        if len(np.unique(y[idx]))==2: molecule.append(roc_auc_score(y[idx],fusion[idx])-roc_auc_score(y[idx],br[idx]))
        chosen=rng.choice(sources,len(sources),True)
        ci=np.concatenate([np.flatnonzero(d.source_component.to_numpy()==s) for s in chosen])
        if len(np.unique(y[ci]))==2: cluster.append(roc_auc_score(y[ci],fusion[ci])-roc_auc_score(y[ci],br[ci]))
    def block(v):
        v=np.asarray(v); return {"mean":float(v.mean()),"ci95":list(map(float,np.quantile(v,[.025,.975]))),
                                 "probability_delta_gt_zero":float((v>0).mean()),"n":int(len(v))}
    report={"evidence_level":"posthoc_exploratory_fixed_rank_fusion","n":n,"weight":a.weight,
            "adapter_column":a.column,"baseline_auc":base_auc,"fusion_auc":fusion_auc,
            "observed_delta":fusion_auc-base_auc,"molecule_bootstrap_delta":block(molecule),
            "source_cluster_bootstrap_delta":block(cluster),
            "claim_boundary":"The fusion was explored after component results and requires prospective freezing."}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))
if __name__=="__main__": main()
