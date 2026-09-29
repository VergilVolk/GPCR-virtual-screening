#!/usr/bin/env python3
"""Frozen nine-target LIT-PCBA audit of official and GPCR-tuned DrugCLIP."""
from __future__ import annotations
import argparse, copy, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from finetune_drugclip_gpcr_screening import binary_metrics
from finetune_drugclip_muscarinic_triplet import Proj


def score(arc,state):
    mol=Proj(copy.deepcopy(state['mol_project'])).eval();pocket=Proj(copy.deepcopy(state['pocket_project'])).eval()
    with torch.inference_mode():return (mol(torch.as_tensor(arc['molecule_representations'].astype(np.float32)))@pocket(torch.as_tensor(arc['pocket_representations'].astype(np.float32))).T).numpy()


def metrics(table,scores,targets):
    per={}
    for target in targets:
        keep=table.target.eq(target).to_numpy();per[target]=binary_metrics(table.loc[keep,'label'].to_numpy(int),scores[keep])
    names=['roc_auc','pr_auc','bedroc_alpha20','ef1pct','ef5pct']
    return {'macro':{name:float(np.mean([per[t][name] for t in targets])) for name in names},'per_target':per}


def main():
    p=argparse.ArgumentParser();p.add_argument('--representations',type=Path,required=True);p.add_argument('--pairs',type=Path,required=True);p.add_argument('--official-projection',type=Path,required=True);p.add_argument('--checkpoint',action='append',nargs=2,metavar=('NAME','PATH'),required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--bootstrap',type=int,default=1000);a=p.parse_args()
    arc=np.load(a.representations,allow_pickle=False);ids=list(map(str,arc['molecule_ids']));mi={v:i for i,v in enumerate(ids)};targets=list(map(str,arc['pocket_ids']));ti={v:i for i,v in enumerate(targets)};table=pd.read_csv(a.pairs);table=table[table.canonical_smiles.astype(str).isin(mi)&table.target.isin(ti)].copy().reset_index(drop=True);pm=np.asarray([mi[v] for v in table.canonical_smiles.astype(str)]);pt=np.asarray([ti[v] for v in table.target]);table['murcko_scaffold']=[MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(v)) or v for v in table.canonical_smiles]
    states={'official':torch.load(a.official_projection,map_location='cpu')};states.update({name:torch.load(path,map_location='cpu') for name,path in a.checkpoint});pred={name:matrix[pm,pt] for name,state in states.items() for matrix in [score(arc,state)]};reports={name:metrics(table,value,targets) for name,value in pred.items()}
    groups={};
    for (target,scaffold),idx in table.groupby(['target','murcko_scaffold']).groups.items():groups[(target,scaffold)]=np.asarray(list(idx),dtype=int)
    rng=np.random.default_rng(20260926);ci={name:{metric:[] for metric in reports[name]['macro']} for name in pred if name!='official'}
    by_target={target:[key for key in groups if key[0]==target] for target in targets}
    for _ in range(a.bootstrap):
        sample=[]
        for target in targets:
            keys=by_target[target];chosen=rng.integers(0,len(keys),len(keys));sample.extend(np.concatenate([groups[keys[i]] for i in chosen]))
        idx=np.asarray(sample,dtype=int);frame=table.iloc[idx].reset_index(drop=True);base=metrics(frame,pred['official'][idx],targets)['macro']
        for name in ci:
            tuned=metrics(frame,pred[name][idx],targets)['macro']
            for metric in ci[name]:ci[name][metric].append(tuned[metric]-base[metric])
    ci={name:{metric:list(map(float,np.quantile(values,[.025,.5,.975]))) for metric,values in block.items()} for name,block in ci.items()}
    report={'protocol':'Frozen nine-target LIT-PCBA subset; no target label used for tuning or checkpoint selection. Scaffold-cluster bootstrap.','n_pairs':len(table),'n_molecules':len(np.unique(pm)),'targets':targets,'metrics':reports,'delta_vs_official_bootstrap_95ci':ci,'claim_boundary':'External non-GPCR screening preservation/transfer audit; no affinity, efficacy, or PAM claim.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');out=table[['pair_id','target','label','murcko_scaffold']].copy();[out.__setitem__(name,value) for name,value in pred.items()];out.to_csv(a.output.with_suffix('.predictions.csv'),index=False);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
