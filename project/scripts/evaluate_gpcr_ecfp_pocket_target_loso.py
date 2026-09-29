#!/usr/bin/env python3
"""Asymmetric ECFP3-to-frozen-DrugCLIP-pocket transfer across GPCR targets."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F

from evaluate_muscarinic_ecfp_pocket_baseline import fps
from evaluate_gpcr_drugclip_target_loso import TARGETS, metrics, scaffold_bootstrap_delta


class Adapter(nn.Module):
    def __init__(self,seed):
        super().__init__();torch.manual_seed(seed);self.project=nn.Linear(2048,128,bias=False);nn.init.normal_(self.project.weight,std=.01)
    def forward(self,x):return F.normalize(self.project(x),dim=-1)


def fit(x,pocket,labels,seen,seed,args,random_labels=False):
    model=Adapter(seed);local={v:i for i,v in enumerate(seen)};y=torch.as_tensor([local[int(v)] for v in labels]);
    if random_labels:
        g=torch.Generator().manual_seed(seed);y=y[torch.randperm(len(y),generator=g)]
    count=torch.bincount(y,minlength=len(seen)).float();weight=count.sum()/(len(seen)*count.clamp_min(1));seen_t=torch.as_tensor(seen);opt=torch.optim.AdamW(model.parameters(),lr=args.lr,weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        opt.zero_grad();score=(model(x)@pocket.T)[:,seen_t];loss=F.cross_entropy(score/args.temperature,y,weight=weight);loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);opt.step()
    return model


def main():
    p=argparse.ArgumentParser();p.add_argument('--representations',type=Path,required=True);p.add_argument('--benchmark',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seeds',default='20260925,20260926,20260927');p.add_argument('--epochs',type=int,default=100);p.add_argument('--lr',type=float,default=1e-3);p.add_argument('--weight-decay',type=float,default=1e-2);p.add_argument('--temperature',type=float,default=.07);p.add_argument('--bootstrap',type=int,default=1000);p.add_argument('--with-random-control',action='store_true');a=p.parse_args();torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(a.representations,allow_pickle=False);table0=pd.read_csv(a.benchmark);idx={v:i for i,v in enumerate(map(str,arc['molecule_ids']))};table=table0[table0.canonical_smiles.astype(str).isin(idx)].copy().reset_index(drop=True);order=np.asarray([idx[v] for v in table.canonical_smiles.astype(str)]);x=fps(table.canonical_smiles.astype(str).tolist(),radius=3);pocket=F.normalize(torch.as_tensor(arc['pocket_embeddings'].astype(np.float32)),dim=-1);labels=np.asarray([TARGETS.index(v) for v in table.target]);official=arc['scores'].astype(np.float32).T[order]
    seeds=list(map(int,a.seeds.split(',')));pred=[np.full_like(official,np.nan) for _ in seeds];rnd=[np.full_like(official,np.nan) for _ in seeds];audit={};scaffolds=table.murcko_scaffold.astype(str).to_numpy()
    for held in range(len(TARGETS)):
        test=labels==held;held_scaf=set(scaffolds[test]);train=(labels!=held)&~np.asarray([s in held_scaf for s in scaffolds]);seen=[j for j in range(len(TARGETS)) if j!=held];audit[TARGETS[held]]={'test_n':int(test.sum()),'train_n':int(train.sum()),'held_scaffold_overlap':0}
        for k,seed in enumerate(seeds):
            model=fit(x[train],pocket,labels[train],seen,seed+1009*held,a,False)
            with torch.inference_mode():pred[k][test]=(model(x[test])@pocket.T).numpy()
            if a.with_random_control:
                control=fit(x[train],pocket,labels[train],seen,seed+1009*held,a,True)
                with torch.inference_mode():rnd[k][test]=(control(x[test])@pocket.T).numpy()
    tuned=np.mean(pred,axis=0);report={'protocol':'leave-one-entire-GPCR-target-out; held-target scaffolds purged; held pocket absent from training loss','training_audit':audit,'metrics':{'official_drugclip':metrics(labels,official),'ecfp3_to_frozen_pocket':metrics(labels,tuned)},'ecfp3_minus_official_macro_recall1_scaffold_bootstrap_95ci':scaffold_bootstrap_delta(table,labels,tuned,official,a.bootstrap,20260926),'hyperparameters':{'epochs':a.epochs,'lr':a.lr,'temperature':a.temperature,'seeds':seeds},'claim_boundary':'Four-target retrospective GPCR target retrieval; not binding, potency, PAM efficacy, or broad GPCR SOTA.'}
    if a.with_random_control:
        random=np.mean(rnd,axis=0);report['metrics']['random_label_control']=metrics(labels,random)
    out=table[['canonical_molecule_id','target','murcko_scaffold']].copy()
    for name,val in [('official',official),('ecfp3',tuned)]:
        for j,target in enumerate(TARGETS):out[f'{name}_{target}']=val[:,j]
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');out.to_csv(a.output.with_suffix('.csv'),index=False);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
