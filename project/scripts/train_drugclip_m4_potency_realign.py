#!/usr/bin/env python3
"""M4-specific DrugCLIP pocket re-alignment for functional pEC50 ranking."""
from __future__ import annotations
import argparse, json, random
from pathlib import Path
import numpy as np, pandas as pd, torch
from rdkit import Chem
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error
from torch.nn import functional as F
from train_drugclip_m4_pocket_realign import PocketRealign


def canonical(s):
    m=Chem.MolFromSmiles(str(s))
    if m is None: raise ValueError(s)
    return Chem.MolToSmiles(m,canonical=True)


def pair_indices(y, source, rng, max_pairs, delta):
    hi,lo=np.where((y[:,None]-y[None,:]>=delta) & (source[:,None]!=source[None,:]))
    if not len(hi): return np.empty((0,2),dtype=int)
    pick=rng.choice(len(hi),min(max_pairs,len(hi)),replace=False)
    return np.c_[hi[pick],lo[pick]]


def fit_predict(xtr,pockets,ytr,strain,xte,mode,seed,args):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    x=torch.tensor(xtr,dtype=torch.float32); p=torch.tensor(pockets,dtype=torch.float32)
    xv=torch.tensor(xte,dtype=torch.float32); s=torch.tensor(strain,dtype=torch.long)
    mean,std=float(ytr.mean()),float(ytr.std()); yn=(ytr-mean)/max(std,1e-6)
    y=torch.tensor(yn,dtype=torch.float32)
    model=PocketRealign(x.shape[1],args.rank,int(strain.max())+1)
    opt=torch.optim.AdamW(model.parameters(),lr=args.lr,weight_decay=args.weight_decay)
    ox,op=F.normalize(x,dim=1),F.normalize(p,dim=1); rng=np.random.default_rng(seed)
    use_rank=mode in {'reg_rank','reg_rank_domain'}; use_domain=mode=='reg_rank_domain'
    for _ in range(args.epochs):
        model.train(); opt.zero_grad()
        zm,zp,_,_,pred,source_logit=model(x,p,args.domain_weight if use_domain else 0.0)
        loss=F.smooth_l1_loss(pred,y,beta=0.5)
        if use_rank:
            pairs=pair_indices(ytr,strain,rng,args.max_pairs,args.pair_delta)
            if len(pairs):
                hi,lo=map(torch.tensor,pairs.T)
                loss=loss+args.rank_weight*F.relu(args.margin-pred[hi]+pred[lo]).mean()
        preserve=(1-(zm*ox).sum(1)).mean()+(1-(zp*op).sum(1)).mean()
        loss=loss+args.preserve_weight*preserve
        if use_domain: loss=loss+args.domain_loss_weight*F.cross_entropy(source_logit,s)
        loss.backward(); opt.step()
    model.eval()
    with torch.no_grad(): prediction=model(xv,p)[4].numpy()*std+mean
    return prediction


def metric(y,p):
    return {'spearman':float(spearmanr(y,p).statistic),'mae':float(mean_absolute_error(y,p)),
            'rmse':float(mean_squared_error(y,p)**0.5),'n':int(len(y))}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--embeddings',type=Path,required=True)
    ap.add_argument('--benchmark',type=Path,required=True); ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--protocol',choices=['source','scaffold'],required=True)
    ap.add_argument('--seeds',default='20260924,20260925,20260926'); ap.add_argument('--epochs',type=int,default=300)
    ap.add_argument('--rank',type=int,default=16); ap.add_argument('--lr',type=float,default=1e-3)
    ap.add_argument('--weight-decay',type=float,default=1e-3); ap.add_argument('--pair-delta',type=float,default=0.5)
    ap.add_argument('--margin',type=float,default=0.2); ap.add_argument('--rank-weight',type=float,default=0.5)
    ap.add_argument('--preserve-weight',type=float,default=0.2); ap.add_argument('--domain-weight',type=float,default=0.2)
    ap.add_argument('--domain-loss-weight',type=float,default=0.2); ap.add_argument('--max-pairs',type=int,default=1024)
    args=ap.parse_args(); torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(args.embeddings,allow_pickle=False); ids=list(map(str,arc['molecule_ids'])); lookup={v:i for i,v in enumerate(ids)}
    d=pd.read_csv(args.benchmark); cs=d.canonical_smiles.map(canonical)
    if missing:=[v for v in cs if v not in lookup]: raise ValueError(missing[:5])
    x=arc['molecule_embeddings'].astype(np.float32)[[lookup[v] for v in cs]]; pockets=arc['pocket_embeddings'].astype(np.float32)
    y=d.pEC50.to_numpy(float); sources=sorted(d.source_component.unique()); smap={v:i for i,v in enumerate(sources)}
    source=d.source_component.map(smap).to_numpy(int); folds=(d.source_fold if args.protocol=='source' else d.scaffold_fold).to_numpy(int)
    seeds=[int(v) for v in args.seeds.split(',')]; modes=['regression','reg_rank','reg_rank_domain']; records=[]; reports=[]
    for seed in seeds:
        pred={m:np.full(len(y),np.nan) for m in modes}
        for fold in sorted(np.unique(folds)):
            tr,te=folds!=fold,folds==fold; unique=sorted(np.unique(source[tr])); dense={v:i for i,v in enumerate(unique)}
            ds=np.asarray([dense[v] for v in source[tr]])
            for mode in modes: pred[mode][te]=fit_predict(x[tr],pockets,y[tr],ds,x[te],mode,seed+100*int(fold),args)
        reports.append({'seed':seed,'metrics':{m:metric(y,pred[m]) for m in modes}})
        for i in range(len(y)): records.append({'canonical_molecule_id':d.iloc[i].canonical_molecule_id,'pEC50':y[i],
            'fold':int(folds[i]),'seed':seed,**{m:float(pred[m][i]) for m in modes}})
    agg={}
    for m in modes:
        vals=[r['metrics'][m]['spearman'] for r in reports]
        agg[m]={'spearman_mean':float(np.mean(vals)),'spearman_min':float(np.min(vals)),'spearman_max':float(np.max(vals)),
                'per_seed':[r['metrics'][m] for r in reports]}
    out={'evidence_level':'retrospective_m4_specific_potency_realign','protocol':args.protocol,'n':len(y),
         'aggregate':agg,'hyperparameters':{k:v for k,v in vars(args).items() if k not in {'embeddings','benchmark','output'}},
         'claim_boundary':'Functional potency ranking; not prospective PAM confirmation.'}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(out,indent=2),encoding='utf-8')
    pd.DataFrame(records).to_csv(args.output.with_suffix('.predictions.csv'),index=False); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
