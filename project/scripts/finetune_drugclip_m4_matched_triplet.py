#!/usr/bin/env python3
"""M4 functional fine-tuning with source- and chemistry-matched triplets.

The GPCR-specialized pocket projection is frozen. For every experimental
inactive in each training fold, functional PAM positives are selected from the
same source component by ECFP similarity. A same-source random-pair arm is the
sampling control.
"""
from __future__ import annotations
import argparse, copy, json, random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import average_precision_score, roc_auc_score
from torch import nn
from torch.nn import functional as F

class Proj(nn.Module):
    def __init__(self,s):
        super().__init__(); self.linear1=nn.Linear(s['linear1.weight'].shape[1],s['linear1.weight'].shape[0]); self.linear2=nn.Linear(s['linear2.weight'].shape[1],s['linear2.weight'].shape[0]); self.load_state_dict(s)
    def forward(self,x): return F.normalize(self.linear2(F.relu(self.linear1(x))),dim=-1)
class Model(nn.Module):
    def __init__(self,s): super().__init__(); self.mol_project=Proj(s['mol_project']); self.pocket_project=Proj(s['pocket_project'])
    def retrieval(self,m,p):
        zm,zp=self.mol_project(m),self.pocket_project(p); states=zm@zp.T; return zm,zp,.1*torch.logsumexp(states/.1,dim=1)

def make_pairs(frame,mode,k,seed):
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048); fps=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in frame.canonical_smiles]; rng=np.random.default_rng(seed); pairs=[]
    for ni in np.flatnonzero(frame.target.to_numpy()==0):
        candidates=np.flatnonzero((frame.target.to_numpy()==1)&(frame.source_component.to_numpy()==frame.iloc[ni].source_component))
        if not len(candidates): candidates=np.flatnonzero(frame.target.to_numpy()==1)
        if mode=='matched':
            sim=np.asarray(DataStructs.BulkTanimotoSimilarity(fps[ni],[fps[j] for j in candidates])); chosen=candidates[np.argsort(-sim)[:min(k,len(candidates))]]
        else: chosen=rng.choice(candidates,min(k,len(candidates)),replace=False)
        pairs.extend((int(pi),int(ni)) for pi in chosen)
    return np.asarray(pairs,dtype=np.int64)

def fit(x,p,frame,initial,mode,seed,args):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed); x=torch.as_tensor(x,dtype=torch.float32); p=torch.as_tensor(p,dtype=torch.float32); model=Model(copy.deepcopy(initial)); frozen=Model(copy.deepcopy(initial)).eval()
    for q in model.pocket_project.parameters(): q.requires_grad_(False)
    for q in frozen.parameters(): q.requires_grad_(False)
    with torch.no_grad(): z0m,z0p,_=frozen.retrieval(x,p)
    pairs=make_pairs(frame,mode,args.positive_k,seed); pair_t=torch.as_tensor(pairs,dtype=torch.long); opt=torch.optim.AdamW(model.mol_project.parameters(),lr=args.lr,weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        opt.zero_grad(); zm,zp,r=model.retrieval(x,p); delta=r[pair_t[:,0]]-r[pair_t[:,1]]; rank=F.relu(args.margin-delta).mean(); preserve=(1-(zm*z0m).sum(1)).mean(); loss=rank+args.preserve_weight*preserve; loss.backward(); torch.nn.utils.clip_grad_norm_(model.mol_project.parameters(),1.); opt.step()
    return model,pairs

def metric(y,s): return {'roc_auc':float(roc_auc_score(y,s)),'average_precision':float(average_precision_score(y,s))}
def bootstrap(y,a,b,n=3000,seed=20260924):
    rng=np.random.default_rng(seed); v=[]
    for _ in range(n):
        i=rng.integers(0,len(y),len(y))
        if np.unique(y[i]).size==2: v.append(roc_auc_score(y[i],a[i])-roc_auc_score(y[i],b[i]))
    return list(map(float,np.quantile(v,[.025,.5,.975])))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--representations',type=Path,required=True); ap.add_argument('--projection',type=Path,required=True); ap.add_argument('--benchmark',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--protocol',choices=['series','source'],default='series'); ap.add_argument('--seeds',default='20260924,20260925,20260926'); ap.add_argument('--epochs',type=int,default=300); ap.add_argument('--lr',type=float,default=1e-4); ap.add_argument('--weight-decay',type=float,default=1e-3); ap.add_argument('--margin',type=float,default=.05); ap.add_argument('--preserve-weight',type=float,default=.5); ap.add_argument('--positive-k',type=int,default=3); ap.add_argument('--bootstrap',type=int,default=3000); args=ap.parse_args(); torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(args.representations,allow_pickle=False); initial=torch.load(args.projection,map_location='cpu'); table=pd.read_csv(args.benchmark); lookup={str(v):i for i,v in enumerate(arc['molecule_ids'])}; order=np.asarray([lookup[str(v)] for v in table.canonical_smiles]); x=arc['molecule_representations'].astype(np.float32)[order]; p=arc['pocket_representations'].astype(np.float32); y=table.target.to_numpy(int); fold_column='series_holdout_fold' if args.protocol=='series' else 'source_fold'; folds=table[fold_column].to_numpy(int); assigned=folds>=0; modes=['frozen','random','matched']; rows=[]
    frozen=Model(copy.deepcopy(initial)).eval()
    for seed in map(int,args.seeds.split(',')):
        pred={m:np.full(len(y),np.nan) for m in modes}
        for fold in sorted(np.unique(folds[assigned])):
            tr=assigned&(folds!=fold); te=assigned&(folds==fold); train=table[tr].reset_index(drop=True); xt=x[tr]
            with torch.no_grad(): pred['frozen'][te]=frozen.retrieval(torch.as_tensor(x[te]),torch.as_tensor(p))[2].numpy()
            for mode in ['random','matched']:
                model,_=fit(xt,p,train,initial,mode,seed+101*int(fold),args)
                with torch.no_grad(): pred[mode][te]=model.retrieval(torch.as_tensor(x[te]),torch.as_tensor(p))[2].numpy()
        for i in np.flatnonzero(assigned): rows.append({'canonical_molecule_id':table.iloc[i].canonical_molecule_id,'target':int(y[i]),'fold':int(folds[i]),'seed':seed,**{m:float(pred[m][i]) for m in modes}})
    frame=pd.DataFrame(rows); ens=frame.groupby('canonical_molecule_id',sort=False)[modes].mean(); truth=table.set_index('canonical_molecule_id').loc[ens.index,'target'].to_numpy(int)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    # Fit deployable models only after the OOF protocol and hyperparameters are
    # frozen.  Each checkpoint uses every molecule eligible for this protocol;
    # external datasets must never be used to select among these seeds.
    full_table=table.loc[assigned].reset_index(drop=True); full_x=x[assigned]; checkpoints=[]
    for seed in map(int,args.seeds.split(',')):
        model,pairs=fit(full_x,p,full_table,initial,'matched',seed,args)
        checkpoint_path=args.output.with_name(f'{args.output.stem}.seed{seed}.projection.pt')
        torch.save({
            'mol_project':model.mol_project.state_dict(),
            'pocket_project':model.pocket_project.state_dict(),
            'logit_scale':initial.get('logit_scale'),
            'source_projection':str(args.projection),
            'protocol':args.protocol,
            'training_mode':'source_chemotype_matched_functional_triplet',
            'seed':seed,
            'n_training_molecules':int(len(full_table)),
            'n_training_pairs':int(len(pairs)),
        },checkpoint_path)
        checkpoints.append(str(checkpoint_path))
    report={'method':'GPCR-initialized DrugCLIP M4 source-chemotype matched functional triplet','protocol':args.protocol+' holdout','n':len(truth),'metrics':{m:metric(truth,ens[m].to_numpy()) for m in modes},'matched_minus_frozen_auc_bootstrap_95ci':bootstrap(truth,ens.matched.to_numpy(),ens.frozen.to_numpy(),args.bootstrap),'matched_minus_random_auc_bootstrap_95ci':bootstrap(truth,ens.matched.to_numpy(),ens.random.to_numpy(),args.bootstrap),'deployable_checkpoints':checkpoints,'checkpoint_selection':'Unweighted three-seed ensemble fixed before external evaluation.','claim_boundary':'Retrospective functional PAM-vs-inactive ranking; not wet-lab confirmation.'}
    args.output.write_text(json.dumps(report,indent=2),encoding='utf-8'); frame.to_csv(args.output.with_suffix('.predictions.csv'),index=False); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
