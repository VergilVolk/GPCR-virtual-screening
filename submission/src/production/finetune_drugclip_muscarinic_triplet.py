#!/usr/bin/env python3
"""True DrugCLIP projection fine-tuning on explicit muscarinic triplets.

Triplet: (active-subtype pocket, same molecule, inactive-subtype pocket).
Because the molecule is identical on both sides, molecular-size and chemotype
effects cancel within every training example. Randomly oriented triplets are an
equal-size negative control.
"""
from __future__ import annotations
import argparse, copy, json, random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from torch import nn
from torch.nn import functional as F

def scaffold(s):
    m=Chem.MolFromSmiles(s); return MurckoScaffold.MurckoScaffoldSmiles(mol=m) or Chem.MolToSmiles(m,canonical=True)
def folds_for(smiles,n=5):
    groups={}
    for s in smiles: groups.setdefault(scaffold(s),[]).append(s)
    bins=[[] for _ in range(n)]; sizes=[0]*n
    for _,members in sorted(groups.items(),key=lambda x:(-len(x[1]),x[0])):
        j=min(range(n),key=lambda k:(sizes[k],k)); bins[j]+=members; sizes[j]+=len(members)
    return {s:i for i,b in enumerate(bins) for s in b}

class Proj(nn.Module):
    def __init__(self,s):
        super().__init__(); self.linear1=nn.Linear(s['linear1.weight'].shape[1],s['linear1.weight'].shape[0]); self.linear2=nn.Linear(s['linear2.weight'].shape[1],s['linear2.weight'].shape[0]); self.load_state_dict(s)
    def forward(self,x): return F.normalize(self.linear2(F.relu(self.linear1(x))),dim=-1)
class Model(nn.Module):
    def __init__(self,s): super().__init__(); self.mol_project=Proj(s['mol_project']); self.pocket_project=Proj(s['pocket_project'])
    def matrix(self,m,p): return self.mol_project(m)@self.pocket_project(p).T

def fit(m,p,d,mi,pi,initial,seed,args,random_orientation=False):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed); model=Model(copy.deepcopy(initial)); frozen=Model(copy.deepcopy(initial)).eval()
    for q in frozen.parameters(): q.requires_grad_(False)
    with torch.no_grad(): z0m=frozen.mol_project(m); z0p=frozen.pocket_project(p)
    mol=torch.tensor([mi[s] for s in d.canonical_smiles]); pos=torch.tensor([pi[s] for s in d.positive_subtype]); neg=torch.tensor([pi[s] for s in d.negative_subtype])
    if random_orientation:
        swap=torch.rand(len(d))<.5; original=pos.clone(); pos[swap]=neg[swap]; neg[swap]=original[swap]
    if args.train_scope == 'pocket':
        for parameter in model.mol_project.parameters(): parameter.requires_grad_(False)
        trainable=model.pocket_project.parameters()
    else:
        trainable=model.parameters()
    opt=torch.optim.AdamW(trainable,lr=args.lr,weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        opt.zero_grad(); score=model.matrix(m,p); delta=score[mol,pos]-score[mol,neg]; loss=F.relu(args.margin-delta).mean(); zm,zp=model.mol_project(m),model.pocket_project(p); loss=loss+args.preserve_weight*((1-(zm*z0m).sum(1)).mean()+(1-(zp*z0p).sum(1)).mean()); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),1.); opt.step()
    return model

def clustered_bootstrap(d,base,tuned,n=5000,seed=20260924):
    rng=np.random.default_rng(seed); molecules=np.asarray(sorted(d.canonical_smiles.unique())); vals=[]
    for _ in range(n):
        sample=rng.choice(molecules,len(molecules),replace=True); idx=np.concatenate([np.flatnonzero(d.canonical_smiles.to_numpy()==s) for s in sample]); vals.append((tuned[idx]>0).mean()-(base[idx]>0).mean())
    return list(map(float,np.quantile(vals,[.025,.5,.975])))

def main():
    a=argparse.ArgumentParser(); a.add_argument('--representations',type=Path,required=True); a.add_argument('--projection',type=Path,required=True); a.add_argument('--contrasts',type=Path,required=True); a.add_argument('--output',type=Path,required=True); a.add_argument('--seeds',default='20260924,20260925,20260926'); a.add_argument('--epochs',type=int,default=300); a.add_argument('--lr',type=float,default=5e-5); a.add_argument('--weight-decay',type=float,default=1e-3); a.add_argument('--margin',type=float,default=.1); a.add_argument('--preserve-weight',type=float,default=.5); a.add_argument('--train-scope',choices=['both','pocket'],default='both'); args=a.parse_args(); torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(args.representations,allow_pickle=False); initial=torch.load(args.projection,map_location='cpu'); mids=list(map(str,arc['molecule_ids'])); pids=list(map(str,arc['pocket_ids'])); mi={v:i for i,v in enumerate(mids)}; pi={v.split('_')[0]:i for i,v in enumerate(pids)}
    d=pd.read_csv(args.contrasts); d=d[~d.overlaps_pacer_m4.astype(bool)].copy().reset_index(drop=True); fmap=folds_for(d.canonical_smiles.unique()); d['fold']=d.canonical_smiles.map(fmap); m=torch.tensor(arc['molecule_representations'],dtype=torch.float32); p=torch.tensor(arc['pocket_representations'],dtype=torch.float32)
    frozen=Model(copy.deepcopy(initial)).eval(); s0=frozen.matrix(m,p).detach().numpy(); base=np.asarray([s0[mi[r.canonical_smiles],pi[r.positive_subtype]]-s0[mi[r.canonical_smiles],pi[r.negative_subtype]] for r in d.itertuples()]); rows=[]; reports=[]
    for seed in map(int,args.seeds.split(',')):
        out={mode:np.full(len(d),np.nan) for mode in ['random_control','targeted']}
        for fold in sorted(d.fold.unique()):
            tr=d[d.fold!=fold]; te=np.flatnonzero(d.fold.to_numpy()==fold)
            for mode in out:
                model=fit(m,p,tr,mi,pi,initial,seed+101*int(fold),args,mode=='random_control'); score=model.matrix(m,p).detach().numpy()
                for idx in te:
                    r=d.iloc[idx]; out[mode][idx]=score[mi[r.canonical_smiles],pi[r.positive_subtype]]-score[mi[r.canonical_smiles],pi[r.negative_subtype]]
        reports.append({'seed':seed,'random_control_accuracy':float((out['random_control']>0).mean()),'targeted_accuracy':float((out['targeted']>0).mean())})
        for i,r in d.iterrows(): rows.append({'canonical_smiles':r.canonical_smiles,'positive_subtype':r.positive_subtype,'negative_subtype':r.negative_subtype,'fold':int(r.fold),'seed':seed,'frozen_delta':float(base[i]),'random_control_delta':float(out['random_control'][i]),'targeted_delta':float(out['targeted'][i])})
    frame=pd.DataFrame(rows); ens=frame.groupby(['canonical_smiles','positive_subtype','negative_subtype'],sort=False)[['frozen_delta','random_control_delta','targeted_delta']].mean().reset_index(); tuned=ens.targeted_delta.to_numpy(); control=ens.random_control_delta.to_numpy()
    report={'method':'DrugCLIP internal projection same-molecule muscarinic triplet','train_scope':args.train_scope,'n_pairs':len(d),'n_molecules':d.canonical_smiles.nunique(),'split':'5-fold molecule Murcko scaffold OOF','frozen_accuracy':float((base>0).mean()),'random_control_accuracy':float((control>0).mean()),'targeted_accuracy':float((tuned>0).mean()),'targeted_minus_frozen_molecule_cluster_bootstrap_95ci':clustered_bootstrap(d,base,tuned),'targeted_minus_random_molecule_cluster_bootstrap_95ci':clustered_bootstrap(d,control,tuned),'seed_reports':reports,'claim_boundary':'Small heterogeneous explicit M1-M5 selectivity set; GPCR-family pilot, not M4 PAM efficacy.'}
    full=fit(m,p,d,mi,pi,initial,20260924,args,False); ck={'mol_project':full.mol_project.state_dict(),'pocket_project':full.pocket_project.state_dict(),'logit_scale':initial['logit_scale'],'source_checkpoint':initial.get('source_checkpoint'),'metadata':{'stage':'M1-M5 explicit selectivity triplet','n_pairs':len(d)}}
    args.output.parent.mkdir(parents=True,exist_ok=True); ckpath=args.output.with_suffix('.projection.pt'); torch.save(ck,ckpath); report['checkpoint']=str(ckpath); args.output.write_text(json.dumps(report,indent=2),encoding='utf-8'); frame.to_csv(args.output.with_suffix('.predictions.csv'),index=False); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
