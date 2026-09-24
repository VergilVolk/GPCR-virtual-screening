#!/usr/bin/env python3
"""Scaffold-held-out pair ranking on explicit M1/M2/M3/M5 selectivity contrasts."""
from __future__ import annotations
import argparse,json,hashlib,random
from pathlib import Path
import numpy as np,pandas as pd,torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from torch.nn import functional as F
from train_drugclip_m4_pocket_realign import LowRankResidual

def scaffold(s):
    m=Chem.MolFromSmiles(s); x=MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    return x or Chem.MolToSmiles(m,canonical=True)

def folds_for(smiles,n=5):
    groups={}
    for s in smiles: groups.setdefault(scaffold(s),[]).append(s)
    bins=[[] for _ in range(n)]; sizes=[0]*n
    for _,members in sorted(groups.items(),key=lambda x:(-len(x[1]),x[0])):
        j=min(range(n),key=lambda k:(sizes[k],k)); bins[j]+=members; sizes[j]+=len(members)
    return {s:i for i,b in enumerate(bins) for s in b}

class Model(torch.nn.Module):
    def __init__(self,dim,rank): super().__init__(); self.mol=LowRankResidual(dim,rank); self.pocket=LowRankResidual(dim,rank)
    def forward(self,m,p): return self.mol(m)@self.pocket(p).T

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--embeddings',type=Path,required=True); ap.add_argument('--contrasts',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True); ap.add_argument('--seeds',default='20260924,20260925,20260926')
    ap.add_argument('--epochs',type=int,default=500); ap.add_argument('--rank',type=int,default=16); ap.add_argument('--lr',type=float,default=1e-3)
    ap.add_argument('--weight-decay',type=float,default=1e-3); ap.add_argument('--margin',type=float,default=0.1); ap.add_argument('--preserve-weight',type=float,default=0.2)
    args=ap.parse_args(); torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    a=np.load(args.embeddings,allow_pickle=False); mids=list(map(str,a['molecule_ids'])); pids=list(map(str,a['pocket_ids']))
    mi={v:i for i,v in enumerate(mids)}; pi={v.split('_')[0]:i for i,v in enumerate(pids)}; d=pd.read_csv(args.contrasts)
    d=d[~d.overlaps_pacer_m4.astype(bool)].copy(); fmap=folds_for(d.canonical_smiles.unique()); d['fold']=d.canonical_smiles.map(fmap)
    molecule=torch.tensor(a['molecule_embeddings'],dtype=torch.float32); pocket=torch.tensor(a['pocket_embeddings'],dtype=torch.float32)
    oi,op=F.normalize(molecule,dim=1),F.normalize(pocket,dim=1)
    zero=np.asarray([a['scores'][pi[r.positive_subtype],mi[r.canonical_smiles]]-a['scores'][pi[r.negative_subtype],mi[r.canonical_smiles]] for r in d.itertuples()])
    reports=[]; rows=[]
    for seed in map(int,args.seeds.split(',')):
        pred=np.full(len(d),np.nan)
        for fold in sorted(d.fold.unique()):
            random.seed(seed+fold*100); np.random.seed(seed+fold*100); torch.manual_seed(seed+fold*100)
            tr=d.fold.to_numpy()!=fold; te=~tr; model=Model(molecule.shape[1],args.rank)
            opt=torch.optim.AdamW(model.parameters(),lr=args.lr,weight_decay=args.weight_decay)
            train=d.iloc[np.flatnonzero(tr)]
            molidx=torch.tensor([mi[s] for s in train.canonical_smiles]); pos=torch.tensor([pi[s] for s in train.positive_subtype]); neg=torch.tensor([pi[s] for s in train.negative_subtype])
            for _ in range(args.epochs):
                opt.zero_grad(); score=model(molecule,pocket); delta=score[molidx,pos]-score[molidx,neg]
                loss=F.relu(args.margin-delta).mean(); zm,zp=model.mol(molecule),model.pocket(pocket)
                loss=loss+args.preserve_weight*((1-(zm*oi).sum(1)).mean()+(1-(zp*op).sum(1)).mean()); loss.backward(); opt.step()
            model.eval(); score=model(molecule,pocket).detach().numpy()
            for idx in np.flatnonzero(te):
                r=d.iloc[idx]; pred[idx]=score[mi[r.canonical_smiles],pi[r.positive_subtype]]-score[mi[r.canonical_smiles],pi[r.negative_subtype]]
        reports.append({'seed':seed,'pair_accuracy':float((pred>0).mean()),'mean_margin':float(pred.mean()),'median_margin':float(np.median(pred))})
        for i,r in d.iterrows(): rows.append({'canonical_smiles':r.canonical_smiles,'positive_subtype':r.positive_subtype,'negative_subtype':r.negative_subtype,'fold':int(r.fold),'seed':seed,'delta':float(pred[d.index.get_loc(i)])})
    # Full-data artefact is saved only for downstream initialization; CV metrics remain primary.
    torch.manual_seed(20260924); full=Model(molecule.shape[1],args.rank); opt=torch.optim.AdamW(full.parameters(),lr=args.lr,weight_decay=args.weight_decay)
    molidx=torch.tensor([mi[s] for s in d.canonical_smiles]); pos=torch.tensor([pi[s] for s in d.positive_subtype]); neg=torch.tensor([pi[s] for s in d.negative_subtype])
    for _ in range(args.epochs):
        opt.zero_grad(); score=full(molecule,pocket); loss=F.relu(args.margin-score[molidx,pos]+score[molidx,neg]).mean(); zm,zp=full.mol(molecule),full.pocket(pocket)
        loss=loss+args.preserve_weight*((1-(zm*oi).sum(1)).mean()+(1-(zp*op).sum(1)).mean()); loss.backward(); opt.step()
    args.output.parent.mkdir(parents=True,exist_ok=True); torch.save({'model':full.state_dict(),'rank':args.rank,'pocket_ids':pids},args.output.with_suffix('.pt'))
    report={'evidence_level':'retrospective_muscarinic_family_selectivity_pretraining','n_pairs':len(d),'n_molecules':d.canonical_smiles.nunique(),'n_scaffold_folds':len(set(fmap.values())),
            'zero_shot_pair_accuracy':float((zero>0).mean()),'seed_reports':reports,'mean_oof_accuracy':float(np.mean([r['pair_accuracy'] for r in reports])),
            'claim_boundary':'Small heterogeneous ChEMBL contrast set; family pretraining pilot only.'}
    args.output.write_text(json.dumps(report,indent=2),encoding='utf-8'); pd.DataFrame(rows).to_csv(args.output.with_suffix('.predictions.csv'),index=False); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
