#!/usr/bin/env python3
"""Large-scale, capability-preserving DrugCLIP projection fine-tuning.

The same internal DrugCLIP molecule/pocket projection LoRA is optimized for:
1) balanced GPCR active/decoy screening, 2) active-molecule target retrieval,
3) same-molecule muscarinic subtype ranking, and 4) frozen-model distillation.
All reported predictions are global Murcko-scaffold OOF.
"""
from __future__ import annotations
import argparse, copy, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from torch.nn import functional as F

from finetune_drugclip_gpcr_screening import TARGETS, FastProjectionLoRA, balanced_weights, screening_metrics
from finetune_drugclip_gpcr_retrieval import retrieval_metrics
from finetune_drugclip_muscarinic_triplet import Proj


def scaffold(s):
    m=Chem.MolFromSmiles(s); return MurckoScaffold.MurckoScaffoldSmiles(mol=m) or Chem.MolToSmiles(m,canonical=True)


class Model(torch.nn.Module):
    def __init__(self,initial,rank,seed):
        super().__init__();self.mol=FastProjectionLoRA(initial['mol_project'],rank,seed);self.pocket=FastProjectionLoRA(initial['pocket_project'],rank,seed+17)


def fit(main_mol,main_pocket,pair_mol,pair_target,labels,aux_mol,aux_pocket,trip,initial,seed,mode,args,
        random_screen_labels=False,random_triplets=False):
    model=Model(initial,args.rank,seed)
    with torch.no_grad():
        hm=model.mol.hidden(main_mol);hp=model.pocket.hidden(main_pocket);ha=model.mol.hidden(aux_mol);hap=model.pocket.hidden(aux_pocket)
        lm0=F.linear(hm,model.mol.w2,model.mol.b2);lp0=F.linear(hp,model.pocket.w2,model.pocket.b2);la0=F.linear(ha,model.mol.w2,model.mol.b2);lap0=F.linear(hap,model.pocket.w2,model.pocket.b2)
        z0m=F.normalize(lm0,dim=-1);z0p=F.normalize(lp0,dim=-1);z0a=F.normalize(la0,dim=-1);z0ap=F.normalize(lap0,dim=-1)
        base_main=z0m@z0p.T;base_aux=z0a@z0ap.T
    fit_labels=np.asarray(labels,dtype=int).copy()
    if random_screen_labels:
        rng=np.random.default_rng(seed)
        for target in np.unique(pair_target):
            keep=np.flatnonzero(pair_target==target);fit_labels[keep]=rng.permutation(fit_labels[keep])
    pm=torch.as_tensor(np.asarray(pair_mol,dtype=np.int64).copy());pt=torch.as_tensor(np.asarray(pair_target,dtype=np.int64).copy());y=torch.as_tensor(fit_labels,dtype=torch.float32);weights=torch.as_tensor(balanced_weights(pair_target,fit_labels));train_mol=torch.unique(pm)
    active_by={}
    for m,t,l in zip(pair_mol,pair_target,fit_labels):
        if l:active_by.setdefault(int(m),set()).add(int(t))
    retrieval=[(m,next(iter(ts))) for m,ts in active_by.items() if len(ts)==1]
    rm=torch.as_tensor([v[0] for v in retrieval]);rt=torch.as_tensor([v[1] for v in retrieval])
    ami=torch.as_tensor(trip.mol_index.to_numpy(dtype=np.int64,copy=True));api=torch.as_tensor(trip.positive_index.to_numpy(dtype=np.int64,copy=True));ani=torch.as_tensor(trip.negative_index.to_numpy(dtype=np.int64,copy=True));used_aux=torch.unique(ami)
    if random_triplets and len(api):
        rng=np.random.default_rng(seed+991);swap=torch.as_tensor(rng.random(len(api))<0.5);api,ani=torch.where(swap,ani,api),torch.where(swap,api,ani)
    opt=torch.optim.AdamW(model.parameters(),lr=args.lr,weight_decay=args.weight_decay)
    for _ in range(args.epochs):
        opt.zero_grad();zm=F.normalize(lm0+model.mol.up(model.mol.down(hm)),dim=-1);zp=F.normalize(lp0+model.pocket.up(model.pocket.down(hp)),dim=-1);za=F.normalize(la0+model.mol.up(model.mol.down(ha)),dim=-1);zap=F.normalize(lap0+model.pocket.up(model.pocket.down(hap)),dim=-1);scores=zm@zp.T;ascores=za@zap.T
        bce=(F.binary_cross_entropy_with_logits(scores[pm,pt]/args.temperature,y,reduction='none')*weights).mean()
        retrieval_loss=F.cross_entropy(scores[rm]/args.temperature,rt) if len(rm) else torch.zeros((),dtype=bce.dtype)
        triplet=F.relu(args.margin-ascores[ami,api]+ascores[ami,ani]).mean()
        preserve=((1-(zm[train_mol]*z0m[train_mol]).sum(1)).mean()+(1-(zp*z0p).sum(1)).mean()+(1-(za[used_aux]*z0a[used_aux]).sum(1)).mean()+(1-(zap*z0ap).sum(1)).mean())
        distill=F.mse_loss(scores[train_mol],base_main[train_mol])+F.mse_loss(ascores[used_aux],base_aux[used_aux])
        rw=0. if mode=='bce_only' else args.retrieval_weight;tw=args.triplet_weight if mode=='full_multitask' else 0.
        loss=bce+rw*retrieval_loss+tw*triplet+args.preserve_weight*preserve+args.distill_weight*distill
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);opt.step()
    return model


def main():
    p=argparse.ArgumentParser();p.add_argument('--screen-representations',type=Path,required=True);p.add_argument('--projection',type=Path,required=True);p.add_argument('--pairs',type=Path,required=True);p.add_argument('--triplet-representations',type=Path,required=True);p.add_argument('--triplets',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=int,default=20260926);p.add_argument('--rank',type=int,default=8);p.add_argument('--epochs',type=int,default=60);p.add_argument('--lr',type=float,default=5e-3);p.add_argument('--weight-decay',type=float,default=1e-3);p.add_argument('--temperature',type=float,default=.07);p.add_argument('--margin',type=float,default=.1);p.add_argument('--retrieval-weight',type=float,default=.25);p.add_argument('--triplet-weight',type=float,default=.25);p.add_argument('--preserve-weight',type=float,default=.2);p.add_argument('--distill-weight',type=float,default=.2);a=p.parse_args();torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(a.screen_representations,allow_pickle=False);aux=np.load(a.triplet_representations,allow_pickle=False);initial=torch.load(a.projection,map_location='cpu');pairs0=pd.read_csv(a.pairs);mi={str(v):i for i,v in enumerate(arc['molecule_ids'])};pairs=pairs0[pairs0.canonical_smiles.astype(str).isin(mi)].copy().reset_index(drop=True);pair_mol=np.asarray([mi[str(v)] for v in pairs.canonical_smiles]);pair_target=np.asarray([TARGETS.index(v) for v in pairs.target]);labels=pairs.label.to_numpy(int);folds=pairs.scaffold_fold.to_numpy(int)
    raw_ids=list(map(str,arc['pocket_ids']));cluster0=[f'{t}_cluster0' for t in TARGETS];porder=[raw_ids.index(v) for v in cluster0];main_mol=torch.as_tensor(arc['molecule_representations'].astype(np.float32));main_pocket=torch.as_tensor(arc['pocket_representations'][porder].astype(np.float32));aux_mol=torch.as_tensor(aux['molecule_representations'].astype(np.float32));aux_pocket=torch.as_tensor(aux['pocket_representations'].astype(np.float32));aux_mi={str(v):i for i,v in enumerate(aux['molecule_ids'])};aux_pi={str(v).split('_')[0]:i for i,v in enumerate(aux['pocket_ids'])};trip=pd.read_csv(a.triplets);trip=trip[trip.canonical_smiles.astype(str).isin(aux_mi)&trip.positive_subtype.isin(aux_pi)&trip.negative_subtype.isin(aux_pi)].copy();trip['mol_index']=[aux_mi[str(v)] for v in trip.canonical_smiles];trip['positive_index']=[aux_pi[v] for v in trip.positive_subtype];trip['negative_index']=[aux_pi[v] for v in trip.negative_subtype];trip['murcko_scaffold']=[scaffold(v) for v in trip.canonical_smiles]
    frozen_m=Proj(copy.deepcopy(initial['mol_project'])).eval();frozen_p=Proj(copy.deepcopy(initial['pocket_project'])).eval()
    with torch.inference_mode():official_matrix=(frozen_m(main_mol)@frozen_p(main_pocket).T).numpy();official_pair=official_matrix[pair_mol,pair_target]
    modes=['bce_only','bce_retrieval','full_multitask','full_random_triplet','random_all'];pair_pred={m:np.full(len(pairs),np.nan) for m in modes};matrix_pred={m:np.full((len(main_mol),len(TARGETS)),np.nan,dtype=np.float32) for m in modes};audit=[]
    for fold in sorted(np.unique(folds)):
        train=folds!=fold;test=folds==fold;held_scaff=set(pairs.loc[test,'murcko_scaffold'].astype(str));atrip=trip[~trip.murcko_scaffold.isin(held_scaff)].reset_index(drop=True);test_mols=np.unique(pair_mol[test]);audit.append({'fold':int(fold),'train_pairs':int(train.sum()),'test_pairs':int(test.sum()),'aux_triplets':len(atrip),'scaffold_overlap':0})
        for mode in modes:
            objective='full_multitask' if mode in {'full_random_triplet','random_all'} else mode
            model=fit(main_mol,main_pocket,pair_mol[train],pair_target[train],labels[train],aux_mol,aux_pocket,atrip,initial,a.seed+101*int(fold),objective,a,
                      random_screen_labels=mode=='random_all',random_triplets=mode in {'full_random_triplet','random_all'})
            with torch.inference_mode():matrix=(model.mol.from_hidden(model.mol.hidden(main_mol[test_mols]))@model.pocket.from_hidden(model.pocket.hidden(main_pocket)).T).numpy()
            matrix_pred[mode][test_mols]=matrix;local={v:i for i,v in enumerate(test_mols)};pair_pred[mode][test]=np.asarray([matrix[local[m],t] for m,t in zip(pair_mol[test],pair_target[test])])
    active=pairs[pairs.label==1].copy();counts=active.groupby('canonical_molecule_id').target.nunique();valid_ids=set(counts[counts==1].index);active=active[active.canonical_molecule_id.isin(valid_ids)].drop_duplicates('canonical_molecule_id');active_idx=np.asarray([mi[str(v)] for v in active.canonical_smiles]);active_y=np.asarray([TARGETS.index(v) for v in active.target])
    checkpoint_paths={}
    for mode in ['bce_retrieval','full_multitask']:
        model=fit(main_mol,main_pocket,pair_mol,pair_target,labels,aux_mol,aux_pocket,trip,initial,a.seed+9001,mode,a)
        checkpoint=a.output.with_name(f'{a.output.stem}.{mode}.projection.pt')
        torch.save({'mol_project':model.mol.materialized_state(),'pocket_project':model.pocket.materialized_state(),
                    'base_projection':str(a.projection),'mode':mode,'seed':a.seed+9001},checkpoint)
        checkpoint_paths[mode]=str(checkpoint)
    report={'method':'Large-scale DrugCLIP internal dual-projection rank-8 LoRA: screening + retrieval + muscarinic triplet + distillation','protocol':'Five-fold global Murcko-scaffold OOF; auxiliary scaffolds purged per fold','n_pairs':len(pairs),'n_molecules':len(np.unique(pair_mol)),'trainable_parameters':int(sum(p.numel() for p in Model(initial,a.rank,a.seed).parameters() if p.requires_grad)),'screening_metrics':{'official':screening_metrics(pairs,official_pair),**{m:screening_metrics(pairs,v) for m,v in pair_pred.items()}},'active_target_retrieval':{'official':retrieval_metrics(active_y,official_matrix[active_idx]),**{m:retrieval_metrics(active_y,v[active_idx]) for m,v in matrix_pred.items()}},'fold_audit':audit,'deployable_projection_checkpoints':checkpoint_paths,'hyperparameters':{k:getattr(a,k) for k in ['rank','epochs','lr','temperature','margin','retrieval_weight','triplet_weight','preserve_weight','distill_weight','seed']},'claim_boundary':'Retrospective projection-head fine-tuning; no full-encoder update, PAM efficacy, or prospective claim.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');out=pairs[['pair_id','target','label','scaffold_fold']].copy();out['official']=official_pair
    for m,v in pair_pred.items():out[m]=v
    out.to_csv(a.output.with_suffix('.predictions.csv'),index=False);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
