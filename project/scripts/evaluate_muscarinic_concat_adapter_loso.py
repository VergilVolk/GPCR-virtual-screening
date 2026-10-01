#!/usr/bin/env python3
"""Joint ECFP/DrugCLIP molecule adapter into frozen DrugCLIP pocket space."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F

from evaluate_muscarinic_ecfp_pocket_baseline import fps
from evaluate_muscarinic_hybrid_target_loso import TARGETS, metrics


class ConcatAdapter(nn.Module):
    def __init__(self, seed: int, drug_init: float):
        super().__init__(); torch.manual_seed(seed)
        self.project = nn.Linear(2048 + 128, 128, bias=False)
        with torch.no_grad():
            nn.init.normal_(self.project.weight[:, :2048], std=0.01)
            self.project.weight[:, 2048:].zero_()
            self.project.weight[:, 2048:].add_(drug_init * torch.eye(128))

    def forward(self, ecfp, drug):
        return F.normalize(self.project(torch.cat([ecfp, drug], dim=-1)), dim=-1)


def fit_model(ecfp, drug, pocket, frame, mi, pi, seed, args):
    model = ConcatAdapter(seed, args.drug_init)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    mol = torch.as_tensor([mi[str(v)] for v in frame.canonical_smiles], dtype=torch.long)
    pos = torch.as_tensor([pi[str(v)] for v in frame.positive_subtype], dtype=torch.long)
    neg = torch.as_tensor([pi[str(v)] for v in frame.negative_subtype], dtype=torch.long)
    for _ in range(args.epochs):
        opt.zero_grad(); score = model(ecfp, drug) @ pocket.T
        loss = F.relu(args.margin - score[mol, pos] + score[mol, neg]).mean()
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); opt.step()
    return model


def pair_scores(model, ecfp, drug, pocket, frame, mi, pi):
    with torch.inference_mode(): matrix = (model(ecfp, drug) @ pocket.T).numpy()
    return np.asarray([matrix[mi[str(r.canonical_smiles)], pi[str(r.positive_subtype)]] -
                       matrix[mi[str(r.canonical_smiles)], pi[str(r.negative_subtype)]]
                       for r in frame.itertuples()])


def main():
    p=argparse.ArgumentParser(); p.add_argument('--representations',type=Path,required=True); p.add_argument('--benchmark-dir',type=Path,required=True); p.add_argument('--output',type=Path,required=True); p.add_argument('--seeds',default='20260924,20260925,20260926'); p.add_argument('--epochs',type=int,default=100); p.add_argument('--lr',type=float,default=1e-3); p.add_argument('--weight-decay',type=float,default=1e-2); p.add_argument('--margin',type=float,default=.1); p.add_argument('--drug-init',type=float,default=.1); a=p.parse_args(); torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(a.representations,allow_pickle=False); ids=list(map(str,arc['molecule_ids'])); mi={v:i for i,v in enumerate(ids)}; pids=[str(v).split('_')[0] for v in arc['pocket_ids']]; pi={v:i for i,v in enumerate(pids)}
    ecfp=fps(ids); drug=F.normalize(torch.as_tensor(arc['molecule_embeddings'].astype(np.float32)),dim=-1); pocket=F.normalize(torch.as_tensor(arc['pocket_embeddings'].astype(np.float32)),dim=-1)
    reports={}; rows=[]
    for held in TARGETS:
        tr0=pd.read_csv(a.benchmark_dir/f'{held}_train.csv'); te0=pd.read_csv(a.benchmark_dir/f'{held}_test.csv'); tr=tr0[tr0.canonical_smiles.astype(str).isin(mi)].reset_index(drop=True); te=te0[te0.canonical_smiles.astype(str).isin(mi)].reset_index(drop=True)
        values=[]
        for seed in map(int,a.seeds.split(',')):
            model=fit_model(ecfp,drug,pocket,tr,mi,pi,seed,a); value=pair_scores(model,ecfp,drug,pocket,te,mi,pi); values.append(value)
            with torch.no_grad():
                ecfp_norm=float(model.project.weight[:,:2048].norm()); drug_norm=float(model.project.weight[:,2048:].norm())
            for i,row in te.iterrows(): rows.append({'held_subtype':held,'canonical_smiles':row.canonical_smiles,'positive_subtype':row.positive_subtype,'negative_subtype':row.negative_subtype,'seed':seed,'concat_delta':float(value[i]),'ecfp_block_norm':ecfp_norm,'drug_block_norm':drug_norm})
        mean=np.mean(values,axis=0); reports[held]={'train_pairs':len(tr),'test_pairs':len(te),'test_molecules':te.canonical_smiles.nunique(),'concat_adapter':metrics(te,mean,held),'per_seed_accuracy':[float((v>0).mean()) for v in values]}
    macro={k:float(np.mean([reports[h]['concat_adapter'][k] for h in TARGETS])) for k in ['accuracy','direction_balanced_accuracy']}
    report={'method':'Trainable concatenated ECFP plus frozen DrugCLIP molecule adapter into frozen DrugCLIP pocket space','split':'strict leave-one-muscarinic-target-out with molecule-disjoint folds','folds':reports,'macro_average':macro,'hyperparameters':{'epochs':a.epochs,'lr':a.lr,'margin':a.margin,'drug_init':a.drug_init,'seeds':a.seeds},'claim_boundary':'Muscarinic subtype activity transfer; not PAM efficacy or binding affinity.'}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(report,indent=2),encoding='utf-8'); pd.DataFrame(rows).to_csv(a.output.with_suffix('.predictions.csv'),index=False); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
