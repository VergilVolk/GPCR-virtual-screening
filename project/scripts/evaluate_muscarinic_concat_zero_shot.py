#!/usr/bin/env python3
"""Train the joint molecule adapter on non-M4 triplets and test frozen M4."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F

from evaluate_muscarinic_ecfp_pocket_baseline import fps, summary
from evaluate_muscarinic_concat_adapter_loso import fit_model, pair_scores


def score_archive(model, arc, frame):
    ids=list(map(str,arc['molecule_ids'])); mi={v:i for i,v in enumerate(ids)}
    pids=[str(v).split('_')[0] for v in arc['pocket_ids']]; pi={v:i for i,v in enumerate(pids)}
    x=fps(ids); drug=F.normalize(torch.as_tensor(arc['molecule_embeddings'].astype(np.float32)),dim=-1); pocket=F.normalize(torch.as_tensor(arc['pocket_embeddings'].astype(np.float32)),dim=-1)
    return pair_scores(model,x,drug,pocket,frame,mi,pi)


def main():
    p=argparse.ArgumentParser(); p.add_argument('--train-representations',type=Path,required=True); p.add_argument('--train-triplets',type=Path,required=True); p.add_argument('--strict-test',type=Path,required=True); p.add_argument('--zero-representations',type=Path,required=True); p.add_argument('--zero-test',type=Path,required=True); p.add_argument('--output',type=Path,required=True); p.add_argument('--seeds',default='20260924,20260925,20260926'); p.add_argument('--epochs',type=int,default=100); p.add_argument('--lr',type=float,default=1e-3); p.add_argument('--weight-decay',type=float,default=1e-2); p.add_argument('--margin',type=float,default=.1); p.add_argument('--drug-init',type=float,default=.1); a=p.parse_args(); torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    train_arc=np.load(a.train_representations,allow_pickle=False); zero_arc=np.load(a.zero_representations,allow_pickle=False); trip=pd.read_csv(a.train_triplets); strict=pd.read_csv(a.strict_test); zero=pd.read_csv(a.zero_test)
    ids=list(map(str,train_arc['molecule_ids'])); mi={v:i for i,v in enumerate(ids)}; pids=[str(v).split('_')[0] for v in train_arc['pocket_ids']]; pi={v:i for i,v in enumerate(pids)}
    trip=trip[trip.canonical_smiles.astype(str).isin(mi)].reset_index(drop=True); x=fps(ids); drug=F.normalize(torch.as_tensor(train_arc['molecule_embeddings'].astype(np.float32)),dim=-1); pocket=F.normalize(torch.as_tensor(train_arc['pocket_embeddings'].astype(np.float32)),dim=-1)
    strict_values=[]; zero_values=[]; rows=[]
    for seed in map(int,a.seeds.split(',')):
        model=fit_model(x,drug,pocket,trip,mi,pi,seed,a); sv=score_archive(model,train_arc,strict); zv=score_archive(model,zero_arc,zero); strict_values.append(sv); zero_values.append(zv)
        rows.append({'seed':seed,'strict_accuracy':float((sv>0).mean()),'zero_m4_accuracy':float((zv>0).mean())})
    sm=np.mean(strict_values,axis=0); zm=np.mean(zero_values,axis=0)
    report={'method':'Joint ECFP plus frozen DrugCLIP molecule adapter into frozen DrugCLIP pocket space','training_pairs':len(trip),'strict_seen_subtypes':summary(strict,sm),'unseen_m4':summary(zero,zm),'per_seed':rows,'hyperparameters':{'epochs':a.epochs,'lr':a.lr,'margin':a.margin,'drug_init':a.drug_init,'seeds':a.seeds},'claim_boundary':'M4-vs-muscarinic subtype activity transfer; not PAM efficacy or binding affinity.'}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(report,indent=2),encoding='utf-8'); so=strict.copy();so['concat_delta']=sm;zo=zero.copy();zo['concat_delta']=zm;so.to_csv(a.output.with_suffix('.strict_predictions.csv'),index=False);zo.to_csv(a.output.with_suffix('.m4_predictions.csv'),index=False);print(json.dumps(report,indent=2))
if __name__=='__main__': main()
