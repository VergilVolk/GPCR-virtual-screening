#!/usr/bin/env python3
"""Select a molecule adapter on non-M4 data, then evaluate frozen M4 once."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from torch import nn
from torch.nn import functional as F

from evaluate_muscarinic_ecfp_pocket_baseline import summary


VARIANTS = {"ecfp2_bit": (2, False), "ecfp3_bit": (3, False), "ecfp2_count": (2, True)}


def feature(smiles, variant):
    if variant == "ecfp2_ecfp3_bit":
        return torch.cat([feature(smiles, "ecfp2_bit"), feature(smiles, "ecfp3_bit")], dim=1)
    radius, count = VARIANTS[variant]
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=2048)
    out = np.zeros((len(smiles), 2048), np.float32)
    for i, value in enumerate(smiles):
        mol = Chem.MolFromSmiles(value)
        fp = gen.GetCountFingerprint(mol) if count else gen.GetFingerprint(mol)
        DataStructs.ConvertToNumpyArray(fp, out[i])
    if count:
        out = np.log1p(out)
    return torch.as_tensor(out)


class Adapter(nn.Module):
    def __init__(self, width, seed):
        super().__init__(); torch.manual_seed(seed); self.project=nn.Linear(width,128,bias=False); nn.init.normal_(self.project.weight,std=.01)
    def forward(self,x): return F.normalize(self.project(x),dim=-1)


def train(x,pocket,trip,mi,pi,seed,args):
    model=Adapter(x.shape[1],seed); opt=torch.optim.AdamW(model.parameters(),lr=args.lr,weight_decay=args.weight_decay)
    mol=torch.as_tensor([mi[str(v)] for v in trip.canonical_smiles]); pos=torch.as_tensor([pi[str(v)] for v in trip.positive_subtype]); neg=torch.as_tensor([pi[str(v)] for v in trip.negative_subtype])
    for _ in range(args.epochs):
        opt.zero_grad(); score=model(x)@pocket.T; loss=F.relu(args.margin-score[mol,pos]+score[mol,neg]).mean(); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),1.); opt.step()
    return model


def score(model,x,pocket,table,mi,pi):
    with torch.inference_mode(): matrix=(model(x)@pocket.T).numpy()
    return np.asarray([matrix[mi[str(r.canonical_smiles)],pi[str(r.positive_subtype)]]-
                       matrix[mi[str(r.canonical_smiles)],pi[str(r.negative_subtype)]]
                       for r in table.itertuples()])


def main():
    p=argparse.ArgumentParser();p.add_argument('--train-representations',type=Path,required=True);p.add_argument('--train-triplets',type=Path,required=True);p.add_argument('--strict-test',type=Path,required=True);p.add_argument('--zero-representations',type=Path,required=True);p.add_argument('--zero-test',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seeds',default='20260924,20260925,20260926');p.add_argument('--epochs',type=int,default=100);p.add_argument('--lr',type=float,default=1e-3);p.add_argument('--weight-decay',type=float,default=1e-2);p.add_argument('--margin',type=float,default=.1);a=p.parse_args();torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(a.train_representations,allow_pickle=False);zero_arc=np.load(a.zero_representations,allow_pickle=False);ids=list(map(str,arc['molecule_ids']));mi={v:i for i,v in enumerate(ids)};pids=[str(v).split('_')[0] for v in arc['pocket_ids']];pi={v:i for i,v in enumerate(pids)};pocket=F.normalize(torch.as_tensor(arc['pocket_embeddings'].astype(np.float32)),dim=-1)
    trip=pd.read_csv(a.train_triplets);trip=trip[trip.canonical_smiles.astype(str).isin(mi)].reset_index(drop=True);strict=pd.read_csv(a.strict_test);zero=pd.read_csv(a.zero_test);seeds=list(map(int,a.seeds.split(',')))
    candidate={}; trained={}
    for variant in ['ecfp2_bit','ecfp3_bit','ecfp2_count','ecfp2_ecfp3_bit']:
        x=feature(ids,variant); vals=[]; models=[]
        for seed in seeds:
            model=train(x,pocket,trip,mi,pi,seed,a); models.append(model); vals.append(score(model,x,pocket,strict,mi,pi))
        mean=np.mean(vals,axis=0);candidate[variant]={'strict_accuracy':float((mean>0).mean()),'per_seed':[float((v>0).mean()) for v in vals],'feature_width':int(x.shape[1])};trained[variant]=(models,x)
    selected=sorted(candidate,key=lambda v:(-candidate[v]['strict_accuracy'],candidate[v]['feature_width'],v))[0]
    zids=list(map(str,zero_arc['molecule_ids']));zmi={v:i for i,v in enumerate(zids)};zpids=[str(v).split('_')[0] for v in zero_arc['pocket_ids']];zpi={v:i for i,v in enumerate(zpids)};zpocket=F.normalize(torch.as_tensor(zero_arc['pocket_embeddings'].astype(np.float32)),dim=-1);zx=feature(zids,selected)
    zvals=[score(model,zx,zpocket,zero,zmi,zpi) for model in trained[selected][0]];zmean=np.mean(zvals,axis=0)
    a.output.parent.mkdir(parents=True,exist_ok=True);checkpoints=[]
    for seed,model in zip(seeds,trained[selected][0]):
        path=a.output.with_name(f'{a.output.stem}.{selected}.seed{seed}.pt');torch.save({'state_dict':model.state_dict(),'variant':selected,'fingerprint_radius':3,'fingerprint_size':2048,'seed':seed,'training_pairs':len(trip),'pocket_ids':pids,'claim_boundary':'Subtype activity adapter; not PAM efficacy.'},path);checkpoints.append(str(path))
    report={'selection_rule':'Highest ensemble accuracy on molecule-disjoint non-M4 strict panel; ties prefer smaller feature width; frozen before M4 evaluation','non_m4_candidates':candidate,'selected':selected,'frozen_m4':summary(zero,zmean),'m4_per_seed':[float((v>0).mean()) for v in zvals],'checkpoints':checkpoints,'deployment_rule':'Unweighted mean of the three fixed-seed score differences.','hyperparameters':{'epochs':a.epochs,'lr':a.lr,'margin':a.margin,'seeds':a.seeds},'claim_boundary':'Muscarinic subtype activity transfer; not PAM efficacy or binding affinity.'}
    a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');out=zero.copy();out['selected_adapter_delta']=zmean;out.to_csv(a.output.with_suffix('.m4_predictions.csv'),index=False);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
