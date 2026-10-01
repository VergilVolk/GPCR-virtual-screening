#!/usr/bin/env python3
"""Strong ECFP-to-frozen-pocket baseline for muscarinic triplet transfer."""
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


def fps(smiles, radius=2):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=2048)
    out = np.zeros((len(smiles), 2048), np.float32)
    for i, value in enumerate(smiles):
        DataStructs.ConvertToNumpyArray(gen.GetFingerprint(Chem.MolFromSmiles(value)), out[i])
    return torch.as_tensor(out)


class ECFPPocket(nn.Module):
    def __init__(self, seed):
        super().__init__(); torch.manual_seed(seed)
        self.project = nn.Linear(2048, 128, bias=False)
        nn.init.normal_(self.project.weight, std=0.01)
    def forward(self, x): return F.normalize(self.project(x), dim=-1)


def fit(smiles, triplets, pocket, mol_index, pocket_index, seed, epochs, lr, margin):
    x = fps(smiles); model = ECFPPocket(seed); opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-2)
    mi = torch.as_tensor([mol_index[str(v)] for v in triplets.canonical_smiles], dtype=torch.long)
    pi = torch.as_tensor([pocket_index[str(v)] for v in triplets.positive_subtype], dtype=torch.long)
    ni = torch.as_tensor([pocket_index[str(v)] for v in triplets.negative_subtype], dtype=torch.long)
    for _ in range(epochs):
        opt.zero_grad(); z = model(x); score = z @ pocket.T
        loss = F.relu(margin - score[mi, pi] + score[mi, ni]).mean()
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); opt.step()
    return model


def score(model, archive, table):
    ids = list(map(str, archive["molecule_ids"])); mi = {v:i for i,v in enumerate(ids)}
    pids = [str(v).split("_")[0] for v in archive["pocket_ids"]]; pi = {v:i for i,v in enumerate(pids)}
    with torch.inference_mode():
        z = model(fps(ids)); pocket = F.normalize(torch.as_tensor(archive["pocket_embeddings"].astype(np.float32)), dim=-1)
        matrix = (z @ pocket.T).numpy()
    return np.asarray([matrix[mi[str(r.canonical_smiles)], pi[str(r.positive_subtype)]]
                       - matrix[mi[str(r.canonical_smiles)], pi[str(r.negative_subtype)]]
                       for r in table.itertuples()])


def summary(table, values):
    directions = {name: float((values[group.index] > 0).mean())
                  for name, group in table.groupby("m4_direction")} if "m4_direction" in table else {}
    return {"n": int(len(table)), "molecules": int(table.canonical_smiles.nunique()),
            "pair_accuracy": float((values > 0).mean()),
            "direction_accuracy": directions,
            "direction_balanced_accuracy": float(np.mean(list(directions.values()))) if directions else None}


def main():
    p=argparse.ArgumentParser(); p.add_argument('--train-representations',type=Path,required=True); p.add_argument('--train-triplets',type=Path,required=True); p.add_argument('--strict-test',type=Path,required=True); p.add_argument('--zero-representations',type=Path,required=True); p.add_argument('--zero-test',type=Path,required=True); p.add_argument('--output',type=Path,required=True); p.add_argument('--seeds',default='20260924,20260925,20260926'); p.add_argument('--epochs',type=int,default=300); p.add_argument('--lr',type=float,default=1e-3); p.add_argument('--margin',type=float,default=.1); a=p.parse_args(); torch.set_num_threads(max(1,min(8,torch.get_num_threads())))
    arc=np.load(a.train_representations,allow_pickle=False); zero=np.load(a.zero_representations,allow_pickle=False)
    trip=pd.read_csv(a.train_triplets); strict=pd.read_csv(a.strict_test); ztest=pd.read_csv(a.zero_test)
    ids=list(map(str,arc['molecule_ids'])); mi={v:i for i,v in enumerate(ids)}; pids=[str(v).split('_')[0] for v in arc['pocket_ids']]; pi={v:i for i,v in enumerate(pids)}
    # Keep the benchmark aligned with the frozen DrugCLIP archive.  A handful of
    # chemically valid rows can fail upstream 3D embedding; they must be audited
    # and removed explicitly rather than raising midway through training.
    available = set(ids)
    missing_train = sorted(set(trip.canonical_smiles.astype(str)) - available)
    trip = trip[trip.canonical_smiles.astype(str).isin(available)].reset_index(drop=True)
    missing_strict = sorted(set(strict.canonical_smiles.astype(str)) - available)
    if missing_strict:
        raise ValueError(f"Strict-test molecules absent from frozen archive: {len(missing_strict)}")
    pocket=F.normalize(torch.as_tensor(arc['pocket_embeddings'].astype(np.float32)),dim=-1)
    strict_values=[]; zero_values=[]; seed_rows=[]
    for seed in map(int,a.seeds.split(',')):
        model=fit(ids,trip,pocket,mi,pi,seed,a.epochs,a.lr,a.margin)
        sv=score(model,arc,strict); zv=score(model,zero,ztest); strict_values.append(sv); zero_values.append(zv)
        seed_rows.append({'seed':seed,'strict_accuracy':float((sv>0).mean()),'zero_m4_accuracy':float((zv>0).mean())})
    strict_mean=np.mean(strict_values,axis=0); zero_mean=np.mean(zero_values,axis=0)
    report={'method':'ECFP4 linear projection to frozen DrugCLIP pocket space with triplet loss','training_audit':{'usable_triplets':int(len(trip)),'missing_train_molecules':len(missing_train),'missing_train_smiles':missing_train},'strict_seen_subtypes':summary(strict,strict_mean),'unseen_m4':summary(ztest,zero_mean),'per_seed':seed_rows,'hyperparameters':{'epochs':a.epochs,'lr':a.lr,'margin':a.margin,'seeds':a.seeds},'claim_boundary':'Strong 2D-molecule plus frozen-3D-pocket baseline; not PAM efficacy.'}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(report,indent=2),encoding='utf-8')
    strict_out = strict.copy(); zero_out = ztest.copy()
    strict_out['ecfp_pocket_delta'] = strict_mean
    zero_out['ecfp_pocket_delta'] = zero_mean
    strict_out.to_csv(a.output.with_suffix('.strict_predictions.csv'), index=False)
    zero_out.to_csv(a.output.with_suffix('.m4_predictions.csv'), index=False)
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
