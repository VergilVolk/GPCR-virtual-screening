#!/usr/bin/env python3
"""Strong ligand-only baseline for the exact 13-target LOSO protocol."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.sparse import csr_matrix
from sklearn.linear_model import LogisticRegression
from finetune_drugclip_gpcr_screening import binary_metrics


METRICS=['roc_auc','pr_auc','bedroc_alpha20','ef1pct','ef5pct']


def sparse_ecfp(smiles):
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048);rows=[];cols=[]
    for i,value in enumerate(smiles):
        bits=list(gen.GetFingerprint(Chem.MolFromSmiles(value)).GetOnBits());rows.extend([i]*len(bits));cols.extend(bits)
    return csr_matrix((np.ones(len(rows),dtype=np.float32),(rows,cols)),shape=(len(smiles),2048))


def summary(table,scores,targets):
    per={}
    for target in targets:
        keep=table.target.eq(target).to_numpy();per[target]=binary_metrics(table.loc[keep,'label'].to_numpy(int),scores[keep])
    return {'macro':{m:float(np.mean([per[t][m] for t in targets])) for m in METRICS},'per_target':per}


def main():
    p=argparse.ArgumentParser();p.add_argument('--predictions',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();table=pd.read_csv(a.predictions);targets=list(table.target.drop_duplicates());unique=list(table.canonical_smiles.drop_duplicates());mi={v:i for i,v in enumerate(unique)};x=sparse_ecfp(unique);row=np.asarray([mi[v] for v in table.canonical_smiles]);labels=table.label.to_numpy(int);scaff=table.murcko_scaffold.astype(str).to_numpy();pred=np.full(len(table),np.nan);audit={}
    for held,target in enumerate(targets):
        test=table.target.eq(target).to_numpy();held_scaff=set(scaff[test]);train=(~test)&~np.asarray([v in held_scaff for v in scaff]);model=LogisticRegression(C=1.,class_weight='balanced',solver='liblinear',max_iter=2000,random_state=20260926);model.fit(x[row[train]],labels[train]);pred[test]=model.predict_proba(x[row[test]])[:,1];audit[target]={'train_pairs':int(train.sum()),'test_pairs':int(test.sum()),'scaffold_overlap':0}
    report={'protocol':'Ligand-only ECFP4 logistic under exact 13-target LOSO and global held-scaffold purge','metrics':summary(table,pred,targets),'training_audit':audit,'claim_boundary':'Retrospective ligand-only baseline; no efficacy claim.'};a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');out=table[['target','label','canonical_smiles','murcko_scaffold']].copy();out['ecfp4_logistic']=pred;out.to_csv(a.output.with_suffix('.predictions.csv'),index=False);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
