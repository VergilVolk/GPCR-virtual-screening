#!/usr/bin/env python3
"""Audit whether the M4R active/decoy benchmark is separable by 2D chemistry alone."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

def ef(y, s, f=.01):
    n=max(1,int(np.ceil(len(y)*f))); return float(y[np.argsort(-s)[:n]].mean()/y.mean())

def main():
    p=argparse.ArgumentParser(); p.add_argument('--benchmark',type=Path,required=True); p.add_argument('--output',type=Path,required=True); a=p.parse_args()
    d=pd.read_csv(a.benchmark); gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    x=np.asarray([gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in d.canonical_smiles],dtype=np.float32)
    y=d.target.to_numpy(int); folds=d.scaffold_fold.to_numpy(int); pred=np.full(len(y),np.nan)
    for fold in sorted(np.unique(folds)):
        tr,te=folds!=fold,folds==fold
        m=LogisticRegression(C=.1,class_weight='balanced',max_iter=1000,solver='liblinear',random_state=20260924)
        m.fit(x[tr],y[tr]); pred[te]=m.predict_proba(x[te])[:,1]
    report={'method':'ECFP4 logistic regression','split':'same 5-fold Murcko scaffold OOF','n':len(y),
            'roc_auc':float(roc_auc_score(y,pred)),'average_precision':float(average_precision_score(y,pred)),
            'ef1pct':ef(y,pred),'interpretation':'High performance indicates active-decoy chemical bias; it is not evidence of learned M4 binding.'}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(report,indent=2),encoding='utf-8'); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
