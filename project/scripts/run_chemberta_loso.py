# -*- coding: utf-8 -*-
"""Frozen ChemBERTa strong baseline for M4 PAM potency LOSO.

Uses the public `seyonec/ChemBERTa-zinc-base-v1` encoder (44M parameters),
mean-pools token embeddings, caches them, and evaluates Ridge/ExtraTrees on the
same leave-one-source-component-out benchmark as PACER-Rank.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.preprocessing import StandardScaler
from transformers import AutoModel, AutoTokenizer

PROJECT=Path(__file__).resolve().parents[1]
DATA=PROJECT/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv"
OUT=PROJECT/"results"/"chemberta_loso_v1"
CACHE=OUT/"chemberta_embeddings.npy"
MODEL_NAME="seyonec/ChemBERTa-zinc-base-v1"

def embed(smiles):
    if CACHE.exists(): return np.load(CACHE)
    tokenizer=AutoTokenizer.from_pretrained(MODEL_NAME); model=AutoModel.from_pretrained(MODEL_NAME)
    model.eval(); torch.set_num_threads(8); chunks=[]
    for start in range(0,len(smiles),32):
        batch=tokenizer(smiles[start:start+32],padding=True,truncation=True,max_length=256,return_tensors="pt")
        with torch.no_grad(): hidden=model(**batch).last_hidden_state
        mask=batch["attention_mask"].unsqueeze(-1); pooled=(hidden*mask).sum(1)/mask.sum(1).clamp(min=1)
        chunks.append(pooled.numpy()); print("embedded",min(start+32,len(smiles)),"/",len(smiles),flush=True)
    value=np.vstack(chunks).astype(np.float32); np.save(CACHE,value); return value

def main():
    OUT.mkdir(parents=True,exist_ok=True); d=pd.read_csv(DATA).reset_index(drop=True); X=embed(d.canonical_smiles.tolist())
    y=d.pEC50.to_numpy(float); groups=d.source_component.to_numpy(str)
    test_groups=[g for g,n in d.source_component.value_counts().items() if n>=8]
    methods=["ChemBERTaRidge","ChemBERTaExtraTrees"]; predictions={m:np.full(len(d),np.nan) for m in methods}; folds=[]
    for group in test_groups:
        te=groups==group; tr=~te; scaler=StandardScaler(); xtr=scaler.fit_transform(X[tr]); xte=scaler.transform(X[te])
        ridge=Ridge(alpha=100.0).fit(xtr,y[tr]); predictions[methods[0]][te]=ridge.predict(xte)
        tree=ExtraTreesRegressor(n_estimators=600,min_samples_leaf=2,max_features=.5,random_state=42,n_jobs=-1).fit(X[tr],y[tr]); predictions[methods[1]][te]=tree.predict(X[te])
        fold={"group":group,"n":int(te.sum()),"methods":{}}
        for m in methods:
            p=predictions[m][te]; fold["methods"][m]={"Spearman":float(spearmanr(y[te],p).statistic),"MAE":float(mean_absolute_error(y[te],p))}
        folds.append(fold); print(group,fold["methods"],flush=True)
    report={"model":MODEL_NAME,"frozen_encoder":True,"folds":folds,"aggregate":{}}
    for m in methods:
        r=np.array([f["methods"][m]["Spearman"] for f in folds]); a=np.array([f["methods"][m]["MAE"] for f in folds])
        report["aggregate"][m]={"macro_Spearman":float(np.mean(r)),"median_Spearman":float(np.median(r)),"worst_Spearman":float(np.min(r)),"positive_series":int((r>0).sum()),"macro_MAE":float(np.mean(a))}
    out=d[["canonical_molecule_id","pEC50","source_component"]].copy()
    for m,p in predictions.items():out[m]=p
    out.to_csv(OUT/"loso_predictions.csv",index=False)
    with open(OUT/"loso_metrics.json","w",encoding="utf-8") as f:json.dump(report,f,ensure_ascii=False,indent=2)
    print(json.dumps(report["aggregate"],ensure_ascii=False,indent=2))

if __name__=="__main__":main()
