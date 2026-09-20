# -*- coding: utf-8 -*-
"""PACER-Rank v0.1: leave-one-medicinal-chemistry-series-out potency ranking.

Primary endpoint: within-series Spearman correlation on an entirely unseen
source component. This avoids pooled cross-series score artifacts.

Models:
  - Tanimoto kNN
  - Random Forest / Extra Trees absolute pEC50 regression
  - NeuralERM: absolute pEC50 regression
  - PACERRank: absolute regression + within-training-series RankNet loss
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.metrics import mean_absolute_error

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
OUT = PROJECT / "results" / "pacer_potency_loso_v01"
SEEDS = [42, 43, 44]
MIN_TEST_SIZE = 8

DESC = [
    Descriptors.MolWt, Descriptors.MolLogP, Descriptors.NumHDonors,
    Descriptors.NumHAcceptors, Descriptors.TPSA,
    Descriptors.NumRotatableBonds, Descriptors.NumAromaticRings,
    Descriptors.FractionCSP3, Descriptors.HeavyAtomCount,
]


def set_seed(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed); torch.set_num_threads(6)


def features(smiles):
    fps, desc = [], []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros(2048, np.float32); DataStructs.ConvertToNumpyArray(fp, arr)
        fps.append(arr); desc.append([fn(mol) for fn in DESC])
    return np.asarray(fps), np.asarray(desc, np.float32)


def preprocess(train_fp, test_fp, train_desc, test_desc):
    selected = np.argsort(train_fp.var(0))[-512:]
    mean, std = train_desc.mean(0), train_desc.std(0); std[std < 1e-8] = 1
    return (
        np.hstack([train_fp[:, selected], (train_desc-mean)/std]).astype(np.float32),
        np.hstack([test_fp[:, selected], (test_desc-mean)/std]).astype(np.float32),
    )


class Ranker(torch.nn.Module):
    def __init__(self, n):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(n, 96), torch.nn.LayerNorm(96), torch.nn.GELU(),
            torch.nn.Dropout(.15), torch.nn.Linear(96, 48), torch.nn.GELU(),
            torch.nn.Linear(48, 1),
        )
    def forward(self, x): return self.net(x).squeeze(-1)


def within_domain_pairs(groups, y, max_pairs_per_group=1500, seed=42):
    rng = np.random.RandomState(seed); left, right, target = [], [], []
    for group in sorted(set(groups)):
        idx = np.where(groups == group)[0]
        candidates = [(int(i), int(j)) for pos, i in enumerate(idx) for j in idx[pos+1:]
                      if abs(y[i]-y[j]) >= 0.30]
        if len(candidates) > max_pairs_per_group:
            chosen = rng.choice(len(candidates), max_pairs_per_group, replace=False)
            candidates = [candidates[k] for k in chosen]
        for i, j in candidates:
            left.append(i); right.append(j); target.append(1.0 if y[j] > y[i] else 0.0)
    return (torch.tensor(left, dtype=torch.long), torch.tensor(right, dtype=torch.long),
            torch.tensor(target, dtype=torch.float32))


def neural_predict(train_x, test_x, y, groups, seed, rank_weight):
    set_seed(seed); xt = torch.from_numpy(train_x); yt = torch.from_numpy(y.astype(np.float32))
    mean, std = float(y.mean()), max(float(y.std()), 1e-6)
    y_scaled = (yt-mean)/std
    left, right, pair_target = within_domain_pairs(groups, y, seed=seed)
    model = Ranker(train_x.shape[1]); opt = torch.optim.AdamW(model.parameters(),1e-3,weight_decay=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=140)
    for _ in range(140):
        model.train(); opt.zero_grad(); pred = model(xt)
        loss = torch.nn.functional.smooth_l1_loss(pred, y_scaled)
        if rank_weight > 0 and len(left):
            difference = pred[right]-pred[left]
            rank_loss = torch.nn.functional.binary_cross_entropy_with_logits(difference, pair_target)
            loss = loss + rank_weight*rank_loss
        loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),5); opt.step(); sched.step()
    model.eval()
    with torch.no_grad(): return model(torch.from_numpy(test_x)).numpy()*std+mean


def knn_predict(train_fp, test_fp, y, k=5):
    output=[]
    train_bits=[DataStructs.CreateFromBitString(''.join('1' if v else '0' for v in row)) for row in train_fp]
    test_bits=[DataStructs.CreateFromBitString(''.join('1' if v else '0' for v in row)) for row in test_fp]
    for fp in test_bits:
        sims=np.asarray(DataStructs.BulkTanimotoSimilarity(fp,train_bits)); idx=np.argsort(sims)[-k:]
        output.append(np.average(y[idx],weights=np.maximum(sims[idx],1e-6)**2))
    return np.asarray(output)


def main():
    OUT.mkdir(parents=True,exist_ok=True); df=pd.read_csv(DATA).reset_index(drop=True)
    fps,desc=features(df.canonical_smiles); y=df.pEC50.to_numpy(float); groups=df.source_component.to_numpy(str)
    test_groups=[g for g,n in df.source_component.value_counts().items() if n>=MIN_TEST_SIZE]
    methods=["TanimotoKNN","RandomForest","ExtraTrees","NeuralERM","PACERRank"]
    predictions={m:np.full(len(df),np.nan) for m in methods}; fold_reports=[]
    for group in test_groups:
        test=groups==group; train=~test; train_x,test_x=preprocess(fps[train],fps[test],desc[train],desc[test])
        predictions["TanimotoKNN"][test]=knn_predict(fps[train],fps[test],y[train])
        for name,cls in [("RandomForest",RandomForestRegressor),("ExtraTrees",ExtraTreesRegressor)]:
            model=cls(n_estimators=600,min_samples_leaf=2,random_state=42,n_jobs=-1); model.fit(train_x,y[train])
            predictions[name][test]=model.predict(test_x)
        for name,weight in [("NeuralERM",0.0),("PACERRank",0.5)]:
            seed_preds=[neural_predict(train_x,test_x,y[train],groups[train],seed,weight) for seed in SEEDS]
            predictions[name][test]=np.vstack(seed_preds).mean(0)
        fold={"group":group,"n":int(test.sum()),"pEC50_mean":float(y[test].mean()),"methods":{}}
        for method in methods:
            pred=predictions[method][test]; rho=float(spearmanr(y[test],pred).statistic)
            fold["methods"][method]={"Spearman":rho,"MAE":float(mean_absolute_error(y[test],pred))}
        fold_reports.append(fold); print(group,fold["methods"],flush=True)
    report={"primary_metric":"macro within-series Spearman","min_test_size":MIN_TEST_SIZE,"folds":fold_reports,"aggregate":{}}
    for method in methods:
        rhos=np.asarray([f["methods"][method]["Spearman"] for f in fold_reports]); maes=np.asarray([f["methods"][method]["MAE"] for f in fold_reports])
        report["aggregate"][method]={
            "macro_Spearman":float(np.nanmean(rhos)),"median_Spearman":float(np.nanmedian(rhos)),
            "worst_Spearman":float(np.nanmin(rhos)),"positive_rho_series":int(np.sum(rhos>0)),
            "n_series":int(len(rhos)),"macro_MAE":float(np.mean(maes)),
        }
    out=df[["canonical_molecule_id","canonical_smiles","pEC50","source_component"]].copy()
    for method,pred in predictions.items():out[method]=pred
    out.to_csv(OUT/"loso_predictions.csv",index=False)
    with open(OUT/"loso_metrics.json","w",encoding="utf-8") as f:json.dump(report,f,ensure_ascii=False,indent=2)
    print(json.dumps(report["aggregate"],ensure_ascii=False,indent=2))


if __name__=="__main__":main()
