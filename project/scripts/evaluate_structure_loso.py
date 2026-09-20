# -*- coding: utf-8 -*-
"""LOSO ablation for label-free 7TRQ pose features.

Primary endpoint is macro within-series Spearman. Pooled metrics are omitted to
avoid the cross-series Simpson artifact documented in the research charter.
"""
from __future__ import annotations
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from lightgbm import LGBMRegressor

RDLogger.DisableLog("rdApp.*")
PROJECT=Path(__file__).resolve().parents[1]
DATA=PROJECT/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv"
STATE=os.environ.get("M4_STATE_EVAL","7TRQ").upper()
STATE_DIR={"7TRQ":"pacer_structure_loso_v01","7TRP":"pacer_structure_7trp_v01","7TRS":"pacer_structure_7trs_v01"}[STATE]
FEAT=PROJECT/"results"/STATE_DIR/"docking_features.csv"
OUT=PROJECT/"results"/STATE_DIR
MIN_N=8
DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,
      Descriptors.NumHAcceptors,Descriptors.TPSA,Descriptors.NumRotatableBonds,
      Descriptors.NumAromaticRings,Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]


def twod(smiles):
    fps=[]; desc=[]
    for s in smiles:
        m=Chem.MolFromSmiles(s); fp=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048)
        a=np.zeros(2048,np.float32); DataStructs.ConvertToNumpyArray(fp,a)
        fps.append(a); desc.append([f(m) for f in DESC])
    return np.asarray(fps),np.asarray(desc,np.float32)


def aggregate(rows):
    rho=np.asarray([r["Spearman"] for r in rows],float)
    mae=np.asarray([r["MAE"] for r in rows],float)
    return {"macro_Spearman":float(np.nanmean(rho)),"median_Spearman":float(np.nanmedian(rho)),
            "worst_Spearman":float(np.nanmin(rho)),"positive_rho_series":int(np.sum(rho>0)),
            "n_series":len(rows),"macro_MAE":float(np.mean(mae))}


def main():
    d=pd.read_csv(DATA); f=pd.read_csv(FEAT)
    f=f[(f.seed==42)&(f.error.fillna("")=="")].drop_duplicates("canonical_molecule_id")
    d=d.merge(f,on="canonical_molecule_id",how="inner")
    if len(d)<400: raise RuntimeError(f"structural feature run incomplete: n={len(d)}")
    fp,desc=twod(d.canonical_smiles); y=d.pEC50.to_numpy(float); groups=d.source_component.astype(str).to_numpy()
    base_cols=["vina_affinity","n_heavy_atoms"]
    ifp_cols=["crystal_ifp_jaccard","crystal_centroid_distance_A","unique_pocket_residue_contacts",
              "pocket_atom_pair_contacts","ligand_atoms_contacting_pocket"]
    hub_cols=["tyr439_min_distance_A","tyr439_atom_pair_contacts"]
    anchor_cols=[f"res{r}_{suffix}" for r in (89,184,186,439) for suffix in ("min_A","pairs")]
    gate_cols=[f"res435_{suffix}" for suffix in ("min_A","pairs")]
    species_probe_cols=[f"res{r}_{suffix}" for r in (432,433) for suffix in ("min_A","pairs")]
    coupling_cols=anchor_cols+gate_cols+species_probe_cols
    residue_cols=[c for c in d if c.startswith("res") and (c.endswith("_min_A") or c.endswith("_pairs"))]
    struct_cols=base_cols+ifp_cols+hub_cols+residue_cols
    sx=d[struct_cols].to_numpy(float)
    methods={"VinaOnly":np.full(len(d),np.nan),"2D_RF":np.full(len(d),np.nan),
             "2D_LightGBM":np.full(len(d),np.nan),
             "IFP_RF":np.full(len(d),np.nan),"Hub_RF":np.full(len(d),np.nan),
             "Literature_Coupling_RF":np.full(len(d),np.nan),
             "Structure_RF":np.full(len(d),np.nan),"2D_plus_Structure_RF":np.full(len(d),np.nan),
             "2D_plus_Structure_LightGBM":np.full(len(d),np.nan)}
    folds=[]
    test_groups=[g for g,n in d.source_component.value_counts().items() if n>=MIN_N]
    for g in test_groups:
        te=groups==g; tr=~te
        selected=np.argsort(fp[tr].var(0))[-512:]
        mean,std=desc[tr].mean(0),desc[tr].std(0); std[std<1e-8]=1
        x2tr=np.hstack([fp[tr][:,selected],(desc[tr]-mean)/std]); x2te=np.hstack([fp[te][:,selected],(desc[te]-mean)/std])
        methods["VinaOnly"][te]=-d.loc[te,"vina_affinity"].to_numpy(float)
        configs={
            "2D_RF":(x2tr,x2te),
            "IFP_RF":(d.loc[tr,ifp_cols].to_numpy(float),d.loc[te,ifp_cols].to_numpy(float)),
            "Hub_RF":(d.loc[tr,hub_cols].to_numpy(float),d.loc[te,hub_cols].to_numpy(float)),
            "Literature_Coupling_RF":(d.loc[tr,coupling_cols].to_numpy(float),d.loc[te,coupling_cols].to_numpy(float)),
            "Structure_RF":(sx[tr],sx[te]),
            "2D_plus_Structure_RF":(np.hstack([x2tr,sx[tr]]),np.hstack([x2te,sx[te]])),
        }
        for name,(xtr,xte) in configs.items():
            model=make_pipeline(SimpleImputer(strategy="median"),StandardScaler(),
                                RandomForestRegressor(n_estimators=600,min_samples_leaf=2,random_state=42,n_jobs=-1))
            model.fit(xtr,y[tr]); methods[name][te]=model.predict(xte)
        for name,xtr,xte in [
            ("2D_LightGBM",x2tr,x2te),
            ("2D_plus_Structure_LightGBM",np.hstack([x2tr,sx[tr]]),np.hstack([x2te,sx[te]]))]:
            model=LGBMRegressor(n_estimators=500,learning_rate=.025,num_leaves=15,max_depth=5,
                min_child_samples=12,subsample=.8,colsample_bytree=.7,reg_alpha=.1,reg_lambda=1.0,
                random_state=42,n_jobs=4,verbosity=-1)
            model.fit(xtr,y[tr]); methods[name][te]=model.predict(xte)
        fold={"group":g,"n":int(te.sum()),"methods":{}}
        for name,pred in methods.items():
            rho=float(spearmanr(y[te],pred[te]).statistic)
            fold["methods"][name]={"Spearman":rho,"MAE":float(mean_absolute_error(y[te],pred[te]))}
        folds.append(fold); print(g,fold["methods"],flush=True)
    report={"state":STATE,"primary_metric":"macro within-series Spearman","n_molecules":len(d),"folds":folds,
            "aggregate":{name:aggregate([f["methods"][name] for f in folds]) for name in methods},
            "interpretation_rule":"Structure is promoted only if fusion improves macro and worst-series Spearman over its matched 2D baseline, especially LightGBM."}
    pred=d[["canonical_molecule_id","source_component","pEC50"]].copy()
    for k,v in methods.items():pred[k]=v
    pred.to_csv(OUT/"structure_loso_predictions.csv",index=False)
    (OUT/"structure_loso_metrics.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(report["aggregate"],ensure_ascii=False,indent=2))


if __name__=="__main__":main()
