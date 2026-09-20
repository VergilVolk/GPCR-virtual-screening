# -*- coding: utf-8 -*-
"""Frozen LOSO test of 7TRQ/7TRP/7TRS M4 state-ensemble features."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem,Descriptors
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1]; DATA=P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv"
STATE_FILES={"q":P/"results"/"pacer_structure_loso_v01"/"docking_features.csv",
             "p":P/"results"/"pacer_structure_7trp_v01"/"docking_features.csv",
             "s":P/"results"/"pacer_structure_7trs_v01"/"docking_features.csv"}
OUT=P/"results"/"pacer_structure_ensemble_v01"
DESC=[Descriptors.MolWt,Descriptors.MolLogP,Descriptors.NumHDonors,Descriptors.NumHAcceptors,
      Descriptors.TPSA,Descriptors.NumRotatableBonds,Descriptors.NumAromaticRings,
      Descriptors.FractionCSP3,Descriptors.HeavyAtomCount]

def feat(smiles):
    fp=[]; ds=[]
    for s in smiles:
        m=Chem.MolFromSmiles(s); b=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048); a=np.zeros(2048,np.float32)
        DataStructs.ConvertToNumpyArray(b,a); fp.append(a); ds.append([f(m) for f in DESC])
    return np.asarray(fp),np.asarray(ds,np.float32)

def model(): return LGBMRegressor(n_estimators=500,learning_rate=.025,num_leaves=15,max_depth=5,min_child_samples=12,
    subsample=.8,colsample_bytree=.7,reg_alpha=.1,reg_lambda=1,random_state=42,n_jobs=4,verbosity=-1)

def main():
    d=pd.read_csv(DATA)
    state_cols={}
    for tag,path in STATE_FILES.items():
        f=pd.read_csv(path); f=f[(f.seed==42)&(f.error.fillna("")=="")].drop_duplicates("canonical_molecule_id")
        if len(f)<400: raise RuntimeError(f"{tag} incomplete n={len(f)}")
        cols=[c for c in f if c in {"vina_affinity","n_heavy_atoms","pocket_atom_pair_contacts",
            "ligand_atoms_contacting_pocket","unique_pocket_residue_contacts","crystal_ifp_jaccard",
            "crystal_centroid_distance_A","tyr439_min_distance_A","tyr439_atom_pair_contacts"}
            or c.startswith("res")]
        state_cols[tag]=[f"{tag}_{c}" for c in cols]
        f=f[["canonical_molecule_id"]+cols].rename(columns={c:f"{tag}_{c}" for c in cols})
        d=d.merge(f,on="canonical_molecule_id",how="inner")
    fp,ds=feat(d.canonical_smiles); y=d.pEC50.to_numpy(float); groups=d.source_component.astype(str).to_numpy()
    methods={k:np.full(len(d),np.nan) for k in ["2D","2D_7TRQ","2D_7TRP","2D_7TRS","2D_Ensemble"]}; folds=[]
    for g,n in d.source_component.value_counts().items():
        if n<8:continue
        te=groups==g;tr=~te;sel=np.argsort(fp[tr].var(0))[-1024:];mean,std=ds[tr].mean(0),ds[tr].std(0);std[std<1e-8]=1
        x2tr=np.hstack([fp[tr][:,sel],(ds[tr]-mean)/std]);x2te=np.hstack([fp[te][:,sel],(ds[te]-mean)/std])
        configs={"2D":(x2tr,x2te)}
        for tag,name in [("q","2D_7TRQ"),("p","2D_7TRP"),("s","2D_7TRS")]:
            sx=d[state_cols[tag]].to_numpy(float); configs[name]=(np.hstack([x2tr,sx[tr]]),np.hstack([x2te,sx[te]]))
        allcols=sum(state_cols.values(),[]); sx=d[allcols].to_numpy(float)
        configs["2D_Ensemble"]=(np.hstack([x2tr,sx[tr]]),np.hstack([x2te,sx[te]]))
        for name,(a,b) in configs.items():m=model();m.fit(a,y[tr]);methods[name][te]=m.predict(b)
        fold={"group":g,"n":int(te.sum()),"methods":{}}
        for name,pred in methods.items():fold["methods"][name]={"Spearman":float(spearmanr(y[te],pred[te]).statistic),"MAE":float(mean_absolute_error(y[te],pred[te]))}
        folds.append(fold);print(g,fold["methods"],flush=True)
    agg={}
    for name in methods:
        r=np.array([f["methods"][name]["Spearman"] for f in folds]);ma=np.array([f["methods"][name]["MAE"] for f in folds])
        agg[name]={"macro_Spearman":float(r.mean()),"median_Spearman":float(np.median(r)),"worst_Spearman":float(r.min()),
                   "positive_rho_series":int((r>0).sum()),"macro_MAE":float(ma.mean())}
    report={"primary":"macro within-series Spearman","folds":folds,"aggregate":agg,
            "promotion_rule":"Ensemble must improve macro and worst vs 2D and single states without selecting states on test labels."}
    OUT.mkdir(parents=True,exist_ok=True);(OUT/"ensemble_loso_metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(agg,indent=2))

if __name__=="__main__":main()
