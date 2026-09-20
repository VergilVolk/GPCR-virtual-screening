# -*- coding: utf-8 -*-
"""Label-free uncertainty and selective prediction on frozen M4 LOSO outputs."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

RDLogger.DisableLog("rdApp.*")
P=Path(__file__).resolve().parents[1]; OUT=P/"results"/"pacer_selective_v01"; OUT.mkdir(parents=True,exist_ok=True)
FILES={
"base":P/"results"/"pacer_lightgbm_loso_v01"/"loso_predictions.csv",
"classic":P/"results"/"pacer_potency_loso_v01"/"loso_predictions.csv",
"bert":P/"results"/"chemberta_loso_v1"/"loso_predictions.csv",
"chemprop":P/"results"/"pacer_chemprop_loso_v01"/"loso_predictions.csv",
"Q":P/"results"/"pacer_structure_loso_v01"/"structure_loso_predictions.csv",
"P":P/"results"/"pacer_structure_7trp_v01"/"structure_loso_predictions.csv",
"S":P/"results"/"pacer_structure_7trs_v01"/"structure_loso_predictions.csv"}
COVERAGES=(0.5,0.7,0.8,0.9,1.0)

def z(x):
    x=np.asarray(x,float);s=x.std();return (x-x.mean())/(s if s>1e-8 else 1)

def ranks01(x):
    return pd.Series(x).rank(method="average",pct=True).to_numpy(float)

def fps(smiles):return [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s),2,nBits=2048) for s in smiles]

def main():
    d=pd.read_csv(FILES["base"])[["canonical_molecule_id","canonical_smiles","pEC50","source_component","LightGBM"]]
    c=pd.read_csv(FILES["classic"])[["canonical_molecule_id","RandomForest","ExtraTrees","NeuralERM"]];d=d.merge(c,on="canonical_molecule_id")
    d=d.merge(pd.read_csv(FILES["bert"])[["canonical_molecule_id","ChemBERTaRidge"]],on="canonical_molecule_id")
    d=d.merge(pd.read_csv(FILES["chemprop"])[["canonical_molecule_id","Chemprop_DMPNN"]],on="canonical_molecule_id")
    for tag in "QPS":
        x=pd.read_csv(FILES[tag])[["canonical_molecule_id","IFP_RF"]].rename(columns={"IFP_RF":f"{tag}_IFP"});d=d.merge(x,on="canonical_molecule_id")
    d=d[d.LightGBM.notna()].reset_index(drop=True); allfp=fps(d.canonical_smiles)
    fold_rows=[]; mol_rows=[]
    for group,x in d.groupby("source_component"):
        idx=x.index.to_numpy();train=np.where(d.source_component.to_numpy()!=group)[0]
        model_cols=["LightGBM","RandomForest","ExtraTrees","ChemBERTaRidge","Chemprop_DMPNN"]
        u_model=np.std(np.vstack([z(x[c]) for c in model_cols]),axis=0)
        u_state=np.std(np.vstack([z(x[c]) for c in ["Q_IFP","P_IFP","S_IFP"]]),axis=0)
        maxsim=[]
        trfp=[allfp[i] for i in train]
        for i in idx:maxsim.append(max(DataStructs.BulkTanimotoSimilarity(allfp[i],trfp)))
        u_ad=1-np.asarray(maxsim)
        uncertainties={"ModelDisagreement":u_model,"StateDisagreement":u_state,"ApplicabilityDistance":u_ad,
                       "CombinedRank":(ranks01(u_model)+ranks01(u_state)+ranks01(u_ad))/3}
        y=x.pEC50.to_numpy(float);pred=x.LightGBM.to_numpy(float);err=np.abs(y-pred)
        fold={"group":group,"n":len(x),"uncertainties":{}}
        for name,u in uncertainties.items():
            vals={"error_correlation":float(spearmanr(u,err).statistic),"coverage":{}}
            order=np.argsort(u)
            for cov in COVERAGES:
                n=max(5,int(np.ceil(cov*len(x))));sel=order[:n]
                vals["coverage"][str(cov)]={"n":int(n),"Spearman":float(spearmanr(y[sel],pred[sel]).statistic),
                                              "MAE":float(mean_absolute_error(y[sel],pred[sel]))}
            fold["uncertainties"][name]=vals
        fold_rows.append(fold)
        for j,di in enumerate(idx):mol_rows.append({"canonical_molecule_id":d.loc[di,"canonical_molecule_id"],"source_component":group,
            "pEC50":y[j],"prediction":pred[j],"abs_error":err[j],**{k:v[j] for k,v in uncertainties.items()}})
    aggregate={}
    for name in fold_rows[0]["uncertainties"]:
        aggregate[name]={"macro_error_correlation":float(np.nanmean([f["uncertainties"][name]["error_correlation"] for f in fold_rows])),"coverage":{}}
        for cov in COVERAGES:
            vals=[f["uncertainties"][name]["coverage"][str(cov)] for f in fold_rows]
            aggregate[name]["coverage"][str(cov)]={"macro_Spearman":float(np.nanmean([v["Spearman"] for v in vals])),
                "worst_Spearman":float(np.nanmin([v["Spearman"] for v in vals])),"macro_MAE":float(np.mean([v["MAE"] for v in vals]))}
    report={"base_model":"LightGBM","uncertainty_is_label_free":True,"folds":fold_rows,"aggregate":aggregate,
            "promotion_rule":"Uncertainty must show positive error correlation and improve macro/worst risk as coverage decreases."}
    (OUT/"selective_metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8");pd.DataFrame(mol_rows).to_csv(OUT/"selective_predictions.csv",index=False)
    print(json.dumps(aggregate,indent=2))

if __name__=="__main__":main()
