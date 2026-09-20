# -*- coding: utf-8 -*-
"""Test whether M4 PAM activity cliffs coincide with receptor-state IFP shifts."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

P=Path(__file__).resolve().parents[1]; PAIRS=P/"data"/"benchmarks"/"m4_pam_v1"/"activity_cliffs"/"cliff_vs_smooth.csv"
MOLS=P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv";OUT=P/"results"/"pacer_cliff_state_switch_v01";OUT.mkdir(parents=True,exist_ok=True)
FILES={"Q":P/"results"/"pacer_structure_loso_v01"/"docking_features.csv",
"P":P/"results"/"pacer_structure_7trp_v01"/"docking_features.csv",
"S":P/"results"/"pacer_structure_7trs_v01"/"docking_features.csv"}

def main():
    mol=pd.read_csv(MOLS)[["canonical_molecule_id","source_component"]]; state={}
    summary_cols=["crystal_ifp_jaccard","crystal_centroid_distance_A","unique_pocket_residue_contacts"]
    residue_cols=[f"res{r}_min_A" for r in (89,92,93,96,184,186,190,423,432,433,435,436,439)]
    for tag,path in FILES.items():
        x=pd.read_csv(path);x=x[(x.seed==42)&(x.error.fillna("")=="")].drop_duplicates("canonical_molecule_id").set_index("canonical_molecule_id")
        a=x[summary_cols+residue_cols].to_numpy(float);a[:,3:]=np.clip(a[:,3:],0,10);mean,std=a.mean(0),a.std(0);std[std<1e-8]=1
        state[tag]={k:v for k,v in zip(x.index,(a-mean)/std)}
    pairs=pd.read_csv(PAIRS);source=dict(zip(mol.canonical_molecule_id,mol.source_component));rows=[]
    for r in pairs.itertuples():
        feats={}
        allvec=[]
        ifps=[]
        for tag in "QPS":
            a,b=state[tag][r.mol_a],state[tag][r.mol_b];delta=np.abs(a-b);allvec.extend(delta);ifps.append(delta[0])
            feats[f"{tag}_summary_distance"]=float(np.linalg.norm(delta[:3]));feats[f"{tag}_residue_distance"]=float(np.linalg.norm(delta[3:]))
        feats["ensemble_state_distance"]=float(np.linalg.norm(allvec));feats["ifp_shift_distance"]=float(np.linalg.norm(ifps))
        raw_a=np.array([state[t][r.mol_a][0] for t in "QPS"]);raw_b=np.array([state[t][r.mol_b][0] for t in "QPS"])
        feats["preferred_state_switch"]=int(np.argmax(raw_a)!=np.argmax(raw_b))
        rows.append({"pair_id":r.pair_id,"label":int(r.is_cliff),"tanimoto":r.tanimoto,"source_pair_fold":r.source_pair_fold,
                     "same_source_component":r.same_source_component,"source_component":source.get(r.mol_a) if r.same_source_component else "cross",**feats})
    d=pd.DataFrame(rows);feature_cols=[c for c in d if c.endswith("distance") or c=="preferred_state_switch"]
    simple={}
    for subset,mask in {"all":np.ones(len(d),bool),"same_source":d.same_source_component.eq(1).to_numpy()}.items():
        y=d.label.to_numpy()[mask];simple[subset]={"n":int(mask.sum()),"n_cliff":int(y.sum()),"metrics":{}}
        for c in ["tanimoto"]+feature_cols:
            score=-d[c].to_numpy()[mask] if c=="tanimoto" else d[c].to_numpy()[mask]
            simple[subset]["metrics"][c]={"ROC_AUC":float(roc_auc_score(y,score)),"PR_AUC":float(average_precision_score(y,score))}
    cv={}
    valid=d.source_pair_fold>=0;y=d.label.to_numpy();fold=d.source_pair_fold.to_numpy()
    configs={"SimilarityOnly":["tanimoto"],"StateOnly":feature_cols,"SimilarityPlusState":["tanimoto"]+feature_cols}
    for name,cols in configs.items():
        pred=np.full(len(d),np.nan);reports=[]
        for f in sorted(set(fold[valid])):
            te=(fold==f)&valid;tr=(fold!=f)&valid
            model=make_pipeline(StandardScaler(),LogisticRegression(C=1,class_weight="balanced",max_iter=2000))
            model.fit(d.loc[tr,cols],y[tr]);pred[te]=model.predict_proba(d.loc[te,cols])[:,1]
            reports.append({"fold":int(f),"n":int(te.sum()),"ROC_AUC":float(roc_auc_score(y[te],pred[te])),"PR_AUC":float(average_precision_score(y[te],pred[te]))})
        cv[name]={"folds":reports,"macro_ROC_AUC":float(np.mean([r["ROC_AUC"] for r in reports])),
                  "worst_ROC_AUC":float(np.min([r["ROC_AUC"] for r in reports])),"pooled_ROC_AUC_reference_only":float(roc_auc_score(y[valid],pred[valid]))}
        d[f"pred_{name}"]=pred
    report={"question":"Do functional activity cliffs coincide with receptor-state contact rearrangement?","simple":simple,"source_fold_cv":cv,
            "rule":"State switching is supported only if StateOnly exceeds chance and SimilarityPlusState improves macro and worst over SimilarityOnly."}
    d.to_csv(OUT/"pair_state_features.csv",index=False);(OUT/"metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
