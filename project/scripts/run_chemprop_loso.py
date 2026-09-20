# -*- coding: utf-8 -*-
"""Chemprop v2 D-MPNN baseline using the frozen leave-one-series-out split."""
from __future__ import annotations
import json, os, shutil, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

PROJECT=Path(__file__).resolve().parents[1]
DATA=PROJECT/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv"
OUT=PROJECT/"results"/"pacer_chemprop_loso_v01"

def main():
    OUT.mkdir(parents=True,exist_ok=True); d=pd.read_csv(DATA); pred=np.full(len(d),np.nan); folds=[]
    groups=d.source_component.astype(str).to_numpy(); test_groups=[g for g,n in d.source_component.value_counts().items() if n>=8]
    for fold_idx,g in enumerate(test_groups):
        fd=OUT/g; fd.mkdir(exist_ok=True); pred_file=fd/"predictions.csv"
        te=groups==g
        if not pred_file.exists():
            if (fd/"model").exists(): shutil.rmtree(fd/"model")
            work=pd.DataFrame({"smiles":d.canonical_smiles,"pEC50":d.pEC50})
            split=np.full(len(d),"train",object); split[te]="test"
            # Validation is a full scaffold fold from the remaining data.
            val=(d.scaffold_fold.to_numpy()==fold_idx%5)&(~te); split[val]="val"; work["split"]=split
            inp=fd/"data.csv"; work.to_csv(inp,index=False)
            cmd=["chemprop","train","-i",str(inp),"-s","smiles","--target-columns","pEC50",
                 "--splits-column","split","-t","regression","--epochs","40","--patience","7",
                 "--accelerator","cpu","--devices","1","-n","0","-b","64","--message-hidden-dim","128",
                 "--ffn-hidden-dim","128","--depth","3","--dropout","0.1","--pytorch-seed","42",
                 "-o",str(fd/"model"),"-q"]
            env=dict(os.environ); env["PYTHONUTF8"]="1"; env["PYTHONIOENCODING"]="utf-8"
            run=subprocess.run(cmd,capture_output=True,text=True,env=env,encoding="utf-8",errors="replace")
            (fd/"train_stdout.txt").write_text(run.stdout+"\nSTDERR\n"+run.stderr,encoding="utf-8")
            if run.returncode: raise RuntimeError(f"Chemprop fold {g} failed; see {fd/'train_stdout.txt'}")
            ckpts=list((fd/"model").rglob("*.ckpt"))
            if not ckpts: raise RuntimeError(f"No checkpoint for {g}")
            test=work.loc[te,["smiles","pEC50"]]; test.to_csv(fd/"test.csv",index=False)
            run=subprocess.run(["chemprop","predict","-i",str(fd/"test.csv"),"-s","smiles",
                "--model-paths",str(ckpts[0]),"-o",str(pred_file),"--accelerator","cpu","--devices","1","-n","0","-q"],
                capture_output=True,text=True,env=env,encoding="utf-8",errors="replace")
            if run.returncode: raise RuntimeError(run.stderr[-1000:])
        p=pd.read_csv(pred_file)
        pred_cols=[c for c in p if c.startswith("pred_")]
        if not pred_cols: raise RuntimeError(f"No Chemprop prediction column in {pred_file}")
        pred[te]=p[pred_cols[0]].to_numpy(float)
        rho=float(spearmanr(d.loc[te,"pEC50"],pred[te]).statistic); mae=float(mean_absolute_error(d.loc[te,"pEC50"],pred[te]))
        folds.append({"group":g,"n":int(te.sum()),"Spearman":rho,"MAE":mae}); print(g,rho,mae,flush=True)
    r=np.array([x["Spearman"] for x in folds]); m=np.array([x["MAE"] for x in folds])
    report={"model":"Chemprop 2.3.1 D-MPNN","folds":folds,"aggregate":{"macro_Spearman":float(r.mean()),
        "median_Spearman":float(np.median(r)),"worst_Spearman":float(r.min()),"positive_rho_series":int((r>0).sum()),
        "n_series":len(r),"macro_MAE":float(m.mean())}}
    out=d[["canonical_molecule_id","canonical_smiles","pEC50","source_component"]].copy(); out["Chemprop_DMPNN"]=pred
    out.to_csv(OUT/"loso_predictions.csv",index=False); (OUT/"loso_metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report["aggregate"],indent=2))

if __name__=="__main__":main()
