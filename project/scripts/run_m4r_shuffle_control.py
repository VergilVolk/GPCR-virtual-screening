# -*- coding: utf-8 -*-
"""Training-label shuffle control for the scaffold-held-out M4R classifier."""
import json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
from lightgbm import LGBMClassifier
from sklearn.metrics import roc_auc_score
P=Path(__file__).resolve().parents[1];B=P/"data"/"benchmarks"/"m4r_enrichment_v1";O=P/"results"/"pacer_m4r_enrichment_v01";O.mkdir(parents=True,exist_ok=True)
def main():
 d=pd.read_csv(B/"molecules.csv");X=sparse.load_npz(B/"feature_cache"/"ecfp4_2048.npz");D=np.load(B/"feature_cache"/"descriptors.npy");X=sparse.hstack([X,sparse.csr_matrix(D)],format="csr");y=d.is_active.to_numpy(int);fold=d.scaffold_fold.to_numpy(int);rows=[]
 for seed in (1000,2000,3000):
  for f in sorted(np.unique(fold)):
   tr=fold!=f;te=~tr;rng=np.random.default_rng(seed+f);ys=rng.permutation(y[tr]);m=LGBMClassifier(n_estimators=120,learning_rate=.06,num_leaves=31,max_depth=7,min_child_samples=50,class_weight="balanced",random_state=seed+f,n_jobs=6,verbosity=-1).fit(X[tr],ys);s=m.predict_proba(X[te])[:,1];order=np.argsort(s)[::-1];k=int(np.ceil(.005*te.sum()));prev=y[te].mean();ef=(y[te][order[:k]].mean()/prev);rows.append({"seed":seed,"fold":int(f),"ROC_AUC":float(roc_auc_score(y[te],s)),"EF0.5%":float(ef)})
 by=[]
 for seed,x in pd.DataFrame(rows).groupby("seed"):by.append({"seed":int(seed),"macro_ROC_AUC":float(x.ROC_AUC.mean()),"macro_EF0.5%":float(x["EF0.5%"].mean())})
 report={"control":"shuffle training labels independently in every outer scaffold fold","seeds":by,"folds":rows,"mean_macro_ROC_AUC":float(np.mean([x["macro_ROC_AUC"] for x in by])),"sd_macro_ROC_AUC":float(np.std([x["macro_ROC_AUC"] for x in by],ddof=1)),"mean_macro_EF0.5%":float(np.mean([x["macro_EF0.5%"] for x in by])),"sd_macro_EF0.5%":float(np.std([x["macro_EF0.5%"] for x in by],ddof=1))};(O/"shuffle_control.json").write_text(json.dumps(report,indent=2),encoding="utf-8");print(json.dumps({k:v for k,v in report.items() if k!="folds"},indent=2))
if __name__=="__main__":main()
