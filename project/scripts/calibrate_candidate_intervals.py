# -*- coding: utf-8 -*-
"""Attach honest retrospective source-LOSO residual intervals to candidates."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];OUT=P/"results"/"pacer_candidates_v01"/"final"
def main():
 o=pd.read_csv(P/"results"/"pacer_lightgbm_loso_v01"/"loso_predictions.csv").dropna(subset=["LightGBM"]).copy();o["abs_error"]=(o.pEC50-o.LightGBM).abs();rows=[]
 for g,x in o.groupby("source_component"):
  cal=o[o.source_component!=g];q=float(np.quantile(cal.abs_error,.9,method="higher"));cov=float((x.abs_error<=q).mean());rows.append({"group":g,"n":len(x),"q90":q,"coverage":cov})
 q=float(np.quantile(o.abs_error,.9,method="higher"));f=pd.read_csv(OUT/"final_candidate_hypotheses.csv");f["potency_reference_lower90"]=f.global_LGBM_pEC50_ref-q if "global_LGBM_pEC50_ref" in f else np.nan;f["potency_reference_upper90"]=f.global_LGBM_pEC50_ref+q if "global_LGBM_pEC50_ref" in f else np.nan
 # final selector intentionally omits global reference; recover it from audit table.
 if f.potency_reference_lower90.isna().all():
  a=pd.read_csv(OUT/"all_predock_with_structure.csv").set_index("candidate_id");p=np.asarray([a.loc[x,"global_LGBM_pEC50_ref"] for x in f.candidate_id]);f["global_LGBM_pEC50_ref"]=p;f["potency_reference_lower90"]=p-q;f["potency_reference_upper90"]=p+q
 f["interval_status"]="retrospective LOSO residual range; not prospective coverage guarantee";f.to_csv(OUT/"final_candidate_hypotheses.csv",index=False)
 report={"method":"absolute residual 90th percentile from source-LOSO LightGBM","global_q90_pEC50":q,"macro_group_coverage":float(np.mean([x["coverage"] for x in rows])),"worst_group_coverage":float(np.min([x["coverage"] for x in rows])),"folds":rows,"warning":"Source shift violates ordinary exchangeability; interval is a retrospective uncertainty scale, not a calibrated guarantee for generated molecules."};(OUT/"interval_audit.json").write_text(json.dumps(report,indent=2),encoding="utf-8");print(json.dumps({k:v for k,v in report.items() if k!="folds"},indent=2))
if __name__=="__main__":main()
