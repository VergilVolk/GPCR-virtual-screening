#!/usr/bin/env python3
"""Build a leakage-audited human M1--M5 ChEMBL activity matrix."""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import pandas as pd
from rdkit import Chem

TARGETS={"M1":"CHEMBL216","M2":"CHEMBL211","M3":"CHEMBL245","M4":"CHEMBL1821","M5":"CHEMBL2035"}
ALLOWED_TYPES={"Ki","Kd","IC50","EC50"}

def canonical(s):
    m=Chem.MolFromSmiles(str(s)); return Chem.MolToSmiles(m,canonical=True) if m else None

def fetch(target,limit=1000):
    url="https://www.ebi.ac.uk/chembl/api/data/activity.json"; offset=0; rows=[]
    while True:
        params={"target_chembl_id":target,"pchembl_value__isnull":"false","standard_relation":"=",
                "limit":limit,"offset":offset}
        for attempt in range(5):
            try:
                with urlopen(url+"?"+urlencode(params),timeout=90) as response:
                    payload=json.load(response)
                break
            except Exception:
                if attempt==4: raise
                time.sleep(2**attempt)
        rows.extend(payload["activities"])
        if not payload["page_meta"]["next"]: break
        offset+=limit
    return rows

def exclusion_smiles(root:Path):
    files=[root/"pam_vs_inactive.csv",root/"potency_molecules.csv"]+list(root.glob("external_*/*.csv"))
    values=set()
    for path in files:
        try: frame=pd.read_csv(path)
        except Exception: continue
        for col in ["canonical_smiles","smiles"]:
            if col in frame:
                values.update(v for v in frame[col].dropna().map(canonical) if v)
    return values,files

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--benchmark-root",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True); args=ap.parse_args()
    records=[]; raw_counts={}
    for subtype,target in TARGETS.items():
        raw=fetch(target); raw_counts[subtype]=len(raw)
        for a in raw:
            if a.get("standard_type") not in ALLOWED_TYPES or a.get("data_validity_comment") is not None: continue
            smi=canonical(a.get("canonical_smiles")); value=a.get("pchembl_value")
            if smi is None or value is None: continue
            records.append({"subtype":subtype,"target_chembl_id":target,"canonical_smiles":smi,
                "molecule_chembl_id":a.get("molecule_chembl_id"),"pchembl_value":float(value),
                "standard_type":a.get("standard_type"),"assay_type":a.get("assay_type"),
                "document_chembl_id":a.get("document_chembl_id"),"document_year":a.get("document_year")})
    frame=pd.DataFrame(records); excluded,files=exclusion_smiles(args.benchmark_root)
    grouped=(frame.groupby(["subtype","canonical_smiles"],as_index=False)
             .agg(pchembl_median=("pchembl_value","median"),n_measurements=("pchembl_value","size"),
                  molecule_chembl_id=("molecule_chembl_id","first"),
                  first_year=("document_year","min"),last_year=("document_year","max")))
    grouped["overlaps_pacer_m4"]=grouped.canonical_smiles.isin(excluded)
    pivot=grouped.pivot(index="canonical_smiles",columns="subtype",values="pchembl_median")
    contrast=[]
    for smi,row in pivot.iterrows():
        for positive in ["M1","M2","M3","M5"]:  # M4 is held out for transfer
            if pd.isna(row.get(positive)) or row[positive]<6.0: continue
            for negative in ["M1","M2","M3","M5"]:
                if positive==negative or pd.isna(row.get(negative)) or row[negative]>5.0: continue
                contrast.append({"canonical_smiles":smi,"positive_subtype":positive,"negative_subtype":negative,
                                 "positive_pchembl":row[positive],"negative_pchembl":row[negative],
                                 "delta_pchembl":row[positive]-row[negative],"overlaps_pacer_m4":smi in excluded})
    contrast=pd.DataFrame(contrast)
    args.output.mkdir(parents=True,exist_ok=True); grouped.to_csv(args.output/"subtype_activity_matrix_long.csv",index=False)
    contrast.to_csv(args.output/"explicit_selectivity_contrasts.csv",index=False)
    audit={"source":"ChEMBL REST API","targets":TARGETS,"raw_counts":raw_counts,
           "n_curated_subtype_molecule_rows":len(grouped),"n_unique_molecules":grouped.canonical_smiles.nunique(),
           "n_pacer_overlap_rows":int(grouped.overlaps_pacer_m4.sum()),
           "n_explicit_non_m4_contrasts":len(contrast),
           "n_explicit_non_m4_contrasts_after_exclusion":int((~contrast.overlaps_pacer_m4).sum()) if len(contrast) else 0,
           "positive_threshold_pchembl":6.0,"negative_threshold_pchembl":5.0,
           "pacer_exclusion_files":[str(p) for p in files],
           "claim_boundary":"Heterogeneous ChEMBL binding/functional data; family pretraining only, not M4 PAM labels."}
    (args.output/"audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8"); print(json.dumps(audit,indent=2))
if __name__=="__main__": main()
