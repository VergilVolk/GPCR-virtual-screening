# -*- coding: utf-8 -*-
"""Enrich missing publication years from NCBI PubMed ESummary without overwriting frozen inputs."""
from __future__ import annotations
import json,time,subprocess
from urllib.parse import urlencode
from pathlib import Path
import pandas as pd

P=Path(__file__).resolve().parents[1];SRC=P/"data"/"benchmarks"/"m4_pam_v1"/"source_metadata.csv"
MOLS=P/"data"/"benchmarks"/"m4_pam_v1"/"potency_molecules.csv";OUT=P/"data"/"benchmarks"/"m4_pam_v1"/"temporal"

def main():
    OUT.mkdir(parents=True,exist_ok=True);d=pd.read_csv(SRC);pmids=[str(int(x)) for x in d.pmid.dropna().unique()]
    cache=OUT/"pubmed_esummary_cache.json";records=json.loads(cache.read_text(encoding="utf-8")) if cache.exists() else {}
    missing=[p for p in pmids if p not in records]
    for i in range(0,len(missing),10):
        ids=missing[i:i+10];url="https://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        query=url+"?"+urlencode({"db":"pubmed","id":",".join(ids),"retmode":"json","tool":"PACER_M4"})
        run=subprocess.run(["curl.exe","-fsSL","--retry","5","--connect-timeout","20",query],capture_output=True,text=True,encoding="utf-8",errors="replace")
        if run.returncode: raise RuntimeError(run.stderr[-1000:])
        res=json.loads(run.stdout)["result"]
        for pmid in ids:
            x=res.get(pmid,{});date=x.get("pubdate") or x.get("epubdate") or "";year=None
            for token in str(date).split():
                if token[:4].isdigit() and 1900<=int(token[:4])<=2100:year=int(token[:4]);break
            records[pmid]={"year":year,"pubdate":date,"title":x.get("title"),"fulljournalname":x.get("fulljournalname"),"doi":next((a.get("value") for a in x.get("articleids",[]) if a.get("idtype")=="doi"),None)}
        cache.write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding="utf-8")
        time.sleep(.35)
    enriched=d.copy();filled=0
    for idx,row in enriched.iterrows():
        if pd.isna(row.pmid):continue
        rec=records.get(str(int(row.pmid)),{});year=rec.get("year")
        if pd.isna(row.year_benchmark) and year:
            enriched.loc[idx,"year_benchmark"]=year;enriched.loc[idx,"year_provenance"]="NCBI_PubMed_ESummary";filled+=1
        if pd.isna(row.year) and year:enriched.loc[idx,"year"]=year
        if pd.isna(row.doi) and rec.get("doi"):enriched.loc[idx,"doi"]=rec["doi"]
    enriched.to_csv(OUT/"source_metadata_enriched.csv",index=False,encoding="utf-8-sig")
    mol=pd.read_csv(MOLS);year_map=dict(zip(enriched.source_id,enriched.year_benchmark));mol["publication_year_enriched"]=mol.primary_source_id.map(year_map)
    mol["temporal_split_enriched"]=pd.cut(mol.publication_year_enriched,bins=[0,2018,2021,9999],labels=["train_le_2018","validation_2019_2021","test_ge_2022"])
    mol.to_csv(OUT/"potency_molecules_temporal.csv",index=False,encoding="utf-8-sig")
    audit={"n_source_rows":len(d),"n_pmids_queried":len(pmids),"n_years_filled":filled,"remaining_missing_years":int(enriched.year_benchmark.isna().sum()),
           "potency_year_coverage":float(mol.publication_year_enriched.notna().mean()),"temporal_counts":mol.temporal_split_enriched.value_counts(dropna=False).to_dict(),
           "cutoffs_preregistered":{"train":"<=2018","validation":"2019-2021","test":">=2022"}}
    (OUT/"audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2,default=int),encoding="utf-8");print(json.dumps(audit,ensure_ascii=False,indent=2,default=int))

if __name__=="__main__":main()
