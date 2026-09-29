#!/usr/bin/env python3
"""Build a leakage-controlled M1--M5 activity benchmark for DrugCLIP."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def scaffold(s):
    m=Chem.MolFromSmiles(s); v=MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False)
    return v or 'ACYCLIC:'+hashlib.sha1(s.encode()).hexdigest()[:12]

def assign(frame,n):
    g=frame.drop_duplicates('canonical_smiles').groupby('murcko_scaffold').size().sort_values(ascending=False)
    loads=np.zeros(n,int); out={}
    for name,size in g.items():
        f=int(np.argmin(loads)); out[name]=f; loads[f]+=int(size)
    return out

def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',type=Path,required=True); p.add_argument('--output-dir',type=Path,required=True); p.add_argument('--folds',type=int,default=5); a=p.parse_args()
    raw=pd.read_csv(a.input); d=raw[~raw.overlaps_pacer_m4.astype(bool)].copy()
    d=d[(d.pchembl_median>=6)|(d.pchembl_median<=5)].copy(); d['target']=(d.pchembl_median>=6).astype(int)
    d['murcko_scaffold']=d.canonical_smiles.map(scaffold); amap=assign(d,a.folds); d['scaffold_fold']=d.murcko_scaffold.map(amap).astype(int)
    unique=d[['canonical_smiles','murcko_scaffold','scaffold_fold']].drop_duplicates('canonical_smiles').reset_index(drop=True)
    unique.insert(0,'canonical_molecule_id',[f'MUSCDOM{i:05d}' for i in range(len(unique))])
    d=d.merge(unique[['canonical_molecule_id','canonical_smiles']],on='canonical_smiles',how='left')
    cols=['canonical_molecule_id','canonical_smiles','subtype','target','pchembl_median','n_measurements','murcko_scaffold','scaffold_fold']
    a.output_dir.mkdir(parents=True,exist_ok=True); d[cols].to_csv(a.output_dir/'pairs.csv',index=False); unique.to_csv(a.output_dir/'molecules.csv',index=False)
    audit={'raw_rows':len(raw),'pacer_overlap_rows_excluded':int(raw.overlaps_pacer_m4.astype(bool).sum()),'usable_pairs':len(d),'unique_molecules':len(unique),'positives':int(d.target.sum()),'negatives':int((d.target==0).sum()),'by_subtype_label':d.groupby(['subtype','target']).size().rename('n').reset_index().to_dict('records'),'claim_boundary':'Heterogeneous ChEMBL M1-M5 activity domain adaptation; not allosteric or PAM-specific.'}
    (a.output_dir/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8'); print(json.dumps(audit,indent=2))
if __name__=='__main__': main()
