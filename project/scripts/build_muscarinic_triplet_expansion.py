#!/usr/bin/env python3
"""Build expanded same-molecule M1/M2/M3/M5 activity-ranking triplets.

The original strict contrast molecules are frozen as an external test set and
are removed entirely from the expanded training set.
"""
from __future__ import annotations
import argparse, itertools, json
from pathlib import Path
import pandas as pd

def main():
    p=argparse.ArgumentParser(); p.add_argument('--activity',type=Path,required=True); p.add_argument('--strict',type=Path,required=True); p.add_argument('--output-dir',type=Path,required=True); p.add_argument('--min-positive',type=float,default=6.0); p.add_argument('--min-delta',type=float,default=1.0); a=p.parse_args()
    d=pd.read_csv(a.activity); d=d[(~d.overlaps_pacer_m4.astype(bool)) & d.subtype.isin(['M1','M2','M3','M5'])].copy(); rows=[]
    for smi,g in d.groupby('canonical_smiles',sort=False):
        for pos,neg in itertools.permutations(g.itertuples(),2):
            delta=float(pos.pchembl_median-neg.pchembl_median)
            if pos.pchembl_median>=a.min_positive and delta>=a.min_delta:
                rows.append({'canonical_smiles':smi,'positive_subtype':pos.subtype,'negative_subtype':neg.subtype,'positive_pchembl':pos.pchembl_median,'negative_pchembl':neg.pchembl_median,'delta_pchembl':delta})
    expanded=pd.DataFrame(rows); strict=pd.read_csv(a.strict); strict=strict[~strict.overlaps_pacer_m4.astype(bool)].copy(); test_molecules=set(strict.canonical_smiles); train=expanded[~expanded.canonical_smiles.isin(test_molecules)].copy(); molecules=expanded[['canonical_smiles']].drop_duplicates().reset_index(drop=True); molecules.insert(0,'canonical_molecule_id',[f'MUSCTRI{i:05d}' for i in range(len(molecules))])
    a.output_dir.mkdir(parents=True,exist_ok=True); train.to_csv(a.output_dir/'train_triplets.csv',index=False); strict.to_csv(a.output_dir/'strict_external_test.csv',index=False); molecules.to_csv(a.output_dir/'molecules.csv',index=False)
    audit={'expanded_pairs_before_test_exclusion':len(expanded),'expanded_molecules':expanded.canonical_smiles.nunique(),'strict_external_pairs':len(strict),'strict_external_molecules':len(test_molecules),'training_pairs':len(train),'training_molecules':train.canonical_smiles.nunique(),'test_molecule_overlap_with_training':len(test_molecules & set(train.canonical_smiles)),'rules':f'Non-M4 subtypes only; positive pChEMBL >= {a.min_positive}; within-molecule delta >= {a.min_delta}; every strict-test molecule removed from training.','claim_boundary':'Heterogeneous ChEMBL activity ranking; not allosteric or PAM-specific.'}
    (a.output_dir/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8'); print(json.dumps(audit,indent=2))
if __name__=='__main__': main()
