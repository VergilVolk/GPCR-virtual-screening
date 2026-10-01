#!/usr/bin/env python3
"""Frozen external endpoint audit for large-scale DrugCLIP projection tuning."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

from evaluate_drugclip_m4_external_panel import retrieval, ci_bootstrap, paired_delta_ci


def main():
    p=argparse.ArgumentParser();p.add_argument('--representations',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--official-projection',type=Path,required=True);p.add_argument('--checkpoint',action='append',nargs=2,metavar=('NAME','PATH'),required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--bootstrap',type=int,default=3000);a=p.parse_args()
    arc=np.load(a.representations,allow_pickle=False);ids=list(map(str,arc['molecule_ids']));idx={v:i for i,v in enumerate(ids)};table=pd.read_csv(a.manifest)
    keys=[str(mid) if str(mid) in idx else str(smi) for mid,smi in zip(table.external_molecule_id,table.canonical_smiles)];order=np.asarray([idx[v] for v in keys]);mol=arc['molecule_representations'].astype(np.float32)[order];pocket=arc['pocket_representations'].astype(np.float32)
    methods={'official':torch.load(a.official_projection,map_location='cpu')}
    for name,path in a.checkpoint:methods[name]=torch.load(path,map_location='cpu')
    for name,state in methods.items():table[name]=retrieval(mol,pocket,state)
    eligible=table[table.eligible.astype(bool)&~table.exact_training_overlap_computed.astype(bool)].copy();cols=['label_binary','potency_value','ordinal_class',*methods];grouped=eligible.groupby(['dataset','endpoint','external_molecule_id'],as_index=False)[cols].mean();reports={}
    for dataset,frame in grouped.groupby('dataset',sort=False):
        binary=frame.label_binary.notna().all();y=frame.label_binary.to_numpy(float) if binary else (frame.potency_value.to_numpy(float) if frame.potency_value.notna().all() else frame.ordinal_class.to_numpy(float));block={'endpoint':frame.endpoint.iloc[0],'n':len(frame),'outcome':'binary' if binary else 'ordered','methods':{},'paired_vs_official':{}}
        for j,name in enumerate(methods):
            score=frame[name].to_numpy(float)
            if binary:block['methods'][name]={'roc_auc':float(roc_auc_score(y,score)),'average_precision':float(average_precision_score(y,score)),'roc_auc_bootstrap_95ci':ci_bootstrap(y,score,'auc',20260926+j,a.bootstrap)}
            else:block['methods'][name]={'spearman':float(spearmanr(y,score).statistic),'spearman_bootstrap_95ci':ci_bootstrap(y,score,'spearman',20260926+j,a.bootstrap)}
            if name!='official':block['paired_vs_official'][name]=paired_delta_ci(y,score,frame.official.to_numpy(float),'auc' if binary else 'spearman',20261026+j,a.bootstrap)
        reports[dataset]=block
    report={'method_freeze':'Checkpoints selected without external endpoint labels. Exact training overlaps excluded.','n_endpoint_rows':len(table),'n_zero_shot_rows':len(eligible),'n_unique_zero_shot_units':len(grouped),'checkpoints':{name:path for name,path in a.checkpoint},'datasets':reports,'claim_boundary':'Cross-source retrospective stress test; no prospective or wet-lab PAM validation.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');table.to_csv(a.output.with_suffix('.predictions.csv'),index=False);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
