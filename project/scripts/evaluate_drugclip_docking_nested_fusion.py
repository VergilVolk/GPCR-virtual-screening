#!/usr/bin/env python3
"""Nested target-LOSO fusion of tuned DrugCLIP and published docking scores."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from finetune_drugclip_gpcr_screening import TARGETS, screening_metrics, binary_metrics


def percentile(frame, column):
    out=np.zeros(len(frame),dtype=float)
    for target in TARGETS:
        keep=frame.target.eq(target).to_numpy();out[keep]=rankdata(frame.loc[keep,column].to_numpy(),method='average')/keep.sum()
    return out


def present_target_bedroc(frame, scores):
    values=[]
    target_values=frame.target.to_numpy()
    for target in frame.target.unique():
        keep=target_values==target
        values.append(binary_metrics(frame.loc[keep,'label'].to_numpy(int),scores[keep])['bedroc_alpha20'])
    return float(np.mean(values))


def main():
    p=argparse.ArgumentParser();p.add_argument('--drugclip-predictions',type=Path,required=True);p.add_argument('--docking-predictions',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    table=pd.read_csv(a.drugclip_predictions);dock=pd.read_csv(a.docking_predictions).pivot(index='pair_id',columns='method',values='score');table=table.join(dock,on='pair_id');methods=list(dock.columns);weights=[0.,.25,.5,.75,1.]
    table['drugclip_rank']=percentile(table,'tuned')
    for method in methods:table[f'{method}_rank']=percentile(table,method)
    fused=np.full(len(table),np.nan);selected_docking=np.full(len(table),np.nan);selection={}
    for held in TARGETS:
        train=~table.target.eq(held);test=~train;best=None
        for method in methods:
            d=table[f'{method}_rank'].to_numpy();z=table.drugclip_rank.to_numpy()
            for weight in weights:
                score=weight*z+(1-weight)*d
                metric=present_target_bedroc(table.loc[train].reset_index(drop=True),score[train])
                candidate=(metric,method,weight)
                if best is None or candidate>best:best=candidate
        metric,method,weight=best;d=table[f'{method}_rank'].to_numpy();z=table.drugclip_rank.to_numpy();fused[test]=weight*z[test]+(1-weight)*d[test];selected_docking[test]=d[test]
        selection[held]={'selected_on_other_targets':method,'drugclip_weight':weight,'training_macro_bedroc20':metric,'held_labels_used_for_selection':False}
    metrics={'official_drugclip':screening_metrics(table,table.official.to_numpy()),'target_loso_drugclip':screening_metrics(table,table.tuned.to_numpy()),'nested_selected_docking':screening_metrics(table,selected_docking),'nested_fusion':screening_metrics(table,fused)}
    report={'protocol':'Outer leave-one-target-out; docking method and fixed rank-fusion weight selected only on the other three targets by macro BEDROC20.','selection':selection,'metrics':metrics,'claim_boundary':'Retrospective nested fusion on four GPCRs; exploratory until replicated on additional targets. No PAM efficacy claim.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');out=table[['pair_id','target','label']].copy();out['nested_selected_docking']=selected_docking;out['nested_fusion']=fused;out.to_csv(a.output.with_suffix('.predictions.csv'),index=False);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
