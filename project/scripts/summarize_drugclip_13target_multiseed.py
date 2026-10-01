#!/usr/bin/env python3
"""Aggregate 13-target LOSO seeds and compute scaffold-cluster confidence intervals."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from finetune_drugclip_gpcr_screening import binary_metrics


METRICS=['roc_auc','pr_auc','bedroc_alpha20','ef1pct','ef5pct']


def summary(table,scores,targets):
    per={}
    for target in targets:
        keep=table.target.eq(target).to_numpy();per[target]=binary_metrics(table.loc[keep,'label'].to_numpy(int),scores[keep])
    return {'macro':{name:float(np.mean([per[t][name] for t in targets])) for name in METRICS},'per_target':per}


def main():
    p=argparse.ArgumentParser();p.add_argument('--prediction',type=Path,action='append',required=True);p.add_argument('--ecfp-predictions',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--bootstrap',type=int,default=1000);a=p.parse_args();frames=[pd.read_csv(v) for v in a.prediction];base=frames[0][['target','label','canonical_smiles','murcko_scaffold']].copy();targets=list(base.target.drop_duplicates());official=np.mean([v.official.to_numpy(float) for v in frames],axis=0);tuned_by=[v.tuned.to_numpy(float) for v in frames];random_by=[v.random.to_numpy(float) for v in frames];tuned=np.mean(tuned_by,axis=0);random=np.mean(random_by,axis=0);ecfp=None
    for frame in frames[1:]:
        if not frame[['target','label','canonical_smiles']].equals(base[['target','label','canonical_smiles']]):raise ValueError('Prediction row mismatch')
    result={'official':summary(base,official,targets),'tuned':summary(base,tuned,targets),'random':summary(base,random,targets)}
    if a.ecfp_predictions:
        ef=pd.read_csv(a.ecfp_predictions)
        if not ef[['target','label','canonical_smiles']].equals(base[['target','label','canonical_smiles']]):raise ValueError('ECFP row mismatch')
        ecfp=ef.ecfp4_logistic.to_numpy(float);result['ecfp4_logistic']=summary(base,ecfp,targets)
    groups={(t,s):np.asarray(list(idx),dtype=int) for (t,s),idx in base.groupby(['target','murcko_scaffold']).groups.items()};by_target={t:[key for key in groups if key[0]==t] for t in targets};rng=np.random.default_rng(20260926);comparisons=['tuned_minus_official','tuned_minus_random']+(['tuned_minus_ecfp4'] if ecfp is not None else []);delta={comparison:{m:[] for m in METRICS} for comparison in comparisons}
    for _ in range(a.bootstrap):
        idx=[]
        for target in targets:
            keys=by_target[target];idx.extend(np.concatenate([groups[keys[i]] for i in rng.integers(0,len(keys),len(keys))]))
        idx=np.asarray(idx,dtype=int);sample=base.iloc[idx].reset_index(drop=True);values={'official':summary(sample,official[idx],targets)['macro'],'tuned':summary(sample,tuned[idx],targets)['macro'],'random':summary(sample,random[idx],targets)['macro']}
        if ecfp is not None:values['ecfp4']=summary(sample,ecfp[idx],targets)['macro']
        for metric in METRICS:delta['tuned_minus_official'][metric].append(values['tuned'][metric]-values['official'][metric]);delta['tuned_minus_random'][metric].append(values['tuned'][metric]-values['random'][metric])
        if ecfp is not None:
            for metric in METRICS:delta['tuned_minus_ecfp4'][metric].append(values['tuned'][metric]-values['ecfp4'][metric])
    ci={c:{m:list(map(float,np.quantile(v,[.025,.5,.975]))) for m,v in block.items()} for c,block in delta.items()};directions={m:int(sum(result['tuned']['per_target'][t][m]>result['official']['per_target'][t][m] for t in targets)) for m in METRICS};report={'protocol':'Mean of independent 13-target LOSO seeds; target-stratified scaffold-cluster bootstrap','seeds':len(frames),'n_pairs':len(base),'n_targets':len(targets),'metrics':result,'per_seed_tuned_macro':[summary(base,v,targets)['macro'] for v in tuned_by],'improved_targets_out_of_13':directions,'bootstrap_95ci':ci,'claim_boundary':'Retrospective target-transfer validation; no efficacy, prospective, or SOTA claim.'};a.output.write_text(json.dumps(report,indent=2),encoding='utf-8');out=base.copy();out['official']=official;out['tuned']=tuned;out['random']=random;out.to_csv(a.output.with_suffix('.predictions.csv'),index=False);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
