"""Replica-aware analysis of the published iperoxo RMSD source data."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "data" / "papers" / "Wang2022_M4_SourceData.xlsx"
OUT = PROJECT / "results" / "pacer_public_pam_dynamic_endpoint_v01"
BLOCK = 250  # 50 ns at 0.2 ns/sample
RNG = np.random.default_rng(20260901)


def exact_permutation(control, pam):
    values=np.r_[control,pam]; observed=control.mean()-pam.mean(); null=[]
    for idx in itertools.combinations(range(6),3):
        a=np.asarray(idx); b=np.asarray([i for i in range(6) if i not in idx])
        null.append(values[a].mean()-values[b].mean())
    null=np.asarray(null)
    return observed,float(np.sum(null>=observed)/len(null)),null


def hierarchical_bootstrap(control, pam, n=10000):
    arrays=[control,pam]; effects=[]
    for _ in range(n):
        cond=[]
        for group in arrays:
            chosen=RNG.integers(0,3,3); means=[]
            for r in chosen:
                blocks=group[r].reshape(-1,BLOCK)
                draw=blocks[RNG.integers(0,len(blocks),len(blocks))]
                means.append(draw.mean())
            cond.append(np.mean(means))
        effects.append(cond[0]-cond[1])
    return np.asarray(effects)


def markdown_table(df):
    cols=list(df.columns)
    lines=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
    for row in df.itertuples(index=False,name=None):
        lines.append('| '+' | '.join(f'{v:.4f}' if isinstance(v,(float,np.floating)) else str(v) for v in row)+' |')
    return '\n'.join(lines)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    data=pd.read_excel(SOURCE,sheet_name='Sheet1')
    expected=[f'iperoxo-sim{i}' for i in range(1,4)]+[f'iperoxo-LY-sim{i}' for i in range(1,4)]
    if list(data.columns)!=expected or data.shape!=(5000,6) or data.isna().any().any():
        raise RuntimeError((data.shape,list(data.columns),int(data.isna().sum().sum())))
    control=data.iloc[:,:3].to_numpy(float).T; pam=data.iloc[:,3:].to_numpy(float).T
    rows=[]
    for condition,arr in [('iperoxo_only',control),('iperoxo_LY2119620',pam)]:
        for i,x in enumerate(arr,1): rows.append({'condition':condition,'replica':i,'mean_RMSD_A':x.mean(),
            'median_RMSD_A':np.median(x),'sd_RMSD_A':x.std(ddof=1),'post100ns_mean_RMSD_A':x[500:].mean()})
    reps=pd.DataFrame(rows); reps.to_csv(OUT/'replica_summary.csv',index=False)
    cmean=np.array([x.mean() for x in control]); pmean=np.array([x.mean() for x in pam])
    effect,p_exact,_=exact_permutation(cmean,pmean)
    boot=hierarchical_bootstrap(control,pam)
    post_effect=control[:,500:].mean(axis=1).mean()-pam[:,500:].mean(axis=1).mean()
    paired=cmean-pmean
    audit={'source_rows':5000,'replicas_per_condition':3,'trajectory_is_independent_unit':True,
           'mean_RMSD_reduction_A':float(effect),'paired_replica_differences_A':paired.tolist(),
           'all_three_paired_positive':bool(np.all(paired>0)),'exact_one_sided_p':p_exact,
           'hierarchical_block_bootstrap_CI95_A':[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))],
           'post100ns_effect_A':float(post_effect),'support':bool(np.all(paired>0) and p_exact<=.05 and
             np.quantile(boot,.025)>0 and post_effect>0),
           'claim_boundary':'One published LY2119620 system; dynamic orthosteric stabilization, not general PAM efficacy prediction.'}
    (OUT/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    # Flat, auditable processed time series.
    long=data.assign(sample=np.arange(5000),time_ns=np.arange(5000)*.2).melt(
        id_vars=['sample','time_ns'],var_name='trajectory',value_name='iperoxo_RMSD_A')
    long.to_csv(OUT/'processed_rmsd_timeseries.csv',index=False)
    report='# PACER Public PAM Dynamic Endpoint v0.1\n\n'+markdown_table(reps)
    report+='\n\n```json\n'+json.dumps(audit,indent=2)+'\n```\n'
    (OUT/'REPORT.md').write_text(report,encoding='utf-8')
    print(reps.to_string(index=False));print(json.dumps(audit,indent=2))


if __name__=='__main__': main()
