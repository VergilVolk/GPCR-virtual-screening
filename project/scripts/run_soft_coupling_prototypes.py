"""Deterministic ordered soft prototypes for the frozen M4 coupling coordinate."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr


PROJECT=Path(__file__).resolve().parents[1]
OUT=PROJECT/'results'/'pacer_soft_coupling_prototypes_v01'
FRAME=PROJECT/'results'/'pacer_orthosteric_stability_annotation_v01'/'frame_ach_stability.csv'
QUANTILES=np.array([.10,.30,.50,.70,.90])


def membership(values,centers,bandwidth):
    logits=-.5*((values[:,None]-centers[None,:])/bandwidth)**2
    logits-=logits.max(1,keepdims=True); x=np.exp(logits)
    return x/x.sum(1,keepdims=True)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(FRAME); rows=[]; frame_rows=[]
    for heldout in range(1,7):
        train=df.replica!=heldout; test=df.replica==heldout
        centers=np.quantile(df.loc[train,'signature_projection'],QUANTILES)
        bandwidth=float(np.median(np.diff(centers)))
        if bandwidth<=0: raise RuntimeError(centers)
        mtrain=membership(df.loc[train,'signature_projection'].to_numpy(),centers,bandwidth)
        mtest=membership(df.loc[test,'signature_projection'].to_numpy(),centers,bandwidth)
        train_occ=mtrain.mean(0); test_occ=mtest.mean(0)
        entropy=-np.sum(test_occ[test_occ>0]*np.log(test_occ[test_occ>0]))
        continuity=float(np.mean(1-.5*np.abs(np.diff(mtest,axis=0)).sum(1)))
        rmsd=df.loc[test,'ACh_pose_RMSD_global_A'].to_numpy()
        weighted=np.array([(mtest[:,j]*rmsd).sum()/mtest[:,j].sum() for j in range(5)])
        monotonic=float(spearmanr(centers,-weighted).statistic)
        rows.append({'heldout_replica':heldout,'occupancy_JSD':float(jensenshannon(test_occ,train_occ,base=2)**2),
                     'soft_temporal_continuity':continuity,'effective_prototypes':float(np.exp(entropy)),
                     'max_occupancy':float(test_occ.max()),'center_vs_minus_ACh_RMSD_spearman':monotonic,
                     **{f'center_{j}':centers[j] for j in range(5)},
                     **{f'occupancy_{j}':test_occ[j] for j in range(5)},
                     **{f'ACh_RMSD_{j}':weighted[j] for j in range(5)}})
        globals_=df.loc[test,'global_frame'].to_numpy()
        for i,g in enumerate(globals_):
            frame_rows.append({'heldout_replica':heldout,'global_frame':int(g),
                               **{f'prototype_{j}_membership':float(mtest[i,j]) for j in range(5)}})
    folds=pd.DataFrame(rows); folds.to_csv(OUT/'fold_metrics.csv',index=False)
    pd.DataFrame(frame_rows).to_csv(OUT/'heldout_memberships.csv',index=False)
    summary={
        'mean_occupancy_JSD':float(folds.occupancy_JSD.mean()),
        'mean_soft_temporal_continuity':float(folds.soft_temporal_continuity.mean()),
        'mean_effective_prototypes':float(folds.effective_prototypes.mean()),
        'mean_max_occupancy':float(folds.max_occupancy.mean()),
        'positive_ACh_monotonicity_folds':int((folds.center_vs_minus_ACh_RMSD_spearman>0).sum()),
    }
    audit={**summary,'deterministic_head':True,'ACh_coordinates_used_to_fit_prototypes':False,
           'development_go':bool(summary['mean_occupancy_JSD']<=.30 and
             summary['mean_soft_temporal_continuity']>=.90 and summary['mean_effective_prototypes']>=3 and
             summary['mean_max_occupancy']<=.60 and summary['positive_ACh_monotonicity_folds']>=5),
           'functional_prediction_go':False,
           'claim_boundary':'Retrospective soft-state representation; not independent PAM efficacy validation.'}
    (OUT/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    (OUT/'REPORT.md').write_text('# PACER Soft Coupling Prototypes v0.1\n\n'+folds.to_markdown(index=False,floatfmt='.4f')+
       '\n\n```json\n'+json.dumps(audit,indent=2)+'\n```\n',encoding='utf-8')
    print(folds[['heldout_replica','occupancy_JSD','soft_temporal_continuity','effective_prototypes','max_occupancy','center_vs_minus_ACh_RMSD_spearman']].to_string(index=False))
    print(json.dumps(audit,indent=2))


if __name__=='__main__': main()
