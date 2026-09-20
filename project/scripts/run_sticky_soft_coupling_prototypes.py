"""Autocorrelation-calibrated sticky ordered prototypes for M4 trajectories."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr


PROJECT=Path(__file__).resolve().parents[1]
OUT=PROJECT/'results'/'pacer_sticky_soft_coupling_prototypes_v02'
FRAME=PROJECT/'results'/'pacer_orthosteric_stability_annotation_v01'/'frame_ach_stability.csv'
Q=np.array([.10,.30,.50,.70,.90]); K=5


def autocorrelation_tau(x,maxlag=100):
    x=x-x.mean(); denom=float(x@x)
    if denom<=0:return 1
    for lag in range(1,min(maxlag,len(x)-1)+1):
        if float(x[:-lag]@x[lag:]/denom)<=np.exp(-1): return lag
    return maxlag


def transition(tau):
    stay=float(np.exp(-1/max(tau,1))); a=np.zeros((K,K))
    for i in range(K):
        a[i,i]=stay; neigh=[j for j in (i-1,i+1) if 0<=j<K]
        for j in neigh:a[i,j]=(1-stay)/len(neigh)
    return a


def posterior(values,centers,bw,a):
    emit=-.5*((values[:,None]-centers[None,:])/bw)**2
    la=np.log(np.clip(a,1e-300,None)); alpha=np.empty_like(emit); beta=np.zeros_like(emit)
    alpha[0]=emit[0]-np.log(K)
    for t in range(1,len(values)):alpha[t]=emit[t]+logsumexp(alpha[t-1][:,None]+la,axis=0)
    for t in range(len(values)-2,-1,-1):beta[t]=logsumexp(la+emit[t+1][None,:]+beta[t+1][None,:],axis=1)
    gamma=alpha+beta; gamma-=logsumexp(gamma,axis=1,keepdims=True)
    return np.exp(gamma)


def main():
    OUT.mkdir(parents=True,exist_ok=True); df=pd.read_csv(FRAME); rows=[]; memberships=[]
    for heldout in range(1,7):
        train_reps=[r for r in range(1,7) if r!=heldout]
        train_values=df[df.replica!=heldout].signature_projection.to_numpy()
        centers=np.quantile(train_values,Q); bw=float(np.median(np.diff(centers)))
        taus=[autocorrelation_tau(df[df.replica==r].signature_projection.to_numpy()) for r in train_reps]
        tau=float(np.median(taus)); a=transition(tau)
        train_posts=[]
        for r in train_reps:train_posts.append(posterior(df[df.replica==r].signature_projection.to_numpy(),centers,bw,a))
        mtrain=np.vstack(train_posts)
        test=df[df.replica==heldout]; mtest=posterior(test.signature_projection.to_numpy(),centers,bw,a)
        train_occ=mtrain.mean(0); test_occ=mtest.mean(0)
        entropy=-np.sum(test_occ[test_occ>0]*np.log(test_occ[test_occ>0]))
        continuity=float(np.mean(1-.5*np.abs(np.diff(mtest,axis=0)).sum(1)))
        rmsd=test.ACh_pose_RMSD_global_A.to_numpy()
        weighted=np.array([(mtest[:,j]*rmsd).sum()/mtest[:,j].sum() for j in range(K)])
        mono=float(spearmanr(centers,-weighted).statistic)
        rows.append({'heldout_replica':heldout,'train_tau_ns':tau,'p_stay':a[0,0],
                     'occupancy_JSD':float(jensenshannon(test_occ,train_occ,base=2)**2),
                     'posterior_temporal_continuity':continuity,'effective_prototypes':float(np.exp(entropy)),
                     'max_occupancy':float(test_occ.max()),'center_vs_minus_ACh_RMSD_spearman':mono,
                     **{f'occupancy_{j}':test_occ[j] for j in range(K)},
                     **{f'ACh_RMSD_{j}':weighted[j] for j in range(K)}})
        for i,g in enumerate(test.global_frame):memberships.append({'heldout_replica':heldout,'global_frame':int(g),
            **{f'prototype_{j}_posterior':float(mtest[i,j]) for j in range(K)}})
    folds=pd.DataFrame(rows);folds.to_csv(OUT/'fold_metrics.csv',index=False);pd.DataFrame(memberships).to_csv(OUT/'heldout_posteriors.csv',index=False)
    s={'mean_occupancy_JSD':float(folds.occupancy_JSD.mean()),
       'mean_posterior_temporal_continuity':float(folds.posterior_temporal_continuity.mean()),
       'mean_effective_prototypes':float(folds.effective_prototypes.mean()),
       'mean_max_occupancy':float(folds.max_occupancy.mean()),
       'positive_ACh_monotonicity_folds':int((folds.center_vs_minus_ACh_RMSD_spearman>0).sum())}
    audit={**s,'tau_fit_on_training_replicas_only':True,'development_go':bool(s['mean_occupancy_JSD']<=.30 and
       s['mean_posterior_temporal_continuity']>=.90 and s['mean_effective_prototypes']>=3 and
       s['mean_max_occupancy']<=.60 and s['positive_ACh_monotonicity_folds']>=5),
       'functional_prediction_go':False,'claim_boundary':'Retrospective temporally regularized state head; no efficacy prediction.'}
    (OUT/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    (OUT/'REPORT.md').write_text('# PACER Sticky Soft Coupling Prototypes v0.2\n\n'+folds.to_markdown(index=False,floatfmt='.4f')+
      '\n\n```json\n'+json.dumps(audit,indent=2)+'\n```\n',encoding='utf-8')
    print(folds[['heldout_replica','train_tau_ns','occupancy_JSD','posterior_temporal_continuity','effective_prototypes','max_occupancy','center_vs_minus_ACh_RMSD_spearman']].to_string(index=False));print(json.dumps(audit,indent=2))


if __name__=='__main__':main()
