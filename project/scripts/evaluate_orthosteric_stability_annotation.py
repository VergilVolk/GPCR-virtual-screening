"""Independent ACh-pose stability annotation of frozen M4 trajectory states."""

from __future__ import annotations

import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from scipy.stats import rankdata


PROJECT=Path(__file__).resolve().parents[1]
PUBLIC=PROJECT/'data'/'m4_gamd_figshare_33283491'
OUT=PROJECT/'results'/'pacer_orthosteric_stability_annotation_v01'
SIGNATURE=PROJECT/'results'/'pacer_consensus_pam_signature_v01'/'frame_signature_scores.csv'
ASSIGN=PROJECT/'results'/'pacer_anchor_guided_triplet_v03'/'heldout_assignments.csv'
RESIDS=[687,690,691,694,695,696,706,781,782,783,784,785,788,852,856,859,868,869,871,872,875]


def align(mobile,reference,coords):
    mc=mobile.mean(0); rc=reference.mean(0)
    rot,_=Rotation.align_vectors(reference-rc,mobile-mc)
    return rot.apply(coords-mc)+rc


def rmsd(a,b): return float(np.sqrt(np.mean(np.sum((a-b)**2,axis=1))))


def zrank(x):
    r=rankdata(x); return (r-r.mean())/r.std(ddof=0)


def eta_squared(v,labels):
    total=np.sum((v-v.mean())**2)
    return float(sum(np.sum(labels==c)*(v[labels==c].mean()-v.mean())**2 for c in np.unique(labels))/total)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    top=str(PUBLIC/'sys-protein.pdb')
    u0=mda.Universe(top,str(PUBLIC/'run1-GaMD-stride1000ps.nc'))
    prot0=u0.select_atoms('name CA and resid '+' '.join(map(str,RESIDS)))
    ach0=u0.select_atoms('resid 907 and not name H*')
    u0.trajectory[0]; ref_prot=prot0.positions.astype(float).copy(); ref_ach=ach0.positions.astype(float).copy()
    if len(prot0)!=21 or len(ach0)!=10: raise RuntimeError((len(prot0),len(ach0)))
    rows=[]
    for rep in range(1,7):
        u=mda.Universe(top,str(PUBLIC/f'run{rep}-GaMD-stride1000ps.nc'))
        prot=u.select_atoms('name CA and resid '+' '.join(map(str,RESIDS)))
        ach=u.select_atoms('resid 907 and not name H*')
        local_ref=None
        for frame,ts in enumerate(u.trajectory):
            transformed=align(prot.positions.astype(float),ref_prot,ach.positions.astype(float))
            if local_ref is None: local_ref=transformed.copy()
            rows.append({'replica':rep,'frame_ns':frame,'global_frame':(rep-1)*500+frame,
                         'ACh_pose_RMSD_global_A':rmsd(transformed,ref_ach),
                         'ACh_pose_RMSD_local_A':rmsd(transformed,local_ref)})
    frame=pd.DataFrame(rows)
    sig=pd.read_csv(SIGNATURE)[['replica','frame_ns','signature_projection']]
    frame=frame.merge(sig,on=['replica','frame_ns'],validate='one_to_one')
    frame.to_csv(OUT/'frame_ach_stability.csv',index=False)

    rep_rows=[]; zr_sig=[]; zr_rmsd=[]
    for rep in range(1,7):
        x=frame[frame.replica==rep]
        a=zrank(x.signature_projection.to_numpy()); b=zrank(-x.ACh_pose_RMSD_global_A.to_numpy())
        al=zrank(-x.ACh_pose_RMSD_local_A.to_numpy())
        rho=float(np.mean(a*b)); rho_local=float(np.mean(a*al))
        rep_rows.append({'replica':rep,'rho_signature_vs_minus_global_RMSD':rho,
                         'rho_signature_vs_minus_local_RMSD':rho_local,
                         'mean_global_RMSD_A':x.ACh_pose_RMSD_global_A.mean(),
                         'mean_local_RMSD_A':x.ACh_pose_RMSD_local_A.mean()})
        zr_sig.append(a); zr_rmsd.append(b)
    rep_df=pd.DataFrame(rep_rows); rep_df.to_csv(OUT/'replica_correlations.csv',index=False)
    observed=rep_df.rho_signature_vs_minus_global_RMSD.mean()
    rng=np.random.default_rng(20260901); null=[]
    for _ in range(10000):
        null.append(np.mean([np.mean(a*np.roll(b,int(rng.integers(25,475)))) for a,b in zip(zr_sig,zr_rmsd)]))
    null=np.asarray(null); p=float((1+np.sum(null>=observed))/(len(null)+1))

    assign=pd.read_csv(ASSIGN); fold_rows=[]
    for (method,rep,seed),g in assign.groupby(['method','heldout_replica','seed']):
        g=g.sort_values('global_frame'); values=frame.set_index('global_frame').loc[g.global_frame,'ACh_pose_RMSD_global_A'].to_numpy()
        fold_rows.append({'method':method,'heldout_replica':rep,'seed':seed,
                          'eta_squared_ACh_RMSD':eta_squared(values,g.state.to_numpy())})
    folds=pd.DataFrame(fold_rows); folds.to_csv(OUT/'cluster_alignment.csv',index=False)
    pivot=folds.groupby(['method','heldout_replica']).eta_squared_ACh_RMSD.mean().unstack(0)
    delta=pivot['AnchorGuidedTriplet-16']-pivot['PCA-10']
    audit={'triplet_never_used_ligand_coordinates':True,'mean_within_replica_rho':float(observed),
           'positive_rho_replicas':int((rep_df.rho_signature_vs_minus_global_RMSD>0).sum()),
           'circular_shift_p':p,'triplet_minus_pca_eta_squared':float(delta.mean()),
           'triplet_better_folds':int((delta>0).sum()),
           'continuous_coupling_support':bool(observed>0 and (rep_df.rho_signature_vs_minus_global_RMSD>0).sum()>=4 and p<=.05),
           'triplet_cluster_increment_support':bool(delta.mean()>=.05 and (delta>0).sum()>=4),
           'claim_boundary':'ACh pose-stability annotation in one MK-97 system; not PAM efficacy validation.'}
    (OUT/'audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    report='# PACER Orthosteric-Stability Annotation v0.1\n\n'+rep_df.to_markdown(index=False,floatfmt='.4f')
    report+='\n\n## Cluster eta-squared\n\n'+pivot.reset_index().to_markdown(index=False,floatfmt='.4f')
    report+='\n\n```json\n'+json.dumps(audit,indent=2)+'\n```\n'
    (OUT/'REPORT.md').write_text(report,encoding='utf-8')
    print(rep_df.to_string(index=False));print(pivot.to_string());print(json.dumps(audit,indent=2))


if __name__=='__main__': main()
