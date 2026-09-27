#!/usr/bin/env python
"""G2-B supplemental robustness: exact empirical RBF direction vs seeded RFF.

Uses frozen G2-B R2-fitted med/MAD and per-region bandwidth. RFF sweeps
and exact RKHS inner products operate on the *same* balanced 10-frame/window
subset. Also replays the original 512D, original-seed, full-frame direction
and verifies against frozen G2-B output. No inference, no change to old gate.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

import analyze_pacer_fkg_g2b as original

DEFAULT_REGIONS = (
    'compound110_extension', 'cooperativity_mutagenesis', 'distal_control',
    'orthosteric_activation_core', 'pam_contact_consensus',
    'pam_contact_union', 'stable_core_control',
)
CONTEXTS = original.CONTEXTS
SIGNS = original.SIGNS


def exact_gram(x, bandwidth):
    """Shared finite-sample RBF Gram. Operates in bounded chunks of rows."""
    x = np.asarray(x, dtype=np.float32)
    sq = np.einsum('ij,ij->i', x, x).astype(np.float64)
    n = len(x)
    gram = np.empty((n, n), dtype=np.float64)
    for start in range(0, n, 64):
        z = x[start:start+64].astype(np.float64)
        d2 = np.maximum(np.sum(z*z, axis=1)[:, None] + sq[None, :]
                        - 2*z @ x.T.astype(np.float64), 0.0)
        gram[start:start+len(z)] = np.exp(-d2/(2*bandwidth*bandwidth))
    gram = (gram + gram.T)/2
    np.fill_diagonal(gram, 1.0)
    return gram


def signed_weights(context_groups, signs):
    """Groups map (replica, context) to contiguous slices in concatenation."""
    n = sum(s.stop-s.start for s in context_groups.values())
    weights = {2:np.zeros(n, dtype=np.float64), 3:np.zeros(n,dtype=np.float64)}
    for (rep, context), sl in context_groups.items():
        if context in signs:
            weights[rep][sl] = signs[context]/(sl.stop-sl.start)
    return weights[2],weights[3]


def cosine_from_gram(gram, w2, w3):
    cross = float(w2 @ gram @ w3)
    n2 = max(float(w2 @ gram @ w2), 0.0)
    n3 = max(float(w3 @ gram @ w3), 0.0)
    denom = np.sqrt(n2*n3)
    return float(np.clip(cross/denom, -1, 1)) if denom > 1e-12 else None


def rff_cosine(features, w2, w3):
    a = w2 @ features
    b = w3 @ features
    return original.cosine(a,b)


def gram_errors(gram, approx, idx):
    # R3 is excluded from med/MAD and bandwidth fits. Report held-out R3
    # pairwise kernel error separately; all are still same-molecule data.
    inds = np.asarray(idx,dtype=int)
    i,j = np.triu_indices(len(inds),1)
    if len(i)==0:
        raise ValueError('need >=2 held-out samples')
    e = np.abs(gram[np.ix_(inds,inds)][i,j] - approx[np.ix_(inds,inds)][i,j])
    return float(np.mean(e)),float(np.max(e))


def validate(args):
    if args.output_root.exists():
        raise FileExistsError(f'Output already exists, refusing overwrite: {args.output_root}')
    if not args.input_root.is_dir() or not args.graph.is_file() or not args.fkg_audit.is_file() or not args.g2b_audit.is_file():
        raise FileNotFoundError('Missing frozen input root, graph, FKG audit, or G2-B audit')
    if args.sample_frames < 2 or args.sample_frames > 100:
        raise ValueError('--sample-frames must be 2..100')
    if len(set(args.rff_dims))!=len(args.rff_dims) or min(args.rff_dims)<16 or args.seeds<2:
        raise ValueError('Use unique dimensions >=16 and at least 2 seeds')
    base=json.loads(args.g2b_audit.read_text(encoding='utf-8'))
    fkg=json.loads(args.fkg_audit.read_text(encoding='utf-8'))
    graph=json.loads(args.graph.read_text(encoding='utf-8'))
    if base['G2B_STATUS']!='COMPLETED_DESCRIPTIVE_NOT_INFERENTIAL':
        raise ValueError('Unexpected frozen G2-B status')
    expected={'rff_dim':512,'calibration_frames_per_context_window':8,'block_frames':20,'seed':271828,'calibration_replica':2}
    for k,v in expected.items():
        if base['settings'].get(k)!=v:
            raise ValueError(f'Original G2-B {k} differs from expected {v}; do not silently change calibration')
    if len(graph['nodes'])!=270 or len(graph['regions'])!=9:
        raise ValueError('Frozen graph shape mismatch')
    graph_hash=original.sha256(args.graph)
    if graph_hash!=fkg['provenance']['graph_sha256'] or graph_hash!=base['provenance']['graph_sha256']:
        raise ValueError('Frozen graph SHA256 mismatch')
    if original.sha256(args.fkg_audit)!=base['provenance']['fkg_audit_sha256']:
        raise ValueError('G2-B used a different FKG audit file')
    names=sorted(graph['regions'])
    if set(args.regions)-set(names):
        raise ValueError('Unknown region(s): '+str(sorted(set(args.regions)-set(names))))
    old_cal={x['region']:x for x in base['calibration']}
    if set(names)!=set(old_cal):
        raise ValueError('G2-B calibration region mismatch')
    manifest={x['path']:x['sha256'] for x in fkg['provenance']['input_files']}
    if len(manifest)!=40:
        raise ValueError('Expected 40 frozen NPZ SHA256 records')
    for rep in (2,3):
        for win in range(5):
            for context in CONTEXTS:
                path=original.get_path(args.input_root,rep,win,context)
                rel=path.relative_to(args.input_root).as_posix()
                if manifest.get(rel)!=original.sha256(path):
                    raise ValueError(f'Frozen embedding SHA256 mismatch: {path}')
    return base,graph,{'input_sha256_verified':40,'graph_sha256':graph_hash,
        'original_fkg_sha256':original.sha256(args.fkg_audit),
        'original_g2b_sha256':original.sha256(args.g2b_audit),
        'original_g2b_settings':base['settings']}


def write_csv(path,records):
    if not records: raise ValueError('Empty output records')
    with path.open('w',encoding='utf-8-sig',newline='') as stream:
        fields=list(records[0]);w=csv.DictWriter(stream,fieldnames=fields)
        w.writeheader();w.writerows(records)


def analyze(args,base,graph,provenance):
    med=np.asarray(base['channel_median'],dtype=np.float64)
    scale=np.asarray(base['channel_scale'],dtype=np.float64)
    if med.shape!=(128,) or scale.shape!=(128,) or (scale<=0).any():
        raise ValueError('Bad frozen channel calibration')
    sorted_regions=sorted(graph['regions'])
    calib={x['region']:x for x in base['calibration']}
    prior={(x['region'],x['axis']):x for x in base['cross_replica_direction']}
    frame_ids=np.linspace(0,99,args.sample_frames,dtype=int)
    summary=[];seed_rows=[];baseline_rows=[]
    for region in args.regions:
        ix=np.asarray([x['embedding_index'] for x in graph['regions'][region]],dtype=int)
        width=float(calib[region]['bandwidth'])
        if width<=0 or len(ix)!=calib[region]['n_residues']:
            raise ValueError('Frozen bandwidth/region length mismatch')
        # Load region data once, keeping memory bounded by one region.
        full={};parts=[];groups={};offset=0;holdout=[]
        for rep in (2,3):
            for context in CONTEXTS:
                windows=[]
                for win in range(5):
                    with np.load(original.get_path(args.input_root,rep,win,context),allow_pickle=False) as file:
                        raw=file['residue_features'].astype(np.float32)
                    x=original.unit_normalize(raw)
                    x=((x-med)/scale).astype(np.float32)[:,ix,:].reshape(100,-1)
                    windows.append(x)
                joined=np.concatenate(windows)
                full[(rep,context)]=joined
                chosen=np.concatenate([x[frame_ids] for x in windows],axis=0)
                parts.append(chosen)
                sl=slice(offset,offset+len(chosen))
                groups[(rep,context)]=sl
                if rep==3:holdout.extend(range(sl.start,sl.stop))
                offset=sl.stop
        subset=np.concatenate(parts).astype(np.float32)
        gram=exact_gram(subset,width)
        ex={axis:cosine_from_gram(gram,*signed_weights(groups,signs)) for axis,signs in SIGNS.items()}
        # Baseline is full 500 frames/context, exact replay of original 512-D
        # RFF seed and original per-window block aggregation.
        original_seed=base['settings']['seed']+sorted_regions.index(region)*101
        W0,b0=original.make_rff(len(ix)*128,width,512,original_seed)
        mu={}
        for key,x in full.items():
            phi=original.rff(x,W0,b0)
            blocks=phi.reshape(5,5,20,512).mean(axis=2)
            mu[key]=blocks.mean(axis=(0,1)).astype(np.float64)
        for axis,signs in SIGNS.items():
            delta={rep:sum(sign*mu[(rep,c)] for c,sign in signs.items()) for rep in (2,3)}
            replay=original.cosine(delta[2],delta[3])
            frozen=prior[(region,axis)]['rff_pooled_direction_cosine']
            diff=abs(replay-frozen) if replay is not None and frozen is not None else None
            if diff is None or diff>args.replay_tolerance:
                raise ValueError(f'Original G2-B replay failed for {region}/{axis}: replay={replay} frozen={frozen} diff={diff}')
            baseline_rows.append({'region':region,'axis':axis,'frozen_full_512_cosine':frozen,
                'replayed_full_512_cosine':replay,'absolute_replay_error':diff,
                'exact_subsample_cosine':ex[axis], 'subsample_frames_per_context':args.sample_frames*5})
        print(f'G2B_ROBUSTNESS: {region} exact kernel and 512-D replay checked',flush=True)
        # For every RFF setting use *identical subset* as exact Gram.
        for dim in args.rff_dims:
            for seed_index in range(args.seeds):
                seed=original_seed+seed_index*100003
                W,b=original.make_rff(len(ix)*128,width,dim,seed)
                phi=original.rff(subset,W,b).astype(np.float64)
                approx=phi@phi.T
                r3mae,r3max=gram_errors(gram,approx,holdout)
                for axis,signs in SIGNS.items():
                    w2,w3=signed_weights(groups,signs)
                    val=rff_cosine(phi,w2,w3)
                    seed_rows.append({'region':region,'axis':axis,'rff_dim':dim,'seed':seed,
                        'rff_subsample_cosine':val,'exact_subsample_cosine':ex[axis],
                        'absolute_cosine_error':abs(val-ex[axis]) if val is not None and ex[axis] is not None else None,
                        'r3_heldout_kernel_mae':r3mae,'r3_heldout_kernel_maxae':r3max})
        del full,parts,gram,subset
    for region in args.regions:
        for axis in SIGNS:
            rows=[r for r in seed_rows if r['region']==region and r['axis']==axis]
            for dim in args.rff_dims:
                subsetrows=[r for r in rows if r['rff_dim']==dim]
                vals=np.asarray([r['rff_subsample_cosine'] for r in subsetrows if r['rff_subsample_cosine'] is not None],dtype=float)
                errs=np.asarray([r['absolute_cosine_error'] for r in subsetrows if r['absolute_cosine_error'] is not None],dtype=float)
                summary.append({'region':region,'axis':axis,'rff_dim':dim,
                    'exact_subsample_cosine':next(r['exact_subsample_cosine'] for r in subsetrows),
                    'rff_seed_min':float(vals.min()) if len(vals) else None,
                    'rff_seed_median':float(np.median(vals)) if len(vals) else None,
                    'rff_seed_max':float(vals.max()) if len(vals) else None,
                    'max_absolute_cosine_error':float(errs.max()) if len(errs) else None,
                    'max_r3_heldout_kernel_mae':max(r['r3_heldout_kernel_mae'] for r in subsetrows),
                    'sign_stable_across_seeds':bool(len(vals)==args.seeds and (np.all(vals>0) or np.all(vals<0)))})
    return summary,seed_rows,baseline_rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',type=Path,required=True)
    p.add_argument('--graph',type=Path,required=True)
    p.add_argument('--fkg-audit',type=Path,required=True)
    p.add_argument('--g2b-audit',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--regions',nargs='+',default=list(DEFAULT_REGIONS))
    p.add_argument('--sample-frames',type=int,default=10)
    p.add_argument('--rff-dims',nargs='+',type=int,default=[256,512,1024])
    p.add_argument('--seeds',type=int,default=3)
    p.add_argument('--replay-tolerance',type=float,default=0.005)
    p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args()
    base,graph,prov=validate(a)
    if a.preflight_only:
        print(json.dumps({'G2B_ROBUSTNESS_PREFLIGHT':'PASS', 'regions':a.regions,
            'rff_dims':a.rff_dims,'seeds':a.seeds,'sample_frames':a.sample_frames,
            **prov,'output_untouched':str(a.output_root)},indent=2,ensure_ascii=False))
        return
    summaries,seeds,baseline=analyze(a,base,graph,prov)
    a.output_root.mkdir(parents=True,exist_ok=False)
    report={'G2B_ROBUSTNESS_STATUS':'COMPLETE_DESCRIPTIVE','provenance':prov,
        'configuration':{'regions':a.regions,'sample_frames_per_window_context':a.sample_frames,
            'rff_dims':a.rff_dims,'seeds':a.seeds,'replay_tolerance':a.replay_tolerance},
        'summary':summaries,'seed_level':seeds,'baseline_replay':baseline,
        'interpretation_boundary':[
            'Exact kernel is empirical RBF mean-embedding cosine on the same deterministic balanced subsample, not an infinite-sample population kernel.',
            'Subsample estimates are NOT directly interchangeable with frozen G2-B full-500-frame-per-context estimates.',
            'R3 holdout refers only to frozen preprocessing/bandwidth calibration, not to molecule or chemical generalization.',
            'Finite 2-replica data, serially correlated frames and overlapping regions prohibit p-value or efficacy claims.',
            'No gate changes, no additional model fit and no random seed selected post hoc.']}
    (a.output_root/'G2B_ROBUSTNESS_AUDIT.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    write_csv(a.output_root/'G2B_ROBUSTNESS_SUMMARY.csv',summaries)
    write_csv(a.output_root/'G2B_ROBUSTNESS_SEEDS.csv',seeds)
    write_csv(a.output_root/'G2B_BASELINE_REPLAY.csv',baseline)
    reportlines=['# PACER-FKG G2-B supplementary RFF/exact-kernel audit','',
        'Status: descriptive only; original G2-B gate, samples and report remain frozen.',
        f'Frozen input hashes verified: {prov["input_sha256_verified"]}/40; original G2-B SHA256 `{prov["original_g2b_sha256"]}`.',
        '', '| Region | Axis | full frozen 512D cosine | exact subsample cosine | 1024D max absolute error |',
        '|---|---|---:|---:|---:|']
    bm={(x['region'],x['axis']):x for x in baseline}
    sm={(x['region'],x['axis'],x['rff_dim']):x for x in summaries}
    for region in a.regions:
        for axis in SIGNS:
            old=bm[(region,axis)]; strongest=sm[(region,axis,max(a.rff_dims))]
            reportlines.append(f'| {region} | {axis} | {old["frozen_full_512_cosine"]:.3f} | '
                f'{old["exact_subsample_cosine"]:.3f} | {strongest["max_absolute_cosine_error"]:.3f} |')
    reportlines += ['','## Boundaries','',*('- '+x for x in report['interpretation_boundary'])]
    (a.output_root/'G2B_ROBUSTNESS_REPORT.md').write_text('\n'.join(reportlines)+'\n',encoding='utf-8')
    print(json.dumps({'G2B_ROBUSTNESS_STATUS':'COMPLETE_DESCRIPTIVE',
        'regions':len(a.regions),'axis_dim_rows':len(summaries),
        'seed_rows':len(seeds),'baseline_replay_rows':len(baseline),
        'output':str(a.output_root)},indent=2))

if __name__=='__main__':main()
