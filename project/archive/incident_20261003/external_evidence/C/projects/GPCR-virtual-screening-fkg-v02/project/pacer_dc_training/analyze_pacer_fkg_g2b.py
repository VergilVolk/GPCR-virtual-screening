#!/usr/bin/env python
"""PACER-FKG G2-B: frozen common RBF-RFF space and descriptive block bootstrap.

No modification to frozen FKG results and no inference from independent frames.
Calibrate channel median/MAD and per-region RBF bandwidth on R2 only;
apply without refitting to R3. Compute all contrasts in identical RFF coordinates.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

CONTEXTS = ('apo', 'probe_only', 'candidate_no_probe', 'candidate_probe')
SIGNS = {
    'synergy_interaction': {'candidate_probe': 1, 'probe_only': -1, 'candidate_no_probe': -1, 'apo': 1},
    'intrinsic_agonism': {'candidate_no_probe': 1, 'apo': -1},
    'conditional_pam_effect': {'candidate_probe': 1, 'probe_only': -1},
}


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def get_path(root, r, w, c):
    return root / f'replica_{r:02d}' / f'window_{w:03d}' / 'atom14' / f'{c}_w{w:03d}.geom2vec.npz'


def unit_normalize(a):
    return a - a.mean(axis=1, keepdims=True)


def fit_channels(samples):
    flat = samples.reshape(-1, samples.shape[-1]).astype(np.float64)
    med = np.median(flat, axis=0)
    scale = np.maximum(1.4826 * np.median(np.abs(flat - med), axis=0), 1e-6)
    return med, scale


def shared_bandwidth(data):
    x = data.astype(np.float64)
    d2 = np.sum(x*x, axis=1)[:,None] + np.sum(x*x, axis=1)[None,:] - 2*x @ x.T
    upper = np.maximum(d2[np.triu_indices(len(x), 1)], 0)
    pos = upper[upper > 1e-12]
    if not len(pos):
        raise ValueError('degenerate bandwidth calibration')
    return float(np.sqrt(np.median(pos)))


def make_rff(input_dim, width, n_features, seed):
    rng = np.random.default_rng(seed)
    W = rng.normal(0, 1/width, size=(input_dim, n_features)).astype(np.float32)
    bias = rng.uniform(0, 2*np.pi, size=n_features).astype(np.float32)
    return W, bias


def rff(x, W, bias):
    return np.cos(x.astype(np.float32) @ W + bias, dtype=np.float32) * np.float32(np.sqrt(2/W.shape[1]))


def cosine(a,b):
    n = float(np.linalg.norm(a)*np.linalg.norm(b))
    return float(np.dot(a,b)/n) if n > 1e-12 else None


def percentile_or_none(values, pct):
    vals = [x for x in values if x is not None and np.isfinite(x)]
    return float(np.percentile(vals,pct)) if vals else None


def bootstrap_direction(blocks, signs, iterations, seed):
    """Resample complete 100-frame windows, then contiguous 20-frame blocks.
    Independently sample both replicas and contexts. No frame-level bootstrap.
    """
    rng = np.random.default_rng(seed)
    out=[]
    contexts = tuple(signs)
    for _ in range(iterations):
        deltas=[]
        for r in (2,3):
            selected_windows = rng.integers(0,5,5)
            delta = None
            for c in contexts:
                x = blocks[r][c] # [windows,blocks,features]
                blockids = rng.integers(0,x.shape[1],(5,x.shape[1]))
                v=x[selected_windows[:,None],blockids,:].mean(axis=(0,1))
                delta = signs[c]*v if delta is None else delta+signs[c]*v
            deltas.append(delta)
        out.append(cosine(*deltas))
    return {'bootstrap_valid':sum(x is not None for x in out),
            'cosine_ci025':percentile_or_none(out,2.5),
            'cosine_ci50':percentile_or_none(out,50),
            'cosine_ci975':percentile_or_none(out,97.5),
            'fraction_cosine_positive':float(np.mean([x>0 for x in out if x is not None])) if any(x is not None for x in out) else None}


def check_inputs(args):
    if args.output_root.exists():
        raise FileExistsError(f'Output exists: {args.output_root}; results are immutable')
    if args.block_frames < 1 or 100 % args.block_frames:
        raise ValueError('--block-frames must divide 100 exactly')
    if args.rff_dim < 16 or args.bootstrap < 1:
        raise ValueError('invalid --rff-dim or --bootstrap')
    graph=json.loads(args.graph.read_text(encoding='utf-8'))
    old=json.loads(args.fkg_audit.read_text(encoding='utf-8'))
    if len(graph['nodes'])!=270 or len(graph['regions'])!=9:
        raise ValueError('unexpected frozen graph')
    if old['provenance']['graph_sha256']!=sha256(args.graph):
        raise ValueError('graph hash mismatch to frozen FKG audit')
    manifest={r['path']:r['sha256'] for r in old['provenance']['input_files']}
    if len(manifest)!=40:
        raise ValueError('expected exact frozen 40-file manifest')
    count=0
    for rep in (2,3):
        for win in range(5):
            frame_ref=None
            for c in CONTEXTS:
                p=get_path(args.input_root,rep,win,c)
                rel=p.relative_to(args.input_root).as_posix()
                if not p.is_file() or manifest.get(rel)!=sha256(p):
                    raise ValueError(f'frozen NPZ missing or SHA256 mismatch: {p}')
                with np.load(p,allow_pickle=False) as d:
                    feat=d['residue_features']
                    ids=d['frame_ids']
                    seq=str(d['sequence'].item())
                    if feat.shape!=(100,270,128) or not np.isfinite(feat).all():
                        raise ValueError(f'feature QC failed: {p}')
                    if len(seq)!=270 or any(seq[i]!=node['site'][0] for i,node in enumerate(graph['nodes'])):
                        raise ValueError(f'sequence QC failed: {p}')
                    if ids.shape!=(100,) or np.any(np.diff(ids)<=0):
                        raise ValueError(f'frame ID QC failed: {p}')
                    if frame_ref is None:frame_ref=ids
                    elif not np.array_equal(frame_ref,ids):raise ValueError(f'context frames mismatch: {p}')
                count+=1
    g2a={}
    for name,path in [('atom14',args.g2a_atom14),('pbc',args.g2a_pbc)]:
        if path is not None:
            d=json.loads(path.read_text(encoding='utf-8'))
            g2a[name]={'path':str(path),'sha256':sha256(path),'status':d.get('G2A_STATUS',d.get('G2A_PBC_STATUS'))}
    return graph, {'inputs_sha256_verified':count,'graph_sha256':sha256(args.graph),
                   'fkg_audit_sha256':sha256(args.fkg_audit),'g2a':g2a,
                   'fkg_gate_unchanged':True}


def calibrate(args, graph):
    # Fixed equally spaced frame sampling in R2 only, all contexts and windows balanced.
    rows=[]
    for w in range(5):
        for c in CONTEXTS:
            with np.load(get_path(args.input_root,2,w,c),allow_pickle=False) as d:
                a=d['residue_features'][np.linspace(0,99,args.calibration_frames,dtype=int)].astype(np.float32)
                rows.append(unit_normalize(a))
    pooled=np.concatenate(rows,axis=0)
    med,scale=fit_channels(pooled)
    standardized=(pooled-med)/scale
    calibration={}
    for region,members in sorted(graph['regions'].items()):
        ix=[m['embedding_index'] for m in members]
        x=standardized[:,ix,:].reshape(len(standardized),-1)
        width=shared_bandwidth(x)
        calibration[region]={'indices':ix,'bandwidth':width,'n_residues':len(ix)}
    return med,scale,calibration,standardized


def rff_error(x,width,W,bias):
    # Calibration diagnostic, not independent validation; lower errors support approximation.
    x=x.astype(np.float32)
    z=rff(x,W,bias)
    gram=z@z.T
    sq=np.sum(x*x,axis=1)
    dist=np.maximum(sq[:,None]+sq[None,:]-2*x@x.T,0)
    exact=np.exp(-dist/(2*width*width))
    iu=np.triu_indices(len(x),1)
    return {'calibration_kernel_mae':float(np.mean(np.abs(exact[iu]-gram[iu]))),
            'calibration_kernel_maxae':float(np.max(np.abs(exact[iu]-gram[iu])))}


def analyse(args,graph,provenance):
    med,scale,cal,standardized=calibrate(args,graph)
    per_region={}
    diagnostics=[]
    for i,(region,info) in enumerate(sorted(cal.items())):
        ix=info['indices']; width=info['bandwidth']
        seed=args.seed+i*101
        W,bias=make_rff(len(ix)*128,width,args.rff_dim,seed)
        x=standardized[:,ix,:].reshape(len(standardized),-1)
        err=rff_error(x,width,W,bias)
        diagnostics.append({'region':region,'calibration_replica':2,'calibration_samples':len(x),
                            'bandwidth':width,'rff_dimension':args.rff_dim,'rff_seed':seed,
                            'n_residues':len(ix),**err})
        per_region[region]={'rff_W':W,'rff_bias':bias,'indices':ix,'blocks':{2:{},3:{}}}
    # One NPZ read per file. Transform all regions in shared R2-fitted coordinates.
    for rep in (2,3):
        for c in CONTEXTS:
            for region in per_region:per_region[region]['blocks'][rep][c]=[]
            for w in range(5):
                with np.load(get_path(args.input_root,rep,w,c),allow_pickle=False) as d:
                    a=unit_normalize(d['residue_features'].astype(np.float32))
                a=((a-med)/scale).astype(np.float32)
                for region,info in per_region.items():
                    x=a[:,info['indices'],:].reshape(100,-1)
                    z=rff(x,info['rff_W'],info['rff_bias'])
                    meanblocks=z.reshape(100//args.block_frames,args.block_frames,args.rff_dim).mean(axis=1)
                    info['blocks'][rep][c].append(meanblocks)
                print(f'encoded R{rep} W{w} {c}',flush=True)
    for info in per_region.values():
        for rep in (2,3):
            for c in CONTEXTS:info['blocks'][rep][c]=np.stack(info['blocks'][rep][c])
    summary=[];windowrows=[]
    for rid,(region,info) in enumerate(per_region.items()):
        blocks=info['blocks']
        for ai,(axis,signs) in enumerate(SIGNS.items()):
            deltas={rep:sum(signs[c]*blocks[rep][c].mean(axis=1) for c in signs) for rep in (2,3)}
            # [windows, features] per replica, same RFF basis
            matched=[cosine(deltas[2][w],deltas[3][w]) for w in range(5)]
            for w in range(5):
                windowrows.append({'region':region,'axis':axis,'window':w,
                    'direction_cosine':matched[w], 'r2_norm':float(np.linalg.norm(deltas[2][w])),
                    'r3_norm':float(np.linalg.norm(deltas[3][w]))})
            pooled=cosine(deltas[2].mean(axis=0),deltas[3].mean(axis=0))
            boot=bootstrap_direction(blocks,signs,args.bootstrap,args.seed+rid*1001+ai*41)
            summary.append({'region':region,'axis':axis,'rff_pooled_direction_cosine':pooled,
                'matched_window_cosines':matched,'matched_window_median_cosine':percentile_or_none(matched,50),
                'r2_contrast_norm':float(np.linalg.norm(deltas[2].mean(axis=0))),
                'r3_contrast_norm':float(np.linalg.norm(deltas[3].mean(axis=0))),**boot})
    return summary,windowrows,diagnostics,med,scale


def write_csv(path,rows):
    if not rows:raise ValueError('empty CSV')
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',type=Path,required=True)
    p.add_argument('--graph',type=Path,required=True)
    p.add_argument('--fkg-audit',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--g2a-atom14',type=Path)
    p.add_argument('--g2a-pbc',type=Path)
    p.add_argument('--rff-dim',type=int,default=512)
    p.add_argument('--calibration-frames',type=int,default=8)
    p.add_argument('--block-frames',type=int,default=20)
    p.add_argument('--bootstrap',type=int,default=400)
    p.add_argument('--seed',type=int,default=271828)
    p.add_argument('--preflight-only',action='store_true')
    args=p.parse_args()
    if args.calibration_frames<2 or args.calibration_frames>100:raise ValueError('invalid calibration frame count')
    graph,provenance=check_inputs(args)
    if args.preflight_only:
        print(json.dumps({'G2B_PREFLIGHT':'PASS',**provenance,'output_untouched':str(args.output_root)},indent=2,ensure_ascii=False))
        return
    summaries,windowrows,calibration,med,scale=analyse(args,graph,provenance)
    report={'G2B_STATUS':'COMPLETED_DESCRIPTIVE_NOT_INFERENTIAL',
        'method':'R2-only frozen robust channel calibration; common per-region RBF kernel RFF; R2/R3 directional cosine',
        'provenance':provenance,'settings':{'rff_dim':args.rff_dim,'calibration_frames_per_context_window':args.calibration_frames,
          'block_frames':args.block_frames,'bootstrap_draws':args.bootstrap,'seed':args.seed,'calibration_replica':2},
        'calibration':calibration,'channel_median':med.tolist(),'channel_scale':scale.tolist(),
        'cross_replica_direction':summaries,
        'limitations':['Random Fourier features approximate the RBF RKHS; calibration MAE is in-sample.',
          'R2 is calibration AND evaluated replica; R3 is held out for preprocessing and bandwidth, but not an independent compound.',
          '5 contiguous windows and 20-frame blocks are serially correlated; bootstrap is descriptive, not a population CI or p-value.',
          'A positive cross-replica cosine is RKHS contrast direction agreement, not pharmacological sign.',
          'Regions overlap; do not count passing regions as independent evidence.',
          'Shared kernel/preprocessing differs from original per-window estimator; original gate remains unchanged.',
          'Two replicas alone cannot establish generalization or efficacy.']}
    args.output_root.mkdir(parents=True,exist_ok=False)
    (args.output_root/'G2B_DIRECTION_AUDIT.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    write_csv(args.output_root/'G2B_DIRECTION_SUMMARY.csv',
              [{k:(json.dumps(v) if isinstance(v,list) else v) for k,v in x.items()} for x in summaries])
    write_csv(args.output_root/'G2B_WINDOW_DIRECTIONS.csv',windowrows)
    write_csv(args.output_root/'G2B_CALIBRATION.csv',calibration)
    lines=['# PACER-FKG G2-B common-kernel direction audit','',
      'Status: descriptive only. No gate modification, drug efficacy inference, or independent-replica significance.',
      '',f'Frozen input checks: {provenance["inputs_sha256_verified"]}/40; RFF dimension: {args.rff_dim}; bootstrap: {args.bootstrap}.',
      '', '| Region | Axis | pooled R2/R3 cosine | block-bootstrap 2.5–97.5% | calibration RFF MAE |',
      '|---|---|---:|---:|---:|']
    errors={x['region']:x['calibration_kernel_mae'] for x in calibration}
    for x in summaries:
        ci=f'{x["cosine_ci025"]:.3f} to {x["cosine_ci975"]:.3f}' if x['cosine_ci025'] is not None else 'undefined'
        pooled=f'{x["rff_pooled_direction_cosine"]:.3f}' if x['rff_pooled_direction_cosine'] is not None else 'undefined'
        lines.append(f'| {x["region"]} | {x["axis"]} | {pooled} | {ci} | {errors[x["region"]]:.4f} |')
    lines+=['','## Interpretation boundary','',*('- '+x for x in report['limitations'])]
    (args.output_root/'G2B_DIRECTION_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({'G2B_STATUS':report['G2B_STATUS'],'frozen_inputs':provenance['inputs_sha256_verified'],
       'region_axis_rows':len(summaries),'window_rows':len(windowrows),'calibration_rows':len(calibration),
       'output':str(args.output_root)},indent=2))

if __name__=='__main__':main()
