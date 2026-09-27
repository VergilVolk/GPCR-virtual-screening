#!/usr/bin/env python
"""Independent, read-only DCD -> frozen raw coordinate and PBC diagnostic.

Never overwrites prior G2-A audits. Does not unwrap, reimage, or modify data.
PBC findings are diagnostics, not a claim of universal PBC correctness.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import tempfile
from pathlib import Path

import numpy as np

SYSTEMS = ('apo','probe_only','compound110__candidate_no_probe','compound110__candidate_probe')
CTX = {'apo':'apo','probe_only':'probe_only','compound110__candidate_no_probe':'candidate_no_probe','compound110__candidate_probe':'candidate_probe'}
REPS = (2,3)
WINDOWS = range(5)


def hash_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(8<<20), b''):
            h.update(block)
    return h.hexdigest()


def inputs(source: Path, raw_root: Path):
    out=[]
    for system in SYSTEMS:
        pdb=source/'topology'/system/'minimized.pdb'
        for rep in REPS:
            fname='trajectory_corrected.dcd' if (system,rep)==('apo',3) else 'trajectory.dcd'
            dcd=source/'trajectories'/system/f'replica_{rep:02d}'/fname
            for w in WINDOWS:
                unit=raw_root/f'replica_{rep:02d}'/f'window_{w:03d}'
                raw=unit/'raw'/f'{CTX[system]}_w{w:03d}.npz'
                manifest=unit/'trajectory_windows.json'
                out.append(dict(system=system,replica=rep,window=w,pdb=pdb,dcd=dcd,raw=raw,manifest=manifest,context=CTX[system]))
    return out


def check_metadata(rec):
    j=json.loads(rec['manifest'].read_text(encoding='utf-8'))
    rows=[row for row in j['rows'] if row.get('context')==rec['context']]
    if len(rows)!=1: raise ValueError(f'metadata context not unique: {rec["manifest"]}')
    m=rows[0]
    if m.get('status')!='ok': raise ValueError(f'raw metadata status not ok: {rec["raw"]}')
    if int(m['frame_count'])!=100: raise ValueError('not 100 frames')
    start=int(m['start_frame']);end=int(m['end_frame'])
    if end-start+1!=100: raise ValueError('window boundaries not 100')
    if rec['window'] != start//100 or start%100: raise ValueError(f'unexpected window start {start}')
    # The original record must identify the corrected R3 apo trajectory.
    expected='trajectory_corrected.dcd' if (rec['system'],rec['replica'])==('apo',3) else 'trajectory.dcd'
    src=Path(m['trajectory_path'].replace('\\','/'))
    if src.name!=expected: raise ValueError(f'wrong trajectory provenance {src} expected {expected}')
    if rec['system'] not in src.parts and not (rec['system']=='apo' and 'apo' in src.parts):
        raise ValueError(f'wrong source system in metadata: {src}')
    return start,end


def prepare_topology(pdb: Path, tmp: Path) -> Path:
    # Same CONECT-free parsing route as original trajectory extractor.
    dest=tmp/pdb.parent.name/'minimized_noconect.pdb'
    dest.parent.mkdir(parents=True,exist_ok=True)
    with pdb.open('r',encoding='utf-8',errors='replace') as r, dest.open('w',encoding='utf-8') as w:
        for line in r:
            if not line.startswith('CONECT'): w.write(line)
    return dest


def compare_raw(raw, sel, xyz, rep, window, context, tolerance):
    with np.load(raw,allow_pickle=False) as z:
        stored=z['coords'];names=z['atom_names'].astype(str);rnames=z['resnames'].astype(str)
        ridx=z['res_index'];resids=z['resids']
    if stored.shape != xyz.shape: raise ValueError(f'{raw}: xyz shape {stored.shape} vs {xyz.shape}')
    if names.tolist()!=sel.names.tolist(): raise ValueError(f'{raw}: atom_names differ')
    if rnames.tolist()!=sel.residues.resnames.tolist(): raise ValueError(f'{raw}: residue names differ')
    if not np.array_equal(resids.astype(int),sel.residues.resids.astype(int)):
        raise ValueError(f'{raw}: residue ids differ')
    # MDAnalysis atom.resindex is universe-global; map to selection-local indices.
    mapping={int(r.resindex):i for i,r in enumerate(sel.residues)}
    expected=np.fromiter((mapping[int(k)] for k in sel.resindices),dtype=np.int64,count=len(sel))
    if not np.array_equal(ridx,expected):raise ValueError(f'{raw}: residue atom mapping differs')
    diff=np.abs(stored.astype(np.float64)-xyz.astype(np.float64))
    dmax=float(diff.max())
    # DCD coordinates use lossy precision; compare tolerance 0.05 A.
    if not np.isfinite(dmax) or dmax>tolerance:
        raise ValueError(f'{raw}: coordinate difference max={dmax:.6g} A > {tolerance} A')
    return dmax


def minimum_image_lengths(vectors,box):
    from MDAnalysis.lib.distances import minimize_vectors
    return np.linalg.norm(minimize_vectors(np.asarray(vectors,dtype=np.float64),np.asarray(box,dtype=np.float32)),axis=-1)


def valid_box(box):
    if box is None: return False
    b=np.asarray(box,dtype=float)
    return b.shape==(6,) and bool(np.isfinite(b).all()) and bool((b[:3]>0).all()) and bool(((b[3:]>0)&(b[3:]<180)).all())


def analyze_frame(ca,box,frame,prev_ca,prev_box, counters, examples):
    if not valid_box(box):
        counters['invalid_box_frames']+=1
        if len(examples['invalid_box'])<12: examples['invalid_box'].append(frame)
        return None
    counters['valid_box_frames']+=1
    direct=np.linalg.norm(np.diff(ca,axis=0),axis=1)
    mic=minimum_image_lengths(np.diff(ca,axis=0),box)
    # 8 A for contiguous receptor CA is a conservative severe-discontinuity threshold.
    # Mark split when direct is enormous but MIC is plausible (< 5 A).
    split=np.where((direct>8.0)&(mic<5.0))[0]
    severe=np.where((direct>8.0)&(mic>=5.0))[0]
    counters['ca_pbc_split_pairs']+=len(split)
    counters['ca_non_pbc_severe_pairs']+=len(severe)
    counters['max_ca_direct_A']=max(counters['max_ca_direct_A'],float(direct.max(initial=0)))
    if len(split) and len(examples['ca_pbc_split'])<12:
        for i in split[:12-len(examples['ca_pbc_split'])]:
            examples['ca_pbc_split'].append(dict(frame=frame,index=int(i),direct_A=float(direct[i]),mic_A=float(mic[i])))
    if len(severe) and len(examples['ca_non_pbc_severe'])<12:
        for i in severe[:12-len(examples['ca_non_pbc_severe'])]:
            examples['ca_non_pbc_severe'].append(dict(frame=frame,index=int(i),direct_A=float(direct[i]),mic_A=float(mic[i])))
    if prev_ca is not None and valid_box(prev_box):
        jump=ca-prev_ca
        jd=np.linalg.norm(jump,axis=1)
        # Receptor CA periodic image jump: direct > min(12 A, 35% min box)
        # and MIC < 6 A; use present-frame box. Informational, not a hard PBC proof.
        threshold=max(12.0,0.35*min(float(x) for x in box[:3]))
        suspect=np.where(jd>threshold)[0]
        if len(suspect):
            jm=minimum_image_lengths(jump[suspect],box)
            wrapped=suspect[jm<6.0]
            counters['temporal_pbc_jump_atoms']+=len(wrapped)
            if len(wrapped) and len(examples['temporal_pbc_jump'])<12:
                for k in wrapped[:12-len(examples['temporal_pbc_jump'])]:
                    examples['temporal_pbc_jump'].append(dict(frame=frame,index=int(k),direct_A=float(jd[k]),mic_A=float(jm[np.where(suspect==k)[0][0]])))
    return np.asarray(box,dtype=float)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-root',type=Path,required=True)
    ap.add_argument('--raw-root',type=Path,required=True)
    ap.add_argument('--output-root',type=Path,required=True)
    ap.add_argument('--coord-tolerance-A',type=float,default=0.05)
    ap.add_argument('--preflight-only',action='store_true')
    args=ap.parse_args()
    allr=inputs(args.source_root,args.raw_root)
    missing=[]; metadata=[]
    for rec in allr:
        for k in ('dcd','pdb','raw','manifest'):
            p=rec[k]
            if not p.is_file() or p.stat().st_size==0:missing.append(str(p))
    if missing:
        print(json.dumps({'G2A_PBC_PREFLIGHT':'FAIL','missing_count':len(set(missing)),'missing':sorted(set(missing))},indent=2))
        raise SystemExit(2)
    for rec in allr:
        s,e=check_metadata(rec); metadata.append((s,e))
    if args.output_root.exists():
        raise SystemExit(f'refusing to overwrite output path: {args.output_root}')
    print(json.dumps({'G2A_PBC_PREFLIGHT':'PASS','trajectory_files':8,'topology_files':4,'raw_windows':40,'metadata_windows':40,'output_untouched':str(args.output_root)},indent=2))
    if args.preflight_only:return

    import MDAnalysis as mda
    detail=[];per_trajectory=[]
    with tempfile.TemporaryDirectory(prefix='g2a_pbc_') as td:
        tmp=Path(td)
        for system in SYSTEMS:
            pdb=args.source_root/'topology'/system/'minimized.pdb'
            clean=prepare_topology(pdb,tmp)
            for rep in REPS:
                records=[(rec,(s,e)) for rec,(s,e) in zip(allr,metadata) if rec['system']==system and rec['replica']==rep]
                dcd=records[0][0]['dcd']
                u=mda.Universe(str(clean),str(dcd))
                sel=u.select_atoms('chainID E and protein')
                if len(sel.residues)!=270:raise ValueError(f'{system} R{rep}: selected {len(sel.residues)} residues, expected 270')
                if int((sel.names=='CA').sum())!=270:raise ValueError(f'{system} R{rep}: expected 270 CA atoms')
                max_end=max(e for _,(_,e) in records)
                if len(u.trajectory)<=max_end:raise ValueError(f'{system} R{rep}: insufficient DCD frames: {len(u.trajectory)} <= {max_end}')
                byframe={}
                for rec,(s,e) in records:
                    if set(byframe).intersection(range(s,e+1)):
                        raise ValueError('overlapping windows')
                    for idx in range(s,e+1):byframe[idx]=rec['window']
                frame_coords={w:[] for w in WINDOWS}
                ctr=dict(valid_box_frames=0,invalid_box_frames=0,ca_pbc_split_pairs=0,ca_non_pbc_severe_pairs=0,temporal_pbc_jump_atoms=0,max_ca_direct_A=0.0)
                examples={k:[] for k in ('invalid_box','ca_pbc_split','ca_non_pbc_severe','temporal_pbc_jump')}
                prev_ca=None;prev_box=None
                # All raw frames are 0-499 in frozen source metadata; also inspect window joins.
                ca_sel=sel.select_atoms('name CA')
                for ts in u.trajectory[:max_end+1]:
                    fi=int(ts.frame)
                    if fi not in byframe:continue
                    xyz=sel.positions.astype(np.float32,copy=True)
                    ca=ca_sel.positions.astype(np.float64,copy=True)
                    box=np.asarray(ts.dimensions,dtype=np.float64) if ts.dimensions is not None else None
                    analyze_frame(ca,box,fi,prev_ca,prev_box,ctr,examples)
                    prev_ca=ca;prev_box=box
                    frame_coords[byframe[fi]].append(xyz)
                if sum(len(v) for v in frame_coords.values())!=500:raise ValueError('500 window frames not observed')
                per_trajectory.append({'system':system,'replica':rep,'dcd':str(dcd),'dcd_sha256':hash_file(dcd),'pdb':str(pdb),'pdb_sha256':hash_file(pdb),'frames_total':len(u.trajectory),'frames_scanned':500,'receptor_atoms':len(sel),'receptor_residues':len(sel.residues),'pbc':ctr,'examples':examples})
                for rec,(s,e) in records:
                    arr=np.stack(frame_coords[rec['window']])
                    dmax=compare_raw(rec['raw'],sel,arr,rep,rec['window'],rec['context'],args.coord_tolerance_A)
                    detail.append(dict(system=system,context=rec['context'],replica=rep,window=rec['window'],start_frame=s,end_frame=e,raw_path=str(rec['raw']),raw_sha256=hash_file(rec['raw']),max_abs_coordinate_difference_A=dmax,raw_match='PASS'))
                print(f"AUDITED {system} R{rep}: 5 windows raw match; boxes {ctr['valid_box_frames']}/500; CA split {ctr['ca_pbc_split_pairs']}; temporal jumps {ctr['temporal_pbc_jump_atoms']}",flush=True)
    status='PASS_PBC_SCREEN_NOT_FULL_PROOF'
    if any(x['pbc']['invalid_box_frames'] for x in per_trajectory):status='PBC_BOX_INVALID'
    if any(x['pbc']['ca_pbc_split_pairs'] or x['pbc']['temporal_pbc_jump_atoms'] for x in per_trajectory):status='PBC_DISCONTINUITY_DETECTED'
    if any(x['pbc']['ca_non_pbc_severe_pairs'] for x in per_trajectory):status='GEOMETRY_DISCONTINUITY_DETECTED'
    report={'G2A_PBC_STATUS':status,'scope':'40 frozen raw windows, eight source DCDs, four PDBs; receptor CA geometry','source_root':str(args.source_root),'raw_root':str(args.raw_root),'input_count':40,'raw_matched_count':len(detail),'trajectory_count':len(per_trajectory),'trajectory_audit':per_trajectory,'window_audit':detail,'limitations':['No atom-by-atom universal PBC proof: tests adjacent receptor CA geometry and temporal CA jumps.','No reimaging/unwrap; no modification of frozen sources.','Trajectory SHA256 hashes are from copied DCD files; without original hashes cannot independently certify copy provenance.','MDAnalysis may auto-convert DCD Angstrom units; raw matching is an essential identity check.','If periodic cuts never affect adjacent CA atoms, other noncovalent distances still require explicit minimum image when relevant.']}
    args.output_root.mkdir(parents=True,exist_ok=False)
    out=args.output_root
    (out/'G2A_PBC_AUDIT.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    with (out/'G2A_PBC_WINDOWS.csv').open('w',encoding='utf-8',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(detail[0]));wr.writeheader();wr.writerows(detail)
    with (out/'G2A_PBC_TRAJECTORIES.csv').open('w',encoding='utf-8',newline='') as f:
        fields=['system','replica','frames_total','frames_scanned','receptor_atoms','receptor_residues','valid_box_frames','invalid_box_frames','ca_pbc_split_pairs','ca_non_pbc_severe_pairs','temporal_pbc_jump_atoms','max_ca_direct_A','dcd_sha256','pdb_sha256']
        wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader()
        for x in per_trajectory:wr.writerow({**{k:v for k,v in x.items() if k!='pbc'},**x['pbc']})
    (out/'G2A_PBC_REPORT.md').write_text('# G2-A independent PBC audit\n\n'+f'**Status:** {status}\n\n40 frozen raw windows were compared to eight copied original DCD trajectories.\n\n'+ '\n'.join(f"- {x['system']} R{x['replica']}: boxes {x['pbc']['valid_box_frames']}/500; adjacent-CA PBC splits {x['pbc']['ca_pbc_split_pairs']}; non-PBC severe CA pairs {x['pbc']['ca_non_pbc_severe_pairs']}; temporal wrapping events {x['pbc']['temporal_pbc_jump_atoms']}." for x in per_trajectory)+'\n\nSee audit JSON for examples and limitations.\n',encoding='utf-8')
    print(json.dumps({'G2A_PBC_STATUS':status,'raw_matched':len(detail),'trajectories':len(per_trajectory),'output':str(out)},indent=2))

if __name__=='__main__':main()
