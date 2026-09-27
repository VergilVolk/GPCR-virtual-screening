#!/usr/bin/env python3
"""G2-A independent, read-only input audit for PACER-FKG Geom2Vec R2/R3.

The upstream raw NPZ archive does not contain periodic cell dimensions. Geometry
checks therefore SCREEN for broken imaging, but cannot certify a full PBC audit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

# Explicit independent atom14 reference, not imported from the extractor under test.
ATOM14 = {
    'A': 'N CA C O CB', 'R': 'N CA C O CB CG CD NE CZ NH1 NH2',
    'N': 'N CA C O CB CG OD1 ND2', 'D': 'N CA C O CB CG OD1 OD2',
    'C': 'N CA C O CB SG', 'Q': 'N CA C O CB CG CD OE1 NE2',
    'E': 'N CA C O CB CG CD OE1 OE2', 'G': 'N CA C O',
    'H': 'N CA C O CB CG ND1 CD2 CE1 NE2',
    'I': 'N CA C O CB CG1 CG2 CD1', 'L': 'N CA C O CB CG CD1 CD2',
    'K': 'N CA C O CB CG CD CE NZ', 'M': 'N CA C O CB CG SD CE',
    'F': 'N CA C O CB CG CD1 CD2 CE1 CE2 CZ',
    'P': 'N CA C O CB CG CD', 'S': 'N CA C O CB OG',
    'T': 'N CA C O CB OG1 CG2',
    'W': 'N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'Y': 'N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH',
    'V': 'N CA C O CB CG1 CG2',
}
ATOM14 = {k: v.split() for k, v in ATOM14.items()}
THREE = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),
                 'ARNDCQEGHILKMFPSTWYV'))
RENAMES = {'HSE':'HIS','HSD':'HIS','HSP':'HIS'}
CONTEXTS = ('apo','probe_only','candidate_no_probe','candidate_probe')
BOND_LIMITS = {'N_CA':2.5,'CA_C':2.5,'C_O':2.0,'PEPTIDE_C_N':3.0,'CA_NEIGHBOR':8.0}


def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()


def read_sequence(path):
    with path.open(newline='',encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=1 or 'seqres' not in rows[0]: raise ValueError('missing or ambiguous seqres CSV')
    return rows[0]['seqres'].strip()


def inspect_reference_constants(path):
    """Independently load the frozen official mdgen table and compare all 20 AAs.

    This executes the *user-supplied local* residue_constants.py. The file
    belongs to their pinned mdgen checkout, not external/untrusted downloads.
    """
    if path is None:return {'status':'NOT_SUPPLIED'}
    p=Path(path)
    if not p.is_file():return {'status':'MISSING','path':str(p)}
    try:
        spec=importlib.util.spec_from_file_location('g2a_mdgen_reference',p)
        rc=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(rc)
        differences={}
        for one,names in ATOM14.items():
            three=next(k for k,v in THREE.items() if v==one)
            actual=[x for x in rc.restype_name_to_atom14_names[three] if x]
            if names!=actual:differences[one]={'independent':names,'mdgen':actual}
        return {'status':'MATCH' if not differences else 'MISMATCH',
                'path':str(p),'sha256':sha256(p),'differences':differences}
    except Exception as e:
        return {'status':'LOAD_FAILED','path':str(p),'sha256':sha256(p),
                'error':f'{type(e).__name__}: {e}'}


def validate_raw_mapping(atom14, raw, sequence, atol=1e-3):
    """Reconstruct mapping solely from raw atom identities and independent table."""
    coords=np.asarray(raw['coords'])
    names=np.asarray(raw['atom_names']).astype(str)
    resnames=np.asarray(raw['resnames']).astype(str)
    resids=np.asarray(raw['resids']) if 'resids' in raw else None
    ridx=np.asarray(raw['res_index'])
    if coords.ndim!=3 or coords.shape[2]!=3 or len(names)!=coords.shape[1] or len(ridx)!=len(names):
        raise ValueError('bad raw coordinate/name/index shapes')
    if atom14.shape != (coords.shape[0],len(sequence),14,3):
        raise ValueError(f'raw/atom14 dimension mismatch: {coords.shape} {atom14.shape}')
    if len(resnames)!=len(sequence) or not np.array_equal(np.unique(ridx),np.arange(len(sequence))):
        raise ValueError('residue index not consecutive or residue count mismatch')
    if resids is not None and len(resids)!=len(sequence): raise ValueError('bad resids length')
    if not np.isfinite(coords).all() or not np.isfinite(atom14).all():
        raise ValueError('nonfinite coordinates')
    gotseq=''.join(THREE.get(RENAMES.get(n,n),'?') for n in resnames)
    if gotseq!=sequence:raise ValueError('raw residue identity/sequence mismatch')
    valid=np.zeros((len(sequence),14),dtype=bool)
    ca_index=np.full(len(sequence),-1,dtype=int)
    c_index=np.full(len(sequence),-1,dtype=int)
    n_index=np.full(len(sequence),-1,dtype=int)
    o_index=np.full(len(sequence),-1,dtype=int)
    duplicates=[]
    max_error=0.0; mismatch_count=0; wrong_mask_count=0
    for r,aa in enumerate(sequence):
        atomidx=np.where(ridx==r)[0]
        ren={'CD':'CD1'} if aa=='I' else {}
        loc={}
        for i in atomidx:
            nm=ren.get(names[i],names[i])
            if nm in loc:duplicates.append([r,nm])
            loc[nm]=i
        for nm,storage in [('CA',ca_index),('C',c_index),('N',n_index),('O',o_index)]: storage[r]=loc.get(nm,-1)
        for s,nm in enumerate(ATOM14[aa]):
            i=loc.get(nm)
            if i is None:
                # Missing atom is a zero-padded slot in the frozen converter.
                n_bad=int(np.count_nonzero(np.any(np.abs(atom14[:,r,s,:])>atol,axis=-1)))
                wrong_mask_count+=n_bad
            else:
                valid[r,s]=True
                err=np.max(np.abs(atom14[:,r,s,:]-coords[:,i,:]),axis=1)
                mismatch_count+=int(np.count_nonzero(err>atol))
                max_error=max(max_error,float(err.max()))
        for s in range(len(ATOM14[aa]),14):
            wrong_mask_count+=int(np.count_nonzero(np.any(np.abs(atom14[:,r,s,:])>atol,axis=-1)))
    if duplicates:raise ValueError(f'duplicate canonical names, first={duplicates[:4]}')
    if np.any(ca_index<0) or np.any(c_index<0) or np.any(n_index<0) or np.any(o_index<0):
        raise ValueError('one or more backbone atoms missing')
    return {'valid_mask':valid,'ca':ca_index,'c':c_index,'n':n_index,'o':o_index,
            'max_mapping_abs_A':max_error,'mapping_mismatch_slots_frames':mismatch_count,
            'incorrect_missing_slot_frames':wrong_mask_count,'coords':coords}


def geometry_screen(mapping):
    c=mapping['coords']; n=mapping['n']; ca=mapping['ca']; carb=mapping['c']; o=mapping['o']
    def measure(a,b):return np.linalg.norm(a-b,axis=-1)
    b={
      'N_CA':measure(c[:,n,:],c[:,ca,:]),
      'CA_C':measure(c[:,ca,:],c[:,carb,:]),
      'C_O':measure(c[:,carb,:],c[:,o,:]),
      'PEPTIDE_C_N':measure(c[:,carb[:-1],:],c[:,n[1:],:]),
      'CA_NEIGHBOR':measure(c[:,ca[:-1],:],c[:,ca[1:],:]),
    }
    # Per-frame centering removes whole-protein translation; large internal
    # jumps may flag imaging, but flexible loops can also move.
    x=c[:,ca,:]
    xc=x-np.mean(x,axis=1,keepdims=True)
    jumps=measure(xc[1:],xc[:-1]) if len(xc)>1 else np.array([0.])
    out={}
    for k,v in b.items():
        out[k]={'max_A':float(v.max()),'violations':int(np.count_nonzero(v>BOND_LIMITS[k]))}
    out['centered_CA_temporal']={'max_A':float(jumps.max()),'jumps_gt_15A':int(np.count_nonzero(jumps>15.0))}
    out['geometry_flags']=sum(x['violations'] for k,x in out.items() if 'violations' in x)+out['centered_CA_temporal']['jumps_gt_15A']
    return out


def process_row(root,rawroot,featroot,row,expect_seq):
    rel=Path(row['input']); path=root/rel
    if not path.is_file():return {'status':'MISSING_ATOM14','input':row['input']}
    actual_hash=sha256(path)
    result={'input':row['input'],'atom14_sha256':actual_hash,'hash_ok':actual_hash==row['input_sha256']}
    if not result['hash_ok']:
        result['status']='FAIL_HASH';return result
    a=np.load(path,mmap_mode='r')
    seq=read_sequence(path.with_suffix('.csv'))
    if seq!=expect_seq or a.shape!=(100,len(expect_seq),14,3):
        result.update(status='FAIL_SEQUENCE_OR_SHAPE',shape=list(a.shape),seq_length=len(seq));return result
    feat=featroot/Path(row['output'])
    if not feat.is_file(): result.update(status='MISSING_FEATURES');return result
    with np.load(feat) as z:
        f=np.asarray(z['frame_ids'])
        fseq=str(np.asarray(z['sequence']).item())
        feature_shape=list(z['residue_features'].shape)
    if not np.array_equal(f,np.arange(100)) or fseq!=seq or feature_shape!=[100,len(seq),128]:
        result.update(status='FAIL_FEATURE_ALIGNMENT',feature_shape=feature_shape);return result
    rawpath=rawroot/rel.parent.parent/'raw'/rel.name.replace('.npy','.npz')
    if not rawpath.is_file():
        result.update(status='RAW_MISSING',raw_expected=str(rawpath));return result
    try:
        with np.load(rawpath) as z:
            m=validate_raw_mapping(a,z,seq)
            dim=np.asarray(z['dimensions']) if 'dimensions' in z else None
        geom=geometry_screen(m)
        result.update(mapping_max_error_A=m['max_mapping_abs_A'],
            mismatch_slots_frames=m['mapping_mismatch_slots_frames'],
            missing_slot_errors=m['incorrect_missing_slot_frames'],
            mask_sha256=hashlib.sha256(m['valid_mask'].tobytes()).hexdigest(),
            geometry=geom,raw_sha256=sha256(rawpath),
            box_status='DIMENSIONS_PRESENT_NOT_INDEPENDENTLY_VERIFIED' if dim is not None else 'BOX_MISSING',
            status=('PASS_GEOMETRY_SCREEN' if geom['geometry_flags']==0 and m['mapping_mismatch_slots_frames']==0 and m['incorrect_missing_slot_frames']==0 else 'FAIL_GEOMETRY_OR_MAPPING'))
    except Exception as e:
        result.update(status='FAIL_RAW_CHECK',reason=str(e))
    return result


def run(args):
    audit=json.loads((args.features_root/'batch_audit.json').read_text(encoding='utf-8'))
    rows=audit.get('rows',[])
    if len(rows)!=40 or audit.get('files_ok')!=40:raise ValueError('expected exactly 40 manifest rows')
    report=json.loads(args.fkg_report.read_text(encoding='utf-8'))
    seq=report['provenance']['sequence']
    parsed=[]
    for row in rows:
        m=re.fullmatch(r'replica_(02|03)/window_00([0-4])/atom14/(apo|probe_only|candidate_no_probe|candidate_probe)_w00([0-4])\.npy',row['input'].replace('\\','/'))
        if not m or m[2]!=m[4]:raise ValueError('noncanonical or mismatched manifest path: '+row['input'])
        parsed.append((int(m[1]),int(m[2]),m[3]))
    expected={(r,w,c) for r in (2,3) for w in range(5) for c in CONTEXTS}
    if set(parsed)!=expected or len(parsed)!=len(expected):raise ValueError('duplicate/missing unit in manifest')
    # Preflight is strictly read-only and does not create output directories.
    expected_paths=[args.atom14_root/Path(row['input']) for row in rows]
    expected_raw=[args.raw_root/Path(row['input']).parent.parent/'raw'/Path(row['input']).name.replace('.npy','.npz') for row in rows]
    available={'atom14':sum(p.is_file() for p in expected_paths),'raw':sum(p.is_file() for p in expected_raw),
               'features':sum((args.features_root/Path(row['output'])).is_file() for row in rows)}
    if args.preflight_only:
        print(json.dumps({'G2A_PREFLIGHT':'READY' if all(v==40 for v in available.values()) else 'INPUTS_INCOMPLETE',
                          'available':available,'expected':40,'output_untouched':str(args.output_root)},ensure_ascii=False,indent=2))
        return
    if args.output_root.exists():raise FileExistsError('refusing to overwrite existing G2-A output')
    results=[]; mask_group=defaultdict(set)
    for row,(r,w,c) in zip(rows,parsed):
        v=process_row(args.atom14_root,args.raw_root,args.features_root,row,seq)
        v.update(replica=r,window=w,context=c)
        results.append(v)
        if 'mask_sha256' in v:mask_group[(r,w)].add(v['mask_sha256'])
    masks_mismatch=[f'R{r}_W{w:03d}' for (r,w),v in mask_group.items() if len(v)!=1]
    statuses=defaultdict(int)
    for row in results: statuses[row['status']]+=1
    mdgen_reference=inspect_reference_constants(args.mdgen_constants)
    all_masks={x['mask_sha256'] for x in results if 'mask_sha256' in x}
    masks_cross_context_and_windows=(len(all_masks)==1)
    geometry_valid=(statuses['PASS_GEOMETRY_SCREEN']==40 and not masks_mismatch and masks_cross_context_and_windows)
    if not geometry_valid: overall='FAIL_OR_INCOMPLETE'
    elif mdgen_reference['status']=='MATCH': overall='PARTIAL_PBC_UNVERIFIED'
    elif mdgen_reference['status']=='MISMATCH': overall='FAIL_MDGEN_ATOM_ORDER'
    else: overall='REFERENCE_UNVERIFIED_PBC_UNVERIFIED'
    output={'method':'G2-A independent atom14/raw identity + PBC-geometry screen',
      'status':overall,
      'atom14_geometry_screen':'PASS' if geometry_valid else 'FAIL_OR_INCOMPLETE',
      'full_pbc_validation':'NOT_ESTABLISHED_NO_CELL_TRAJECTORY_RECONSTRUCTION',
      'count':len(results),'status_counts':dict(statuses),'same_masks_within_units':not masks_mismatch,
      'mask_mismatch_units':masks_mismatch,'cross_window_context_mask_match':masks_cross_context_and_windows,'upstream_batch_audit_sha256':sha256(args.features_root/'batch_audit.json'),
      'fkg_report_sha256':sha256(args.fkg_report),
      'independent_atom14_reference':'hardcoded audit reference; official mdgen source required for full geometry-screen pass',
      'mdgen_reference_source':mdgen_reference,
      'cutoffs_A':BOND_LIMITS,'rows':results,
      'claim_boundary':'Geometry screens can reject suspect imaging, not prove proper PBC treatment without original unit-cell trajectory, ligand and protein image review.'}
    args.output_root.mkdir(parents=True,exist_ok=False)
    (args.output_root/'G2A_AUDIT.json').write_text(json.dumps(output,indent=2,ensure_ascii=False),encoding='utf-8')
    with (args.output_root/'G2A_PER_INPUT.csv').open('w',encoding='utf-8-sig',newline='') as f:
        keys=['replica','window','context','input','status','hash_ok','mapping_max_error_A','mismatch_slots_frames','missing_slot_errors','mask_sha256','box_status']
        wr=csv.DictWriter(f,keys,extrasaction='ignore');wr.writeheader();wr.writerows(results)
    lines=['# PACER-FKG G2-A · atom14 / PBC independent QC','',f"Overall: **{output['status']}**",
           f"Atom14/raw mapping and geometry screen: **{output['atom14_geometry_screen']}**",
           f"Status counts: `{json.dumps(dict(statuses),sort_keys=True)}`",
           f"Four-context matched residue masks: **{'PASS' if not masks_mismatch and masks_cross_context_and_windows else 'FAIL'}**",
           f"Official mdgen atom14 mapping match: **{mdgen_reference['status']}**",
           '','**PBC boundary:** Raw windows currently lack original cell/imaging provenance. Passing this screen is not proof of fully correct periodic handling.',
           '','To complete G2-A: independently compare official mdgen atom14 table and validate original box-aware DCD + topology, including receptor–ligand imaging.']
    (args.output_root/'G2A_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({'G2A_STATUS':output['status'],'counts':dict(statuses),
         'mask_mismatches':masks_mismatch,'output':str(args.output_root)},ensure_ascii=False,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--atom14-root',type=Path,required=True,help='root with replica_02/window_000/atom14/...')
    p.add_argument('--raw-root',type=Path,help='root with replica_02/window_000/raw/...; defaults to atom14-root')
    p.add_argument('--features-root',type=Path,default=Path('project/results/pacer_dc_geom2vec_R2R3_full_v01'))
    p.add_argument('--fkg-report',type=Path,default=Path('project/results/pacer_dc_fkg_R2R3_v01/PACER_FKG_R2R3_AUDIT.json'))
    p.add_argument('--mdgen-constants',type=Path,help='optional frozen mdgen residue_constants.py provenance')
    p.add_argument('--output-root',type=Path,default=Path('project/results/pacer_dc_fkg_G2A_v01'))
    p.add_argument('--preflight-only',action='store_true')
    a=p.parse_args();a.raw_root=a.raw_root or a.atom14_root
    run(a)

if __name__=='__main__':main()
