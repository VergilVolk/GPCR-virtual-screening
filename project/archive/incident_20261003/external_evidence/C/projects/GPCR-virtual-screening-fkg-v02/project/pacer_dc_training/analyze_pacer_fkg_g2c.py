#!/usr/bin/env python
"""PACER-FKG G2-C: frozen graph ablation of *signed* node-level RFF contrasts.

G2-B regional RFF maps have different dimensions/bandwidths, so are NOT graph
compatible. This experiment defines a new R2-calibrated, shared atom-level RBF
kernel over 128-dimensional residue descriptors; the graph propagates SIGNED
shared-RFF contrast vectors, never nonnegative kernel distances. Descriptive.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix

CONTEXTS = ('apo','probe_only','candidate_no_probe','candidate_probe')
SIGNS = {
 'synergy_interaction': {'candidate_probe':1,'probe_only':-1,'candidate_no_probe':-1,'apo':1},
 'intrinsic_agonism': {'candidate_no_probe':1,'apo':-1},
 'conditional_pam_effect': {'candidate_probe':1,'probe_only':-1},
}

def sha256(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''): h.update(block)
 return h.hexdigest()

def path_for(root,r,w,c):
 return root/f'replica_{r:02d}'/f'window_{w:03d}'/'atom14'/f'{c}_w{w:03d}.geom2vec.npz'

def norm(a,b):
 n=float(np.linalg.norm(a)*np.linalg.norm(b))
 return float(np.clip(np.sum(a*b)/n,-1,1)) if n>1.e-10 else None

def diff_graph(v,transition,alpha=.65,steps=20):
 """v: (n_nodes,n_rff) signed contrasts; row-stochastic propagation."""
 result=np.array(v,copy=True)
 for _ in range(steps): result=(1-alpha)*v+alpha*(transition@result)
 return result

def graph_edges(graph):
 edges={}
 for e in graph['edges']:
  i,j=int(e['source']),int(e['target'])
  if i==j: continue
  pair=tuple(sorted((i,j)))
  weight=1. if e['edge_type']=='backbone' else float(e['support_fraction'])
  if pair in edges:
   old=edges[pair]
   if old[0]!='backbone' and e['edge_type']=='backbone':edges[pair]=('backbone',1.)
   elif old[0]!='backbone':edges[pair]=('contact',max(old[1],weight))
  else:edges[pair]=('backbone' if e['edge_type']=='backbone' else 'contact',weight)
 return edges

def permute_contacts(edges,rng,swaps=20_000):
 """Degree-preserving contact-edge rewiring; backbone is completely fixed.

Each contact edge carries its source weight. The degree sequence is fixed,
not the per-node weighted strength; this is explicitly disclosed in report.
"""
 backbone={x for x,(kind,_) in edges.items() if kind=='backbone'}
 contacts=[(a,b,w) for (a,b),(kind,w) in edges.items() if kind=='contact']
 taken=backbone|{(a,b) for a,b,_ in contacts}
 if len(contacts)<2: raise ValueError('insufficient contact edges')
 accepted=0
 for _ in range(swaps):
  i,j=rng.choice(len(contacts),size=2,replace=False)
  a,b,wi=contacts[i]; c,d,wj=contacts[j]
  if rng.integers(2): c,d=d,c
  if len({a,b,c,d})<4: continue
  x,y=tuple(sorted((a,d))),tuple(sorted((c,b)))
  if x==y or x in taken or y in taken:continue
  taken.remove(tuple(sorted((a,b))));taken.remove(tuple(sorted((c,d))))
  taken.add(x);taken.add(y)
  contacts[i]=(x[0],x[1],wi);contacts[j]=(y[0],y[1],wj)
  accepted+=1
 if accepted==0:raise ValueError('contact rewiring failed')
 out={e:('backbone',1.) for e in backbone}
 out.update({(a,b):('contact',w) for a,b,w in contacts})
 return out,accepted

def transition(edges,n=270):
 row=[];col=[];dat=[]
 for (i,j),(_,w) in edges.items():
  row.extend((i,j));col.extend((j,i));dat.extend((w,w))
 matrix=csr_matrix((dat,(row,col)),shape=(n,n),dtype=np.float64)
 degree=np.asarray(matrix.sum(axis=1)).ravel()
 if (degree<=0).any(): raise ValueError('graph has isolated nodes')
 return csr_matrix((1/degree,(np.arange(n),np.arange(n))),shape=(n,n))@matrix

def get_graph_region(graph,region):
 return np.asarray([int(m['embedding_index']) for m in graph['regions'][region]],dtype=int)

def preflight(args):
 if args.output_root.exists():raise FileExistsError(f'refusing existing output directory: {args.output_root}')
 if args.rff_dim<32 or args.n_null<2 or args.block_frames<1 or 100%args.block_frames or args.bootstrap<1:
  raise ValueError('invalid RFF, null count, block length or bootstrap settings')
 for p in (args.graph,args.fkg_audit,args.g2b_audit):
  if not p.is_file():raise FileNotFoundError(p)
 graph=json.loads(args.graph.read_text(encoding='utf-8'))
 fkg=json.loads(args.fkg_audit.read_text(encoding='utf-8'))
 g2b=json.loads(args.g2b_audit.read_text(encoding='utf-8'))
 gh=sha256(args.graph); fh=sha256(args.fkg_audit)
 if len(graph['nodes'])!=270 or len(graph['regions'])!=9 or gh!=fkg['provenance']['graph_sha256']:
  raise ValueError('frozen graph mismatch')
 if g2b.get('G2B_STATUS')!='COMPLETED_DESCRIPTIVE_NOT_INFERENTIAL' or g2b['provenance']['graph_sha256']!=gh or g2b['provenance']['fkg_audit_sha256']!=fh:
  raise ValueError('frozen G2B provenance mismatch')
 med=np.asarray(g2b['channel_median'],dtype=np.float64)
 scale=np.asarray(g2b['channel_scale'],dtype=np.float64)
 if med.shape!=(128,) or scale.shape!=(128,) or (scale<=0).any():raise ValueError('bad frozen preprocessing')
 manifest={m['path']:m['sha256'] for m in fkg['provenance']['input_files']}
 if len(manifest)!=40:raise ValueError('expected 40 frozen inputs')
 for r in (2,3):
  for w in range(5):
   matched=None
   for c in CONTEXTS:
    p=path_for(args.input_root,r,w,c);rel=p.relative_to(args.input_root).as_posix()
    if not p.is_file() or sha256(p)!=manifest.get(rel):raise ValueError(f'input hash mismatch: {p}')
    with np.load(p,allow_pickle=False) as z:
     a=z['residue_features'];ids=z['frame_ids'];seq=str(z['sequence'].item())
     if a.shape!=(100,270,128) or not np.isfinite(a).all():raise ValueError(f'bad embedding: {p}')
     if seq!=fkg['provenance']['sequence']:raise ValueError(f'sequence mismatch: {p}')
     if ids.shape!=(100,) or np.any(np.diff(ids)<=0):raise ValueError(f'bad frames: {p}')
     if matched is None:matched=ids
     elif not np.array_equal(matched,ids):raise ValueError(f'context frames mismatch: {p}')
 es=graph_edges(graph)
 return graph,g2b,med,scale,es,{'input_hashes_verified':40,'graph_sha256':gh,'fkg_audit_sha256':fh,'g2b_audit_sha256':sha256(args.g2b_audit),'graph_edges_unique':len(es),'g2b_calibration_replica':g2b['settings']['calibration_replica']}

def fit_shared_width(args,med,scale):
 """R2 only, four contexts/five windows, deterministic balanced frame/node subsample."""
 rng=np.random.default_rng(args.seed)
 samples=[]
 nodeids=np.linspace(0,269,32,dtype=int)
 frameids=np.linspace(0,99,8,dtype=int)
 for w in range(5):
  for c in CONTEXTS:
   with np.load(path_for(args.input_root,2,w,c),allow_pickle=False) as z:
    a=z['residue_features'][frameids].astype(np.float64)
   a-=a.mean(axis=1,keepdims=True)
   samples.append(((a-med)/scale)[:,nodeids,:].reshape(-1,128))
 x=np.concatenate(samples,axis=0)
 # Fit a single 128D bandwidth; no per-node/per-region refit and no R3 leakage.
 idx=rng.choice(len(x),size=min(len(x),1000),replace=False)
 v=x[idx];sq=np.sum(v*v,axis=1)
 d2=np.maximum(sq[:,None]+sq[None,:]-2*v@v.T,0)
 upper=d2[np.triu_indices(len(v),1)];positive=upper[upper>1e-12]
 if not len(positive):raise ValueError('degenerate node bandwidth')
 return float(np.sqrt(np.median(positive))),len(v)

def rff_map(a,weights,bias):
 x=np.asarray(a,dtype=np.float32).reshape(-1,128)
 z=np.sqrt(2./len(bias))*np.cos(x@weights+bias)
 return z.reshape(100,270,len(bias)).astype(np.float32)

def load_contrasts(args,med,scale,width):
 rng=np.random.default_rng(args.seed+771)
 weights=rng.normal(0,1/width,size=(128,args.rff_dim)).astype(np.float32)
 bias=rng.uniform(0,2*np.pi,size=args.rff_dim).astype(np.float32)
 # Each context tensor: five windows x five contiguous blocks x 270 x RFF.
 blocks={r:{c:[] for c in CONTEXTS} for r in (2,3)}
 for r in (2,3):
  for w in range(5):
   for c in CONTEXTS:
    with np.load(path_for(args.input_root,r,w,c),allow_pickle=False) as z:a=z['residue_features'].astype(np.float32)
    a-=a.mean(axis=1,keepdims=True)
    a=(a-med)/scale
    mapped=rff_map(a,weights,bias)
    blocks[r][c].append(mapped.reshape(100//args.block_frames,args.block_frames,270,args.rff_dim).mean(axis=1))
   print(f'G2C encoded R{r} W{w} four contexts',flush=True)
 for r in (2,3):
  for c in CONTEXTS:blocks[r][c]=np.stack(blocks[r][c])
 return {axis:{r:sum(sign*blocks[r][c] for c,sign in signs.items()) for r in (2,3)} for axis,signs in SIGNS.items()}

def regional_cosines(v,graph):
 return {name:norm(v[2][get_graph_region(graph,name)].mean(axis=0),v[3][get_graph_region(graph,name)].mean(axis=0)) for name in sorted(graph['regions'])}

def region_vectors(node_contrast,indices):return node_contrast[indices].mean(axis=0)

def percentile(a,p):return float(np.percentile(a,p)) if len(a) else None

def bootstrap_pair(no,real,idx,n_draws,seed):
 """Paired shared-window/block resampling for true-minus-no-graph contrast; descriptive."""
 rng=np.random.default_rng(seed);deltas=[];base=[];true=[]
 # Slice the region once; selecting all 270 residues in every draw is wasteful.
 no={r:no[r][:,:,idx,:] for r in (2,3)}
 real={r:real[r][:,:,idx,:] for r in (2,3)}
 for _ in range(n_draws):
  res={}
  for name,collection in [('no',no),('real',real)]:
   res[name]={}
  # Common sampled windows/blocks across models and across both replicas.
  for r in (2,3):
   wins=rng.integers(0,5,size=5)
   bi=rng.integers(0,no[r].shape[1],size=(5,no[r].shape[1]))
   for name,collection in [('no',no),('real',real)]:
    selected=collection[r][wins[:,None],bi]
    res[name][r]=selected.mean(axis=(0,1,2))
  u=norm(res['no'][2],res['no'][3]);v=norm(res['real'][2],res['real'][3])
  if u is not None and v is not None:base.append(u);true.append(v);deltas.append(v-u)
 return {'draws_valid':len(deltas),'delta_cosine_ci025':percentile(deltas,2.5),
         'delta_cosine_ci975':percentile(deltas,97.5),'delta_cosine_median':percentile(deltas,50),
         'true_cosine_ci025':percentile(true,2.5),'true_cosine_ci975':percentile(true,97.5)}

def analyze(args,graph,med,scale,edges,provenance):
 width,ncal=fit_shared_width(args,med,scale)
 groups=load_contrasts(args,med,scale,width)
 real=transition(edges)
 rng=np.random.default_rng(args.seed+1007)
 nulls=[];swapcounts=[]
 for j in range(args.n_null):
  rewired,accepted=permute_contacts(edges,rng,swaps=args.null_swaps)
  nulls.append(transition(rewired));swapcounts.append(accepted)
 nodes_by_region={name:get_graph_region(graph,name) for name in sorted(graph['regions'])}
 rows=[];matched=[];nrows=[]
 for ai,(axis,repblocks) in enumerate(groups.items()):
  node={r:repblocks[r].mean(axis=(0,1)).astype(np.float64) for r in (2,3)}
  graphnode={r:diff_graph(node[r],real,args.alpha,args.steps) for r in (2,3)}
  nullnode=[{r:diff_graph(node[r],m,args.alpha,args.steps) for r in (2,3)} for m in nulls]
  # Precompute block-level graph maps for paired sensitivity, not null graphs.
  graphblocks={r:np.stack([diff_graph(z,real,args.alpha,args.steps) for window in repblocks[r] for z in window]).reshape(repblocks[r].shape) for r in (2,3)}
  for ri,(name,idx) in enumerate(nodes_by_region.items()):
   none=[region_vectors(node[r],idx) for r in (2,3)]
   true=[region_vectors(graphnode[r],idx) for r in (2,3)]
   ncos=norm(*none);tcos=norm(*true)
   shuffled=[norm(region_vectors(n[2],idx),region_vectors(n[3],idx)) for n in nullnode]
   valid=[x for x in shuffled if x is not None]
   delta=tcos-ncos if tcos is not None and ncos is not None else None
   relative=[tcos-x for x in valid] if tcos is not None else []
   boot=bootstrap_pair(repblocks,graphblocks,idx,args.bootstrap,args.seed+ai*100+ri)
   row={'region':name,'axis':axis,'n_residues':len(idx),'no_graph_cosine':ncos,'true_graph_cosine':tcos,
        'true_minus_no_cosine':delta,'r2_no_norm':float(np.linalg.norm(none[0])),'r3_no_norm':float(np.linalg.norm(none[1])),
        'r2_true_norm':float(np.linalg.norm(true[0])),'r3_true_norm':float(np.linalg.norm(true[1])),
        'null_median_cosine':percentile(valid,50),'null_q025_cosine':percentile(valid,2.5),
        'null_q975_cosine':percentile(valid,97.5),'true_minus_null_median':percentile(relative,50),
        'null_count_valid':len(valid),**boot}
   rows.append(row)
   for j,value in enumerate(shuffled):nrows.append({'region':name,'axis':axis,'null_id':j,'cosine':value,'accepted_edge_swaps':swapcounts[j]})
   for w in range(5):
    a=[repblocks[r][w].mean(axis=0)[idx].mean(axis=0) for r in (2,3)]
    b=[graphblocks[r][w].mean(axis=0)[idx].mean(axis=0) for r in (2,3)]
    matched.append({'region':name,'axis':axis,'window':w,'no_graph_cosine':norm(*a),'true_graph_cosine':norm(*b)})
 return rows,matched,nrows,{'shared_node_bandwidth':width,'node_calibration_samples':ncal,'rff_dim':args.rff_dim,
       'rff_seed':args.seed+771,'n_null':args.n_null,'swaps_accepted':swapcounts,'alpha':args.alpha,'steps':args.steps,
       'block_frames':args.block_frames,'bootstrap':args.bootstrap},provenance

def write_csv(path,rows):
 if not rows:raise ValueError('empty CSV')
 fields=list(rows[0])
 with path.open('w',encoding='utf-8-sig',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(rows)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--input-root',type=Path,required=True)
 p.add_argument('--graph',type=Path,required=True)
 p.add_argument('--fkg-audit',type=Path,required=True)
 p.add_argument('--g2b-audit',type=Path,required=True)
 p.add_argument('--output-root',type=Path,required=True)
 p.add_argument('--rff-dim',type=int,default=256)
 p.add_argument('--seed',type=int,default=271828)
 p.add_argument('--n-null',type=int,default=12)
 p.add_argument('--null-swaps',type=int,default=12000)
 p.add_argument('--alpha',type=float,default=.65)
 p.add_argument('--steps',type=int,default=20)
 p.add_argument('--block-frames',type=int,default=20)
 p.add_argument('--bootstrap',type=int,default=200)
 p.add_argument('--preflight-only',action='store_true')
 args=p.parse_args()
 if not (0<args.alpha<1) or args.steps<1 or args.null_swaps<100:raise ValueError('invalid diffusion settings')
 graph,old,med,scale,edges,provenance=preflight(args)
 if args.preflight_only:
  print(json.dumps({'G2C_PREFLIGHT':'PASS',**provenance,'n_regions':len(graph['regions']),
    'n_axes':len(SIGNS),'n_null':args.n_null,'rff_dim':args.rff_dim,'output_untouched':str(args.output_root)},indent=2));return
 rows,windows,nullrows,cal,_=analyze(args,graph,med,scale,edges,provenance)
 if len(rows)!=27 or len(windows)!=135 or len(nullrows)!=27*args.n_null:raise ValueError('unexpected row counts')
 report={'G2C_STATUS':'COMPLETED_DESCRIPTIVE_NOT_INFERENTIAL','method':'signed shared-node RFF contrast diffusion',
  'provenance':provenance,'calibration':cal,'region_axis_summary':rows,
  'limitations':['New shared 128D node kernel, NOT the per-region G2-B RFF estimand; across-stage cosine magnitudes must not be compared.',
  'Signed RFF node deltas are graph-diffused; not the nonnegative FKG sqrt-U2 scores.',
  'Contact edges rewired with unweighted degree preserved and backbone fixed; node weighted strengths are NOT preserved.',
  'Same R2-fitted channel median/MAD as frozen G2-B; shared node bandwidth fitted on R2 only.',
  'Null controls are descriptive conditional on one frozen graph; no p-values or randomization inference.',
  'Only two replica trajectories, five serially correlated windows and overlapping regions; resampling is sensitivity only.',
  'G2-C does not certify PAM efficacy or graph causal mechanism; G2-B remains frozen.']}
 args.output_root.mkdir(parents=True,exist_ok=False)
 (args.output_root/'G2C_GRAPH_ABLATION_AUDIT.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
 write_csv(args.output_root/'G2C_REGION_AXIS.csv',rows)
 write_csv(args.output_root/'G2C_MATCHED_WINDOWS.csv',windows)
 write_csv(args.output_root/'G2C_SHUFFLED_GRAPHS.csv',nullrows)
 lines=['# PACER-FKG G2-C signed graph ablation','','Descriptive only; no gate update, efficacy or causal graph claim.','',
        f'Inputs verified: {provenance["input_hashes_verified"]}/40; shared-node RFF dimension: {args.rff_dim}; shuffled graphs: {args.n_null}.','',
        '| Region | Axis | no graph | frozen graph | null median | true-minus-no | descriptive delta block interval |',
        '|---|---|---:|---:|---:|---:|---:|']
 for r in rows:
  fmt=lambda x:f'{x:.3f}' if x is not None else 'N/A'
  ci=fmt(r['delta_cosine_ci025'])+' to '+fmt(r['delta_cosine_ci975'])
  lines.append(f'| {r["region"]} | {r["axis"]} | {fmt(r["no_graph_cosine"])} | {fmt(r["true_graph_cosine"])} | {fmt(r["null_median_cosine"])} | {fmt(r["true_minus_no_cosine"])} | {ci} |')
 lines+=['','## Interpretation boundaries','',*('- '+v for v in report['limitations'])]
 (args.output_root/'G2C_GRAPH_ABLATION_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 print(json.dumps({'G2C_STATUS':report['G2C_STATUS'],'inputs':40,'region_axis_rows':len(rows),
  'matched_windows':len(windows),'null_records':len(nullrows),'outputs':5,'output':str(args.output_root)},indent=2))
if __name__=='__main__':main()
