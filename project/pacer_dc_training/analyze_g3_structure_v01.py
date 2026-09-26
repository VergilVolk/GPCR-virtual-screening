import json
import numpy as np
from pathlib import Path

p=Path('project/results/pacer_dc_four_context_v01/compound110')
regions=json.loads((p/'G2_REGION_MAP_v02.json').read_text(encoding='utf-8'))['regions']
ids=[v['embedding_index'] for v in regions['ACh_pocket']]
contexts=('candidate_probe','probe_only','candidate_no_probe','apo')
print('CONTEXT REP W2_MEAN_A W4_MEAN_A DELTA_A')
for c in contexts:
 for r in (2,3):
  frames=[]
  for w in (2,4):
   f=p/('replica_%02d'%r)/('window_%03d'%w)/'raw'/('%s_w%03d.npz'%(c,w))
   with np.load(f,allow_pickle=False) as z:
    ca=z['coords'][:,z['atom_names']=='CA',:].astype(np.float64)
    assert ca.shape==(100,270,3)
    frames.append(ca)
  x=np.concatenate(frames)
  ref=x[0]
  y=x-x.mean(axis=1,keepdims=True)
  u,_,vt=np.linalg.svd(np.einsum('fai,aj->fij',y,ref-ref.mean(axis=0)))
  u[:,:,-1]*=np.linalg.det(u@vt)[:,None]
  aligned=y@(u@vt)+ref.mean(axis=0)
  rms=np.sqrt(np.mean(np.sum((aligned[:,ids]-ref[ids])**2,axis=2),axis=1))
  assert len(rms)==200 and np.isfinite(rms).all() and abs(rms[0])<1e-6
  a=float(rms[:100].mean())
  b=float(rms[100:].mean())
  print(c,'R%d'%r,round(a,6),round(b,6),round(b-a,6))
print('QC PASS')
