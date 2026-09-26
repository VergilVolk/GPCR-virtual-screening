import json,itertools
import numpy as np
from pathlib import Path

p=Path('project/results/pacer_dc_four_context_v01/compound110')
regions=json.loads((p/'G2_REGION_MAP_v02.json').read_text(encoding='utf-8'))['regions']
ids=[v['embedding_index'] for v in regions['ACh_pocket']]
ii,jj=np.triu_indices(len(ids),1)
assert len(ids)==12 and len(ii)==66
cs=('candidate_probe','probe_only','candidate_no_probe','apo')
d={}
for r in (2,3):
 for w in range(5):
  for c in cs:
   f=p/('replica_%02d'%r)/('window_%03d'%w)/'raw'/('%s_w%03d.npz'%(c,w))
   with np.load(f,allow_pickle=False) as z:
    ca=z['coords'][:,z['atom_names']=='CA',:][:,ids,:].astype(np.float64)
   assert ca.shape==(100,12,3)
   d[r,w,c]=np.linalg.norm(ca[:,ii]-ca[:,jj],axis=2).mean(axis=0)
g={(r,w):d[r,w,cs[0]]-d[r,w,cs[1]]-d[r,w,cs[2]]+d[r,w,cs[3]] for r in (2,3) for w in range(5)}
n=np.linalg.norm
pairs=list(itertools.combinations(range(5),2))
result={}
print('FROM TO GEOMETRY_COS R2_SHIFT_RMS_A R3_SHIFT_RMS_A')
for i,j in pairs:
 a=g[2,j]-g[2,i]
 b=g[3,j]-g[3,i]
 assert n(a)>0 and n(b)>0
 cosine=float(np.dot(a,b)/(n(a)*n(b)))
 result[i,j]=cosine
 print('W%d'%i,'W%d'%j,round(cosine,6),round(float(n(a)/np.sqrt(66)),6),round(float(n(b)/np.sqrt(66)),6))
assert abs(result[2,4]-0.295906)<0.00001
assert all(np.isfinite(v) for v in result.values())
print('QC PASS')
