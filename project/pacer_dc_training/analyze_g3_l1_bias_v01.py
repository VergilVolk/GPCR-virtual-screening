import json
from pathlib import Path
import numpy as np
p=Path('project/results/pacer_dc_four_context_v01/compound110')
regions=json.loads((p/'G2_REGION_MAP_v02.json').read_text(encoding='utf-8'))['regions']
ids=[v['embedding_index'] for v in regions['ACh_pocket']]
cs={'CA':'candidate_probe','A':'probe_only','C':'candidate_no_probe','0':'apo'}
norm=np.linalg.norm
def cos(x,y):
 return float(np.dot(x,y)/(norm(x)*norm(y)))
print('WINDOW CA_C_COS A_APO_COS BIAS_COS CANCEL_RATIO DINT_COS DINT_BIAS_NORM R2_NORM R3_NORM REL_DISTANCE')
for w in range(5):
 z={}
 for r in (2,3):
  for k,c in cs.items():
   f=p/f'five_layer_diagnostic_R{r}_W{w}_v01'/f'L1_{c}_R{r}_W{w}.npy'
   assert f.is_file(),f'Missing: {f}'
   arr=np.load(f,mmap_mode='r')
   assert arr.shape==(100,270,21),f'Unexpected shape: {f}'
   z[r,k]=np.asarray(arr[:,ids,:],dtype=np.float64).mean(axis=(0,1))
 d={k:z[2,k]-z[3,k] for k in cs}
 P=d['CA']-d['A']
 A=d['C']-d['0']
 I=P-A
 t={r:z[r,'CA']-z[r,'A']-z[r,'C']+z[r,'0'] for r in (2,3)}
 assert np.allclose(I,t[2]-t[3],rtol=1e-10,atol=1e-10)
 values=[cos(d['CA'],d['C']),cos(d['A'],d['0']),cos(P,A),float(norm(I)/(norm(P)+norm(A))),cos(t[2],t[3]),float(norm(I))]
 assert np.all(np.isfinite(values))
 print(f'W{w}',*[round(v,6) for v in values],*[round(float(v),6) for v in (norm(t[2]),norm(t[3]),2*norm(I)/(norm(t[2])+norm(t[3])))])
print('QC PASS')
