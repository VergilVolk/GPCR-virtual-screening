"""Independent descriptive analysis of frozen regional vectors. No project imports.
Run from any directory with Python + numpy + matplotlib. Writes only beside this file.
No MD/encoder/RFF/graph evaluation on trajectory data, fitting, thresholds or labels.
"""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['MPLCONFIGDIR']=str(Path(__file__).resolve().parent/'_mplconfig')
import json, csv, hashlib, sys, platform, subprocess, itertools
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[2]
ROOT = REPO/'project/results/pacer_stage4_pbcfix01_fkg_v02_v01'
SRC = ROOT/'phase2b_graph_region'
CS = ['PACERGEN02432','PACERGEN00123','PACERGEN00462']
BS = ['STATE_MOTION','SIGNED_DRIFT']
KS = ['Delta_AGO','Delta_PAM','Delta_INT']
RS = ['pam_contact_consensus','pam_contact_union','cooperativity_mutagenesis',
      'compound110_extension','orthosteric_activation_core','orthosteric_contact_union',
      'intracellular_microswitches','stable_core_control','distal_control']
LABELS = ['PAM consensus','PAM union','Cooperativity','110 extension','Orth. activation',
          'Orth. union','Intracell. switches','Stable control','Distal control']
FORM = {'Delta_AGO':{'C':1.,'A':-1.},'Delta_PAM':{'CP':1.,'P':-1.},
        'Delta_INT':{'CP':1.,'C':-1.,'P':-1.,'A':1.}}
PAIRS = [(1,2),(1,3),(2,3)]
manifest = {}
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def record(p, role, expected=None):
    p=Path(p).resolve(); h=sha(p)
    if expected is not None: assert h==expected, ('SHA256 mismatch',str(p),h,expected)
    manifest[str(p)]={'path':str(p),'bytes':p.stat().st_size,'sha256':h,'role':role,
                      'expected_sha256':expected,'verified':expected is not None}
    return p
def load(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def art(a,role):
    p=record(a['path'],role,a['sha256']); assert p.stat().st_size==a['bytes']; return p
def dump(name,x): (OUT/name).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def csvout(name,rows):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def cos(a,b):
    den=np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.clip(np.dot(a,b)/den,-1,1)) if den>1e-12 else np.nan
def close(a,b): assert np.allclose(a,b,rtol=1e-10,atol=1e-12,equal_nan=True), (a,b)

receipt=load(record(SRC/'PHASE2B_FREEZE_RECEIPT_v01.json','final Phase2b receipt'))
dp=art(receipt['result'],'frozen scientific result'); d=load(dp)
assert d['schema']=='pacer.stage4.pbcfix01.frozen.graph_region.v1'
assert d['status']==receipt['status']=='STAGE4_FULL_FROZEN_EVALUATION_COMPLETE'
assert set(d['results'])==set(CS) and set(d['regions'])==set(RS)
assert not any(d['forbidden_operations_performed'].values())
assert d['temporal_contract']['blocks_per_trajectory']==10
for a in receipt['implementation']: art(a,'frozen implementation')
p2=load(art(d['phase2a_receipt'],'Phase2a receipt'))
if 'phase1_receipt' in p2: art(p2['phase1_receipt'],'Phase1 receipt')
anchor=d['historical_anchor']; fp=art(anchor['freeze'],'historical freeze'); freeze=load(fp)
FROOT=fp.parents[4]
for a in freeze['definitions'].values():
    a=dict(a); a['path']=str(FROOT/a['path']); art(a,'frozen definition')
for a in anchor.get('consumed_assets',[]): art(a,'historical consumed numerical asset')
graph=load(art(d['graph'],'frozen graph'))
for r in RS: assert sorted(d['regions'][r])==sorted(x['embedding_index'] for x in graph['regions'][r])
inputp=record(ROOT/'input_audit/INPUT_AUTHENTICATION_v01.json','full input authentication')
inp=load(inputp); assert inp['status']=='36_OF_36_ACCEPTED' and len(inp['jobs'])==36
assert len({(j['candidate'],j['context'],j['replica']) for j in inp['jobs']})==36
assert all(len(j['segments'])==100 for j in inp['jobs'])
for p in [ROOT/'local_validation/LOCAL_VALIDATION_v01.json',ROOT/'deployment_acceptance/DEPLOYMENT_ACCEPTANCE_v01.json',
          REPO/'docs/PACER_STAGE4_PBCFIX01_FKG_V02_HANDOFF_v01.md',REPO/'docs/PBCFIX01_DOWNLOAD_AND_FKG_ACCEPTANCE_v01.md']:
    record(p,'historical protocol or earlier acceptance snapshot')
arrays={}; errors=[]; checked=0; max_formula=0.
assert len(d['vector_artifacts'])==18
for item in d['vector_artifacts']:
    c,b,rep=item['candidate'],item['branch'],item['replica']
    assert c in CS and b in BS and rep in [1,2,3]
    assert item['contexts']==d['results'][c]['contexts'] and item['cluster']==d['results'][c]['cluster']
    assert item['seed']==d['prospective_replica_policy'][f'R{rep}']['seed']
    p=art(item['artifact'],'frozen NPZ')
    assert p== (SRC/'full'/c/b/f'replica_{rep:02d}'/'BLOCK_REGIONAL_VECTORS.npz').resolve()
    for a in item['upstream']: art(a,'authenticated upstream manifest')
    with np.load(p,allow_pickle=False) as z: a={k:z[k].astype(np.float64) for k in z.files}
    expected={f'context__{x}__{r}' for x in ['A','P','C','CP'] for r in RS}
    branch=d['results'][c]['branches'][b]
    expected|={f'contrast__{k}__{r}' for k in branch for r in RS}
    assert set(a)==expected and all(v.shape==(10,512) and np.isfinite(v).all() for v in a.values())
    assert (c,b,rep) not in arrays; arrays[c,b,rep]=a
    for k,itemk in branch.items():
        if k in FORM: assert itemk['formula']==FORM[k]
        for r in RS:
            v=a[f'contrast__{k}__{r}']; m=v.mean(0); norms=np.linalg.norm(v,axis=1)
            s=itemk['regions'][r]['replicas'][f'R{rep}']
            close(m,s['signed_direction_vector']); close(np.linalg.norm(m),s['mean_vector_norm'])
            for key,val in [('block_magnitude_mean',norms.mean()),('block_magnitude_q10',np.quantile(norms,.1)),
                            ('block_magnitude_median',np.median(norms)),('block_magnitude_q90',np.quantile(norms,.9))]: close(val,s[key])
            assert s['n_correlated_time_blocks']==10
            reconstructed=sum(sign*a[f'context__{x}__{r}'] for x,sign in itemk['formula'].items())
            err=float(np.max(np.abs(v-reconstructed))); max_formula=max(max_formula,err)
            # Contrast and contexts were separately diffused with float32 rounding.
            assert np.allclose(v,reconstructed,rtol=2e-5,atol=2e-7),err
            checked+=1
assert set(arrays)==set(itertools.product(CS,BS,[1,2,3]))
assert set(SRC.rglob('BLOCK_REGIONAL_VECTORS.npz'))=={Path(x['artifact']['path']) for x in d['vector_artifacts']}
def vec(c,b,k,r,rep): return arrays[c,b,rep][f'contrast__{k}__{r}']
rows=[]; blocks=[]; relation=[]; spatial=[]
for c,b,k,r in itertools.product(CS,BS,KS,RS):
    row={'candidate':c,'branch':b,'contrast':k,'region':r,'n_region_residues':len(d['regions'][r])}
    means=[vec(c,b,k,r,i).mean(0) for i in [1,2,3]]
    norms=[np.linalg.norm(v) for v in means]
    for rep,v,m in zip([1,2,3],[vec(c,b,k,r,i) for i in [1,2,3]],means):
        mags=np.linalg.norm(v,axis=1); pre=f'R{rep}_'
        row.update({pre+'mean_vector_norm':np.linalg.norm(m),pre+'block_norm_mean':mags.mean(),
                    pre+'block_norm_min':mags.min(),pre+'block_norm_q10':np.quantile(mags,.1),
                    pre+'block_norm_median':np.median(mags),pre+'block_norm_q90':np.quantile(mags,.9),
                    pre+'block_norm_max':mags.max(),pre+'temporal_coherence':np.linalg.norm(m)/mags.mean(),
                    pre+'early_late_cosine':cos(v[:5].mean(0),v[5:].mean(0))})
        for ctrl in RS[-2:]: row[pre+'ratio_to_'+ctrl]=np.linalg.norm(m)/np.linalg.norm(vec(c,b,k,ctrl,rep).mean(0))
        for i,n in enumerate(mags):
            blocks.append({'candidate':c,'branch':b,'contrast':k,'region':r,'replica':rep,'block':i+1,
                           'block_first_ps':10+i*1000,'block_last_ps':960+i*1000,'norm':n,
                           'npz_key':f'contrast__{k}__{r}'})
    pcs=[]
    for i,j in PAIRS:
        value=cos(means[i-1],means[j-1]); pcs.append(value); row[f'R{i}_R{j}_cosine']=value
        close(value,d['results'][c]['branches'][b][k]['regions'][r]['all_pair_direction_cosines'][f'R{i}_R{j}'])
    row.update({'pair_cosine_mean':np.mean(pcs),'pair_cosine_min':min(pcs),'pair_cosine_max':max(pcs),
                'replica_norm_mean':np.mean(norms),'replica_norm_sd':np.std(norms,ddof=1),
                'replica_norm_cv':np.std(norms,ddof=1)/np.mean(norms),'replica_norm_max_min':max(norms)/min(norms)})
    rows.append(row)
for c,b,r,rep in itertools.product(CS,BS,RS,[1,2,3]):
    ago,pam,inter=[vec(c,b,k,r,rep).mean(0) for k in KS]
    relation.append({'candidate':c,'branch':b,'region':r,'replica':rep,'AGO_PAM_cosine':cos(ago,pam),
                     'INT_norm_over_AGO_plus_PAM_norm':np.linalg.norm(inter)/(np.linalg.norm(ago)+np.linalg.norm(pam)),
                     'INT_equals_PAM_minus_AGO_max_abs':np.max(np.abs(inter-(pam-ago)))})
for c,b,k,rep in itertools.product(CS,BS,KS,[1,2,3]):
    for r,s in itertools.combinations(RS,2):
        spatial.append({'candidate':c,'branch':b,'contrast':k,'replica':rep,'region_a':r,'region_b':s,
                        'cosine':cos(vec(c,b,k,r,rep).mean(0),vec(c,b,k,s,rep).mean(0))})
csvout('STAGE4_PBCFIX01_REGION_COMPARISON_v01.csv',rows)
csvout('BLOCK_MAGNITUDES_v01.csv',blocks); csvout('CONTRAST_RELATIONS_v01.csv',relation)
csvout('CROSS_REGION_DIRECTION_v01.csv',spatial)
lookup={(r['candidate'],r['branch'],r['contrast'],r['region']):r for r in rows}
crows=[]
for c,b,k in itertools.product(CS,BS,KS):
    rr=[lookup[c,b,k,r] for r in RS]
    cr={'candidate':c,'branch':b,'contrast':k,'retained_replicas':'R1;R2;R3','n_regions':9,
        'all_pair_cosine_min':min(r['pair_cosine_min'] for r in rr),
        'regions_all_three_pairs_positive':sum(r['pair_cosine_min']>0 for r in rr),
        'regions_with_negative_pair':sum(r['pair_cosine_min']<0 for r in rr),
        'most_directionally_consistent_region':max(rr,key=lambda r:r['pair_cosine_min'])['region'],
        'least_directionally_consistent_region':min(rr,key=lambda r:r['pair_cosine_min'])['region']}
    for r in RS:
        item=lookup[c,b,k,r]
        for field in ['pair_cosine_min','pair_cosine_mean','replica_norm_mean','replica_norm_cv']:
            cr[r+'__'+field]=item[field]
    crows.append(cr)
csvout('STAGE4_PBCFIX01_CANDIDATE_COMPARISON_v01.csv',crows)
print('Frozen identity and descriptive replay passed; CSV complete',flush=True)

# Structural-only diffusion influence: never apply this operator to MD/model vectors.
# Exact frozen edge construction, to quantify shared regional information without refitting.
n=len(graph['nodes']); edges={}
for e in graph['edges']:
    i,j=int(e['source']),int(e['target'])
    if i==j: continue
    pair=tuple(sorted((i,j))); kind='backbone' if e['edge_type']=='backbone' else 'contact'
    w=1. if kind=='backbone' else float(e['support_fraction']); old=edges.get(pair)
    if old is None or (old[0]!='backbone' and kind=='backbone'): edges[pair]=(kind,w)
    elif old[0]!='backbone': edges[pair]=(kind,max(old[1],w))
adj=np.zeros((n,n))
for (i,j),(_,w) in edges.items(): adj[i,j]=adj[j,i]=w
transition=adj/adj.sum(1,keepdims=True); influence=np.eye(n)
for _ in range(d['graph_steps']):
    influence=(1-d['graph_alpha'])*np.eye(n)+d['graph_alpha']*np.einsum('ij,jk->ik',transition,influence,optimize=False)
weights={r:influence[d['regions'][r]].mean(0) for r in RS}
overlap=[]
for r,s in itertools.product(RS,RS):
    a,z=set(d['regions'][r]),set(d['regions'][s]); common=a&z
    overlap.append({'region_a':r,'region_b':s,'n_a':len(a),'n_b':len(z),'intersection_n':len(common),
                    'intersection_indices':';'.join(map(str,sorted(common))),'jaccard':len(common)/len(a|z),
                    'diffusion_weight_cosine':cos(weights[r],weights[s]),
                    'diffusion_shared_probability_mass':np.minimum(weights[r],weights[s]).sum()})
csvout('REGION_OVERLAP_DIFFUSION_v01.csv',overlap)
csvout('REGION_MEMBERSHIP_v01.csv',[{'region':r,**m} for r in RS for m in graph['regions'][r]])
print('Structural region dependence complete',flush=True)

plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig_sources=[]
def save(fig,name,source,description):
    fig.savefig(OUT/(name+'.png'),dpi=220,bbox_inches='tight'); fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig); fig_sources.append({'figure':name,'source_csv':source,'description':description,
                                      'frozen_result_sha256':sha(dp)})
for b in BS:
    fig,axes=plt.subplots(3,3,figsize=(15,12),layout='constrained',sharex=True,sharey=True)
    for ci,c in enumerate(CS):
        for ki,k in enumerate(KS):
            ax=axes[ci,ki]; m=np.array([[lookup[c,b,k,r][f'R{i}_R{j}_cosine'] for i,j in PAIRS] for r in RS])
            im=ax.imshow(m,vmin=-1,vmax=1,cmap='RdBu',aspect='auto')
            for y,x in itertools.product(range(9),range(3)): ax.text(x,y,f'{m[y,x]:.3f}',ha='center',va='center',fontsize=8,color='white' if abs(m[y,x])>.6 else 'black')
            ax.set_title(f'{c} | {k}'); ax.set_xticks(range(3),['R1-R2','R1-R3','R2-R3']); ax.tick_params(labelbottom=True)
            ax.set_yticks(range(9),LABELS); ax.tick_params(labelleft=ki==0)
    fig.colorbar(im,ax=axes,label='Cosine of replica mean vectors',shrink=.6)
    fig.suptitle(b+' | All replica pairs; all frozen regions (descriptive)')
    save(fig,b+'_replica_direction','STAGE4_PBCFIX01_REGION_COMPARISON_v01.csv','Annotated R1-R2/R1-R3/R2-R3 cosines; fixed -1 to 1 scale')
    fig,axes=plt.subplots(3,3,figsize=(16,12),layout='constrained',sharey=True)
    colors=['#0072B2','#D55E00','#009E73']
    for ci,c in enumerate(CS):
        for ki,k in enumerate(KS):
            ax=axes[ci,ki]
            for rep,color in zip([1,2,3],colors):
                vals=[vec(c,b,k,r,rep) for r in RS]; yy=np.arange(9)+(rep-2)*.22
                qs=np.array([np.quantile(np.linalg.norm(v,axis=1),[.1,.5,.9]) for v in vals])
                ax.hlines(yy,qs[:,0],qs[:,2],color=color,lw=1.5)
                ax.scatter(qs[:,1],yy,color=color,s=14,label=f'R{rep} block median [q10,q90]')
                ax.scatter([np.linalg.norm(v.mean(0)) for v in vals],yy,marker='x',color=color,s=22,label=f'R{rep} norm(mean vector)')
            ax.set_yticks(range(9),LABELS); ax.invert_yaxis(); ax.set_title(c+' | '+k); ax.set_xlabel('Frozen RFF vector L2 norm'); ax.grid(axis='x',alpha=.2)
    # A common branch-only amplitude scale across all panels.
    lim=max(np.quantile(np.linalg.norm(vec(c,b,k,r,i),axis=1),.9) for c,k,r,i in itertools.product(CS,KS,RS,[1,2,3]))*1.06
    for ax in axes.flat: ax.set_xlim(0,lim)
    axes[0,0].legend(fontsize=6,loc='lower right'); fig.suptitle(b+' | Ten correlated blocks per replica; intervals are descriptive, not CIs')
    save(fig,b+'_magnitude_distribution','STAGE4_PBCFIX01_REGION_COMPARISON_v01.csv','Per-replica block q10/median/q90 and norm of temporal mean; branch-specific scales')
    fig,axes=plt.subplots(3,3,figsize=(15,11),layout='constrained',sharey=True)
    for ci,c in enumerate(CS):
        for ki,k in enumerate(KS):
            ax=axes[ci,ki]
            for rep,color in zip([1,2,3],colors):
                ys=[lookup[c,b,k,r][f'R{rep}_ratio_to_stable_core_control'] for r in RS]
                ax.plot(ys,np.arange(9)+(rep-2)*.15,'o',ms=4,color=color,label=f'R{rep}')
            ax.axvline(1,color='gray',lw=1); ax.set_xscale('log'); ax.set_yticks(range(9),LABELS); ax.invert_yaxis()
            ax.set_title(c+' | '+k); ax.set_xlabel('norm(mean vector) / stable control (log)'); ax.grid(axis='x',alpha=.2)
    ratios=[lookup[c,b,k,r][f'R{i}_ratio_to_stable_core_control'] for c,k,r,i in itertools.product(CS,KS,RS,[1,2,3])]
    for ax in axes.flat: ax.set_xlim(min(ratios)*.9,max(ratios)*1.1)
    axes[0,0].legend(); fig.suptitle(b+' | Regional specificity diagnostic; ratio 1 is equality, not an activity threshold')
    save(fig,b+'_region_specificity','STAGE4_PBCFIX01_REGION_COMPARISON_v01.csv','Per-replica region/stable-control norm ratios; distal control retained')
fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
for ax,field,title in zip(axes,['jaccard','diffusion_shared_probability_mass'],['Direct residue Jaccard','Shared structural diffusion mass']):
    m=np.array([[next(x[field] for x in overlap if x['region_a']==r and x['region_b']==s) for s in RS] for r in RS])
    im=ax.imshow(m,vmin=0,vmax=1,cmap='viridis'); ax.set_xticks(range(9),LABELS,rotation=70,ha='right'); ax.set_yticks(range(9),LABELS); ax.set_title(title)
    for y,x in itertools.product(range(9),range(9)): ax.text(x,y,f'{m[y,x]:.2f}',ha='center',va='center',fontsize=6,color='white' if m[y,x]<.6 else 'black')
fig.colorbar(im,ax=axes,shrink=.6); fig.suptitle('Frozen region dependence | Structural operator only; no model rerun')
save(fig,'region_dependence','REGION_OVERLAP_DIFFUSION_v01.csv','Jaccard and sum(min regional influence weights); not independent evidence')
csvout('FIGURE_PROVENANCE_v01.csv',fig_sources)
csvout('INPUT_SHA256_v01.csv',list(manifest.values()))
audit={'status':'PASS','json_npz_identity':'18/18 SHA256 and complete candidate/branch/replica grid',
       'n_arrays':sum(len(a) for a in arrays.values()),'n_json_replica_region_summaries_replayed':checked,
       'n_main_contrast_pair_cosines_replayed':len(rows)*3,'max_context_contrast_formula_absolute_error':max_formula,
       'formula_tolerance':{'rtol':2e-5,'atol':2e-7},'json_summary_tolerance':{'rtol':1e-10,'atol':1e-12},
       'input_authentication_status':inp['status'],'authenticated_trajectories':36,'authenticated_segments':3600,
       'scope':'Frozen artifact identity and descriptive replay; MD files and model pipeline not re-executed',
       'input_hashes_unchanged_after_analysis':all(sha(p)==v['sha256'] for p,v in manifest.items()),
       'retained_replicas':[1,2,3],'excluded_replicas':[],'region_rows':len(rows),'block_rows':len(blocks),
       'candidate_rows':len(crows),'figures_png_pdf':len(fig_sources)}
assert audit['input_hashes_unchanged_after_analysis']; dump('ANALYSIS_VALIDATION_v01.json',audit)
version={'analysis_version':'v01','utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(__file__),
         'report_renderer_sha256':sha(OUT/'write_report_v01.py'),
         'python':sys.version,'numpy':np.__version__,'matplotlib':matplotlib.__version__,
         'platform':platform.platform(),'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
         'command':sys.executable+' '+str(Path(__file__).relative_to(REPO)),
         'method':'L2 norm, cosine, empirical correlated-block quantiles; no inferential tests or fitted scores'}
dump('ANALYSIS_CODE_VERSION_v01.json',version)
print(json.dumps(audit,indent=2)); print('OUTPUT',OUT)
