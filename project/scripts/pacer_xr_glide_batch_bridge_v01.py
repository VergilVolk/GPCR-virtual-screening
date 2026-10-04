"""Schrodinger bridge for audited caches, deterministic libraries and batch parsing."""
import json,sys,math
from pathlib import Path
from schrodinger import structure

def require(ok,message):
    if not ok: raise ValueError(message)

def identity(st,expected):
    variant=st.property.get('s_pacer_prepared_variant_id',st.property.get('s_smoke_ligand_state_id'))
    cid=st.property.get('s_pacer_candidate_id',st.property.get('s_smoke_candidate_id'))
    parent=st.property.get('s_pacer_canonical_smiles',st.property.get('s_sd_parent\\_smiles'))
    require(variant in expected,'Unknown variant '+str(variant))
    meta=expected[variant]
    require(cid==meta['candidate_id'] and parent==meta['canonical_smiles'] and st.title==variant,'Parent/state/variant identity mismatch')
    require(all(math.isfinite(x) for a in st.atom for x in a.xyz),'Nonfinite coordinates')
    return variant,meta

def parse(path,expected,receptor):
    groups={}; seen=False
    for entry,st in enumerate(structure.StructureReader(str(path)),1):
        if st.property.get('b_glide_receptor'):
            require(not seen,'Duplicate receptor'); seen=True
            for axis,c,inner,outer in zip('xyz',receptor['center'],receptor['inner'],receptor['outer']):
                require(abs(st.property['r_glide_gridbox_'+axis+'cent']-c)<1e-5,'Grid center mismatch')
                require(abs(st.property['r_glide_gridbox_lig'+axis+'range']-inner)<1e-5,'Grid inner mismatch')
                require(abs(st.property['r_glide_gridbox_'+axis+'range']-outer)<1e-5,'Grid outer mismatch')
            continue
        variant,meta=identity(st,expected); score=float(st.property['r_i_glide_gscore'])
        require(math.isfinite(score),'Nonfinite gscore')
        cid=meta['candidate_id']; group=groups.setdefault(cid,[])
        group.append(dict(candidate_id=cid,canonical_smiles=meta['canonical_smiles'],prepared_variant_id=variant,molscrub_state_id=meta['molscrub_state_id'],entry_index=entry,candidate_entry_index=len(group)+2,glide_pose_number=st.property.get('i_i_glide_posenum'),glide_gscore=score,receptor_id=receptor['receptor_id'],grid_sha256=receptor['grid_sha256']))
    require(seen and bool(groups),'No valid batch poses')
    for group in groups.values():
        group.sort(key=lambda p:(p['glide_gscore'],p['prepared_variant_id'],p['entry_index']))
        for rank,p in enumerate(group,1): p['pose_rank']=rank
    return groups

def main():
    mode,path=sys.argv[1:]; meta=json.loads(Path(path).read_text()); work=Path(path).parent
    if mode in ['audit','pack']:
        counts={}; expected={v['prepared_variant_id']:v for p in meta['prepared'] for v in p['variants']}
        require(len(expected)==sum(len(p['variants']) for p in meta['prepared']),'Duplicate prepared variant identity')
        writer=structure.StructureWriter(str(work/'library.maegz')) if mode=='pack' else None
        try:
            for parent in meta['prepared']:
                variants=[]
                for st in structure.StructureReader(parent['ligand']):
                    variant,v=identity(st,expected); require(v['candidate_id']==parent['candidate_id'],'Ligand file contains another parent')
                    variants.append(variant)
                    if writer:
                        st.property['s_pacer_candidate_id']=v['candidate_id']; st.property['s_pacer_canonical_smiles']=v['canonical_smiles']
                        st.property['i_pacer_rerun_final_rank']=v['rerun_final_rank']; st.property['s_pacer_molscrub_state_id']=v['molscrub_state_id']; st.property['s_pacer_prepared_variant_id']=variant
                        writer.append(st)
                require(set(variants)=={v['prepared_variant_id'] for v in parent['variants']} and len(variants)==len(parent['variants']),'Prepared variants missing/duplicated')
                counts[parent['candidate_id']]=len(variants)
        finally:
            if writer: writer.close()
        results={}
        if mode=='audit':
            for job in meta.get('completed',[]):
                expected_one={v['prepared_variant_id']:v for v in job['variants']}
                groups=parse(job['poseviewer'],expected_one,job['receptor'])
                require(set(groups)=={job['candidate_id']},'Old poseviewer identity mismatch')
                results[job['candidate_id']+'_'+job['receptor']['receptor_id']]=groups[job['candidate_id']]
        (work/(mode+'_result.json')).write_text(json.dumps(dict(counts=counts,completed_poses=results),indent=2))
    elif mode=='parse':
        expected={v['prepared_variant_id']:v for p in meta['prepared'] for v in p['variants']}
        groups=parse(meta['poseviewer'],expected,meta['receptor'])
        wanted={p['candidate_id'] for p in meta['prepared']}
        require(set(groups)==wanted,'Batch missing candidate poses: '+str(sorted(wanted-set(groups))))
        (work/'batch_poses.json').write_text(json.dumps(groups,indent=2))
    else: raise ValueError('Unknown mode')
    print(mode+' PASS')

if __name__=='__main__': main()
