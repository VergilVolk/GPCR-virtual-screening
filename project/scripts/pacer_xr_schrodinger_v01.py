"""Schrodinger 2023-1 Python bridge: identity stamping and finite pose extraction."""
import json,sys,math
from collections import Counter
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import analyze

def require(condition,message):
    if not condition: raise ValueError(message)

def main():
    mode,metadata_path=sys.argv[1:]
    meta=json.loads(Path(metadata_path).read_text()); work=Path(metadata_path).parent
    if mode=='stamp':
        counts=Counter(); variants=[]
        with structure.StructureWriter(str(work/'docking_ready.maegz')) as writer:
            for st in structure.StructureReader(str(work/'prepared.maegz')):
                state=st.title
                require(state in meta['states'],'LigPrep returned unknown source state')
                counts[state]+=1; variant=state+'__v%02d'%counts[state]
                st.title=variant
                st.property['s_pacer_candidate_id']=meta['candidate_id']
                st.property['s_pacer_canonical_smiles']=meta['canonical_smiles']
                st.property['i_pacer_rerun_final_rank']=meta['rerun_final_rank']
                st.property['s_pacer_molscrub_state_id']=state
                st.property['s_pacer_prepared_variant_id']=variant
                writer.append(st)
                variants.append(dict(candidate_id=meta['candidate_id'],canonical_smiles=meta['canonical_smiles'],rerun_final_rank=meta['rerun_final_rank'],molscrub_state_id=state,prepared_variant_id=variant,prepared_title=variant,prepared_smiles=analyze.generate_smiles(st)))
        require(bool(variants) and set(counts)==set(meta['states']),'LigPrep silently dropped a Molscrub state')
        (work/'variants.json').write_text(json.dumps(variants,indent=2))
    elif mode=='parse':
        accepted={v['prepared_variant_id']:v for v in meta['variants']}; poses=[]; receptor_seen=False
        for entry,st in enumerate(structure.StructureReader(meta['poseviewer']),1):
            if st.property.get('b_glide_receptor'):
                require(not receptor_seen,'Duplicate receptor entry'); receptor_seen=True
                for axis,value in zip('xyz',meta['center']):
                    require(abs(st.property['r_glide_gridbox_'+axis+'cent']-value)<1e-5,'Returned grid center differs')
                continue
            variant=st.property.get('s_pacer_prepared_variant_id',st.property.get('s_smoke_ligand_state_id'))
            candidate=st.property.get('s_pacer_candidate_id',st.property.get('s_smoke_candidate_id'))
            parent=st.property.get('s_pacer_canonical_smiles',st.property.get('s_sd_parent\\_smiles'))
            require(candidate==meta['candidate_id'] and parent==meta['canonical_smiles'],'Pose parent identity mismatch')
            require(variant in accepted and st.title==variant,'Unrecognized prepared variant')
            score=float(st.property['r_i_glide_gscore'])
            require(math.isfinite(score) and all(math.isfinite(x) for a in st.atom for x in a.xyz),'Nonfinite pose or GlideScore')
            poses.append(dict(candidate_id=candidate,canonical_smiles=parent,receptor_id=meta['receptor_id'],grid_sha256=meta['grid_sha256'],molscrub_state_id=accepted[variant]['molscrub_state_id'],prepared_variant_id=variant,entry_index=entry,glide_pose_number=st.property.get('i_i_glide_posenum'),glide_gscore=score))
        require(receptor_seen and bool(poses),'Missing receptor or ligand poses')
        poses.sort(key=lambda p:(p['glide_gscore'],p['prepared_variant_id'],p['entry_index']))
        for rank,pose in enumerate(poses,1): pose['pose_rank']=rank
        (work/'poses.json').write_text(json.dumps(poses,indent=2))
    else: raise ValueError('Unknown bridge operation')
    print(mode+' PASS')

if __name__=='__main__': main()
