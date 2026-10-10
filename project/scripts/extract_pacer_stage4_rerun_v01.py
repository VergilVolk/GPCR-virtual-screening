"""Schrodinger-only extraction; no docking, coordinate generation or minimization."""
import hashlib
import json
import sys
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import analyze


def main():
    root = Path(sys.argv[1]).resolve()
    request_path=root/'manifests/extraction_requests.json'
    if request_path.exists():
        requests=json.loads(request_path.read_text())
    else:
        requests=json.loads((root/'manifests/STAGE4_RERUN_POSE_AUDIT.json').read_text())['candidate_poses']
        for r in requests: r['source_poseviewer']=str(root/r['source_poseviewer'])
    audits = []
    by_source = {}
    for r in requests:
        by_source.setdefault(r['source_poseviewer'], []).append(r)
    for source, rows in by_source.items():
        wanted = {r['entry_index']: r for r in rows}
        for index, st in enumerate(structure.StructureReader(source), 1):
            if st.property.get('b_glide_receptor'):
                st.write(str(root/'inputs/receptors/glide_cluster_01.pdb'))
                st.write(str(root/'inputs/receptors/glide_cluster_01.mae'))
            if index not in wanted:
                continue
            r = wanted.pop(index)
            cid = r['candidate_id']
            assert st.property['s_pacer_candidate_id'] == cid
            assert st.property['s_pacer_canonical_smiles'] == r['canonical_smiles']
            assert st.property['s_pacer_prepared_variant_id'] == r['selected_variant_id']
            assert abs(st.property['r_i_glide_gscore']-r['glide_gscore']) < 1e-8
            st.write(str(root/f'inputs/poses/{cid}_glide_cluster_01.mae'))
            st.write(str(root/f'inputs/ligands/{cid}_glide_cluster_01.sdf'))
            atoms = [dict(index=a.index, atomic_number=a.atomic_number, formal_charge=a.formal_charge,
                          name=a.pdbname, xyz=list(a.xyz)) for a in st.atom]
            bonds = [dict(a=b.atom1.index, b=b.atom2.index, order=b.order) for b in st.bond]
            audits.append(dict(r, title=st.title, state_smiles=analyze.generate_smiles(st),
                               formal_charge=sum(a.formal_charge for a in st.atom),
                               atoms=atoms, bonds=bonds,
                               properties={k: v for k, v in st.property.items()}))
        assert not wanted, 'Missing unique requested pose entries'
    assert len(audits) == 3
    (root/'manifests/extracted_pose_records.json').write_text(json.dumps(audits, indent=2)+'\n')
    print('Extracted three real Glide poses; all original coordinates retained')


if __name__ == '__main__':
    main()
