"""Linux build adapter: reuse validated scientific functions, isolate all outputs."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/legacy'))


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()


def table(name):
    with (ROOT/'manifests'/name).open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))


def verify_inputs():
    # Verify immutable scripts, config, manifests and all inputs on every workflow entry.
    for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
        digest,name=line.split('  ',1)
        path=(ROOT/name).resolve()
        assert ROOT in path.parents and path.is_file(), 'Invalid/missing package member: '+name
        assert sha(path)==digest, 'Package integrity mismatch: '+name
    manifest=json.loads((ROOT/'manifests/STAGE4_RERUN_INPUT_MANIFEST.json').read_text())
    assert manifest['protocol_id']=='pacer_stage4_rerun_v01'
    assert sha(ROOT/'manifests/frozen_ranking.csv')==manifest['frozen_ranking_sha256']
    for name,digest in manifest['assets'].items(): assert sha(ROOT/name)==digest, 'Input hash mismatch: '+name
    return manifest


def setup():
    import pacer_stage4_prospective as s4
    s4.OUT=ROOT
    s4.PROJECT=ROOT
    manifest=verify_inputs()
    systems=table('STAGE4_RERUN_SYSTEM_MANIFEST.csv')
    rows=[dict(candidate_id=r['candidate_id'],canonical_smiles=r['canonical_smiles'],source_state_smiles=r['state_smiles'],
               source_cluster='1',pose_file=r['raw_pose'],pose_file_sha256=r['raw_pose_sha256'],
               receptor_file=str(ROOT/'inputs/receptors/cluster_01_receptor.pdb')) for r in manifest['candidates']]
    s4.frozen_rows=lambda:rows
    s4.systems_for=lambda row:{r['system_id']: ([row['candidate_id']] if r['has_candidate']=='True' else [])+
                          (['ACH'] if r['has_ACH']=='True' else []) for r in systems if r['candidate_id']==row['candidate_id']}
    return s4,rows


def main():
    p=argparse.ArgumentParser(); p.add_argument('--all',action='store_true',required=True); p.parse_args()
    import numpy as np
    from rdkit import Chem
    from openff.units import unit as off_unit
    import build_pacer_stage4_prospective_membrane as legacy
    import build_pacer_dc_openmm_reference_systems as qc
    s4,rows=setup()
    orientation=json.loads((ROOT/'manifests/orientation_pacer_stage4_rerun_v01.json').read_text())
    def orient(row,args,history,qc):
        cid=row['candidate_id']
        candidate=qc.openff_molecule(Chem.SDMolSupplier(str(ROOT/f'inputs/ligands/{cid}_OPM.sdf'),removeHs=False)[0],cid)
        probe=qc.openff_molecule(Chem.SDMolSupplier(str(ROOT/'inputs/probe/ACH_cluster_01_OPM.sdf'),removeHs=False)[0],'ACH')
        for m in [candidate,probe]: m.generate_unique_atom_names()
        audit=dict(orientation,mobile_sha256=sha(ROOT/'inputs/receptors/cluster_01_receptor.pdb'),
                   candidate_pose_sha256=row['pose_file_sha256'])
        # A fresh copy for each candidate is needed because the legacy builder normalizes alongside it.
        d=ROOT/'orientation'/cid; d.mkdir(parents=True,exist_ok=True)
        receptor=d/'oriented_receptor.pdb'
        receptor.write_bytes((ROOT/'inputs/receptors/oriented_cluster_01.pdb').read_bytes())
        return receptor,candidate,probe,audit
    legacy.orient=orient
    legacy.assets=lambda args:{'OPM':ROOT/'inputs/receptors/7TRS_OPM.pdb','ACH':ROOT/'inputs/probe/ACH_cluster_01_OPM.sdf',
                               'receptor':ROOT/'inputs/receptors/cluster_01_receptor.pdb'}
    # The function writes per-candidate base before four systems, so each quartet shares one exact box.
    # Contexts have unique candidate-prefixed IDs even for A/P; no simulation results are aliased.
    for row in rows:
        names=list(s4.systems_for(row))
        existing=[ROOT/'systems'/n for n in names if (ROOT/'systems'/n).exists()]
        if existing:
            assert len(existing)==4, 'Partial candidate build exists; retain evidence and use a new deployment to rebuild'
            for d in existing:
                audit=json.loads((d/'build_audit.json').read_text()); assert audit['status']=='COMPLETE'
                assert audit['alignment']['candidate_pose_sha256']==row['pose_file_sha256']
                for name,digest in audit['outputs'].items(): assert sha(d/name)==digest
            continue
        legacy.build(row,argparse.Namespace(prepare_only=False))
        shared=ROOT/'shared_bases/cluster1'
        saved=ROOT/'shared_bases'/row['candidate_id']
        assert not saved.exists()
        shared.rename(saved)
        for name in names:
            p=ROOT/'systems'/name/'build_audit.json'
            audit=json.loads(p.read_text())
            assert sha(saved/'shared_OPM_POPC_water.pdb')==audit['shared_base_sha256']
            audit['shared_base_path']=str((saved/'shared_OPM_POPC_water.pdb').relative_to(ROOT))
            audit['rerun_config_sha256']=sha(ROOT/'config/pacer_stage4_rerun_v01.json')
            p.write_text(json.dumps(audit,indent=2)+'\n')
    # Full parameterization, residue/disulfide and context audit supplement.
    from openmm import app,unit,XmlSerializer,NonbondedForce
    from scipy.spatial import cKDTree
    for row in table('STAGE4_RERUN_SYSTEM_MANIFEST.csv'):
        d=ROOT/row['system_path']; audit=json.loads((d/'build_audit.json').read_text())
        assert audit['rerun_config_sha256']==sha(ROOT/'config/pacer_stage4_rerun_v01.json')
        assert sha(ROOT/audit['shared_base_path'])==audit['shared_base_sha256']
        pdb=app.PDBFile(str(d/'input.pdb')); system=XmlSerializer.deserialize((d/'system.xml').read_text())
        assert pdb.topology.getNumAtoms()==system.getNumParticles()
        params=json.loads((d/'parameters.json').read_text())
        assert len(params)==int(row['has_candidate']=='True')+int(row['has_ACH']=='True')
        for name,param in params.items(): assert abs(sum(param['partial_charges'])-param['formal_charge'])<1e-5
        atoms=list(pdb.topology.atoms()); xyz=np.asarray(pdb.positions.value_in_unit(unit.angstrom))
        bonded={frozenset([a.index,b.index]) for a,b in pdb.topology.bonds()}
        heavy=[a.index for a in atoms if a.element and a.element.symbol!='H']
        clashes=[]
        for x,y in cKDTree(xyz[heavy]).query_pairs(1.0):
            a,b=atoms[heavy[x]],atoms[heavy[y]]
            if a.residue!=b.residue and frozenset([a.index,b.index]) not in bonded:
                clashes.append([a.index,b.index,a.residue.name,b.residue.name])
        assert not clashes, 'Severe inter-residue heavy atom overlap: '+str(clashes[:20])
        expected_cyx={l[22:26].strip() for l in (ROOT/'inputs/receptors/cluster_01_receptor.pdb').read_text().splitlines()
                      if l.startswith('ATOM') and l[17:20]=='CYX'}
        # PDB output may renumber residues. Compare to the inherited preserved pre-output IDs via the input baseline sequence.
        from openmm import app as omm_app
        baseline=omm_app.PDBFile(str(ROOT/'inputs/receptors/oriented_cluster_01.pdb'))
        original_res=list(baseline.topology.residues())
        topology_res=list(pdb.topology.residues())
        assert [r.name for r in topology_res[:len(original_res)]]==[r.name for r in original_res]
        expected_indices={i for i,r in enumerate(original_res) if r.id in expected_cyx}
        paired_indices={atom.residue.index for a,b in pdb.topology.bonds() if a.name==b.name=='SG' for atom in (a,b)}
        assert expected_indices.issubset(paired_indices), 'Frozen CYX disulfide connectivity missing'
        assert audit['probe_atoms']==(params['ACH']['n_atoms'] if row['has_ACH']=='True' else 0)
        assert not any(r.name=='ACH' for r in original_res), 'Native ACh retained in receptor baseline'
        audit['rerun_QC']=dict(topology_particles_match=True,parameter_charge_sums=True,
            severe_inter_residue_heavy_clashes=clashes,threshold_A=1.0,
            disulfides=[[a.residue.id,b.residue.id] for a,b in pdb.topology.bonds() if a.name==b.name=='SG'])
        nonbonded=[f for f in system.getForces() if isinstance(f,NonbondedForce)]
        assert len(nonbonded)==1
        audit['rerun_QC']['system_net_charge_e']=sum(nonbonded[0].getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
        audit['rerun_QC']['ions']={name:sum(r.name.upper()==name for r in pdb.topology.residues()) for name in ['NA','CL']}
        audit['rerun_QC']['ion_policy']='Inherited baseline neutralization before ligand insertion; no extra counterions silently added'
        (d/'build_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    receipt=dict(status='PASS',systems=12,config_sha256=sha(ROOT/'config/pacer_stage4_rerun_v01.json'),
                 builds={r['system_id']:sha(ROOT/r['system_path']/'build_audit.json') for r in table('STAGE4_RERUN_SYSTEM_MANIFEST.csv')})
    (ROOT/'logs/build_pacer_stage4_rerun_v01.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__': main()
