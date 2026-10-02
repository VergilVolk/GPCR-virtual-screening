"""Contract tests: identity corruption, rigid transform, nonfinite smoke inputs."""
import copy
from pathlib import Path
import sys
import unittest
import tempfile
import json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pacer_stage4_prospective as s4
from run_pacer_stage4_membrane_smoke import finite_metrics

class ContractTests(unittest.TestCase):
    def test_dry_run_writes_manifest_without_creating_smoke_receipts(self):
        previous=s4.OUT
        with tempfile.TemporaryDirectory() as directory:
            try:
                s4.OUT=Path(directory)
                (s4.OUT/'master_manifest.json').write_text('{}')
                rows=[dict(candidate_id=cid,source_cluster=str(cluster),pose_file_sha256=digest)
                      for cid,(cluster,digest) in s4.EXPECTED.items()]
                s4.update_manifest(rows)
                m=json.loads((s4.OUT/'master_manifest.json').read_text())
                self.assertEqual(m['expected_job_count'],36)
                self.assertFalse((s4.OUT/'smoke').exists())
                self.assertFalse(m['production_started'])
            finally: s4.OUT=previous
    def test_rigid_transform_preserves_internal_geometry(self):
        rng=np.random.default_rng(42); xyz=rng.normal(size=(274,3))
        R=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        q=xyz[:270]@R+np.array([3,2,-1])
        mobile=[dict(residue='ALA' if i<270 else 'GLY',resid=str(i),xyz=p.tolist()) for i,p in enumerate(xyz)]
        reference=[dict(residue='ALA',resid=str(i+100),xyz=p.tolist()) for i,p in enumerate(q)]
        audit=s4.alignment(mobile,reference)
        mapped=s4.transform_xyz(xyz,audit)
        self.assertLess(audit['rmsd_A'],1e-12)
        np.testing.assert_allclose(mapped[:270],q,atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(mapped[1:]-mapped[:-1],axis=1),
                                   np.linalg.norm(xyz[1:]-xyz[:-1],axis=1),atol=1e-12)
        reference[10]['residue']='TRP'
        with self.assertRaisesRegex(ValueError,'correspondence'): s4.alignment(mobile,reference)

    @unittest.skipUnless((s4.FROZEN/'md_shortlist_pose_manifest.csv').exists(),'Frozen deployment assets unavailable')
    def test_frozen_identity_and_microstate_charge(self):
        rows=s4.frozen_rows()
        charges={r['candidate_id']:s4.Chem.GetFormalCharge(s4.candidate_molecule(r)) for r in rows}
        self.assertEqual(charges,{'PACER0010':0,'PACER0073':0,'PACER0027':2})
        corrupt=copy.deepcopy(rows[0]); corrupt['source_state_smiles']=rows[1]['source_state_smiles']
        with self.assertRaisesRegex(ValueError,'graph mismatch'): s4.candidate_molecule(corrupt)
        self.assertEqual(sum(len(s4.systems_for(r)) for r in rows),12)

    def test_nonfinite_smoke_rejected(self):
        from openmm import System,unit
        class State:
            def getPotentialEnergy(self): return 0*unit.kilojoule_per_mole
        system=System(); system.addParticle(1)
        xyz=np.zeros((1,3))*unit.nanometer
        forces=np.array([[np.inf,0,0]])*unit.kilojoule_per_mole/unit.nanometer
        with self.assertRaisesRegex(ValueError,'Nonfinite'): finite_metrics(State(),system,xyz,forces)

    @unittest.skipUnless((s4.FROZEN/'md_shortlist_pose_manifest.csv').exists(),'Frozen deployment assets unavailable')
    def test_capped_cluster_topology_and_native_probe_exclusion(self):
        from openmm import app
        import io
        for cluster in [0,4,9]:
            source=s4.PROJECT/f'results/m4_gamd_ensemble/receptors/cluster_{cluster:02d}_receptor.pdb'
            original=source.read_text().splitlines(); normalized=s4.normalized_receptor_lines(original)
            self.assertFalse(any(x[17:20]=='ACH' for x in normalized if x.startswith('ATOM')))
            # Cap/chain repair must not move any retained protein atom.
            self.assertEqual([x[30:54] for x in normalized if x.startswith('ATOM')],
                             [x[30:54] for x in original if x.startswith('ATOM') and x[17:20]!='ACH'])
            pdb=app.PDBFile(io.StringIO('\n'.join(normalized)+'\n'))
            system=app.ForceField('amber14/protein.ff14SB.xml').createSystem(
                pdb.topology,nonbondedMethod=app.NoCutoff,constraints=app.HBonds)
            self.assertEqual(system.getNumParticles(),pdb.topology.getNumAtoms())

if __name__=='__main__': unittest.main()
