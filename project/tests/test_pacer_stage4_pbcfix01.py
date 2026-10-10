import sys
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pacer_stage4_pbcfix01 as fix


class IndependentDistances(unittest.TestCase):
    def test_orthogonal_and_multiple_wraps(self):
        self.assertAlmostEqual(fix.periodic_squared([20.1,-9.8,0.3], [[2,0,0],[0,2,0],[0,0,2]]), .14)

    def test_triclinic(self):
        box = [[3,0,0],[1,3,0],[.5,.7,3]]
        self.assertAlmostEqual(fix.periodic_squared([8.1,-4.4,6.3],box), .14)

    def test_boundary(self):
        self.assertAlmostEqual(fix.periodic_squared([1,0,0],[[2,0,0],[0,2,0],[0,0,2]]),1)


class OpenMMRegression(unittest.TestCase):
    def system(self):
        from openmm import System, CustomExternalForce, Vec3, unit
        system = System()
        system.addParticle(12)
        system.setDefaultPeriodicBoxVectors(Vec3(2,0,0),Vec3(0,2,0),Vec3(0,0,2))
        force = CustomExternalForce(fix.BAD)
        force.addGlobalParameter('k',1000)
        for n in ['x0','y0','z0']: force.addPerParticleParameter(n)
        force.addParticle(0,[0,0,0])
        system.addForce(force)
        return system, [Vec3(20.1,-9.8,.3)]*unit.nanometer

    def test_energy_before_after_and_independent_reference(self):
        system,positions = self.system()
        i = fix.fix(system)
        report = fix.energy_comparison(system,i,positions)
        self.assertAlmostEqual(report['openmm_kj_mol']['before'],250070.,places=5)
        self.assertAlmostEqual(report['openmm_kj_mol']['after'],70.,places=5)
        self.assertGreater(report['openmm_kj_mol']['before'],report['openmm_kj_mol']['after']*1000)

    def test_no_wrap_equality(self):
        from openmm import Vec3, unit
        system,_ = self.system()
        i=fix.fix(system)
        report=fix.energy_comparison(system,i,[Vec3(.1,.2,.3)]*unit.nanometer)
        self.assertAlmostEqual(report['openmm_kj_mol']['before'],report['openmm_kj_mol']['after'])

    def test_unknown_or_double_patch_rejected(self):
        system,_ = self.system()
        fix.fix(system)
        with self.assertRaises(RuntimeError): fix.fix(system)

    def test_lattice_invariance(self):
        from openmm import Vec3, unit
        system,_=self.system()
        i=fix.fix(system)
        for xyz in [(2.1,.2,.3),(.1,-3.8,.3),(.1,.2,6.3)]:
            r=fix.energy_comparison(system,i,[Vec3(*xyz)]*unit.nanometer)
            self.assertAlmostEqual(r['openmm_kj_mol']['after'],70.,places=5)

    def test_triclinic_energy(self):
        from openmm import Vec3,unit
        system,_=self.system()
        system.setDefaultPeriodicBoxVectors(Vec3(3,0,0),Vec3(1,3,0),Vec3(.5,.7,3))
        i=fix.fix(system)
        report=fix.energy_comparison(system,i,[Vec3(8.1,-4.4,6.3)]*unit.nanometer)
        self.assertAlmostEqual(report['openmm_kj_mol']['after'],70.,places=5)


class RecoveryAndAuthorization(unittest.TestCase):
    def test_raw_recovery_ignores_invalid_legacy_states(self):
        from openmm import XmlSerializer,app,Vec3,unit,Platform,NonbondedForce,MonteCarloMembraneBarostat
        system,positions=OpenMMRegression().system()
        nb=NonbondedForce(); nb.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
        nb.setCutoffDistance(.9); nb.addParticle(0,.3,0); system.addForce(nb)
        system.addForce(MonteCarloMembraneBarostat(1*unit.bar,0*unit.bar*unit.nanometer,300*unit.kelvin,
                        MonteCarloMembraneBarostat.XYIsotropic,MonteCarloMembraneBarostat.ZFree,25))
        topology=app.Topology(); chain=topology.addChain(); residue=topology.addResidue('LIG',chain)
        topology.addAtom('C',app.element.carbon,residue)
        topology.setPeriodicBoxVectors([Vec3(2,0,0),Vec3(0,2,0),Vec3(0,0,2)]*unit.nanometer)
        actual_simulation=app.Simulation
        actual_compare=fix.energy_comparison
        actual_platform=Platform.getPlatformByName
        # This fixture runs Reference only. It tests orchestration, not CUDA acceptance.
        def reference_simulation(top,sys,integ,platform,properties):
            return actual_simulation(top,sys,integ,Platform.getPlatformByName('Reference'))
        def reference_compare(sys,index,xyz,platform,device):
            return actual_compare(sys,index,xyz,'Reference')
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)/'new'; old=Path(temp)/'old'; src=old/'systems/old_system'
            src.mkdir(parents=True); (root/'evidence').mkdir(parents=True)
            with (src/'input.pdb').open('w') as f:
                app.PDBFile.writeFile(topology,[Vec3(.1,.2,.3)]*unit.nanometer,f)
            (src/'system.xml').write_text(XmlSerializer.serialize(system))
            (src/'parameters.json').write_text('{}')
            (src/'state.xml').write_text('INVALID OLD STATE MUST NEVER DESERIALIZE')
            (src/'refined_state.xml').write_text('INVALID FAILED REFINEMENT')
            hashes={p.name:fix.sha(p) for p in src.iterdir()}
            fix.dump(root/'evidence/source_snapshot.json',{'systems':{'old_system':hashes}})
            row=dict(system_id='new_system',source_system_id='old_system',source_system_path='systems/old_system')
            with patch.object(fix,'ROOT',root),patch.object(fix,'source',return_value=old),\
                 patch.object(app,'Simulation',side_effect=reference_simulation),\
                 patch.object(Platform,'getPlatformByName',side_effect=lambda name:actual_platform('Reference')),\
                 patch.object(fix,'energy_comparison',side_effect=reference_compare):
                fix.recover(row,'0')
                audit=json.loads((root/'systems/new_system/build_audit.json').read_text())
                self.assertFalse(audit['reused_old_state'])
                self.assertEqual(fix.sha(src/'input.pdb'),fix.sha(root/'systems/new_system/input.pdb'))
                self.assertTrue(audit['minimized']['finite'])
                for name,digest in hashes.items(): self.assertEqual(fix.sha(src/name),digest)
                with patch.object(fix,'identity',return_value={'fixture_platform':'Reference'}):
                    fix.smoke('new_system','0',True)
                    self.assertTrue(fix.check_smoke('new_system',True))
                    smoke=json.loads((root/'single_apo/new_system/audit.json').read_text())
                    self.assertEqual(smoke['steps'],500)
                    self.assertTrue(smoke['checkpoint_roundtrip'])
                with self.assertRaises(FileExistsError): fix.recover(row,'0')

    def test_authorization_and_stale_gate(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); (root/'logs').mkdir()
            ident={'build_audit_sha256':'build','config_sha256':'config'}
            gate=dict(status='PASS',builds={'new':'build'},config_sha256='config',smokes={'new':'smoke'},
                      single_apo_sha256='smoke',storage_budget_sha256='storage')
            path=root/'logs'/f'server_preflight_{fix.VERSION}.json'
            path.write_text(json.dumps(gate))
            def fake_sha(p):
                return 'storage' if Path(p).name=='storage_budget.json' else 'hash'
            with patch.object(fix,'ROOT',root),patch.object(fix,'identity',return_value=ident),\
                 patch.object(fix,'table',return_value=[dict(system_id='new',context='apo')]),\
                 patch.object(fix,'check_smoke',return_value='smoke'),patch.object(fix,'sha',side_effect=fake_sha),\
                 patch.object(fix,'storage_budget',return_value={'sufficient':True}):
                with self.assertRaises(FileNotFoundError): fix.authorized_gate({'system_id':'new'})
                approval=dict(protocol_id=fix.VERSION,preflight_sha256='hash',scientific_review_approved=True,
                              user_explicit_authorization=True,reviewer='fixture',user_authorization_reference='fixture')
                (root/'authorization.json').write_text(json.dumps(approval))
                self.assertEqual(fix.authorized_gate({'system_id':'new'})['authorization_sha256'],'hash')
                gate['builds']['new']='old_build'; path.write_text(json.dumps(gate))
                with self.assertRaises(RuntimeError): fix.authorized_gate({'system_id':'new'})


if __name__=='__main__': unittest.main()
