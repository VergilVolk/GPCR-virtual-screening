import sys
import unittest
from pathlib import Path

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


if __name__=='__main__': unittest.main()
