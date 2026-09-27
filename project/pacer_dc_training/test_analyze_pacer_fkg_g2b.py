import importlib.util
import pathlib
import unittest
import numpy as np

spec=importlib.util.spec_from_file_location('g2b',pathlib.Path(__file__).with_name('analyze_pacer_fkg_g2b.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class G2BTests(unittest.TestCase):
    def test_global_shift_cancelled(self):
        a=np.random.default_rng(1).normal(size=(4,3,2))
        self.assertTrue(np.allclose(m.unit_normalize(a+12),m.unit_normalize(a)))
    def test_robust_channel(self):
        a=np.array([[[1.,2.]],[[3.,6.]],[[5.,10.]]])
        med,scale=m.fit_channels(a)
        self.assertTrue(np.allclose(med,[3,6]))
        self.assertTrue(np.all(scale>0))
    def test_bandwidth(self):
        x=np.array([[0.,0.],[1.,0.],[0.,1.]])
        self.assertAlmostEqual(m.shared_bandwidth(x),1.)
    def test_rff_same_inputs(self):
        W,b=m.make_rff(2,1.,1024,3)
        a=m.rff(np.array([[1.,2.],[1.,2.]]),W,b)
        self.assertTrue(np.allclose(a[0],a[1]))
        self.assertAlmostEqual(float(a[0]@a[1]),1,delta=.12)
    def test_rff_deterministic(self):
        w,b=m.make_rff(3,2.,32,4)
        v,c=m.make_rff(3,2.,32,4)
        self.assertTrue(np.array_equal(w,v));self.assertTrue(np.array_equal(b,c))
    def test_cosine(self):
        self.assertAlmostEqual(m.cosine(np.array([1.,0]),np.array([-1.,0])),-1.)
        self.assertIsNone(m.cosine(np.zeros(2),np.ones(2)))
    def test_bootstrap_agreement(self):
        blocks={r:{c:np.ones((5,5,2))*s for c,s in [('apo',-1),('candidate_no_probe',1)]} for r in (2,3)}
        out=m.bootstrap_direction(blocks,{'candidate_no_probe':1,'apo':-1},50,5)
        self.assertAlmostEqual(out['cosine_ci50'],1.)
        self.assertEqual(out['bootstrap_valid'],50)
    def test_bootstrap_opposite(self):
        blocks={2:{'apo':np.zeros((5,5,2)),'candidate_no_probe':np.ones((5,5,2))},
                3:{'apo':np.ones((5,5,2)),'candidate_no_probe':np.zeros((5,5,2))}}
        out=m.bootstrap_direction(blocks,{'candidate_no_probe':1,'apo':-1},20,3)
        self.assertAlmostEqual(out['cosine_ci50'],-1.)
    def test_csv_handles_all_fields(self):
        import tempfile,csv
        with tempfile.TemporaryDirectory() as t:
            p=pathlib.Path(t)/'t.csv';m.write_csv(p,[{'a':1,'b':2}])
            with p.open(encoding='utf-8-sig') as f:
                self.assertEqual(list(csv.DictReader(f))[0]['b'],'2')

if __name__=='__main__':unittest.main()
