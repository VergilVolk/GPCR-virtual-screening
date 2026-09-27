import importlib.util
from pathlib import Path
import unittest
import numpy as np

spec=importlib.util.spec_from_file_location('rob',Path(__file__).with_name('audit_pacer_fkg_g2b_robustness.py'))
rob=importlib.util.module_from_spec(spec);spec.loader.exec_module(rob)

class TestG2BRobustness(unittest.TestCase):
    def test_exact_gram_matches_explicit_rbf(self):
        x=np.array([[0.,0.],[1.,0.],[0.,2.]])
        g=rob.exact_gram(x,1.)
        self.assertTrue(np.allclose(g,g.T))
        self.assertTrue(np.allclose(np.diag(g),1.))
        self.assertAlmostEqual(g[0,1],np.exp(-.5))
    def test_identical_replicas_direction(self):
        rng=np.random.default_rng(2)
        x=rng.normal(size=(8,2));g=rob.exact_gram(x,2.)
        grp={(2,'apo'):slice(0,2),(2,'candidate_no_probe'):slice(2,4),
             (3,'apo'):slice(4,6),(3,'candidate_no_probe'):slice(6,8)}
        # same context coordinates in R2 and R3
        y=np.concatenate([x[:4],x[:4]])
        g=rob.exact_gram(y,2.)
        w2,w3=rob.signed_weights(grp,{'candidate_no_probe':1,'apo':-1})
        self.assertAlmostEqual(rob.cosine_from_gram(g,w2,w3),1.,places=7)
    def test_opposite_replicas_direction(self):
        x=np.array([[0.],[0.1],[1.],[1.1],[1.],[1.1],[0.],[0.1]])
        g=rob.exact_gram(x,0.8)
        grp={(2,'apo'):slice(0,2),(2,'candidate_no_probe'):slice(2,4),
             (3,'apo'):slice(4,6),(3,'candidate_no_probe'):slice(6,8)}
        w2,w3=rob.signed_weights(grp,{'candidate_no_probe':1,'apo':-1})
        self.assertAlmostEqual(rob.cosine_from_gram(g,w2,w3),-1.,places=7)
    def test_rff_converges_to_exact(self):
        rng=np.random.default_rng(9)
        x=rng.normal(size=(16,2)).astype(np.float32);g=rob.exact_gram(x,1.5)
        grp={(2,'apo'):slice(0,4),(2,'candidate_no_probe'):slice(4,8),
             (3,'apo'):slice(8,12),(3,'candidate_no_probe'):slice(12,16)}
        weights=rob.signed_weights(grp,{'candidate_no_probe':1,'apo':-1})
        W,b=rob.original.make_rff(2,1.5,20000,42)
        approx=rob.rff_cosine(rob.original.rff(x,W,b).astype(float),*weights)
        self.assertAlmostEqual(approx,rob.cosine_from_gram(g,*weights),delta=.06)
    def test_heldout_pairwise_error(self):
        g=np.eye(4); a=np.eye(4)
        self.assertEqual(rob.gram_errors(g,a,[2,3]),(0.,0.))
    def test_zero_norm_returns_undefined(self):
        g=np.ones((4,4)); w2=np.array([-.5,.5,0,0]);w3=np.array([0,0,-.5,.5])
        self.assertIsNone(rob.cosine_from_gram(g,w2,w3))
    def test_different_dims_distinct_deterministic(self):
        for dim in (256,512,1024):
            a,b=rob.original.make_rff(4,2,dim,18)
            c,d=rob.original.make_rff(4,2,dim,18)
            self.assertTrue(np.array_equal(a,c));self.assertTrue(np.array_equal(b,d))

if __name__=='__main__':unittest.main()
