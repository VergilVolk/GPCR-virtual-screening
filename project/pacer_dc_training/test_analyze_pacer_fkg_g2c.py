import unittest
import numpy as np
from analyze_pacer_fkg_g2c import norm, diff_graph, permute_contacts, transition, graph_edges, bootstrap_pair, rff_map

class GraphAblationTests(unittest.TestCase):
 def setUp(self):
  self.edges={(0,1):('backbone',1.),(1,2):('backbone',1.),(2,3):('backbone',1.),(3,4):('backbone',1.),(4,5):('backbone',1.),
              (0,2):('contact',.8),(0,3):('contact',.4),(0,4):('contact',.6),(1,3):('contact',.3),(1,4):('contact',.9),(1,5):('contact',.7),(2,4):('contact',.8),(2,5):('contact',.5),(3,5):('contact',.6)}
 def test_cosine_direction(self):
  self.assertAlmostEqual(norm(np.array([1.,0]),np.array([-1.,0])),-1.)
  self.assertIsNone(norm(np.zeros(2),np.ones(2)))
 def test_signed_diffusion_linearity(self):
  t=transition(self.edges,n=6)
  x=np.arange(12.).reshape(6,2);y=np.ones_like(x)
  np.testing.assert_allclose(diff_graph(x+y,t),diff_graph(x,t)+diff_graph(y,t),atol=1e-12)
 def test_constant_field_preserved(self):
  t=transition(self.edges,n=6)
  a=np.ones((6,5));np.testing.assert_allclose(diff_graph(a,t),a,atol=1e-12)
 def test_weighted_transition_stochastic(self):
  t=transition(self.edges,n=6)
  np.testing.assert_allclose(np.asarray(t.sum(axis=1)).ravel(),np.ones(6))
 def test_rewire_preserves_degree_and_backbone(self):
  rng=np.random.default_rng(7)
  # Dense 6-node graph can have no legal swap. Use a spacious 12-node sparse contact graph.
  edges={(i,i+1):('backbone',1.) for i in range(11)}
  contacts=[(0,4),(1,6),(2,8),(3,10),(4,9),(0,7),(5,11),(2,10)]
  edges.update({x:('contact',.5) for x in contacts})
  result,n=permute_contacts(edges,rng,swaps=3000)
  self.assertGreater(n,0)
  for idx in range(12):
   self.assertEqual(sum(idx in k for k in edges),sum(idx in k for k in result))
  for e in edges:
   if edges[e][0]=='backbone':self.assertEqual(result[e],('backbone',1.))
  self.assertEqual(sorted(v[1] for v in edges.values() if v[0]=='contact'),sorted(v[1] for v in result.values() if v[0]=='contact'))
 def test_reproducible_shuffle(self):
  edges={(i,i+1):('backbone',1.) for i in range(11)}
  edges.update({x:('contact',.5) for x in [(0,4),(1,6),(2,8),(3,10),(4,9),(0,7),(5,11),(2,10)]})
  a,n1=permute_contacts(edges,np.random.default_rng(41),swaps=2000)
  b,n2=permute_contacts(edges,np.random.default_rng(41),swaps=2000)
  self.assertEqual(a,b);self.assertEqual(n1,n2)
 def test_rff_output_shape(self):
  a=np.zeros((100,270,128),dtype=np.float32);w=np.ones((128,32),dtype=np.float32)*.1;b=np.zeros(32,dtype=np.float32)
  z=rff_map(a,w,b)
  self.assertEqual(z.shape,(100,270,32));self.assertTrue(np.isfinite(z).all())
 def test_block_resample_identical_models(self):
  rng=np.random.default_rng(7)
  a={2:rng.normal(size=(5,5,6,8)),3:rng.normal(size=(5,5,6,8))}
  b=bootstrap_pair(a,a,np.array([0,2]),30,12)
  self.assertEqual(b['draws_valid'],30)
  self.assertAlmostEqual(b['delta_cosine_ci025'],0.)
  self.assertAlmostEqual(b['delta_cosine_ci975'],0.)

if __name__=='__main__':unittest.main()
