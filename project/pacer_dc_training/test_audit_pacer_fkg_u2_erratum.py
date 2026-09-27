"""Synthetic mathematical and provenance-oriented tests. No user MD inputs needed."""
from __future__ import annotations
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from audit_pacer_fkg_u2_erratum import cosine, spearman, compute, verify_g2c

class ErratumTests(unittest.TestCase):
    def test_cosine_direction_and_zero(self):
        self.assertAlmostEqual(cosine(np.array([2.,0.]), np.array([4.,0.])), 1.)
        self.assertAlmostEqual(cosine(np.array([2.,0.]), np.array([-4.,0.])), -1.)
        self.assertIsNone(cosine(np.zeros(2), np.ones(2)))

    def test_spearman_with_ties(self):
        self.assertAlmostEqual(spearman([1.,1.,2.], [3.,3.,5.]), 1.)
        self.assertIsNone(spearman([1.,1.,1.], [1.,2.,3.]))

    def fixture(self):
        graph={'regions':{
            'stable_core_control':[{'embedding_index':0}],
            'region_a':[{'embedding_index':1}],
        }}
        old={'cross_replica_summary':[
            {'region':n,'axis':'synergy_interaction','legacy_magnitude_screen':n=='region_a'}
            for n in graph['regions']]}
        # Constant within each replica/window. Region_a larger than stable core
        # in R2 but lower in R3, so corrected diagnostic detects disagreement.
        blocks={r:np.zeros((5,2,2,3), dtype=np.float32) for r in (2,3)}
        blocks[2][:,:,0,0]=1.;blocks[2][:,:,1,0]=3.
        blocks[3][:,:,0,0]=2.;blocks[3][:,:,1,0]=-1.
        return graph,old,{'synergy_interaction':blocks}

    def test_common_kernel_norm_gap_and_direction(self):
        graph,old,groups=self.fixture()
        rows,windows=compute(groups,graph,old,draws=10,seed=4)
        self.assertEqual(len(rows),2)
        self.assertEqual(len(windows),20)
        r=next(x for x in rows if x['region']=='region_a')
        self.assertTrue(r['legacy_magnitude_screen_superseded'])
        self.assertAlmostEqual(r['r2_median_shared_norm_minus_stable'],2.)
        self.assertAlmostEqual(r['r3_median_shared_norm_minus_stable'],-1.)
        self.assertAlmostEqual(r['pooled_cross_replica_direction_cosine'],-1.)
        self.assertEqual(r['new_qualification_gate'],'NOT_DEFINED')
        self.assertEqual(r['block_draws'],10)

    def test_required_g2c_crosscheck_hard_fails(self):
        # verify_g2c requires 27 region-axis keys, then checks each exact cosine.
        items=[{'region':f'r{i}','axis':'synergy_interaction','no_graph_cosine':.5} for i in range(27)]
        prior={'G2C_STATUS':'COMPLETED_DESCRIPTIVE_NOT_INFERENTIAL',
               'region_axis_summary':items}
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'reference.json';p.write_text(json.dumps(prior))
            rows=[{'region':f'r{i}','axis':'synergy_interaction',
                   'pooled_cross_replica_direction_cosine':.5} for i in range(27)]
            self.assertEqual(verify_g2c(rows,p,1e-5)['status'],'PASS')
            rows[5]['pooled_cross_replica_direction_cosine']=.5001
            with self.assertRaises(AssertionError):verify_g2c(rows,p,1e-5)

if __name__=='__main__':unittest.main()
