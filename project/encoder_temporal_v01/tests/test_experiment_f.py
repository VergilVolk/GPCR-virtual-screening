import unittest

import numpy as np

from project.encoder_temporal_v01.experiment_f import (
    contiguous_blocks,
    derive_representations,
    pair_target,
    temporal_descriptors,
    zero_baseline_metrics,
)


class TemporalContractTests(unittest.TestCase):
    def setUp(self):
        self.z = np.arange(100 * 3 * 4, dtype=np.float32).reshape(100, 3, 4)

    def test_delta_and_reverse(self):
        delta, _ = temporal_descriptors(self.z)
        np.testing.assert_array_equal(delta, self.z[1:] - self.z[:-1])
        np.testing.assert_array_equal(self.z[:-1] - self.z[1:], -delta)

    def test_target_reverse(self):
        angle = np.linspace(-1, 1, 300).reshape(100, 3)
        dy, mask = pair_target(angle, np.ones((100, 3), dtype=bool))
        encoded = np.stack((np.sin(angle), np.cos(angle)), axis=-1)
        np.testing.assert_allclose(encoded[:-1] - encoded[1:], -dy, rtol=0, atol=0)
        self.assertTrue(mask.all())

    def test_mc_exact_embedding(self):
        reps = derive_representations(self.z, self.z + 1, self.z + 2)
        dm, am = temporal_descriptors(reps["M"]); dmc, amc = temporal_descriptors(reps["MC"])
        np.testing.assert_array_equal(dmc[..., :4], dm); np.testing.assert_array_equal(amc[..., :4], am)
        self.assertEqual(np.count_nonzero(dmc[..., 4:]), 0); self.assertEqual(np.count_nonzero(amc[..., 4:]), 0)

    def test_blocks_exact_contiguous_exhaustive(self):
        blocks = contiguous_blocks(100)
        self.assertEqual(len(blocks), 5)
        np.testing.assert_array_equal(np.concatenate(blocks), np.arange(100))
        for block in blocks:
            self.assertEqual(len(block), 20); np.testing.assert_array_equal(np.diff(block), np.ones(19))

    def test_no_cross_window_pair(self):
        self.assertEqual(temporal_descriptors(np.zeros((100, 2, 3)))[0].shape[0], 99)
        self.assertEqual(temporal_descriptors(np.zeros((100, 2, 3)))[0].shape[0] * 2, 198)

    def test_zero_baseline(self):
        dy = np.asarray([[3.0, 4.0], [0.0, 0.0]])
        metric = zero_baseline_metrics(dy)
        self.assertAlmostEqual(metric["delta_sin_rmse"], np.sqrt(4.5)); self.assertAlmostEqual(metric["delta_cos_rmse"], np.sqrt(8.0))
        self.assertAlmostEqual(metric["combined_component_rmse"], 2.5)


if __name__ == "__main__":
    unittest.main()
