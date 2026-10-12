"""Regression checks for the sample adapter against the actual Stage4 formula."""
import unittest
import numpy as np
from run_module4_example import sample_descriptors, verify_vendor
from project.pacer_fkg_v02.run_stage4_prospective_phase2a_frozen_apply_v01 import construct


class SampleAdapterTests(unittest.TestCase):
    def test_one_block_equals_formal_first_block(self):
        frames = np.random.default_rng(27101).normal(size=(20, 270, 256)).astype(np.float32)
        # Synthetic feature tensors test only mathematical equivalence; they
        # are never written as MD or used as demonstration inference inputs.
        raw = np.zeros((1000, 270, 256), dtype=np.float32)
        raw[np.arange(0, 100, 5)] = frames
        expected = construct(raw)
        actual = sample_descriptors(frames)
        for branch in actual:
            np.testing.assert_array_equal(actual[branch], expected[branch][:1])

    def test_reject_incomplete_block(self):
        with self.assertRaises(ValueError):
            sample_descriptors(np.zeros((4, 270, 256), dtype=np.float32))

    def test_reject_nonfinite_features(self):
        frames = np.zeros((20, 270, 256), dtype=np.float32)
        frames[0, 0, 0] = np.nan
        with self.assertRaises(ValueError):
            sample_descriptors(frames)

    def test_constant_trajectory_has_no_motion_or_drift(self):
        frames = np.full((20, 270, 256), 2, dtype=np.float32)
        actual = sample_descriptors(frames)
        self.assertEqual(np.count_nonzero(actual["SIGNED_DRIFT"]), 0)
        self.assertEqual(np.count_nonzero(actual["STATE_MOTION"][..., 256:]), 0)
        np.testing.assert_array_equal(actual["STATE_MOTION"][..., :256], 2)

    def test_upstream_source_identity(self):
        self.assertEqual(verify_vendor()["commit"], "371d642ec1061664f16e49fcac702d07fc8d0b51")


if __name__ == "__main__":
    unittest.main()
