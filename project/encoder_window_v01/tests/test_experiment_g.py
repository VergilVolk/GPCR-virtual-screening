import unittest

import numpy as np

from project.encoder_window_v01.experiment_g import block_branches, block_indices, derive_reps, target_branches


class WindowContractTests(unittest.TestCase):
    def setUp(self):
        self.z = np.arange(20 * 3 * 4, dtype=np.float32).reshape(20, 3, 4)

    def test_five_contiguous_blocks(self):
        blocks = block_indices(); self.assertEqual(len(blocks), 5); np.testing.assert_array_equal(np.concatenate(blocks), np.arange(100)); self.assertTrue(all(len(x) == 20 and np.all(np.diff(x) == 1) for x in blocks))

    def test_no_delta_crosses_block_or_window(self):
        self.assertEqual(sum(len(x) - 1 for x in block_indices()), 95); self.assertNotEqual(block_indices()[0][-1] + 1, block_indices()[1][0] - 1)

    def test_static_direct_mean(self):
        np.testing.assert_array_equal(block_branches(self.z)["static"], self.z.mean(0))

    def test_endpoint_identity(self):
        np.testing.assert_array_equal(block_branches(self.z)["signed"], (self.z[-1] - self.z[0]) / 19)

    def test_reversal_identities(self):
        f, r = block_branches(self.z), block_branches(self.z[::-1]); np.testing.assert_array_equal(f["static"], r["static"]); np.testing.assert_array_equal(f["signed"], -r["signed"]); np.testing.assert_array_equal(f["rms"], r["rms"])

    def test_mc_exact_embedding(self):
        reps = derive_reps(self.z, self.z + 1, self.z + 2)
        for branch in ("static", "signed", "rms"):
            m = block_branches(self.z)[branch]; mc = block_branches(reps["MC"])[branch]; np.testing.assert_array_equal(mc[..., :4], m); self.assertEqual(np.count_nonzero(mc[..., 4:]), 0)

    def test_target_identities(self):
        angle = np.linspace(-1, 1, 60).reshape(20, 3); valid = np.ones_like(angle, dtype=bool); f, fm = target_branches(angle, valid); r, rm = target_branches(angle[::-1], valid[::-1]); self.assertTrue(np.array_equal(fm, rm)); np.testing.assert_allclose(f["static"], r["static"], atol=1e-15); np.testing.assert_allclose(f["signed"], -r["signed"], atol=1e-15); np.testing.assert_allclose(f["rms"], r["rms"], atol=1e-15)


if __name__ == "__main__": unittest.main()
