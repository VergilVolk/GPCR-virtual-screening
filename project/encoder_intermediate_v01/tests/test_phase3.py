from __future__ import annotations

import unittest

import numpy as np

from project.encoder_intermediate_v01.phase3 import torsion_targets, wrapped_angle_error


class TorsionTests(unittest.TestCase):
    def test_terminal_masks_and_finite_components(self):
        coords = np.zeros((2, 4, 14, 3), dtype=np.float32)
        for frame in range(2):
            for residue in range(4):
                offset = residue * 3.0
                coords[frame, residue, 0] = [offset, 0.2 * ((residue + frame) % 2), 0.1]
                coords[frame, residue, 1] = [offset + 1.0, 0.8, 0.3]
                coords[frame, residue, 2] = [offset + 2.0, 0.1, 1.0]
        target = torsion_targets(coords)
        self.assertFalse(target["phi_mask"][:, 0].any())
        self.assertFalse(target["psi_mask"][:, -1].any())
        self.assertTrue(target["phi_mask"][:, 1:].all())
        self.assertTrue(target["psi_mask"][:, :-1].all())

    def test_wrapped_angle_error(self):
        observed = np.deg2rad(np.asarray([179.0, -179.0]))
        predicted = np.deg2rad(np.asarray([-179.0, 179.0]))
        np.testing.assert_allclose(wrapped_angle_error(predicted, observed), [2.0, 2.0], atol=1e-12)


if __name__ == "__main__":
    unittest.main()
