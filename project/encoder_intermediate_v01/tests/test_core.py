from __future__ import annotations

import unittest

import torch

from project.encoder_intermediate_v01.core import invariant_residue_features


class ReadoutTests(unittest.TestCase):
    def test_equal_atom_mean_and_vector_norm(self):
        x = torch.tensor([[1.0, 3.0], [3.0, 5.0], [7.0, 11.0]])
        v = torch.zeros((3, 3, 2))
        v[0, 0] = torch.tensor([3.0, 4.0])
        v[1, 1] = torch.tensor([4.0, 3.0])
        v[2, 2] = torch.tensor([5.0, 12.0])
        observed = invariant_residue_features(x, v, torch.tensor([0, 0, 1]), 2)
        expected = torch.tensor([[2.0, 4.0, 3.5, 3.5], [7.0, 11.0, 5.0, 12.0]])
        torch.testing.assert_close(observed, expected)


if __name__ == "__main__":
    unittest.main()

