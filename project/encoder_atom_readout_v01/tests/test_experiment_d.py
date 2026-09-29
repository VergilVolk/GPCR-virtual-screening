from __future__ import annotations

import unittest

import numpy as np
import torch

from project.encoder_atom_readout_v01.experiment_d import (
    BACKBONE_ATOMS,
    CHI1_QUADRUPLETS,
    atom_metadata,
    derive_readouts,
    grouped_mean,
)


class AtomReadoutTests(unittest.TestCase):
    def test_frozen_chi1_mapping_excludes_ala_and_gly(self):
        self.assertNotIn("A", CHI1_QUADRUPLETS)
        self.assertNotIn("G", CHI1_QUADRUPLETS)
        self.assertEqual(CHI1_QUADRUPLETS["I"], ("N", "CA", "CB", "CG1"))
        self.assertEqual(CHI1_QUADRUPLETS["C"], ("N", "CA", "CB", "SG"))

    def test_group_membership_and_zero_sidechain(self):
        residue_index, backbone, sidechain, records = atom_metadata("GA")
        self.assertEqual(records[0]["sidechain_count"], 0)
        self.assertFalse(records[0]["sidechain_present"])
        self.assertEqual(records[1]["sidechain_count"], 1)
        self.assertEqual(records[0]["backbone_atoms"], list(BACKBONE_ATOMS))
        descriptor = torch.arange(len(residue_index) * 2, dtype=torch.float32).reshape(len(residue_index), 2)
        pooled, counts = grouped_mean(
            descriptor,
            torch.as_tensor(residue_index),
            torch.as_tensor(sidechain),
            2,
        )
        torch.testing.assert_close(pooled[0], torch.zeros(2))
        self.assertEqual(float(counts[0]), 0.0)
        self.assertEqual(float(counts[1]), 1.0)

    def test_bs_and_mc_exact_blocks(self):
        stored = {
            "M": np.arange(12, dtype=np.float32).reshape(1, 2, 6),
            "B": np.ones((1, 2, 6), dtype=np.float32),
            "S": np.full((1, 2, 6), 2.0, dtype=np.float32),
        }
        result = derive_readouts(stored)
        width = stored["M"].shape[-1]
        self.assertTrue(np.array_equal(result["C1-BS"][..., :width], stored["B"]))
        self.assertTrue(np.array_equal(result["C1-BS"][..., width:], stored["S"]))
        self.assertTrue(np.array_equal(result["C1-MC"][..., :width], stored["M"]))
        self.assertEqual(np.count_nonzero(result["C1-MC"][..., width:]), 0)


if __name__ == "__main__":
    unittest.main()
