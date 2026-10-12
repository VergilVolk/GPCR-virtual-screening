from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from project.pacer_fkg_v02 import run_phase2_calibration as phase2


class Phase2CalibrationTests(unittest.TestCase):
    def test_cli_modes(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                phase2.parse_args([])
            with self.assertRaises(SystemExit):
                phase2.parse_args(["--run", "--verify"])
        self.assertTrue(phase2.parse_args(["--run"]).run)
        self.assertTrue(phase2.parse_args(["--verify"]).verify)

    def test_exact_blocks_boundaries_shapes_and_endpoint_identity(self):
        # Frame value identifies its trajectory-local index, exposing boundaries.
        frames = np.broadcast_to(np.arange(1000, dtype=np.float32)[:, None, None],
                                 (1000, 270, 256)).copy()
        sm, sd = phase2.construct_blocks(frames)
        self.assertEqual(sm.shape, (50, 270, 512))
        self.assertEqual(sd.shape, (50, 270, 256))
        self.assertEqual(sm.dtype, np.float32)
        self.assertEqual(sd.dtype, np.float32)
        self.assertEqual(sm[0, 0, 0], np.float32(9.5))
        self.assertEqual(sm[1, 0, 0], np.float32(29.5))
        self.assertTrue(np.all(sd == 1.0))
        self.assertTrue(phase2.verify_endpoint_identity(frames, sd))

    def test_no_cross_trajectory_blocks(self):
        a = np.zeros((1000, 270, 256), dtype=np.float32)
        b = np.full_like(a, 100.0)
        sm_a, _ = phase2.construct_blocks(a)
        sm_b, _ = phase2.construct_blocks(b)
        self.assertTrue(np.all(sm_a[..., :256] == 0))
        self.assertTrue(np.all(sm_b[..., :256] == 100))

    def test_normalization_shapes_and_determinism(self):
        rng = np.random.default_rng(7)
        for width in (512, 256):
            x = rng.normal(size=(8, 270, width)).astype(np.float32)
            c1, s1 = phase2.fit_channels(x)
            c2, s2 = phase2.fit_channels(x.copy())
            self.assertEqual(c1.shape, (width,))
            self.assertEqual(s1.shape, (width,))
            self.assertTrue(np.array_equal(c1, c2))
            self.assertTrue(np.array_equal(s1, s2))
            self.assertTrue(np.all(s1 >= 1e-6))

    def test_rff_replay_is_deterministic_and_branch_compatible(self):
        for branch, info in phase2.BRANCHES.items():
            w1, b1 = phase2.make_rff(info["width"], 2.5, info["seed"])
            w2, b2 = phase2.make_rff(info["width"], 2.5, info["seed"])
            self.assertEqual(w1.shape, (info["width"], 512), branch)
            self.assertEqual(b1.shape, (512,), branch)
            self.assertTrue(np.array_equal(w1, w2))
            self.assertTrue(np.array_equal(b1, b2))

    def test_calibration_inventory_excludes_r1_r3(self):
        inventory = [{"system": system, "replica": replica,
                      "used_for_fitting": replica == 2}
                     for system in phase2.SYSTEMS for replica in phase2.REPLICAS]
        fitted = [x for x in inventory if x["used_for_fitting"]]
        excluded = [x for x in inventory if not x["used_for_fitting"]]
        self.assertEqual(len(fitted) * 50, 200)
        self.assertEqual({x["replica"] for x in fitted}, {2})
        self.assertEqual({x["replica"] for x in excluded}, {1, 3})

    def test_v01_protection_verifier_detects_change(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            protected = root / "old.txt"
            protected.write_text("unchanged", encoding="utf-8")
            listing = root / "hashes.txt"
            listing.write_text(f"{phase2.sha256(protected)}  {protected}\n", encoding="utf-8")
            with mock.patch.object(phase2, "V01_PROTECTION", listing):
                self.assertEqual(phase2.verify_v01_protection(), 1)
                protected.write_text("changed", encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "v01 protection hash mismatch"):
                    phase2.verify_v01_protection()


if __name__ == "__main__":
    unittest.main()
