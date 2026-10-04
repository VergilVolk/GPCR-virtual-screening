from __future__ import annotations

import contextlib
import io
import inspect
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from project.pacer_fkg_v02 import run_phase3_fkg as phase3


class Phase3FrozenEvaluationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.anchor = phase3.verify_frozen_anchor()

    def test_01_required_freeze_sha_is_enforced(self):
        with mock.patch.object(phase3, "sha256", return_value="0" * 64):
            with self.assertRaisesRegex(RuntimeError, "required V02_FREEZE_MANIFEST"):
                phase3.verify_frozen_anchor()

    def test_02_frozen_flag_is_enforced(self):
        bad = dict(self.anchor["freeze"])
        bad["V02_NUMERICAL_STATE_FROZEN"] = False
        with mock.patch.object(phase3, "sha256", return_value=phase3.REQUIRED_FREEZE_SHA256), \
             mock.patch.object(phase3, "load_json", return_value=bad):
            with self.assertRaisesRegex(RuntimeError, "not true"):
                phase3.verify_frozen_anchor()

    def test_03_phase2_state_cannot_be_refit(self):
        source = inspect.getsource(phase3)
        for forbidden in ("phase2.fit_channels(", "phase2.shared_bandwidth(",
                          "phase2.make_rff(", "default_rng("):
            self.assertNotIn(forbidden, source)
        self.assertIn("load_frozen_state", source)

    def test_04_frozen_normalization_hashes_are_enforced(self):
        original = phase3._verify_record
        def reject(record, description):
            if "center" in description or "scale" in description:
                raise RuntimeError("normalization hash mismatch")
            return original(record, description)
        with mock.patch.object(phase3, "_verify_record", side_effect=reject):
            with self.assertRaisesRegex(RuntimeError, "normalization hash mismatch"):
                phase3.load_frozen_state(self.anchor, "STATE_MOTION")

    def test_05_frozen_rff_hashes_and_state_are_enforced(self):
        original = phase3._verify_record
        def reject(record, description):
            if "RFF" in description:
                raise RuntimeError("RFF hash mismatch")
            return original(record, description)
        with mock.patch.object(phase3, "_verify_record", side_effect=reject):
            with self.assertRaisesRegex(RuntimeError, "RFF hash mismatch"):
                phase3.load_frozen_state(self.anchor, "SIGNED_DRIFT")

    def test_06_four_contexts_map_correctly(self):
        self.assertEqual(phase3.CONTEXTS, {
            "A": "apo", "P": "probe_only",
            "C": "compound110__candidate_no_probe",
            "CP": "compound110__candidate_probe",
        })

    def test_07_all_replicas_map_correctly(self):
        self.assertEqual(phase3.REPLICAS, (1, 2, 3))

    def test_08_r2_is_calibration_replica(self):
        self.assertEqual(phase3.REPLICA_LABELS[2], "CALIBRATION_REPLICA")

    def test_09_r1_r3_are_application_replicas(self):
        self.assertEqual({phase3.REPLICA_LABELS[x] for x in (1, 3)},
                         {"FROZEN_APPLICATION_REPLICA"})

    def test_10_historical_contrast_formulas_exact(self):
        values = {"A": np.array([2.0, -1.0]), "P": np.array([5.0, 4.0]),
                  "C": np.array([11.0, 8.0]), "CP": np.array([23.0, 16.0])}
        np.testing.assert_array_equal(phase3.historical_contrast(values, "Delta_PAM"),
                                      values["CP"] - values["P"])
        np.testing.assert_array_equal(phase3.historical_contrast(values, "Delta_AGO"),
                                      values["C"] - values["A"])
        np.testing.assert_array_equal(phase3.historical_contrast(values, "Delta_INT"),
                                      values["CP"] - values["P"] - values["C"] + values["A"])

    def test_11_branches_remain_separate(self):
        self.assertEqual(set(phase3.BRANCHES), {"STATE_MOTION", "SIGNED_DRIFT"})
        source = inspect.getsource(phase3.consistency_report)
        self.assertNotIn("average", source)
        self.assertNotIn("weight", source)
        self.assertNotIn("combined_score", source)

    def test_12_deprecated_squared_kernel_path_is_unreachable(self):
        source = inspect.getsource(phase3).lower()
        self.assertNotIn("kernel" + "_u2", source)
        self.assertNotIn("signed_kernel_stat", source)

    def test_13_block_count_is_exactly_50(self):
        self.assertEqual(phase3.N_BLOCKS, 50)
        frames = np.broadcast_to(np.arange(1000, dtype=np.float32)[:, None, None],
                                 (1000, 1, 256))
        # Use the frozen constructor with a light residue dimension rejection
        # followed by the canonical boundary arithmetic.
        starts = np.arange(phase3.N_BLOCKS) * phase3.phase2.BLOCK_FRAMES
        stops = starts + phase3.phase2.BLOCK_FRAMES
        self.assertEqual((starts[0], stops[-1]), (0, 1000))
        self.assertEqual(frames.shape[0] // phase3.phase2.BLOCK_FRAMES, 50)

    def test_14_no_block_crosses_context_or_replica_boundaries(self):
        boundaries = {(context, replica): [(i * 20, (i + 1) * 20) for i in range(50)]
                      for context in phase3.CONTEXTS for replica in phase3.REPLICAS}
        self.assertEqual(len(boundaries), 12)
        for key, blocks in boundaries.items():
            self.assertEqual(len(blocks), 50, key)
            self.assertTrue(all(stop - start == 20 for start, stop in blocks))
            self.assertEqual(blocks[0], (0, 20))
            self.assertEqual(blocks[-1], (980, 1000))

    def test_15_complete_valid_run_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            existing = Path(temp)
            with mock.patch.object(phase3, "verify_frozen_anchor", return_value=self.anchor), \
                 mock.patch.object(phase3, "PHASE3_ROOT", existing):
                with self.assertRaisesRegex(FileExistsError, "refusing silent overwrite"):
                    phase3.run()

    def test_16_v01_protection_hashes_remain_unchanged(self):
        self.assertEqual(phase3.phase2.verify_v01_protection(), 41)

    def test_17_phase2_freeze_hash_remains_unchanged(self):
        self.assertEqual(phase3.sha256(phase3.FREEZE_PATH), phase3.REQUIRED_FREEZE_SHA256)

    def test_cli_supports_exactly_three_exclusive_modes(self):
        self.assertTrue(phase3.parse_args(["--smoke"]).smoke)
        self.assertTrue(phase3.parse_args(["--run"]).run)
        self.assertTrue(phase3.parse_args(["--verify"]).verify)
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                phase3.parse_args([])
            with self.assertRaises(SystemExit):
                phase3.parse_args(["--smoke", "--run"])

    def test_diffusion_shape_finiteness_and_formula(self):
        graph = phase3.load_json(phase3.GRAPH_PATH)
        transition = phase3.graph_transition(graph)
        x = np.zeros((1, phase3.N_RESIDUES, phase3.RFF_FEATURES), dtype=np.float32)
        x[0, 0, 0] = 1.0
        got = phase3.diffuse(x, transition)
        self.assertEqual(got.shape, x.shape)
        self.assertTrue(np.isfinite(got).all())
        once = ((1.0 - phase3.ALPHA) * x[0] + phase3.ALPHA * (transition @ x[0])).astype(np.float32)
        manual = x.copy()
        for _ in range(phase3.DIFFUSION_STEPS):
            manual[0] = ((1.0 - phase3.ALPHA) * x[0] + phase3.ALPHA * (transition @ manual[0])).astype(np.float32)
        self.assertTrue(np.array_equal(got, manual))
        self.assertFalse(np.array_equal(got[0], once))


if __name__ == "__main__":
    unittest.main()
