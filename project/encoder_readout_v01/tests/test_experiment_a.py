"""Synthetic-only numerical and rejection tests; scratch files stay in this worktree."""
import copy
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import experiment_a as exp


class DecompositionTests(unittest.TestCase):
    def setUp(self):
        self.regions = {"pair": np.array([3, 1]), "single": np.array([0]),
                        "overlap": np.array([1, 2, 3])}
        rng = np.random.default_rng(941)
        self.h = (rng.normal(size=(7, 4, 5)) + np.arange(7)[:, None, None] * 100).astype(np.float32)

    def test_analytic_values_and_shapes(self):
        h = np.array([[[1., 10.], [3., 14.], [8., 21.]]], dtype=np.float32)
        regions = {"first_last": np.array([0, 2])}
        out = exp.decompose(h, regions)
        np.testing.assert_array_equal(out["global_mean"], [[4., 15.]])
        np.testing.assert_array_equal(out["local_residual"], [[[-3., -5.], [-1., -1.], [4., 6.]]])
        np.testing.assert_array_equal(out["region_local_mean"], [[[.5, .5]]])
        np.testing.assert_array_equal(out["region_original_mean"], [[[4.5, 15.5]]])
        self.assertTrue(all(a.dtype == np.float64 for a in out.values()))

    def test_reconstruction_order_overlap_and_repeat(self):
        out = exp.decompose(self.h, self.regions)
        report = exp.verify(self.h, out, self.regions)
        self.assertTrue(all(report["repeat_bitwise_equal"].values()))
        np.testing.assert_allclose(out["local_residual"] + out["global_mean"][:, None], self.h, rtol=0, atol=1e-12)
        # No sorting of frames or nodes; a consistent permutation permutes residuals.
        frame_order = [6, 0, 2, 1, 5, 3, 4]
        shuffled = exp.decompose(self.h[frame_order], self.regions)
        np.testing.assert_array_equal(shuffled["local_residual"], out["local_residual"][frame_order])

    def test_global_shift_is_retained_in_global_branch(self):
        h = np.arange(60, dtype=np.float64).reshape(3, 4, 5)
        shifts = np.array([1., 2., 4.])[:, None, None] * 1024
        a, b = exp.decompose(h, self.regions), exp.decompose(h + shifts, self.regions)
        np.testing.assert_array_equal(a["local_residual"], b["local_residual"])
        np.testing.assert_array_equal(b["global_mean"] - a["global_mean"], np.broadcast_to(shifts[:, 0], (3, 5)))

    def test_local_residue_perturbation_is_retained(self):
        h = np.zeros((2, 4, 3), dtype=np.float32)
        perturbed = h.copy()
        perturbed[1, 1, 2] = 8
        regions = {"includes_target": np.array([1, 3]), "excludes_target": np.array([0, 2])}
        out = exp.decompose(perturbed, regions)
        expected_global = np.zeros((2, 3), dtype=np.float64)
        expected_global[1, 2] = 2
        expected_residual = np.zeros_like(perturbed, dtype=np.float64)
        expected_residual[1, :, 2] = [-2, 6, -2, -2]
        np.testing.assert_array_equal(out["global_mean"], expected_global)
        np.testing.assert_array_equal(out["local_residual"], expected_residual)
        np.testing.assert_array_equal(out["local_residual"] + out["global_mean"][:, None], perturbed)
        expected_regions = np.zeros((2, 2, 3), dtype=np.float64)
        expected_regions[1, 0, 2] = 4
        np.testing.assert_array_equal(out["region_original_mean"], expected_regions)
        np.testing.assert_array_equal(out["region_local_mean"] + out["global_mean"][:, None], expected_regions)
        exp.verify(perturbed, out, regions)

    def test_residue_permutation_preserves_identity_and_region_readouts(self):
        # Integer-valued inputs keep the mean exact, isolating ordering from rounding.
        h = np.arange(60, dtype=np.float32).reshape(3, 4, 5)
        order = np.array([2, 0, 3, 1])  # new index -> original residue index
        inverse = np.argsort(order)
        remapped_regions = {name: inverse[ix] for name, ix in self.regions.items()}
        original = exp.decompose(h, self.regions)
        permuted = exp.decompose(h[:, order, :], remapped_regions)
        np.testing.assert_array_equal(permuted["global_mean"], original["global_mean"])
        np.testing.assert_array_equal(permuted["local_residual"], original["local_residual"][:, order, :])
        for name in ("region_local_mean", "region_original_mean"):
            np.testing.assert_array_equal(permuted[name], original[name])
        exp.verify(h[:, order, :], permuted, remapped_regions)

    def test_layout_independent_and_zero_tensor(self):
        a = exp.decompose(self.h, self.regions)
        b = exp.decompose(np.asfortranarray(self.h), self.regions)
        for k in a:
            np.testing.assert_array_equal(a[k], b[k])
        zero = np.zeros_like(self.h)
        report = exp.verify(zero, exp.decompose(zero, self.regions), self.regions)
        self.assertEqual(report["checks"]["residue_reconstruction"]["max_abs"], 0)
        self.assertIsNone(report["checks"]["residue_reconstruction"]["relative_l2"])

    def test_bad_features_regions_and_corrupted_readout_rejected(self):
        for h in (np.full_like(self.h, np.nan), self.h[0], self.h.astype(int), self.h[:0]):
            with self.assertRaises(ValueError):
                exp.decompose(h, self.regions)
        for indices in ([1, 1], [-1], [4], [], [1.5]):
            with self.assertRaises(ValueError):
                exp.decompose(self.h, {"bad": np.asarray(indices)})
        out = exp.decompose(self.h, self.regions)
        out["region_local_mean"][0, 0, 0] += 1
        with self.assertRaises(ValueError):
            exp.verify(self.h, out, self.regions)


class ArchiveContractTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=Path(__file__).parent, prefix="scratch_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.h = np.arange(6 * 5 * 4, dtype=np.float32).reshape(6, 5, 4)
        self.data = {"residue_features": self.h, "global_features": self.h.mean(axis=1),
                     "sequence": np.asarray("ACDEF"), "frame_ids": np.arange(6, dtype=np.int64)}

    def write(self, data=None):
        path = self.root / "input.npz"
        np.savez(path, **(self.data if data is None else data))
        return {"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

    def read(self, row):
        return exp.read_input(self.root, row, "ACDEF", shape=(6, 5, 4))

    def test_valid_input_and_hash_and_missing_rejections(self):
        row = self.write()
        np.testing.assert_array_equal(self.read(row)["residue_features"], self.h)
        with self.assertRaises(ValueError):
            self.read({**row, "sha256": "0" * 64})
        with self.assertRaises(FileNotFoundError):
            self.read({**row, "path": "missing.npz"})
        with self.assertRaises(ValueError):
            self.read({**row, "path": "../escape.npz"})

    def test_dtype_shape_sequence_frame_and_finiteness_rejections(self):
        cases = [
            ("residue_features", self.h.astype(np.float64)),
            ("residue_features", self.h[:, :-1]),
            ("residue_features", np.full_like(self.h, np.nan)),
            ("global_features", np.full((6, 4), np.inf, dtype=np.float32)),
            ("global_features", np.zeros((6, 4), dtype=np.float32)),
            ("global_features", self.data["global_features"].astype(np.float64)),
            ("sequence", np.asarray("WRONG")),
            ("sequence", np.asarray(["ACDEF"])),
            ("frame_ids", np.arange(6, dtype=np.int32)),
            ("frame_ids", np.arange(6, dtype=np.int64)[::-1]),
            ("frame_ids", np.arange(6, dtype=np.int64) + 100),
            ("residue_index", np.arange(5, dtype=np.int64)[::-1]),
            ("optional_numeric", np.asarray([np.inf])),
        ]
        for key, value in cases:
            with self.subTest(key=key, dtype=value.dtype, shape=value.shape):
                row = self.write({**self.data, key: value})
                with self.assertRaises(ValueError):
                    self.read(row)
        data = dict(self.data)
        del data["frame_ids"]
        with self.assertRaises(ValueError):
            self.read(self.write(data))

    def test_no_output_on_failed_preflight_and_existing_path_rejected(self):
        with patch.object(exp, "OUTPUT_BASE", self.root):
            with self.assertRaises(FileExistsError):
                self._existing()
            with self.assertRaises(ValueError):
                exp.check_output_root(self.root.parent / "outside", self.root / "input")
            target = self.root / "run"
            with patch.object(exp, "preflight", side_effect=ValueError("bad input")):
                with self.assertRaises(ValueError):
                    exp.run(self.root / "input", target)
            self.assertFalse(target.exists())
            with self.assertRaises(ValueError):
                exp.check_output_root(self.root / "input" / "run", self.root / "input")

    def _existing(self):
        path = self.root / "existing"
        path.mkdir()
        exp.check_output_root(path, self.root / "input")

    def test_region_identity_and_duplicate_rejections(self):
        graph = {"nodes": [{"embedding_index": i, "site": f"A{i}"} for i in range(3)],
                 "regions": {"a": [{"embedding_index": 2, "site": "A2"},
                                   {"embedding_index": 0, "site": "A0"}]}}
        np.testing.assert_array_equal(exp.region_indices(graph)["a"], [2, 0])
        bad = copy.deepcopy(graph)
        bad["regions"]["a"][0]["site"] = "X2"
        with self.assertRaises(ValueError):
            exp.region_indices(bad)
        bad = copy.deepcopy(graph)
        bad["regions"]["a"].append(bad["regions"]["a"][0])
        with self.assertRaises(ValueError):
            exp.region_indices(bad)


if __name__ == "__main__":
    unittest.main()
