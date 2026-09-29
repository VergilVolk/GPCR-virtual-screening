from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import run_phase1_bs256 as runner

from run_phase1_bs256 import (
    BS_WIDTH,
    N_RESIDUES,
    Job,
    compare_manifest_identity,
    expected_manifest_identity,
    finite_array,
    parse_args,
    phase0_jobs,
    replay_metrics,
    sha256,
    validate_cache,
)


class Phase1BS256Tests(unittest.TestCase):
    def test_modes_are_mutually_exclusive_and_required(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parse_args([])
            with self.assertRaises(SystemExit):
                parse_args(["--smoke", "--full"])
        self.assertTrue(parse_args(["--smoke"]).smoke)
        self.assertTrue(parse_args(["--full"]).full)
        self.assertTrue(parse_args(["--verify"]).verify)

    def test_phase0_job_contracts(self) -> None:
        smoke = phase0_jobs("smoke")
        full = phase0_jobs("full")
        self.assertEqual(len(smoke), 1)
        self.assertEqual((smoke[0].system, smoke[0].replica, smoke[0].frame_ids), ("apo", 2, (0, 1, 2, 3)))
        self.assertEqual(smoke[0].output_shape, (4, N_RESIDUES, BS_WIDTH))
        self.assertEqual(len(full), 12)
        self.assertTrue(all(job.output_shape == (1000, N_RESIDUES, BS_WIDTH) for job in full))

    def test_manifest_identity_detects_any_contract_change(self) -> None:
        job = phase0_jobs("smoke")[0]
        expected = expected_manifest_identity(job)
        self.assertEqual(compare_manifest_identity(dict(expected), expected), [])
        changed = dict(expected)
        changed["selection"] = "protein"
        self.assertEqual(compare_manifest_identity(changed, expected), ["selection"])

    def test_array_validation_requires_exact_shape_dtype_and_finite(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "cache.npy"
            shape = (4, N_RESIDUES, BS_WIDTH)
            np.save(path, np.zeros(shape, dtype=np.float32))
            array, valid = finite_array(path, shape)
            self.assertTrue(valid)
            del array
            np.save(path, np.zeros(shape, dtype=np.float64))
            array, valid = finite_array(path, shape)
            self.assertFalse(valid)
            del array
            bad = np.zeros(shape, dtype=np.float32)
            bad[0, 0, 0] = np.nan
            np.save(path, bad)
            array, valid = finite_array(path, shape)
            self.assertFalse(valid)
            del array

    def test_valid_cache_is_reusable_and_corruption_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with mock.patch.object(runner, "CACHE_ROOT", root / "cache"), mock.patch.object(runner, "REPORT_ROOT", root / "reports"):
                job = phase0_jobs("smoke")[0]
                job.cache_path.parent.mkdir(parents=True)
                job.manifest_path.parent.mkdir(parents=True)
                np.save(job.cache_path, np.zeros(job.output_shape, dtype=np.float32))
                manifest = expected_manifest_identity(job)
                manifest["cache"] = {"path": str(job.cache_path), "sha256": sha256(job.cache_path), "bytes": job.cache_path.stat().st_size}
                manifest["replay"] = {"passed": True, "frame_ids": [0, 1, 2, 3]}
                job.manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                status, source = validate_cache(job, verify_source_hashes=False)
                self.assertEqual(status["status"], "valid")
                self.assertIsNone(source)
                with job.cache_path.open("r+b") as handle:
                    handle.seek(-1, 2)
                    handle.write(b"X")
                with self.assertRaisesRegex(RuntimeError, "SHA256 mismatch"):
                    validate_cache(job, verify_source_hashes=False)

    def test_replay_metric_accepts_identical_and_rejects_large_change(self) -> None:
        reference = np.ones((4, N_RESIDUES, BS_WIDTH), dtype=np.float32)
        exact = replay_metrics(reference, reference.copy())
        self.assertTrue(exact["passed"])
        self.assertTrue(exact["bitwise_identical"])
        changed = reference.copy()
        changed[..., 0] += 1.0
        self.assertFalse(replay_metrics(reference, changed)["passed"])


if __name__ == "__main__":
    unittest.main()
