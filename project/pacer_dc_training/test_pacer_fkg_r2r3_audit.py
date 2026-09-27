"""Engineering-only tests for PACER-FKG R2/R3 preflight; no formal MD analysis."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from run_pacer_fkg_r2r3_audit import preflight, spearman


GRAPH = Path(__file__).resolve().parents[1] / "results" / "pacer_dc_geom2vec_pilot_v01" / "M4_MULTISTRUCTURE_GRAPH_v01.json"
CONTEXTS = ("apo", "probe_only", "candidate_no_probe", "candidate_probe")


class FkgPreflightTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        base = Path(self.tmp.name)
        self.input_root = base / "embeddings"
        self.output_root = base / "new_results"
        self.graph = json.loads(GRAPH.read_text(encoding="utf-8"))
        self.sequence = "".join(node["site"][0] for node in self.graph["nodes"])
        self.paths = {}
        for replica in (2, 3):
            for context in CONTEXTS:
                root = self.input_root / f"replica_{replica:02d}" / "window_000" / "atom14"
                root.mkdir(parents=True, exist_ok=True)
                path = root / f"{context}_w000.geom2vec.npz"
                self.write_npz(path)
                self.paths[(replica, context)] = path
        self.args = SimpleNamespace(
            input_root=self.input_root,
            output_root=self.output_root,
            graph=GRAPH,
            replicas=[2, 3],
            windows=[0],
            expected_frames=4,
        )

    def write_npz(self, path, sequence=None, nonfinite=False):
        features = np.zeros((4, 270, 128), dtype=np.float32)
        if nonfinite:
            features[0, 0, 0] = np.nan
        np.savez_compressed(
            path,
            residue_features=features,
            frame_ids=np.arange(4, dtype=np.int64),
            sequence=np.asarray(self.sequence if sequence is None else sequence),
        )

    def test_spearman_ties_and_constant_input(self):
        self.assertAlmostEqual(spearman([1, 1, 2], [2, 2, 3]), 1.0)
        self.assertIsNone(spearman([1, 1, 1], [2, 3, 4]))

    def test_preflight_accepts_complete_fixture_without_writing(self):
        result = preflight(self.args)
        self.assertEqual(result["n_inputs"], 8)
        self.assertEqual(len(result["graph_sha256"]), 64)
        self.assertFalse(self.output_root.exists())

    def test_preflight_refuses_existing_results(self):
        self.output_root.mkdir()
        with self.assertRaises(FileExistsError):
            preflight(self.args)

    def test_preflight_refuses_missing_context(self):
        self.paths[(3, "candidate_probe")].unlink()
        with self.assertRaises(FileNotFoundError):
            preflight(self.args)

    def test_preflight_refuses_sequence_mismatch(self):
        seq = "A" + self.sequence[1:] if self.sequence[0] != "A" else "M" + self.sequence[1:]
        self.write_npz(self.paths[(2, "apo")], sequence=seq)
        with self.assertRaises(ValueError):
            preflight(self.args)

    def test_preflight_refuses_nonfinite_input(self):
        self.write_npz(self.paths[(2, "apo")], nonfinite=True)
        with self.assertRaises(ValueError):
            preflight(self.args)


if __name__ == "__main__":
    unittest.main()
