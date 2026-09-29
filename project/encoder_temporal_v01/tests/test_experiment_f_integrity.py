import json
import unittest
from pathlib import Path

from project.encoder_atom_readout_v01.experiment_d import sha256


ROOT = Path(__file__).resolve().parents[2] / "results/encoder_temporal_F_v01/run_001"
D_ROOT = Path(__file__).resolve().parents[2] / "results/encoder_atom_readout_D_v01/run_001"


@unittest.skipUnless((ROOT / "EXPERIMENT_F_COMPLETE.json").is_file(), "completed Experiment F artifacts required")
class CompletedAuditTests(unittest.TestCase):
    def test_source_and_cache_hashes_unchanged(self):
        extraction = json.loads((D_ROOT / "extraction_audit.json").read_text())
        self.assertTrue(all(sha256(Path(row["cache_path"])) == row["cache_sha256"] and sha256(Path(row["source_path"])) == row["source_sha256"] for row in extraction["rows"]))

    def test_static_continuity_matches_d(self):
        continuity = json.loads((ROOT / "static_continuity.json").read_text())
        self.assertTrue(continuity["passed"]); self.assertLessEqual(continuity["max_abs_metric_difference"], continuity["tolerance"])

    def test_r2_only_scaler_and_probe_counts(self):
        inventory = json.loads((ROOT / "pair_inventory.json").read_text())
        probe = json.loads((ROOT / "phi_dynamic_probe.json").read_text())
        expected = inventory["splits"]["R2"]["valid_samples"]["phi"]
        self.assertTrue(all(value["count"] == expected for value in probe["scalers"].values()))
        self.assertTrue(all(value["samples"] == expected for value in probe["models"].values()))

    def test_identical_masks_across_representations(self):
        for target in ("phi", "psi", "chi1"):
            probe = json.loads((ROOT / f"{target}_dynamic_probe.json").read_text())
            for split in ("R2", "R3"):
                counts = {value["samples"] for value in probe["metrics"].values() if value["split"] == split}
                self.assertEqual(len(counts), 1)

    def test_no_cross_window_pairs_and_exact_counts(self):
        inventory = json.loads((ROOT / "pair_inventory.json").read_text())
        self.assertEqual(inventory["cross_window_pairs"], 0)
        self.assertTrue(all(group["stored_frame_pairs"] == 99 for group in inventory["groups"]))

    def test_block_partition_and_completion(self):
        blocks = json.loads((ROOT / "block_averaging_audit.json").read_text())
        self.assertTrue(all(item["blocks_per_window"] == 5 and item["block_size"] == 20 and item["decomposition_relative_error"] < 1e-12 for item in blocks.values()))
        integrity = json.loads((ROOT / "integrity_completion_audit.json").read_text())
        self.assertEqual(integrity["status"], "EXPERIMENT_F_COMPLETE")
