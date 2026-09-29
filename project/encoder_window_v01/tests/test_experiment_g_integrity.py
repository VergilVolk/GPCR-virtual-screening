import json
import unittest
from pathlib import Path

import numpy as np

from project.encoder_atom_readout_v01.experiment_d import sha256
from project.encoder_window_v01.experiment_g import block_branches, block_indices, derive_reps


PROJECT = Path(__file__).resolve().parents[2]
ROOT = PROJECT / "results/encoder_window_G_v01/run_001"
D_ROOT = PROJECT / "results/encoder_atom_readout_D_v01/run_001"


@unittest.skipUnless((ROOT / "EXPERIMENT_G_COMPLETE.json").is_file(), "completed Experiment G required")
class CompletedWindowAuditTests(unittest.TestCase):
    def test_predecessor_and_c1_hashes(self):
        provenance = json.loads((ROOT / "provenance.json").read_text())
        for key in ("experiment_d_completion", "experiment_d_integrity", "experiment_f_completion", "experiment_f_integrity", "experiment_f_report", "c1_implementation", "checkpoint"):
            self.assertEqual(sha256(Path(provenance[key]["path"])), provenance[key]["sha256"])

    def test_inventory_has_no_boundary_deltas(self):
        inventory = json.loads((ROOT / "block_inventory.json").read_text())
        self.assertEqual(inventory["blocks_total"], 200); self.assertEqual(inventory["cross_block_deltas"], 0); self.assertEqual(inventory["cross_window_deltas"], 0)
        self.assertTrue(all(row["frames"] == 20 and row["within_block_deltas"] == 19 for row in inventory["rows"]))

    def test_identity_reversal_and_static_continuity(self):
        self.assertTrue(json.loads((ROOT / "endpoint_identity_audit.json").read_text())["passed"])
        self.assertTrue(json.loads((ROOT / "reversal_audit.json").read_text())["passed"])
        self.assertTrue(json.loads((ROOT / "static_branch_continuity.json").read_text())["passed"])

    def test_real_cache_mc_exact_embedding(self):
        extraction = json.loads((D_ROOT / "extraction_audit.json").read_text()); path = Path(extraction["rows"][0]["cache_path"])
        with np.load(path) as archive: reps = derive_reps(archive["M"], archive["B"], archive["S"]); m = archive["M"]
        idx = block_indices()[0]
        for branch in ("static", "signed", "rms"):
            observed, expected = block_branches(reps["MC"][idx])[branch], block_branches(m[idx])[branch]
            np.testing.assert_array_equal(observed[..., :128], expected); self.assertEqual(np.count_nonzero(observed[..., 128:]), 0)

    def test_r2_only_fit_counts_and_identical_masks(self):
        inventory = json.loads((ROOT / "block_inventory.json").read_text()); cross = json.loads((ROOT / "cross_branch_probe.json").read_text())
        for key, scaler in cross["scalers"].items(): self.assertEqual(scaler["count"], inventory["valid_samples"]["R2"][key.split("|")[2]])
        for key, model in cross["models"].items(): self.assertEqual(model["samples"], inventory["valid_samples"]["R2"][key.split("|")[2]])
        for target in ("static_branch", "signed_dynamic_branch", "dynamic_rms_branch"):
            metrics = json.loads((ROOT / f"{target}.json").read_text())["metrics"]
            for split in ("R2", "R3"):
                by_torsion = {}
                for value in metrics.values():
                    if value["split"] == split: by_torsion.setdefault(value["torsion"], set()).add(value["samples"])
                self.assertTrue(all(len(counts) == 1 for counts in by_torsion.values()))

    def test_sources_caches_and_completion(self):
        extraction = json.loads((D_ROOT / "extraction_audit.json").read_text())
        self.assertTrue(all(sha256(Path(row["cache_path"])) == row["cache_sha256"] and sha256(Path(row["source_path"])) == row["source_sha256"] for row in extraction["rows"]))
        marker = json.loads((ROOT / "EXPERIMENT_G_COMPLETE.json").read_text()); integrity = json.loads((ROOT / "integrity_completion_audit.json").read_text())
        self.assertEqual(marker["status"], "EXPERIMENT_G_COMPLETE"); self.assertEqual(integrity["status"], "EXPERIMENT_G_COMPLETE")
