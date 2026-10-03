"""Stage4 contracts, temporal compatibility and fail-closed application tests.

Synthetic arrays here exercise math only. No MD/topology files are fabricated.
"""
import ast
import csv
import importlib
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from project.pacer_fkg_v02 import stage4_prospective_common_v01 as c
from project.pacer_fkg_v02 import run_stage4_prospective_phase1_bs256_v01 as p1
from project.pacer_fkg_v02 import run_stage4_prospective_phase2a_frozen_apply_v01 as p2
from project.pacer_fkg_v02 import run_stage4_prospective_phase2b_graph_region_v01 as p3


class Stage4ContractTests(unittest.TestCase):
    def job(self, candidate="PACER0010", context="CP", replica=3):
        return next(j for j in c.expected_jobs() if (j.candidate, j.context, j.replica) == (candidate, context, replica))

    def progress(self, job, status="complete", ns=10.0):
        return {"system": job.system, "replica": job.replica, "seed": job.seed,
                "status": status, "completed_ns": ns, "target_ns": 10.0,
                "timestep_fs": 2.0, "restraint_k_kj_mol_nm2": 0.0}

    def test_exact_design(self):
        self.assertEqual(c.CANDIDATE_CLUSTERS, {"PACER0073": 0, "PACER0027": 4, "PACER0010": 9})
        jobs = c.expected_jobs()
        self.assertEqual(len(jobs), 36)
        self.assertEqual(len({j.system for j in jobs}), 12)
        self.assertEqual(len({j.key for j in jobs}), 36)
        self.assertEqual(c.SEEDS, {1: 27101, 2: 38201, 3: 49301})
        for candidate, cluster in c.CANDIDATE_CLUSTERS.items():
            self.assertEqual(c.context_map(candidate), {"A": f"cluster{cluster}__apo", "P": f"cluster{cluster}__probe_only",
                                                      "C": f"{candidate}__candidate_no_probe", "CP": f"{candidate}__candidate_probe"})
            subset = [j for j in jobs if j.candidate == candidate]
            self.assertEqual(len(subset), 12)
            for j in subset:
                self.assertEqual(j.cluster, cluster)
                self.assertEqual(j.seed, c.SEEDS[j.replica])

    def test_physical_time_contract(self):
        self.assertEqual((c.RAW_FRAMES, c.RAW_FRAME_SPACING_PS, c.TEMPORAL_STRIDE), (1000, 10, 5))
        self.assertEqual((c.ANALYSIS_FRAMES, c.ANALYSIS_FRAME_SPACING_PS), (200, 50))
        self.assertEqual((c.BLOCK_FRAMES, c.BLOCKS_PER_TRAJECTORY, c.BLOCK_DURATION_NS), (20, 10, 1.0))
        self.assertEqual(c.TEMPORAL_CONTRACT["raw_frame_ids"], list(range(0, 1000, 5)))

    def test_stride_retains_historical_lag_and_endpoint(self):
        frames = np.broadcast_to(np.arange(1000, dtype=np.float32)[:, None, None], (1000, 270, 256))
        branches = p2.construct(frames)
        self.assertEqual(branches["STATE_MOTION"].shape, (10, 270, 512))
        self.assertTrue(np.all(branches["SIGNED_DRIFT"] == 5))
        self.assertTrue(np.all(branches["STATE_MOTION"][..., 256:] == 5))
        self.assertEqual(branches["STATE_MOTION"][0, 0, 0], 47.5)

    def test_exact_stage_b_formula_compatibility(self):
        from project.pacer_fkg_v02.run_cm00734_stage_b_phase2a_frozen_apply_v01 import construct as historical
        values = np.sin(np.arange(1000, dtype=np.float32) / 11)
        frames = np.broadcast_to(values[:, None, None], (1000, 270, 256))
        current = p2.construct(frames)
        sampled = np.tile(frames[::5], (2, 1, 1))
        old = historical(sampled)
        for branch in c.BRANCHES:
            np.testing.assert_array_equal(current[branch], old[branch][:10])

    def test_no_cross_block_differences(self):
        values = (np.arange(1000) // 100).astype(np.float32)
        branches = p2.construct(np.broadcast_to(values[:, None, None], (1000, 270, 256)))
        self.assertTrue(np.all(branches["SIGNED_DRIFT"] == 0))
        self.assertTrue(np.all(branches["STATE_MOTION"][..., 256:] == 0))

    def test_bad_raw_shape_and_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            p2.construct(np.zeros((200, 270, 256), np.float32))
        with self.assertRaises(ValueError):
            p2.construct(np.broadcast_to(np.float32(np.nan), (1000, 270, 256)))

    def test_known_stale_progress_only_physical_complete(self):
        job = self.job()
        p = self.progress(job, "running", 9.9)
        self.assertIn("KNOWN_STALE", c.completion_acceptance(job, p, 5_000_000, 1000, True, True))
        for args in ((4_995_000, 1000, True, True), (5_000_000, 999, True, True),
                     (5_000_000, 1000, False, True), (5_000_000, 1000, True, False)):
            with self.subTest(args=args), self.assertRaises(RuntimeError):
                c.completion_acceptance(job, p, *args)

    def test_stale_not_generalised_to_other_jobs(self):
        for job in c.expected_jobs():
            if (job.candidate, job.context, job.replica) == ("PACER0010", "CP", 3):
                continue
            with self.subTest(key=job.key), self.assertRaises(RuntimeError):
                c.completion_acceptance(job, self.progress(job, "running", 9.9), 5_000_000, 1000, True, True)

    def test_wrong_seed_and_unexpected_bookkeeping_rejected(self):
        job = self.job()
        for field, value in (("seed", 27101), ("system", "cluster0__apo"), ("timestep_fs", 4),
                             ("target_ns", 20), ("completed_ns", 9.8)):
            p = self.progress(job)
            p[field] = value
            with self.subTest(field=field), self.assertRaises(RuntimeError):
                c.completion_acceptance(job, p, 5_000_000, 1000, True, True)

    def test_all_normal_jobs_accept(self):
        for job in c.expected_jobs():
            self.assertEqual(c.completion_acceptance(job, self.progress(job), 5_000_000, 1000, True, True),
                             "physical_complete_and_progress_complete")

    def test_unique_ordered_sequence_projection(self):
        self.assertEqual(c.sequence_projection("XABCDEY", "ABCDE"), (1, 2, 3, 4, 5))
        self.assertEqual(c.sequence_projection("ABCDE", "ABCDE"), (0, 1, 2, 3, 4))
        for source, target in (("AABCDE", "ABCDE"), ("ABCD", "ABCDE"), ("ABCED", "ABCDE")):
            with self.assertRaises(RuntimeError):
                c.sequence_projection(source, target)

    def test_actual_frozen_cluster_construct_projection(self):
        import MDAnalysis as mda
        mapping = c.FROZEN_ROOT / c.MAPPING_REL
        if not mapping.exists():
            self.skipTest("optional historical mapping assets not available")
        with mapping.open(encoding="utf-8-sig", newline="") as f:
            sequence = "".join(c.AA[r["resname"]] for r in csv.DictReader(f))
        for cluster in (0, 4, 9):
            p = c.REPO / f"project/results/m4_gamd_ensemble/receptors/cluster_{cluster:02d}_receptor.pdb"
            if not p.exists():
                self.skipTest("optional actual cluster PDB not available")
            u = mda.Universe(str(p))
            selection, receipt = c.topology_selection(u, sequence)
            self.assertEqual(receipt["construct_standard_residues"], 274)
            self.assertEqual(receipt["receptor_residue_count"], 270)
            self.assertEqual(len(receipt["excluded_construct_resindices"]), 4)
            self.assertEqual(receipt["sequence"], sequence)
            self.assertEqual(c.validate_selected_atoms(u, selection, sequence), 2139)

    def test_duplicate_or_missing_topology_rejected_without_md(self):
        with tempfile.TemporaryDirectory() as d:
            # Directory-only fixture; no molecular files are invented.
            with self.assertRaises(RuntimeError):
                c.discover_inputs(self.job(), Path(d))

    def stage4_topology_fixture(self):
        import MDAnalysis as mda
        mapping = c.FROZEN_ROOT / c.MAPPING_REL
        path = c.SOURCE_ROOT / "systems/cluster0__apo/minimized.pdb"
        if not mapping.exists() or not path.exists():
            self.skipTest("actual Stage4 topology / historical mapping unavailable")
        with mapping.open(encoding="utf-8-sig", newline="") as f:
            sequence = "".join(c.AA[r["resname"]] for r in csv.DictReader(f))
        return mda.Universe(str(path)), sequence

    def test_actual_stage4_all_twelve_split_capped_topologies(self):
        import MDAnalysis as mda
        _, sequence = self.stage4_topology_fixture()
        for system in sorted({j.system for j in c.expected_jobs()}):
            with self.subTest(system=system):
                path = c.discover_inputs(next(j for j in c.expected_jobs() if j.system == system))["topology"]
                before = c.sha256(path)
                # Match the adapter's existing parser-only CONECT exclusion.
                # Keep original atom/residue records, order and source bytes.
                parser_text = "".join(line for line in path.read_text().splitlines(keepends=True)
                                      if not line.startswith("CONECT"))
                u = mda.Universe(io.StringIO(parser_text), format="PDB")
                protein = u.select_atoms("protein")
                self.assertEqual(len(protein.residues), 278)
                self.assertEqual([len(protein.select_atoms(f"chainID {ch}").residues) for ch in ("A", "B")], [198, 80])
                selection, receipt = c.topology_selection(u, sequence)
                self.assertEqual(receipt["construct_standard_residues"], 274)
                self.assertEqual(receipt["receptor_residue_count"], 270)
                # Frozen order excludes the construct's loop-edge HIS and tail LLL.
                self.assertEqual(receipt["excluded_construct_resindices"], [196, 274, 275, 276])
                self.assertEqual([r.ix for r in u.select_atoms(selection).residues],
                                 list(range(1, 196)) + list(range(199, 274)))
                self.assertEqual([x["source_chain"] for x in receipt["mapping"]], ["A"] * 195 + ["B"] * 75)
                self.assertEqual(c.validate_selected_atoms(u, selection, sequence), 2139)
                self.assertEqual(c.sha256(path), before)

    def test_stage4_projection_independent_of_chain_and_segid_labels(self):
        u, sequence = self.stage4_topology_fixture()
        original, _ = c.topology_selection(u, sequence)
        # In-memory labels only: never alter the real PDB or fabricate MD files.
        u.atoms.chainIDs = ""
        u.segments.segids = "renamed"
        changed, receipt = c.topology_selection(u, sequence)
        self.assertEqual(changed, original)
        self.assertEqual({r["source_segid"] for r in receipt["mapping"]}, {"renamed"})

    def test_stage4_duplicate_receptor_projection_rejected(self):
        import MDAnalysis as mda
        u, sequence = self.stage4_topology_fixture()
        protein = u.select_atoms("protein")
        duplicate = mda.Merge(protein, protein)
        with self.assertRaisesRegex(RuntimeError, "missing or ambiguous"):
            c.topology_selection(duplicate, sequence)

    def test_stage4_reversed_fragments_and_changed_contract_rejected(self):
        import MDAnalysis as mda
        u, sequence = self.stage4_topology_fixture()
        reversed_fragments = mda.Merge(u.select_atoms("protein and chainID B"),
                                       u.select_atoms("protein and chainID A"))
        with self.assertRaisesRegex(RuntimeError, "missing or ambiguous"):
            c.topology_selection(reversed_fragments, sequence)
        with self.assertRaisesRegex(RuntimeError, "270-residue contract"):
            c.topology_selection(u, sequence[:-1])

    def test_incomplete_phase1_and_phase2_cache_pairs_fail_closed(self):
        from unittest.mock import MagicMock
        cache, manifest = MagicMock(), MagicMock()
        cache.exists.return_value = True
        cache.is_file.return_value = True
        manifest.exists.return_value = False
        manifest.is_file.return_value = False
        with patch.object(p1, "cache_paths", return_value=(cache, manifest)), self.assertRaises(RuntimeError):
            p1.validate_pair(self.job(), "full")
        dp, rp, mp = MagicMock(), MagicMock(), MagicMock()
        dp.exists.return_value = dp.is_file.return_value = True
        rp.exists.return_value = rp.is_file.return_value = False
        mp.exists.return_value = mp.is_file.return_value = False
        with patch.object(p2, "paths", return_value=(dp, rp, mp)), self.assertRaises(RuntimeError):
            p2.validate_pair(self.job(), "STATE_MOTION", "full")

    def test_valid_existing_cache_skips_inference(self):
        cached = {"cache": {"sha256": "already_authenticated"}}
        with patch.object(p1, "validate_pair", return_value=cached), patch.object(p1, "mapped_topology") as mapping:
            self.assertIs(p1.extract_one(self.job(), "full", None, None, "cuda", {}), cached)
            mapping.assert_not_called()

    def test_private_encoder_selection_does_not_change_historical_module(self):
        from project.pacer_fkg_v02 import run_phase1_bs256 as historical
        selection, resnames = historical.SELECTION, dict(historical.RESNAME_MAP)
        engine = p1.load_phase1_module()
        engine.SELECTION = "resindex 1 2 3"
        self.assertEqual(historical.SELECTION, selection)
        self.assertEqual(historical.RESNAME_MAP, resnames)
        self.assertEqual(engine.RESNAME_MAP["CYX"], "CYS")

    def test_historical_outputs_blocked(self):
        for relative in ("project/results/pacer_fkg_v02_longmd_v01/x.json",
                         "project/cache/pacer_fkg_v02_longmd_v01/x.npy",
                         "project/results/pacer_dc_close_loop_20ns_v01/x.json",
                         "project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01/x.json",
                         "project/pacer_fkg_v02/run_phase3_fkg.py"):
            with self.assertRaises(ValueError):
                c.guard_output(c.REPO / relative)
        self.assertEqual(c.guard_output(c.RESULT_ROOT / "test.json"), (c.RESULT_ROOT / "test.json").resolve())
        with self.assertRaises(ValueError):
            c.guard_output(c.RESULT_ROOT / ".." / "pacer_fkg_v02_longmd_v01" / "test.json")

    def test_missing_source_has_no_import_side_effects(self):
        with tempfile.TemporaryDirectory() as d:
            missing = Path(d) / "not_downloaded"
            availability = c.input_availability(missing)
            self.assertEqual(availability["file_sets_present"], 0)
            self.assertEqual(len(availability["jobs"]), 36)
            with self.assertRaises(RuntimeError):
                c.discover_inputs(self.job(), missing)
        for name in ("stage4_prospective_common_v01", "run_stage4_prospective_phase1_bs256_v01",
                     "run_stage4_prospective_phase2a_frozen_apply_v01", "run_stage4_prospective_phase2b_graph_region_v01"):
            self.assertIsNotNone(importlib.import_module("project.pacer_fkg_v02." + name))

    def test_cli_help_without_download(self):
        for module in (p1, p2, p3):
            result = subprocess.run([sys.executable, module.__file__, "--help"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            for flag in ("--smoke", "--run", "--verify", "--preflight", "--source-root"):
                self.assertIn(flag, result.stdout)

    def test_five_contrasts_match_formulas(self):
        from project.pacer_fkg_v02 import run_phase3_fkg as historical
        contexts = {"A": np.array([1], np.float32), "P": np.array([3], np.float32),
                    "C": np.array([5], np.float32), "CP": np.array([11], np.float32)}
        expected = {"Delta_PAM": 8, "Delta_AGO": 4, "Delta_INT": 4,
                    "Probe_effect": 2, "Probe_background_effect": 6}
        self.assertEqual(set(p3.all_contrasts(historical)), set(expected))
        for name, value in expected.items():
            self.assertEqual(p3.contrast(contexts, name, historical)[0], value)
        with self.assertRaises(ValueError):
            p3.contrast({"C": contexts["C"], "CP": contexts["CP"]}, "Delta_PAM", historical)

    def test_no_refit_call_or_random_seed_selection(self):
        forbidden_calls = {"train", "fit", "calibrate", "replay_state", "seed", "manual_seed",
                           "default_rng", "normal", "run_phase2_calibration"}
        for module in (c, p1, p2, p3):
            tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
            calls = {n.func.id if isinstance(n.func, ast.Name) else n.func.attr
                     for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, (ast.Name, ast.Attribute))}
            self.assertFalse(calls & forbidden_calls, calls & forbidden_calls)
        self.assertFalse(any(c.FORBIDDEN.values()))

    def test_anchor_hash_is_mandatory_and_private_binding(self):
        self.assertEqual(c.REQUIRED_FREEZE_SHA256, "b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd")
        with tempfile.TemporaryDirectory() as d, patch.object(c, "FROZEN_ROOT", Path(d)):
            with self.assertRaises(RuntimeError):
                c.verify_frozen_anchor()
        from project.pacer_fkg_v02 import run_phase3_fkg as historical
        old_root = historical.REPO
        engine = c.private_engine()
        self.assertEqual(historical.REPO, old_root)
        self.assertIsNot(engine, historical)

    def test_no_auto_labels_and_r2_is_evaluation(self):
        self.assertNotIn("CALIBRATION_REPLICA", p3.__dict__)
        result = p3.summarize(np.ones((10, 512), np.float32))
        self.assertEqual(result["replica_role"], "FROZEN_STATE_EVALUATION_ONLY")
        self.assertEqual(result["n_correlated_time_blocks"], 10)
        self.assertNotIn("PAM_probability", result)


if __name__ == "__main__":
    unittest.main()
