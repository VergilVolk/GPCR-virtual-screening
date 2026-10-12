"""Synthetic segment tests and optional read-only real-data tests.

Synthetic trajectories are explicitly fixtures, never production evidence.
Set PACER_STAGE4_PBCFIX01_TEST_ROOT to enable the real two-segment test.
"""
import copy
import csv
import json
import os
import subprocess
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from project.pacer_fkg_v02 import stage4_pbcfix01_common_v01 as c
from project.pacer_fkg_v02 import run_stage4_pbcfix01_phase1_bs256_v01 as p1
from project.pacer_fkg_v02 import run_stage4_pbcfix01_phase2a_frozen_apply_v01 as p2
from project.pacer_fkg_v02 import run_stage4_pbcfix01_phase2b_graph_region_v01 as p3

BUNDLE = c.REPO / f"project/results/{c.PROTOCOL}/server_bundle/{c.PROTOCOL}"


class Contracts(unittest.TestCase):
    def test_actual_36_job_package(self):
        self.assertEqual(len(c.validate_matrix(BUNDLE)), 36)
        self.assertGreater(len(c.verify_sums(BUNDLE)), 5)
        jobs = c.expected_jobs()
        self.assertEqual(len({j.key for j in jobs}), 36)
        self.assertEqual(len({j.system for j in jobs}), 12)
        self.assertEqual(len(c.SEGMENT_ENDS), 100)
        self.assertEqual(c.SEGMENT_ENDS[-1], 5_000_000)
        for candidate in c.CANDIDATE_CLUSTERS:
            self.assertEqual({j.context for j in jobs if j.candidate == candidate}, set(c.ALIASES))
            self.assertEqual({j.replica for j in jobs if j.candidate == candidate}, {1, 2, 3})

    def test_all_matrix_identity_fields_fail_closed(self):
        rows = c.csv_rows(BUNDLE / "manifests/STAGE4_RERUN_JOB_MATRIX.csv")
        systems = c.csv_rows(BUNDLE / "manifests/STAGE4_RERUN_SYSTEM_MANIFEST.csv")
        for field in c.job_row(c.expected_jobs()[0]):
            changed = copy.deepcopy(rows)
            changed[0][field] = "WRONG"
            with self.subTest(field=field), patch.object(c, "csv_rows", side_effect=[changed, systems]):
                with self.assertRaises(RuntimeError):
                    c.validate_matrix(BUNDLE)
        for bad in (rows[:-1], rows + rows[:1], [rows[0]] * 36):
            with patch.object(c, "csv_rows", return_value=bad), self.assertRaises(RuntimeError):
                c.validate_matrix(BUNDLE)

    def test_server_absolute_paths_rebase_without_rewriting_identity(self):
        rows = c.csv_rows(BUNDLE / "manifests/STAGE4_RERUN_JOB_MATRIX.csv")
        systems = c.csv_rows(BUNDLE / "manifests/STAGE4_RERUN_SYSTEM_MANIFEST.csv")
        for row in rows:
            for key in ("output_path", "checkpoint_path"):
                row[key] = c.SERVER_ROOT + "/" + row[key]
        for row in systems:
            row["system_path"] = c.SERVER_ROOT + "/" + row["system_path"]
        original = copy.deepcopy(rows)
        with patch.object(c, "csv_rows", side_effect=[rows, systems]):
            verified = c.validate_matrix(BUNDLE)
        self.assertEqual(list(verified.values()), original)
        for path in ("/data/other/production/job", c.SERVER_ROOT + "/../escape", "C:/old/path"):
            with self.subTest(path=path), self.assertRaises(RuntimeError):
                c.campaign_relative(path)

    def test_analysis_lock_blocks_other_process_and_releases(self):
        with c.analysis_lock():
            with c.analysis_lock():
                pass
            code = "from project.pacer_fkg_v02 import stage4_pbcfix01_common_v01 as c; c.analysis_lock().__enter__()"
            result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("another PBCFIX01", result.stderr)
        with c.analysis_lock():
            pass

    def test_recovery_preserves_orphans_and_protects_history(self):
        c.CACHE_ROOT.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=c.CACHE_ROOT) as temporary:
            base = Path(temporary)
            orphan = base / "orphan.npy"
            orphan.write_bytes(b"failed extraction evidence")
            with self.assertRaises(ValueError):
                c.quarantine_outputs([orphan, c.REPO / "project/results/pacer_fkg_v02_longmd_v01/protected"])
            self.assertTrue(orphan.is_file())
            result = c.quarantine_outputs([orphan])
            preserved = Path(result["moves"][0]["preserved"])
            self.assertFalse(orphan.exists())
            self.assertEqual(preserved.read_bytes(), b"failed extraction evidence")
            token = preserved.relative_to(c.CACHE_ROOT).parts[1]
            for root in (c.RESULT_ROOT, c.CACHE_ROOT):
                target = root / "_recovery" / token
                self.assertTrue(target.resolve().is_relative_to(root.resolve()))
                shutil.rmtree(target)

    def test_phase1_recovery_finds_orphan_and_invalidates_downstream(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            job = c.expected_jobs()[0]
            (base / "full").mkdir()
            orphan = base / "full/C1_BS256.npy"
            orphan.write_bytes(b"unpublished")
            def paths(job, mode):
                return base / mode / "C1_BS256.npy", base / mode / "manifest.json"
            with patch.object(p1, "jobs_from_receipt", return_value=[job]), patch.object(p1, "cache_paths", side_effect=paths), \
                 patch.object(c, "quarantine_outputs", return_value={"moves": []}) as quarantine:
                result = p1.recover(base)
            self.assertEqual(len(result["failed_jobs"]), 1)
            targets = quarantine.call_args.args[0]
            self.assertIn(orphan, targets)
            self.assertIn(c.RESULT_ROOT / "phase2a_frozen_apply", targets)
            self.assertIn(c.RESULT_ROOT / "phase2b_graph_region", targets)

    def test_frozen_assets_if_present(self):
        if not (c.FROZEN_ROOT / c.FREEZE_REL).is_file():
            self.skipTest("frozen assets not installed")
        _, anchor = c.verify_frozen_anchor()
        self.assertEqual(anchor["receipt"]["freeze"]["sha256"], c.REQUIRED_FREEZE_SHA256)

    def test_actual_oriented_cluster01_baseline(self):
        import MDAnalysis as mda
        top = c.REPO / "project/results/pacer_stage4_rerun_v01/server_bundle/pacer_stage4_rerun_v01/inputs/receptors/oriented_cluster_01.pdb"
        if not top.is_file() or not (c.FROZEN_ROOT / c.MAPPING_REL).is_file():
            self.skipTest("oriented baseline/frozen mapping not installed")
        u = mda.Universe(str(top))
        try:
            sequence = c.frozen_sequence()
            selection, mapping = c.topology_selection(u, sequence)
            self.assertEqual(mapping["receptor_residue_count"], 270)
            self.assertEqual(c.validate_selected_atoms(u, selection, sequence), 2139)
        finally:
            u.trajectory.close()

    def test_history_and_path_guards(self):
        from project.pacer_fkg_v02 import stage4_prospective_common_v01 as old
        self.assertNotEqual(c.RESULT_ROOT, old.RESULT_ROOT)
        self.assertNotEqual(c.CACHE_ROOT, old.CACHE_ROOT)
        self.assertEqual(old.CANDIDATE_CLUSTERS, {"PACER0073": 0, "PACER0027": 4, "PACER0010": 9})
        with self.assertRaises(ValueError):
            c.guard_output(old.RESULT_ROOT / "anything.json")
        for name in ("../outside", "C:/outside", "/outside", "a\\b", ""):
            with self.subTest(name=name), self.assertRaises(RuntimeError):
                c.member(BUNDLE, name)

    def test_all_cli_help_without_download(self):
        for name in ("inputs", "phase1_bs256", "phase2a_frozen_apply", "phase2b_graph_region"):
            module = f"project.pacer_fkg_v02.run_stage4_pbcfix01_{name}_v01"
            result = subprocess.run([sys.executable, "-m", module, "--help"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--source-root", result.stdout)

    def test_five_matched_contrasts(self):
        engine = c.private_engine()
        contexts = {"A": np.array([1.]), "P": np.array([3.]), "C": np.array([7.]), "CP": np.array([15.])}
        from project.pacer_fkg_v02 import run_stage4_prospective_phase2b_graph_region_v01 as old
        self.assertEqual(p3.all_contrasts(engine), old.all_contrasts(engine))
        for name in p3.all_contrasts(engine):
            np.testing.assert_array_equal(p3.contrast(contexts, name, engine), old.contrast(contexts, name, engine))

    def test_exact_numeric_phase_reuse(self):
        from project.pacer_fkg_v02 import run_stage4_prospective_phase2a_frozen_apply_v01 as old
        x = np.broadcast_to(np.sin(np.arange(1000, dtype=np.float32))[:, None, None], (1000, 270, 256))
        observed, expected = p2.construct(x), old.construct(x)
        for branch in c.BRANCHES:
            np.testing.assert_array_equal(observed[branch], expected[branch])
        with self.assertRaises(ValueError):
            p3.contrast({"A": x}, "Probe_effect", None)


class SyntheticSegments(unittest.TestCase):
    def setUp(self):
        import MDAnalysis as mda
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.job = c.expected_jobs()[0]
        self.identity = {"job": c.job_row(self.job), "synthetic_fixture": True}
        nres = 270
        u = mda.Universe.empty(nres * 5, n_residues=nres,
                               atom_resindex=np.repeat(np.arange(nres), 5), trajectory=True)
        u.add_TopologyAttr("names", ["N", "CA", "C", "O", "CB"] * nres)
        u.add_TopologyAttr("resnames", ["ALA"] * nres)
        u.add_TopologyAttr("resids", np.arange(1, nres + 1))
        u.add_TopologyAttr("chainIDs", ["A"] * (nres * 5))
        u.dimensions = [60, 60, 60, 90, 90, 90]
        u.atoms.positions = np.zeros((nres * 5, 3))
        top = self.root / "minimized.pdb"
        u.atoms.write(str(top))
        progress = self.root / "progress.json"
        progress.write_text(json.dumps({"identity": self.identity, "status": "RUNNING", "completed_steps": 100000, "completed_ns": .2}))
        dirs = []
        for index, end in enumerate(c.SEGMENT_ENDS[:2]):
            d = self.root / f"{end:08d}"
            d.mkdir()
            dirs.append(d)
            with mda.coordinates.DCD.DCDWriter(str(d / "trajectory.dcd"), n_atoms=len(u.atoms), dt=10, istart=1) as writer:
                for frame in range(10):
                    u.atoms.positions = np.full((len(u.atoms), 3), index * 10 + frame, dtype=np.float32)
                    writer.write(u)
            fields = ["Step", "Time (ps)", "Potential Energy (kJ/mole)", "Kinetic Energy (kJ/mole)", "Temperature (K)", "Box Volume (nm^3)"]
            with (d / "state.csv").open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(fields)
                for step in range(end - 45000, end + 1, 5000):
                    writer.writerow([step, step * .002, -1, 1, 300, 216])
            (d / "checkpoint.chk").write_bytes(b"SYNTHETIC_NOT_AN_OPENMM_CHECKPOINT")
            self.rewrite_receipt(d, index)
        self.paths = {"topology": top, "progress": progress, "segments": dirs,
                      "trajectory": tuple(d / "trajectory.dcd" for d in dirs)}
        self.patch = patch.object(c, "sanitize_topology", side_effect=lambda source, system: source)
        self.patch.start()
        self.addCleanup(self.patch.stop)

    def rewrite_receipt(self, d, index):
        receipt = {"identity": self.identity, "start_step": index * 50000, "end_step": (index + 1) * 50000,
                   "outputs": {name: c.sha256(d / name) for name in ("trajectory.dcd", "state.csv", "checkpoint.chk")}}
        (d / "receipt.json").write_text(json.dumps(receipt))

    def test_chain_coordinates_absolute_time_and_frozen_mapping(self):
        rec = c.inspect_segments(self.job, self.paths, self.identity, "A" * 270)
        self.assertEqual(rec["dcd_frames"], 20)
        self.assertEqual(rec["receptor"]["atom14_heavy_atoms"], 1350)
        c.verify_input_record(rec)
        from dataclasses import replace
        job = replace(self.job, topology=self.paths["topology"], trajectory=self.paths["trajectory"],
                      sequence="A" * 270, selection=rec["selection"])
        result = p1.mapped_topology(p1.load_phase1_module(), job)
        u = result[0]
        try:
            for frame in (0, 9, 10, 19, 5, 10, 15):
                ts = u.trajectory[frame]
                self.assertAlmostEqual(ts.time, (frame + 1) * 10)
                np.testing.assert_allclose(ts.positions, frame)
            xyz = p1.load_phase1_module().load_frame_batch(u, result[1], [9, 10])
            self.assertEqual(xyz.shape, (2, 1350, 3))
        finally:
            u.trajectory.close()

    def test_corrupt_checkpoint_and_changed_source_rejected(self):
        rec = c.inspect_segments(self.job, self.paths, self.identity, "A" * 270)
        (self.paths["segments"][0] / "checkpoint.chk").write_bytes(b"changed")
        with self.assertRaises(RuntimeError):
            c.verify_input_record(rec)
        with self.assertRaises(RuntimeError):
            c.inspect_segments(self.job, self.paths, self.identity, "A" * 270)

    def test_100_segments_1000_frames_without_merged_dcd(self):
        import MDAnalysis as mda
        u = mda.Universe(str(self.paths["topology"]))
        try:
            # Repeated small synthetic files isolate ChainReader indexing; this
            # does not claim 100 independent production segments were verified.
            c.attach_segments(u, self.paths["trajectory"] * 50)
            self.assertEqual(len(u.trajectory), 1000)
            for frame in (0, 9, 10, 19, 20, 995, 999):
                self.assertAlmostEqual(u.trajectory[frame].time, (frame + 1) * 10)
            self.assertAlmostEqual(u.trajectory[999].time, 10000)
            self.assertEqual(len(list(self.root.rglob("*.dcd"))), 2)
        finally:
            u.trajectory.close()

    def test_wrong_dcd_sampling_and_count_rejected_after_valid_hash(self):
        import MDAnalysis as mda
        u = mda.Universe(str(self.paths["topology"]))
        d = self.paths["segments"][0]
        try:
            for frames, dt in ((9, 10), (10, 5)):
                with mda.coordinates.DCD.DCDWriter(str(d / "trajectory.dcd"), n_atoms=len(u.atoms), dt=dt, istart=1) as writer:
                    for frame in range(frames):
                        writer.write(u)
                self.rewrite_receipt(d, 0)
                with self.subTest(frames=frames, dt=dt), self.assertRaisesRegex(RuntimeError, "frame/dt mismatch"):
                    c.inspect_segments(self.job, self.paths, self.identity, "A" * 270)
        finally:
            u.trajectory.close()

    def test_time_gap_even_with_valid_sha_rejected(self):
        d = self.paths["segments"][1]
        path = d / "state.csv"
        path.write_text(path.read_text().replace("55000,110.0", "55000,10.0"))
        self.rewrite_receipt(d, 1)
        with self.assertRaisesRegex(RuntimeError, "time mismatch"):
            c.inspect_segments(self.job, self.paths, self.identity, "A" * 270)

    def test_box_volume_mismatch_with_valid_hash_rejected(self):
        d = self.paths["segments"][0]
        path = d / "state.csv"
        path.write_text(path.read_text().replace(",216", ",432"))
        self.rewrite_receipt(d, 0)
        with self.assertRaisesRegex(RuntimeError, "box volume mismatch"):
            c.inspect_segments(self.job, self.paths, self.identity, "A" * 270)

    def test_atom_count_mismatch_rejected(self):
        import MDAnalysis as mda
        d = self.paths["segments"][0]
        u = mda.Universe.empty(1, trajectory=True)
        u.dimensions = [60, 60, 60, 90, 90, 90]
        with mda.coordinates.DCD.DCDWriter(str(d / "trajectory.dcd"), n_atoms=1, dt=10, istart=1) as writer:
            for _ in range(10):
                writer.write(u)
        self.rewrite_receipt(d, 0)
        with self.assertRaises((ValueError, RuntimeError)):
            c.inspect_segments(self.job, self.paths, self.identity, "A" * 270)

    def test_receipt_gap_and_cross_job_rejected(self):
        d = self.paths["segments"][1]
        path = d / "receipt.json"
        original = c.load_json(path)
        for key, value in (("start_step", 0), ("end_step", 50000), ("identity", {"wrong_candidate": True})):
            changed = dict(original, **{key: value})
            path.write_text(json.dumps(changed))
            with self.subTest(key=key), self.assertRaises(RuntimeError):
                c.inspect_segments(self.job, self.paths, self.identity, "A" * 270)

    def test_incompatible_sequence_and_atom14_rejected(self):
        import MDAnalysis as mda
        with self.assertRaises(RuntimeError):
            c.inspect_segments(self.job, self.paths, self.identity, "G" * 270)
        u = mda.Universe(str(self.paths["topology"]))
        self.addCleanup(u.trajectory.close)
        u.atoms[0].name = "MISSING_N"
        with self.assertRaisesRegex(RuntimeError, "missing atom14"):
            c.validate_selected_atoms(u, "all", "A" * 270)
        with self.assertRaises(RuntimeError):
            c.sequence_projection("A" * 271, "A" * 270)


class ProductionReceipts(unittest.TestCase):
    def test_preparation_and_authorization_hash_chain(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            job = c.expected_jobs()[0]
            sid = job.system
            row = c.validate_matrix(BUNDLE)[c.job_row(job)["job_id"]]
            def write(relative, obj):
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(obj))
                return c.sha256(path)
            def sha(relative):
                return c.sha256(root / relative)
            build_outputs = {}
            for name in ("minimized.pdb", "system.xml", "state.xml"):
                build_outputs[name] = write(f"systems/{sid}/{name}", {"synthetic": name})
            build_sha = write(f"systems/{sid}/build_audit.json", {"protocol_id": c.PROTOCOL, "status": "COMPLETE", "outputs": build_outputs})
            config_sha = write(f"config/{c.PROTOCOL}.json", {"synthetic": True})
            package_sha = write("SHA256SUMS_PBCFIX01.txt", {"synthetic": True})
            snapshot_sha = write("evidence/source_snapshot.json", {"synthetic": True})
            smoke_sha = write(f"smoke/{sid}/audit.json", {"synthetic": True})
            gate_sha = write(f"logs/server_preflight_{c.PROTOCOL}.json", {"status": "PASS", "config_sha256": config_sha,
                              "builds": {sid: build_sha}, "smokes": {sid: smoke_sha}})
            approval_sha = write("authorization.json", {"protocol_id": c.PROTOCOL, "preflight_sha256": gate_sha,
                                 "scientific_review_approved": True, "user_explicit_authorization": True})
            authorization = dict(protocol_id=c.PROTOCOL, system_id=sid, build_audit_sha256=build_sha,
                                 config_sha256=config_sha, package_manifest_sha256=package_sha,
                                 system_sha256=build_outputs["system.xml"], topology_sha256=build_outputs["minimized.pdb"],
                                 state_sha256=build_outputs["state.xml"], source_snapshot_sha256=snapshot_sha,
                                 preflight_sha256=gate_sha, authorization_sha256=approval_sha)
            expected = dict(authorization_identity=authorization, job=row, config_sha256=config_sha,
                            system_sha256=build_outputs["system.xml"], topology_sha256=build_outputs["minimized.pdb"],
                            smoke_sha256=smoke_sha)
            prep_identity = dict(expected, initial_state_sha256=build_outputs["state.xml"])
            previous = build_outputs["state.xml"]
            stage_receipts = []
            base = f"equilibration/{sid}/replica_01"
            for label in ("nvt", "npt", "release_01", "release_02", "release_03", "release_04"):
                outputs = {name: write(f"{base}/{label}/{name}", {"synthetic": label + name})
                           for name in ("state.xml", "state.csv", "trajectory.dcd")}
                stage_receipts.append(write(f"{base}/{label}/receipt.json", {"passed": True,
                                      "identity": dict(prep_identity, stage=label, previous_state_sha256=previous), "outputs": outputs}))
                previous = outputs["state.xml"]
            ready_sha = write(f"{base}/ready_{c.PROTOCOL}.json", {"passed": True, "identity": prep_identity,
                              "stage_receipts": stage_receipts, "start_state": f"{base}/release_04/state.xml",
                              "start_state_sha256": previous})
            expected.update(start_state_sha256=previous, preparation_receipt_sha256=ready_sha)
            observed, assets = c.production_identity(job, root, row)
            self.assertEqual(observed, expected)
            for record in assets:
                c.verify_artifact(record)
            write(f"{base}/release_02/state.xml", {"changed": True})
            with self.assertRaisesRegex(RuntimeError, "SHA256 mismatch"):
                c.production_identity(job, root, row)


class OptionalRealSegments(unittest.TestCase):
    def test_real_two_committed_segments(self):
        source = os.environ.get("PACER_STAGE4_PBCFIX01_TEST_ROOT")
        if not source:
            self.skipTest("new production segments not available; set PACER_STAGE4_PBCFIX01_TEST_ROOT")
        from project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 import sample
        job_id = os.environ.get("PACER_STAGE4_PBCFIX01_TEST_JOB", c.job_row(c.expected_jobs()[0])["job_id"])
        result = sample(Path(source), job_id, 2)
        self.assertEqual(result["frames"], 20)
        self.assertEqual(result["status"], "SAMPLE_ONLY_NOT_FULL_ACCEPTANCE")


if __name__ == "__main__":
    unittest.main()
