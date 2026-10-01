"""Replay the frozen Phase-1 encoder on existing trajectories and compare.

This is verification/replay only: it runs the frozen C1-BS256 extractor on the
first frames of an EXISTING CM00734 trajectory and compares the output with the
frozen Phase-1 cache.  No new MD, no new molecule, no calibration.

    python project/integration_freeze/tools/verify_phase1_replay.py [--device cuda]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from project.pacer_fkg_v02.frozen_runner import engines  # noqa: E402
from project.pacer_fkg_v02.frozen_runner.fkg_phase1 import (  # noqa: E402
    DEFAULT_CHECKPOINT_RELPATH, CandidateJob, frozen_sequence,
)
from project.pacer_fkg_v02.frozen_runner.io_utils import load_npy_copy  # noqa: E402

ANALYSIS = "project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01"
TOPOLOGY_ROOT = f"{ANALYSIS}/topology_sanitized"
PRODUCTION_ROOT = "project/results/pacer_dc_cm00734_stage_b_20ns_v01/production"
CACHE_ROOT = "project/cache/pacer_dc_cm00734_stage_b_20ns_v01/track_b_fkg_v02/bs256"
SYSTEMS = ("CM00734__candidate_no_probe", "CM00734__candidate_probe")
REPLICAS = (1, 2, 3)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--max-jobs", type=int, default=1, help="how many (system, replica) pairs to replay")
    args = parser.parse_args()

    phase1 = engines.load_phase1_module()
    engines.assert_frozen_numerics()
    sequence = frozen_sequence(REPO, phase1)
    checkpoint = REPO / DEFAULT_CHECKPOINT_RELPATH
    geom2vec_source = Path(r"C:\projects\geom2vec-source")
    model, model_record = phase1.load_model(checkpoint, args.device, geom2vec_source)

    frame_ids = list(phase1.REPLAY_FRAME_IDS)
    report = {
        "schema": "pacer.final.frozen_runner.phase1_replay.v1",
        "device": args.device,
        "replay_frame_ids": frame_ids,
        "model_provenance": {
            "checkpoint_sha256": model_record["checkpoint"]["actual_sha256"],
            "strict_load_ok": model_record["checkpoint"].get("strict_load_ok"),
            "geom2vec_source_equivalent": model_record["geom2vec_source"].get("source_equivalent"),
        },
        "jobs": [],
    }

    done = 0
    for system in SYSTEMS:
        if done >= args.max_jobs:
            break
        for replica in REPLICAS:
            if done >= args.max_jobs:
                break
            topology = REPO / TOPOLOGY_ROOT / f"{system}_minimized_no_conect.pdb"
            trajectory = REPO / PRODUCTION_ROOT / system / f"replica_{replica:02d}" / "trajectory.dcd"
            cache = REPO / CACHE_ROOT / system / f"replica_{replica:02d}" / "C1_BS256.npy"
            if not (topology.is_file() and trajectory.is_file() and cache.is_file()):
                report["jobs"].append({"system": system, "replica": replica, "status": "SKIPPED_MISSING_INPUT"})
                continue
            job = CandidateJob(system, replica, topology, trajectory, "", "", sequence, 400)
            universe, atom_indices, atom_z, residue_index, bb, sc, _ = phase1.topology_mapping(job)
            try:
                replayed = phase1.infer_bs256(model, universe, atom_indices, atom_z, residue_index, bb, sc,
                                              frame_ids, args.device)
            finally:
                universe.trajectory.close()
            frozen = load_npy_copy(cache)[: len(frame_ids)]
            metrics = phase1.replay_metrics(frozen, replayed)
            report["jobs"].append({
                "system": system, "replica": replica, "status": "REPLAYED",
                "frozen_cache": str(cache), "shape": list(replayed.shape),
                "bitwise_identical": metrics["bitwise_identical"],
                "max_abs": metrics["max_abs"], "relative_rmse": metrics["relative_rmse"],
                "limit": metrics["limit"], "passed": metrics["passed"],
            })
            done += 1

    replayed_jobs = [j for j in report["jobs"] if j.get("status") == "REPLAYED"]
    report["n_replayed"] = len(replayed_jobs)
    report["all_passed"] = bool(replayed_jobs) and all(j["passed"] for j in replayed_jobs)
    report["all_bitwise_identical"] = bool(replayed_jobs) and all(j["bitwise_identical"] for j in replayed_jobs)

    out = REPO / "project/integration_freeze/receipts/PHASE1_ENCODER_REPLAY_v01.json"
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
