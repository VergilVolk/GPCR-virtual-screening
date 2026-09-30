#!/usr/bin/env python
"""Freeze PACER-DC 20 ns closure Track-B Phase-1 BS256 artifacts."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[2]

CLOSURE = REPO / "project/results/pacer_dc_close_loop_20ns_v01"

CACHE_ROOT = (
    REPO
    / "project/cache/pacer_dc_close_loop_20ns_v01"
    / "track_b_fkg_v02/bs256/full"
)

REPORT_ROOT = (
    CLOSURE
    / "track_b_fkg_v02/phase1_bs256"
)

RUN_FULL = REPORT_ROOT / "RUN_FULL_v01.json"

RUNNER = (
    REPO
    / "project/pacer_fkg_v02"
    / "run_close_loop_phase1_bs256_v01.py"
)

OUTPUT = REPORT_ROOT / "PHASE1_FREEZE_RECEIPT_v01.json"

EXEC_MANIFEST = CLOSURE / "CLOSE_LOOP_EXECUTION_MANIFEST_v01.json"
PREFLIGHT = CLOSURE / "LY2119620_CLOSE_LOOP_PREFLIGHT_v01.json"
TRANSFER = CLOSURE / "TRANSFER_SHA256_v01.txt"
TOPOLOGY_AUDIT = (
    CLOSURE
    / "topology_sanitized"
    / "LY_TOPOLOGY_SANITIZATION_v01.json"
)

EXPECTED_EXTRACTOR_COMMIT = (
    "600a8833d0574566a72ec397650cc98c6dbf0a83"
)

EXPECTED_TRANSFER_SHA = (
    "d0383b9fda1bf6de8d227b4bf26a623c4621218292b602812def9aed41799696"
)

EXPECTED_EXEC_SHA = (
    "d0333145061bc5ff6a4b52df33269cf1893c77472e2c839252c03c84e6a39b77"
)

EXPECTED_PREFLIGHT_SHA = (
    "6b5d6c33ee7405c152ac0890c259da425091927fd87e02d311d9457c20313b3b"
)

SYSTEMS = (
    "LY2119620__candidate_no_probe",
    "LY2119620__candidate_probe",
)

REPLICAS = (1, 2, 3)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(8 * 1024 * 1024),
            b"",
        ):
            h.update(block)
    return h.hexdigest()


def load_json(path: Path):
    return json.loads(
        path.read_text(encoding="utf-8-sig")
    )


def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args],
        cwd=REPO,
        text=True,
    ).strip()


def artifact(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)

    return {
        "path": path.relative_to(REPO).as_posix(),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


def main() -> int:
    if OUTPUT.exists():
        raise FileExistsError(
            f"refusing overwrite: {OUTPUT}"
        )

    # Closure anchors.
    if sha256(EXEC_MANIFEST) != EXPECTED_EXEC_SHA:
        raise RuntimeError(
            "execution manifest SHA mismatch"
        )

    if sha256(PREFLIGHT) != EXPECTED_PREFLIGHT_SHA:
        raise RuntimeError(
            "preflight SHA mismatch"
        )

    if sha256(TRANSFER) != EXPECTED_TRANSFER_SHA:
        raise RuntimeError(
            "transfer manifest SHA mismatch"
        )

    # Code anchor of the extractor itself, independent of later commits.
    runner_rel = RUNNER.relative_to(REPO).as_posix()

    extractor_commit = git(
        "log",
        "-1",
        "--format=%H",
        "--",
        runner_rel,
    )

    if extractor_commit != EXPECTED_EXTRACTOR_COMMIT:
        raise RuntimeError(
            "extractor Git commit changed"
        )

    run = load_json(RUN_FULL)

    if run.get("status") != "PASS":
        raise RuntimeError(
            "RUN_FULL is not PASS"
        )

    if len(run.get("jobs", [])) != 6:
        raise RuntimeError(
            "RUN_FULL does not contain six jobs"
        )

    run_jobs = {
        item["job"]: item
        for item in run["jobs"]
    }

    records = []

    for system in SYSTEMS:
        for replica in REPLICAS:
            key = (
                f"{system}__replica_{replica:02d}"
            )

            cache = (
                CACHE_ROOT
                / system
                / f"replica_{replica:02d}"
                / "C1_BS256.npy"
            )

            manifest = (
                REPORT_ROOT
                / "manifests"
                / "full"
                / f"{key}.json"
            )

            if not cache.is_file():
                raise FileNotFoundError(cache)

            if not manifest.is_file():
                raise FileNotFoundError(manifest)

            arr = np.load(
                cache,
                mmap_mode="r",
                allow_pickle=False,
            )

            if arr.shape != (400, 270, 256):
                raise RuntimeError(
                    f"cache shape mismatch: {key}"
                )

            if arr.dtype != np.float32:
                raise RuntimeError(
                    f"cache dtype mismatch: {key}"
                )

            if not np.isfinite(arr).all():
                raise RuntimeError(
                    f"non-finite cache: {key}"
                )

            cache_sha = sha256(cache)
            item = load_json(manifest)

            if item["cache"]["sha256"] != cache_sha:
                raise RuntimeError(
                    f"manifest/cache SHA mismatch: {key}"
                )

            replay = item["replay"]

            if replay.get("passed") is not True:
                raise RuntimeError(
                    f"replay not PASS: {key}"
                )

            if float(
                replay["relative_rmse"]
            ) >= 1e-4:
                raise RuntimeError(
                    f"replay RMSE exceeds frozen limit: {key}"
                )

            run_item = run_jobs.get(key)

            if run_item is None:
                raise RuntimeError(
                    f"job absent from RUN_FULL: {key}"
                )

            if run_item["cache_sha256"] != cache_sha:
                raise RuntimeError(
                    f"RUN_FULL/cache SHA mismatch: {key}"
                )

            records.append(
                {
                    "system": system,
                    "replica": replica,
                    "shape": [400, 270, 256],
                    "dtype": "float32",
                    "finite": True,
                    "cache": artifact(cache),
                    "manifest": artifact(manifest),
                    "replay": {
                        "passed": True,
                        "frame_ids":
                            replay["frame_ids"],
                        "relative_rmse":
                            replay["relative_rmse"],
                        "limit": 1e-4,
                    },
                }
            )

    receipt = {
        "schema":
            "pacer_dc.close_loop_20ns.track_b.phase1_freeze.v1",

        "status":
            "PHASE1_FROZEN",

        "track":
            "B / frozen C1-BS256 + PACER-FKG-v02",

        "scope":
            "LY2119620 six trajectories only",

        "code": {
            "branch":
                git("branch", "--show-current"),
            "current_head":
                git("rev-parse", "HEAD"),
            "extractor_commit":
                extractor_commit,
            "extractor":
                artifact(RUNNER),
        },

        "closure_anchors": {
            "execution_manifest":
                artifact(EXEC_MANIFEST),
            "preflight":
                artifact(PREFLIGHT),
            "transfer_manifest":
                artifact(TRANSFER),
            "topology_sanitization":
                artifact(TOPOLOGY_AUDIT),
        },

        "run_report":
            artifact(RUN_FULL),

        "caches": records,

        "cache_count":
            len(records),

        "contract": {
            "trajectories": 6,
            "frames_per_trajectory": 400,
            "frame_spacing_ps": 50,
            "trajectory_ns": 20,
            "residues": 270,
            "bs256_width": 256,
            "next_block_frames": 20,
            "next_blocks_per_trajectory": 20,
        },

        "forbidden_operations_performed": {
            "phase2_calibration": False,
            "normalization_refit": False,
            "bandwidth_refit": False,
            "rff_refit": False,
            "graph_refit": False,
            "region_refit": False,
            "outcome_driven_tuning": False,
        },
    }

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUTPUT.write_text(
        json.dumps(receipt, indent=2) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "status":
                    receipt["status"],
                "cache_count":
                    receipt["cache_count"],
                "extractor_commit":
                    extractor_commit,
                "output":
                    str(OUTPUT),
                "output_sha256":
                    sha256(OUTPUT),
            },
            indent=2,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())