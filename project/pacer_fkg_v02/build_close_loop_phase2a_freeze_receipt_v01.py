#!/usr/bin/env python
"""Freeze PACER-DC 20 ns closure Track-B Phase-2A artifacts."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[2]

CLOSURE = REPO / "project/results/pacer_dc_close_loop_20ns_v01"

PHASE1_RECEIPT = (
    CLOSURE
    / "track_b_fkg_v02"
    / "phase1_bs256"
    / "PHASE1_FREEZE_RECEIPT_v01.json"
)

EXPECTED_PHASE1_SHA = (
    "58b0c76459683018b1c850dfe12fd1b4"
    "b616a7fe79727aa818d2aba0d26070ad"
)

EXPECTED_FKG_FREEZE_SHA = (
    "b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd"
)

EXPECTED_RUNNER_COMMIT = (
    "3423bcce5fab452f8d38a5167979d3b6d90261bf"
)

RUNNER = (
    REPO
    / "project/pacer_fkg_v02"
    / "run_close_loop_phase2a_frozen_apply_v01.py"
)

CACHE_ROOT = (
    REPO
    / "project/cache/pacer_dc_close_loop_20ns_v01"
    / "track_b_fkg_v02"
    / "phase2a_frozen_apply"
    / "full"
)

REPORT_ROOT = (
    CLOSURE
    / "track_b_fkg_v02"
    / "phase2a_frozen_apply"
)

RUN_FULL = REPORT_ROOT / "RUN_FULL_v01.json"
OUTPUT = REPORT_ROOT / "PHASE2A_FREEZE_RECEIPT_v01.json"

SYSTEMS = (
    "LY2119620__candidate_no_probe",
    "LY2119620__candidate_probe",
)

REPLICAS = (1, 2, 3)

EXPECTED_OUTPUTS = {
    "state_motion_descriptor": (
        "STATE_MOTION_DESCRIPTOR.npy",
        (20, 270, 512),
    ),
    "signed_drift_descriptor": (
        "SIGNED_DRIFT_DESCRIPTOR.npy",
        (20, 270, 256),
    ),
    "state_motion_rff": (
        "STATE_MOTION_RFF.npy",
        (20, 270, 512),
    ),
    "signed_drift_rff": (
        "SIGNED_DRIFT_RFF.npy",
        (20, 270, 512),
    ),
}


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

    if sha256(PHASE1_RECEIPT) != EXPECTED_PHASE1_SHA:
        raise RuntimeError(
            "Phase-1 receipt SHA mismatch"
        )

    phase1 = load_json(PHASE1_RECEIPT)

    if phase1.get("status") != "PHASE1_FROZEN":
        raise RuntimeError(
            "Phase-1 is not frozen"
        )

    runner_rel = RUNNER.relative_to(REPO).as_posix()

    runner_commit = git(
        "log",
        "-1",
        "--format=%H",
        "--",
        runner_rel,
    )

    if runner_commit != EXPECTED_RUNNER_COMMIT:
        raise RuntimeError(
            "Phase-2A runner commit changed"
        )

    run = load_json(RUN_FULL)

    if run.get("status") != "PHASE2A_PASS":
        raise RuntimeError(
            "RUN_FULL is not PHASE2A_PASS"
        )

    if run.get("job_count") != 6:
        raise RuntimeError(
            "RUN_FULL job count is not six"
        )

    if (
        run.get("historical_fkg_v02_freeze_sha256")
        != EXPECTED_FKG_FREEZE_SHA
    ):
        raise RuntimeError(
            "historical FKG-v02 freeze SHA changed"
        )

    forbidden = run.get(
        "forbidden_operations_performed",
        {},
    )

    if any(bool(value) for value in forbidden.values()):
        raise RuntimeError(
            "RUN_FULL records forbidden operation"
        )

    records = []
    artifact_count = 0

    for system in SYSTEMS:
        for replica in REPLICAS:
            key = f"{system}__replica_{replica:02d}"

            folder = (
                CACHE_ROOT
                / system
                / f"replica_{replica:02d}"
            )

            manifest_path = (
                REPORT_ROOT
                / "manifests"
                / "full"
                / f"{key}.json"
            )

            manifest = load_json(manifest_path)

            if (
                manifest["definitions"]
                ["endpoint_identity_passed"]
                is not True
            ):
                raise RuntimeError(
                    f"endpoint identity failed: {key}"
                )

            frozen = manifest["frozen_numerical_state"]

            for field in (
                "normalization_refit",
                "bandwidth_refit",
                "rff_refit",
                "random_number_generation",
            ):
                if frozen[field] is not False:
                    raise RuntimeError(
                        f"{field} not false: {key}"
                    )

            if (
                frozen["freeze_manifest_sha256"]
                != EXPECTED_FKG_FREEZE_SHA
            ):
                raise RuntimeError(
                    f"freeze SHA mismatch: {key}"
                )

            outputs = {}

            for manifest_key, (filename, shape) in EXPECTED_OUTPUTS.items():
                path = folder / filename

                arr = np.load(
                    path,
                    mmap_mode="r",
                    allow_pickle=False,
                )

                if arr.shape != shape:
                    raise RuntimeError(
                        f"shape mismatch: {key} {filename}"
                    )

                if arr.dtype != np.float32:
                    raise RuntimeError(
                        f"dtype mismatch: {key} {filename}"
                    )

                if not np.isfinite(arr).all():
                    raise RuntimeError(
                        f"non-finite: {key} {filename}"
                    )

                record = artifact(path)

                if (
                    manifest["outputs"][manifest_key]["sha256"]
                    != record["sha256"]
                ):
                    raise RuntimeError(
                        f"manifest SHA mismatch: "
                        f"{key} {filename}"
                    )

                outputs[manifest_key] = record
                artifact_count += 1

            records.append(
                {
                    "system": system,
                    "replica": replica,
                    "endpoint_identity": True,
                    "manifest": artifact(manifest_path),
                    "outputs": outputs,
                }
            )

    if len(records) != 6:
        raise RuntimeError(
            "expected six trajectory records"
        )

    if artifact_count != 24:
        raise RuntimeError(
            "expected 24 numerical artifacts"
        )

    receipt = {
        "schema":
            "pacer_dc.close_loop_20ns.track_b.phase2a_freeze.v1",

        "status":
            "PHASE2A_FROZEN",

        "track":
            "B / frozen C1-BS256 + PACER-FKG-v02",

        "code": {
            "branch":
                git("branch", "--show-current"),
            "current_head":
                git("rev-parse", "HEAD"),
            "runner_commit":
                runner_commit,
            "runner":
                artifact(RUNNER),
        },

        "upstream": {
            "phase1_freeze_receipt":
                artifact(PHASE1_RECEIPT),
            "historical_fkg_v02_freeze_sha256":
                EXPECTED_FKG_FREEZE_SHA,
        },

        "run_report":
            artifact(RUN_FULL),

        "records":
            records,

        "trajectory_count":
            6,

        "artifact_count":
            artifact_count,

        "contract": {
            "frames_per_trajectory":
                400,
            "frame_spacing_ps":
                50,
            "block_frames":
                20,
            "blocks_per_trajectory":
                20,
            "block_duration_ps":
                1000,

            "STATE_MOTION_DESCRIPTOR":
                [20, 270, 512],

            "SIGNED_DRIFT_DESCRIPTOR":
                [20, 270, 256],

            "STATE_MOTION_RFF":
                [20, 270, 512],

            "SIGNED_DRIFT_RFF":
                [20, 270, 512],
        },

        "endpoint_identity":
            "PASS_6_OF_6",

        "forbidden_operations_performed": {
            "phase2_calibration":
                False,
            "normalization_refit":
                False,
            "bandwidth_refit":
                False,
            "rff_refit":
                False,
            "graph_diffusion":
                False,
            "regional_pooling":
                False,
            "context_contrast":
                False,
            "biological_interpretation":
                False,
            "outcome_driven_tuning":
                False,
        },
    }

    OUTPUT.write_text(
        json.dumps(receipt, indent=2) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "status":
                    receipt["status"],
                "trajectory_count":
                    receipt["trajectory_count"],
                "artifact_count":
                    receipt["artifact_count"],
                "runner_commit":
                    runner_commit,
                "historical_fkg_v02_freeze_sha256":
                    EXPECTED_FKG_FREEZE_SHA,
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