#!/usr/bin/env python
"""PACER-DC 20 ns closure Track-B Phase-2A frozen apply-only.

Input:
    six frozen C1-BS256 caches, each [400, 270, 256]

Output per trajectory:
    STATE_MOTION descriptor [20, 270, 512]
    SIGNED_DRIFT descriptor [20, 270, 256]
    STATE_MOTION frozen-RFF [20, 270, 512]
    SIGNED_DRIFT frozen-RFF [20, 270, 512]

No calibration or refitting is permitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


from project.pacer_fkg_v02.run_phase3_fkg import (  # noqa: E402
    REQUIRED_FREEZE_SHA256,
    load_frozen_state,
    normalize_frozen,
    rff_frozen,
    verify_frozen_anchor,
)


VERSION = "PACER_DC_CLOSE_LOOP_PHASE2A_FROZEN_APPLY_v01"

CLOSURE_ROOT = (
    REPO
    / "project/results/pacer_dc_close_loop_20ns_v01"
)

PHASE1_RECEIPT = (
    CLOSURE_ROOT
    / "track_b_fkg_v02"
    / "phase1_bs256"
    / "PHASE1_FREEZE_RECEIPT_v01.json"
)

EXPECTED_PHASE1_RECEIPT_SHA = (
    "58b0c76459683018b1c850dfe12fd1b4"
    "b616a7fe79727aa818d2aba0d26070ad"
)

CACHE_ROOT = (
    REPO
    / "project/cache/pacer_dc_close_loop_20ns_v01"
    / "track_b_fkg_v02"
    / "phase2a_frozen_apply"
)

REPORT_ROOT = (
    CLOSURE_ROOT
    / "track_b_fkg_v02"
    / "phase2a_frozen_apply"
)

SYSTEMS = (
    "LY2119620__candidate_no_probe",
    "LY2119620__candidate_probe",
)

REPLICAS = (1, 2, 3)

N_FRAMES = 400
N_RESIDUES = 270
BS_WIDTH = 256

BLOCK_FRAMES = 20
N_BLOCKS = 20

STATE_MOTION_WIDTH = 512
SIGNED_DRIFT_WIDTH = 256
RFF_WIDTH = 512


@dataclass(frozen=True)
class Job:
    system: str
    replica: int
    input_cache: Path
    input_sha256: str

    @property
    def key(self) -> str:
        return f"{self.system}__replica_{self.replica:02d}"


def utc_now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


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


def atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp = path.with_name(
        f".{path.name}.partial.{os.getpid()}"
    )

    tmp.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )

    os.replace(tmp, path)


def atomic_npy(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp = path.with_name(
        f".{path.name}.partial.{os.getpid()}"
    )

    with tmp.open("wb") as handle:
        np.save(
            handle,
            np.asarray(array),
            allow_pickle=False,
        )

    os.replace(tmp, path)


def artifact(path: Path) -> dict:
    return {
        "path": path.relative_to(REPO).as_posix(),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


def verify_phase1_anchor() -> dict:
    if not PHASE1_RECEIPT.is_file():
        raise FileNotFoundError(PHASE1_RECEIPT)

    actual = sha256(PHASE1_RECEIPT)

    if actual != EXPECTED_PHASE1_RECEIPT_SHA:
        raise RuntimeError(
            "Phase-1 freeze receipt SHA mismatch"
        )

    receipt = load_json(PHASE1_RECEIPT)

    if receipt.get("status") != "PHASE1_FROZEN":
        raise RuntimeError(
            "Phase-1 receipt is not frozen"
        )

    if receipt.get("cache_count") != 6:
        raise RuntimeError(
            "Phase-1 receipt does not contain six caches"
        )

    forbidden = receipt.get(
        "forbidden_operations_performed",
        {},
    )

    if any(bool(value) for value in forbidden.values()):
        raise RuntimeError(
            "Phase-1 receipt records forbidden operation"
        )

    return receipt


def build_jobs(receipt: dict) -> list[Job]:
    records = {
        (item["system"], int(item["replica"])): item
        for item in receipt["caches"]
    }

    jobs = []

    for system in SYSTEMS:
        for replica in REPLICAS:
            record = records.get(
                (system, replica)
            )

            if record is None:
                raise RuntimeError(
                    f"missing Phase-1 receipt record: "
                    f"{system} R{replica}"
                )

            cache_record = record["cache"]

            cache = REPO / cache_record["path"]

            if not cache.is_file():
                raise FileNotFoundError(cache)

            digest = sha256(cache)

            if digest != cache_record["sha256"]:
                raise RuntimeError(
                    f"Phase-1 cache SHA mismatch: "
                    f"{system} R{replica}"
                )

            arr = np.load(
                cache,
                mmap_mode="r",
                allow_pickle=False,
            )

            if arr.shape != (
                N_FRAMES,
                N_RESIDUES,
                BS_WIDTH,
            ):
                raise RuntimeError(
                    f"Phase-1 cache shape mismatch: "
                    f"{system} R{replica}"
                )

            if arr.dtype != np.float32:
                raise RuntimeError(
                    f"Phase-1 cache dtype mismatch: "
                    f"{system} R{replica}"
                )

            if not np.isfinite(arr).all():
                raise RuntimeError(
                    f"non-finite Phase-1 cache: "
                    f"{system} R{replica}"
                )

            jobs.append(
                Job(
                    system=system,
                    replica=replica,
                    input_cache=cache,
                    input_sha256=digest,
                )
            )

    if len(jobs) != 6:
        raise RuntimeError(
            f"expected six jobs, got {len(jobs)}"
        )

    return jobs


def construct_blocks(
    frames: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Exact historical v02 formulas, adapted only for 400→20 blocks."""

    if frames.shape != (
        N_FRAMES,
        N_RESIDUES,
        BS_WIDTH,
    ):
        raise ValueError(
            f"expected "
            f"{(N_FRAMES, N_RESIDUES, BS_WIDTH)}, "
            f"got {frames.shape}"
        )

    blocks = np.asarray(
        frames,
        dtype=np.float32,
    ).reshape(
        N_BLOCKS,
        BLOCK_FRAMES,
        N_RESIDUES,
        BS_WIDTH,
    )

    mu_z = blocks.mean(
        axis=1,
        dtype=np.float64,
    ).astype(np.float32)

    delta = np.diff(
        blocks,
        axis=1,
    )

    signed_drift = delta.mean(
        axis=1,
        dtype=np.float64,
    ).astype(np.float32)

    rms_delta = np.sqrt(
        np.mean(
            np.square(
                delta,
                dtype=np.float64,
            ),
            axis=1,
        )
    ).astype(np.float32)

    state_motion = np.concatenate(
        (
            mu_z,
            rms_delta,
        ),
        axis=-1,
    ).astype(
        np.float32,
        copy=False,
    )

    endpoint = (
        blocks[:, -1]
        - blocks[:, 0]
    ) / np.float32(
        BLOCK_FRAMES - 1
    )

    endpoint_identity = bool(
        np.allclose(
            signed_drift,
            endpoint,
            rtol=2e-5,
            atol=2e-6,
        )
    )

    return (
        state_motion,
        signed_drift,
        endpoint_identity,
    )


def evaluate_job(
    job: Job,
    anchor: dict,
    mode: str,
) -> dict:

    output_dir = (
        CACHE_ROOT
        / mode
        / job.system
        / f"replica_{job.replica:02d}"
    )

    manifest_path = (
        REPORT_ROOT
        / "manifests"
        / mode
        / f"{job.key}.json"
    )

    paths = {
        "state_motion_descriptor":
            output_dir / "STATE_MOTION_DESCRIPTOR.npy",

        "signed_drift_descriptor":
            output_dir / "SIGNED_DRIFT_DESCRIPTOR.npy",

        "state_motion_rff":
            output_dir / "STATE_MOTION_RFF.npy",

        "signed_drift_rff":
            output_dir / "SIGNED_DRIFT_RFF.npy",
    }

    if manifest_path.exists():
        raise FileExistsError(
            f"refusing overwrite: {manifest_path}"
        )

    if any(path.exists() for path in paths.values()):
        raise FileExistsError(
            f"refusing overwrite in {output_dir}"
        )

    frames = np.load(
        job.input_cache,
        mmap_mode="r",
        allow_pickle=False,
    )

    (
        state_motion,
        signed_drift,
        endpoint_identity,
    ) = construct_blocks(frames)

    if state_motion.shape != (
        N_BLOCKS,
        N_RESIDUES,
        STATE_MOTION_WIDTH,
    ):
        raise RuntimeError(
            f"STATE_MOTION shape mismatch: {job.key}"
        )

    if signed_drift.shape != (
        N_BLOCKS,
        N_RESIDUES,
        SIGNED_DRIFT_WIDTH,
    ):
        raise RuntimeError(
            f"SIGNED_DRIFT shape mismatch: {job.key}"
        )

    if not endpoint_identity:
        raise RuntimeError(
            f"SIGNED_DRIFT endpoint identity failed: "
            f"{job.key}"
        )

    if (
        not np.isfinite(state_motion).all()
        or not np.isfinite(signed_drift).all()
    ):
        raise RuntimeError(
            f"non-finite descriptors: {job.key}"
        )

    sm_state = load_frozen_state(
        anchor,
        "STATE_MOTION",
    )

    sd_state = load_frozen_state(
        anchor,
        "SIGNED_DRIFT",
    )

    sm_normalized = normalize_frozen(
        state_motion,
        sm_state,
    )

    sd_normalized = normalize_frozen(
        signed_drift,
        sd_state,
    )

    if (
        not np.isfinite(sm_normalized).all()
        or not np.isfinite(sd_normalized).all()
    ):
        raise RuntimeError(
            f"non-finite normalized descriptors: "
            f"{job.key}"
        )

    sm_rff = rff_frozen(
        sm_normalized,
        sm_state,
    )

    sd_rff = rff_frozen(
        sd_normalized,
        sd_state,
    )

    if sm_rff.shape != (
        N_BLOCKS,
        N_RESIDUES,
        RFF_WIDTH,
    ):
        raise RuntimeError(
            f"STATE_MOTION RFF shape mismatch: "
            f"{job.key}"
        )

    if sd_rff.shape != (
        N_BLOCKS,
        N_RESIDUES,
        RFF_WIDTH,
    ):
        raise RuntimeError(
            f"SIGNED_DRIFT RFF shape mismatch: "
            f"{job.key}"
        )

    if (
        not np.isfinite(sm_rff).all()
        or not np.isfinite(sd_rff).all()
    ):
        raise RuntimeError(
            f"non-finite RFF output: {job.key}"
        )

    atomic_npy(
        paths["state_motion_descriptor"],
        state_motion,
    )

    atomic_npy(
        paths["signed_drift_descriptor"],
        signed_drift,
    )

    atomic_npy(
        paths["state_motion_rff"],
        sm_rff,
    )

    atomic_npy(
        paths["signed_drift_rff"],
        sd_rff,
    )

    manifest = {
        "schema":
            "pacer_dc.close_loop_20ns."
            "phase2a_frozen_apply_manifest.v1",

        "version":
            VERSION,

        "created_at":
            utc_now(),

        "system":
            job.system,

        "replica":
            job.replica,

        "input": {
            "cache":
                artifact(job.input_cache),
            "expected_shape":
                [400, 270, 256],
        },

        "block_contract": {
            "frames":
                N_FRAMES,
            "block_frames":
                BLOCK_FRAMES,
            "blocks":
                N_BLOCKS,
            "frame_spacing_ps":
                50,
            "block_duration_ps":
                1000,
            "partial_blocks":
                False,
        },

        "definitions": {
            "STATE_MOTION":
                "concat(mean(Z), "
                "sqrt(mean(diff(Z)^2)))",

            "SIGNED_DRIFT":
                "mean(diff(Z))",

            "signed_drift_endpoint_identity":
                "(Z_last-Z_first)/(20-1)",

            "endpoint_identity_passed":
                endpoint_identity,
        },

        "frozen_numerical_state": {
            "freeze_manifest_sha256":
                REQUIRED_FREEZE_SHA256,

            "normalization_refit":
                False,

            "bandwidth_refit":
                False,

            "rff_refit":
                False,

            "random_number_generation":
                False,
        },

        "outputs": {
            name: artifact(path)
            for name, path in paths.items()
        },

        "shapes": {
            "STATE_MOTION_DESCRIPTOR":
                [20, 270, 512],

            "SIGNED_DRIFT_DESCRIPTOR":
                [20, 270, 256],

            "STATE_MOTION_RFF":
                [20, 270, 512],

            "SIGNED_DRIFT_RFF":
                [20, 270, 512],
        },

        "finite": True,

        "scope_boundary":
            "Descriptor construction and frozen "
            "normalization/RFF application only; "
            "no graph diffusion, regional pooling, "
            "context contrasts, scoring, biological "
            "interpretation, calibration, or refitting.",
    }

    atomic_json(
        manifest_path,
        manifest,
    )

    return {
        "job":
            job.key,

        "status":
            "PASS",

        "endpoint_identity":
            True,

        "outputs": {
            name: record["sha256"]
            for name, record
            in manifest["outputs"].items()
        },
    }


def run(mode: str) -> dict:
    receipt = verify_phase1_anchor()

    # This performs the historical deterministic frozen-state verification.
    anchor = verify_frozen_anchor()

    if REQUIRED_FREEZE_SHA256 != (
        "b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd"
    ):
        raise RuntimeError(
            "unexpected historical freeze SHA"
        )

    jobs = build_jobs(receipt)

    if mode == "smoke":
        jobs = [jobs[0]]

    results = []

    for index, job in enumerate(jobs, 1):
        print(
            f"=== {mode.upper()} "
            f"{index}/{len(jobs)}: "
            f"{job.key} ===",
            flush=True,
        )

        result = evaluate_job(
            job,
            anchor,
            mode,
        )

        results.append(result)

        print(
            f"[{job.key}] PASS",
            flush=True,
        )

    report = {
        "schema":
            "pacer_dc.close_loop_20ns."
            "phase2a_frozen_apply_run.v1",

        "version":
            VERSION,

        "created_at":
            utc_now(),

        "status":
            "PHASE2A_PASS",

        "mode":
            mode,

        "phase1_freeze_receipt": {
            "path":
                PHASE1_RECEIPT.relative_to(
                    REPO
                ).as_posix(),

            "sha256":
                EXPECTED_PHASE1_RECEIPT_SHA,
        },

        "historical_fkg_v02_freeze_sha256":
            REQUIRED_FREEZE_SHA256,

        "jobs":
            results,

        "job_count":
            len(results),

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
        },
    }

    report_path = (
        REPORT_ROOT
        / f"RUN_{mode.upper()}_v01.json"
    )

    if report_path.exists():
        raise FileExistsError(
            f"refusing overwrite: {report_path}"
        )

    atomic_json(
        report_path,
        report,
    )

    return {
        "status":
            report["status"],

        "mode":
            mode,

        "jobs":
            len(results),

        "freeze_manifest_sha256":
            REQUIRED_FREEZE_SHA256,

        "report":
            str(report_path),
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = parser.add_mutually_exclusive_group(
        required=True
    )

    mode.add_argument(
        "--smoke",
        action="store_true",
    )

    mode.add_argument(
        "--run",
        action="store_true",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    mode = (
        "smoke"
        if args.smoke
        else "full"
    )

    result = run(mode)

    print(
        json.dumps(
            result,
            indent=2,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())