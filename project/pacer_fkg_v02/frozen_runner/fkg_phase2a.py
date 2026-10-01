"""Generic Phase-2A driver: frozen normalization + frozen RFF apply.

No fitting, no RNG, no bandwidth selection.  `verify_frozen_anchor()` must pass
before the first descriptor is written.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from project.pacer_fkg_v02.frozen_runner import blocks as blocklib
from project.pacer_fkg_v02.frozen_runner import engines, paths as pathlib_policy
from project.pacer_fkg_v02.frozen_runner.io_utils import check_npy
from project.pacer_fkg_v02.frozen_runner.spec import RunSpec, sha256_file

VERSION = "PACER_FKG_V02_FROZEN_RUNNER_PHASE2A_v01"


def _now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def _atomic_json(path: Path, payload: Any) -> None:
    tmp = path.with_name(f".{path.name}.partial.{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _atomic_npy(path: Path, array: np.ndarray) -> None:
    tmp = path.with_name(f".{path.name}.partial.{os.getpid()}")
    with tmp.open("wb") as handle:
        np.save(handle, np.asarray(array), allow_pickle=False)
    os.replace(tmp, path)


def phase1_records(spec: RunSpec) -> list[dict[str, Any]]:
    receipt_path = spec.paths.report_root / "phase1" / "PHASE1_FREEZE_RECEIPT_v01.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8-sig"))
    if receipt.get("status") != "PHASE1_FROZEN" or receipt.get("run_id") != spec.run_id:
        raise RuntimeError("Phase-1 receipt is not usable for this run spec")
    expected = len(spec.candidate_systems) * len(spec.replicas)
    if int(receipt.get("cache_count", -1)) != expected:
        raise RuntimeError("Phase-1 receipt cache_count mismatch")
    index = {(x["system"], int(x["replica"])): x for x in receipt["caches"]}
    out = []
    for system in spec.candidate_systems:
        for replica in spec.replicas:
            item = index[(system, replica)]
            cache = Path(item["cache"]["path"])
            if sha256_file(cache) != item["cache"]["sha256"]:
                raise RuntimeError(f"Phase-1 cache hash mismatch for {system} R{replica}")
            out.append({"system": system, "replica": replica, "path": cache, "record": item["cache"]})
    return out


def run(spec: RunSpec, *, repo_root: Path) -> dict[str, Any]:
    root = repo_root.resolve()
    receipt_path = pathlib_policy.guard_new_output_noclobber(
        spec.paths.report_root / "phase2a" / "PHASE2A_FREEZE_RECEIPT_v01.json", root, spec.paths.allowed_output_roots)

    anchor = engines.verify_frozen_anchor()          # fail closed
    engines.assert_frozen_numerics()

    outputs = []
    for item in phase1_records(spec):
        frames = check_npy(item["path"], (spec.n_frames, engines.N_RESIDUES, engines.P2_BS_WIDTH))
        state_motion, signed_drift = blocklib.construct_blocks(
            frames, block_frames=spec.block_frames, n_blocks=spec.blocks_per_trajectory)
        if not blocklib.verify_endpoint_identity(
                frames, signed_drift, block_frames=spec.block_frames,
                n_blocks=spec.blocks_per_trajectory):
            raise RuntimeError(f"SIGNED_DRIFT endpoint identity failed for {item['system']} R{item['replica']}")
        del frames
        descriptors = {"STATE_MOTION": state_motion, "SIGNED_DRIFT": signed_drift}
        for branch in ("STATE_MOTION", "SIGNED_DRIFT"):
            state = engines.load_frozen_state(anchor, branch)
            mapped = engines.rff_frozen(engines.normalize_frozen(descriptors[branch], state), state)
            expected_shape = (spec.blocks_per_trajectory, engines.N_RESIDUES, engines.RFF_FEATURES)
            if mapped.shape != expected_shape or not np.isfinite(mapped).all():
                raise RuntimeError(f"RFF contract failed for {item['system']} R{item['replica']} {branch}")
            directory = spec.paths.phase2a_cache_root / item["system"] / f"replica_{item['replica']:02d}"
            descriptor_path = pathlib_policy.guard_new_output_noclobber(
                directory / f"{branch}.npy", root, spec.paths.allowed_output_roots)
            rff_path = pathlib_policy.guard_new_output_noclobber(
                directory / f"{branch}_RFF.npy", root, spec.paths.allowed_output_roots)
            descriptor_path.parent.mkdir(parents=True, exist_ok=True)
            _atomic_npy(descriptor_path, descriptors[branch])
            _atomic_npy(rff_path, mapped)
            outputs.append({"system": item["system"], "replica": item["replica"], "branch": branch,
                            "input_cache": item["record"],
                            "descriptor": {"path": str(descriptor_path), "sha256": sha256_file(descriptor_path),
                                           "bytes": descriptor_path.stat().st_size},
                            "rff": {"path": str(rff_path), "sha256": sha256_file(rff_path),
                                    "bytes": rff_path.stat().st_size}})

    receipt = {
        "schema": "pacer.final.frozen_runner.phase2a_freeze_receipt.v1", "version": VERSION,
        "created_at": _now(), "status": "PHASE2A_FROZEN", "run_id": spec.run_id,
        "candidate_id": spec.molecule.candidate_id,
        "historical_fkg_v02_freeze_sha256": engines.REQUIRED_FREEZE_SHA256,
        "phase2_verification": anchor["phase2_verification"],
        "geometry": {"blocks_per_trajectory": spec.blocks_per_trajectory, "block_frames": spec.block_frames,
                     "frame_spacing_ps": 50, "block_duration_ns": 1.0},
        "phase1_receipt": {"path": str(spec.paths.report_root / "phase1" / "PHASE1_FREEZE_RECEIPT_v01.json")},
        "outputs": outputs,
        "forbidden_operations_performed": {"calibration": False, "normalization_refit": False,
                                           "bandwidth_refit": False, "rff_refit": False, "graph_refit": False,
                                           "region_refit": False, "threshold_selection": False,
                                           "outcome_driven_tuning": False},
        "run_spec": {"path": str(spec.source_path), "sha256": sha256_file(spec.source_path)},
    }
    _atomic_json(receipt_path, receipt)
    return {"status": "PHASE2A_FROZEN", "receipt": str(receipt_path),
            "receipt_sha256": sha256_file(receipt_path), "outputs": len(outputs)}


def verify(spec: RunSpec, *, repo_root: Path) -> dict[str, Any]:
    root = repo_root.resolve()
    engines.verify_frozen_anchor()
    engines.assert_frozen_numerics()
    receipt_path = spec.paths.report_root / "phase2a" / "PHASE2A_FREEZE_RECEIPT_v01.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8-sig"))
    if receipt.get("status") != "PHASE2A_FROZEN":
        raise RuntimeError("invalid Phase-2A receipt status")
    expected = len(spec.candidate_systems) * len(spec.replicas) * len(engines.BRANCHES)
    if len(receipt.get("outputs", [])) != expected:
        raise RuntimeError(f"Phase-2A receipt output count {len(receipt.get('outputs', []))} != {expected}")
    for item in receipt["outputs"]:
        for key in ("descriptor", "rff"):
            path = Path(item[key]["path"])
            if not path.is_file() or sha256_file(path) != item[key]["sha256"]:
                raise RuntimeError(f"Phase-2A artifact verification failed: {path}")
    return {"status": "PHASE2A_VERIFIED", "receipt_sha256": sha256_file(receipt_path), "outputs": expected}
