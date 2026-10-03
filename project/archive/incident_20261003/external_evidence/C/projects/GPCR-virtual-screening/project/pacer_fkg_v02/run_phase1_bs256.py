#!/usr/bin/env python
"""PACER-FKG v02 Phase 1: frozen C1-BS256 extraction only.

This runner intentionally stops at the frame-level ``[T, 270, 256]`` cache.
It does not calculate temporal branches, normalization, kernels, RFF state,
PACER-FKG scores, contrasts, or biological metrics.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import torch


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from project.encoder_atom_readout_v01.experiment_d import (  # noqa: E402
    BACKBONE_ATOMS,
    atom_metadata,
    extract_readouts,
)
from project.encoder_intermediate_v01.core import (  # noqa: E402
    EXPECTED_CHECKPOINT_SHA256,
    PINNED_GEOM2VEC_COMMIT,
    sha256,
    verify_and_strict_load,
)
from project.encoder_intermediate_v01.run_phase_0_2 import (  # noqa: E402
    DEFAULT_CHECKPOINT,
    source_verification,
)
from project.pacer_dc_training.extract_geom2vec_atom14 import (  # noqa: E402
    ATOM14,
    ATOMIC_NUMBER,
)
from project.pacer_fkg_v02.path_guard import guard_write_path  # noqa: E402


VERSION = "PACER_FKG_V02_PHASE1_BS256_v01"
BASELINE_COMMIT = "0a4b4f238c73add75a09a0875f698250a39bdf33"
BASELINE_TAG = "encoder-candidate-ah-frozen-20260929"
CACHE_ROOT = REPO / "project/cache/pacer_fkg_v02_longmd_v01/bs256"
REPORT_ROOT = REPO / "project/results/pacer_fkg_v02_longmd_v01/encoder"
PROVENANCE_ROOT = REPO / "project/results/pacer_fkg_v02_longmd_v01/provenance"
PHASE0_GATE = REPO / "project/pacer_fkg_v02/PHASE0_GATE.json"
BASELINE_MANIFEST = REPO / "project/pacer_fkg_v02/BASELINE_TRACKED_TREE_SHA256.json"
CANDIDATE_SPEC = REPO / "project/encoder_candidate_v01/CANDIDATE_SPEC.json"
FROZEN_CORE = REPO / "project/encoder_intermediate_v01/core.py"
FROZEN_READOUT = REPO / "project/encoder_atom_readout_v01/experiment_d.py"
RUNNER_PATH = Path(__file__).resolve()
SELECTION = "chainID E and protein"
SYSTEMS = (
    "apo",
    "probe_only",
    "compound110__candidate_no_probe",
    "compound110__candidate_probe",
)
REPLICAS = (1, 2, 3)
N_RESIDUES = 270
BS_WIDTH = 256
BATCH_SIZE = 4  # frozen Experiment-D extraction batch size
SMOKE_FRAME_IDS = (0, 1, 2, 3)
REPLAY_FRAME_IDS = (0, 1, 2, 3)
REPLAY_RELATIVE_RMSE_LIMIT = 1e-4  # frozen Experiment-D numerical replay criterion
RESNAME_MAP = {"HSE": "HIS", "HSD": "HIS", "HSP": "HIS"}
RESTYPE_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


@dataclass(frozen=True)
class Job:
    mode: str
    system: str
    replica: int
    topology: Path
    trajectory: Path
    trajectory_sha256: str
    topology_sha256: str
    sequence: str
    source_frames: int
    frame_ids: tuple[int, ...]

    @property
    def key(self) -> str:
        return f"{self.system}__replica_{self.replica:02d}"

    @property
    def cache_path(self) -> Path:
        leaf = "C1_BS256.frames_0000_0003.npy" if self.mode == "smoke" else "C1_BS256.npy"
        return CACHE_ROOT / self.mode / self.system / f"replica_{self.replica:02d}" / leaf

    @property
    def manifest_path(self) -> Path:
        return REPORT_ROOT / "manifests" / self.mode / f"{self.key}.json"

    @property
    def output_shape(self) -> tuple[int, int, int]:
        return (len(self.frame_ids), N_RESIDUES, BS_WIDTH)


def utc_now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def hash_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def ensure_output_path(path: Path, required_root: Path) -> Path:
    """Apply the shared guard plus the narrower Phase-1 output contract."""
    resolved = guard_write_path(path, REPO)
    root = required_root.resolve(strict=False)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"Phase-1 output is outside required root {root}: {resolved}") from exc
    return resolved


def atomic_json(path: Path, payload: Any) -> None:
    target = ensure_output_path(path, REPORT_ROOT)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.partial.{os.getpid()}")
    ensure_output_path(temporary, REPORT_ROOT)
    try:
        temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def baseline_hashes() -> dict[str, dict]:
    manifest = load_json(BASELINE_MANIFEST)
    if manifest["baseline_commit"] != BASELINE_COMMIT:
        raise RuntimeError("Phase-0 baseline commit differs from the frozen candidate anchor")
    return {item["path"]: item for item in manifest["files"]}


def verify_phase0_and_frozen_candidate() -> dict:
    gate = load_json(PHASE0_GATE)
    checks = gate.get("checks", {})
    required = (
        checks.get("server_manifest_authenticated"),
        checks.get("server_local_checksums_match"),
        checks.get("production_files_complete"),
        checks.get("completion_qc_pass"),
        checks.get("preprocessing_route_compatible"),
        checks.get("atom_residue_compatible"),
        checks.get("repository_protection", {}).get("V01_RESULT_HASHES_UNCHANGED"),
        checks.get("repository_protection", {}).get("BASELINE_TRACKED_FILES_UNCHANGED"),
    )
    if gate.get("decision") not in {"PASS", "PASS_WITH_DOCUMENTED_CAVEAT"} or not all(required):
        raise RuntimeError("Phase-0 gate is not eligible for Phase-1 extraction")

    baseline = baseline_hashes()
    protected = (FROZEN_CORE, FROZEN_READOUT, CANDIDATE_SPEC)
    records = []
    for path in protected:
        relative = path.relative_to(REPO).as_posix()
        expected = baseline.get(relative)
        actual = sha256(path)
        if expected is None or actual != expected["sha256"]:
            raise RuntimeError(f"frozen candidate artifact differs from Phase-0 baseline: {relative}")
        records.append({"path": str(path), "sha256": actual})

    spec = load_json(CANDIDATE_SPEC)
    if spec.get("status") != "FROZEN_NO_TUNABLE_PARAMETERS":
        raise RuntimeError("candidate specification is not frozen")
    if spec.get("checkpoint_sha256") != EXPECTED_CHECKPOINT_SHA256:
        raise RuntimeError("candidate checkpoint hash differs from frozen core")
    if spec.get("geom2vec_source_commit") != PINNED_GEOM2VEC_COMMIT:
        raise RuntimeError("candidate Geom2Vec commit differs from frozen core")
    if spec.get("residue_readout", {}).get("BS_dimension") != BS_WIDTH:
        raise RuntimeError("candidate specification is not BS256")
    return {
        "phase0_gate": {"path": str(PHASE0_GATE), "sha256": sha256(PHASE0_GATE), "decision": gate["decision"]},
        "baseline_commit": BASELINE_COMMIT,
        "baseline_tag": BASELINE_TAG,
        "candidate_spec": {"path": str(CANDIDATE_SPEC), "sha256": sha256(CANDIDATE_SPEC)},
        "frozen_implementation": records,
    }


def phase0_jobs(mode: str) -> list[Job]:
    metadata = load_json(PROVENANCE_ROOT / "TRAJECTORY_METADATA_AUDIT.json")["trajectories"]
    compatibility = load_json(PROVENANCE_ROOT / "LONG_MD_ENCODER_COMPATIBILITY.json")["trajectories"]
    inventory = load_json(PROVENANCE_ROOT / "LONG_MD_FILE_INVENTORY.json")["files"]
    preprocessing = load_json(PROVENANCE_ROOT / "LONG_MD_ENCODER_PREPROCESSING_AUDIT.json")
    compat_by_key = {(x["system"], int(x["replica"])): x for x in compatibility}
    hash_by_path = {x["relative_path"].replace("\\", "/"): x for x in inventory}
    topology_by_system = {
        system: preprocessing["short_reference_authentication"][system]["topology"] for system in SYSTEMS
    }
    jobs = []
    for item in metadata:
        system, replica = item["system"], int(item["replica"])
        if mode == "smoke" and (system != "apo" or replica != 2):
            continue
        if system not in SYSTEMS or replica not in REPLICAS:
            continue
        comp = compat_by_key[(system, replica)]
        if not comp.get("compatible"):
            raise RuntimeError(f"Phase-0 compatibility is not PASS for {system} R{replica}")
        trajectory = Path(item["trajectory_selected"])
        topology = Path(item["topology_selected"])
        relative = f"{system}/replica_{replica:02d}/trajectory.dcd"
        inv = hash_by_path.get(relative)
        topo = topology_by_system[system]
        if inv is None or inv.get("hash_match_status") != "MATCH":
            raise RuntimeError(f"trajectory lacks authenticated Phase-0 checksum: {relative}")
        if Path(topo["path"]).resolve() != topology.resolve() or not topo.get("match"):
            raise RuntimeError(f"topology lacks authenticated Phase-0 checksum: {system}")
        source_frames = int(item["stored_frames"])
        frames = SMOKE_FRAME_IDS if mode == "smoke" else tuple(range(source_frames))
        if mode == "full" and source_frames != 1000:
            raise RuntimeError(f"full extraction requires exactly 1000 frames: {system} R{replica}")
        jobs.append(Job(mode, system, replica, topology, trajectory, inv["local_sha256"],
                        topo["sha256"], comp["sequence"], source_frames, tuple(frames)))
    expected = 1 if mode == "smoke" else 12
    if len(jobs) != expected:
        raise RuntimeError(f"Phase-0 inventory produced {len(jobs)} {mode} jobs; expected {expected}")
    return sorted(jobs, key=lambda x: (SYSTEMS.index(x.system), x.replica))


def authenticate_job_inputs(job: Job) -> dict:
    if not job.trajectory.is_file() or not job.topology.is_file():
        raise FileNotFoundError(f"missing authenticated input for {job.key}")
    started = time.perf_counter()
    trajectory_hash = sha256(job.trajectory)
    if trajectory_hash != job.trajectory_sha256:
        raise RuntimeError(f"trajectory SHA256 mismatch for {job.key}")
    topology_hash = sha256(job.topology)
    if topology_hash != job.topology_sha256:
        raise RuntimeError(f"topology SHA256 mismatch for {job.key}")
    return {
        "trajectory": {"path": str(job.trajectory), "sha256": trajectory_hash, "bytes": job.trajectory.stat().st_size},
        "topology": {"path": str(job.topology), "sha256": topology_hash, "bytes": job.topology.stat().st_size},
        "hash_verification_seconds": time.perf_counter() - started,
    }


def implementation_identity() -> dict:
    return {
        "runner": {"path": str(RUNNER_PATH), "sha256": sha256(RUNNER_PATH)},
        "frozen_core": {"path": str(FROZEN_CORE), "sha256": sha256(FROZEN_CORE)},
        "frozen_readout": {"path": str(FROZEN_READOUT), "sha256": sha256(FROZEN_READOUT)},
        "candidate_spec": {"path": str(CANDIDATE_SPEC), "sha256": sha256(CANDIDATE_SPEC)},
    }


def expected_manifest_identity(job: Job) -> dict:
    return {
        "schema": "pacer_fkg_v02.phase1_bs256_cache_manifest.v1",
        "version": VERSION,
        "mode": job.mode,
        "system": job.system,
        "replica": job.replica,
        "shape": list(job.output_shape),
        "dtype": "float32",
        "frame_ids": list(job.frame_ids),
        "sequence": job.sequence,
        "selection": SELECTION,
        "preprocessing": {"make_whole": False, "center": False, "align": False, "pbc_transform": False},
        "checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "geom2vec_source_commit": PINNED_GEOM2VEC_COMMIT,
        "trajectory_sha256": job.trajectory_sha256,
        "topology_sha256": job.topology_sha256,
        "implementation": implementation_identity(),
    }


def compare_manifest_identity(manifest: dict, expected: dict) -> list[str]:
    mismatches = []
    for key, value in expected.items():
        if manifest.get(key) != value:
            mismatches.append(key)
    return mismatches


def finite_array(path: Path, expected_shape: tuple[int, int, int], chunk: int = 16) -> tuple[np.ndarray, bool]:
    array = np.load(path, mmap_mode="r", allow_pickle=False)
    if array.shape != expected_shape or array.dtype != np.dtype("float32"):
        return array, False
    for start in range(0, len(array), chunk):
        if not np.isfinite(array[start:start + chunk]).all():
            return array, False
    return array, True


def close_memmap(array) -> None:
    mapping = getattr(array, "_mmap", None) if array is not None else None
    if mapping is not None:
        mapping.close()


def validate_cache(job: Job, verify_source_hashes: bool = True) -> tuple[dict, dict | None]:
    cache, manifest_path = job.cache_path, job.manifest_path
    if not cache.exists() and not manifest_path.exists():
        return {"status": "missing", "cache": str(cache), "manifest": str(manifest_path)}, None
    if not cache.is_file() or not manifest_path.is_file():
        raise RuntimeError(f"incomplete cache/manifest pair for {job.key}; refusing overwrite")
    manifest = load_json(manifest_path)
    mismatches = compare_manifest_identity(manifest, expected_manifest_identity(job))
    if mismatches:
        raise RuntimeError(f"cache manifest identity mismatch for {job.key}: {mismatches}; refusing overwrite")
    actual_hash = sha256(cache)
    cache_record = manifest.get("cache", {})
    if cache_record.get("sha256") != actual_hash:
        raise RuntimeError(f"cache SHA256 mismatch for {job.key}; refusing overwrite")
    if Path(cache_record.get("path", "")).resolve() != cache.resolve() or cache_record.get("bytes") != cache.stat().st_size:
        raise RuntimeError(f"cache path/size manifest mismatch for {job.key}; refusing overwrite")
    expected_replay = [frame for frame in REPLAY_FRAME_IDS if frame in job.frame_ids]
    replay = manifest.get("replay", {})
    if replay.get("passed") is not True or replay.get("frame_ids") != expected_replay:
        raise RuntimeError(f"cache replay manifest is not valid for {job.key}; refusing overwrite")
    array, finite = finite_array(cache, job.output_shape)
    if not finite:
        raise RuntimeError(f"cache shape/dtype/finite validation failed for {job.key}; refusing overwrite")
    gly = np.asarray([aa == "G" for aa in job.sequence])
    if np.count_nonzero(array[:, gly, 128:]) != 0:
        raise RuntimeError(f"Gly side-chain block is not exact zero for {job.key}")
    observed_shape, observed_dtype = list(array.shape), str(array.dtype)
    close_memmap(array)
    source = authenticate_job_inputs(job) if verify_source_hashes else None
    return {
        "status": "valid",
        "cache": str(cache),
        "manifest": str(manifest_path),
        "cache_sha256": actual_hash,
        "shape": observed_shape,
        "dtype": observed_dtype,
    }, source


def topology_mapping(job: Job):
    import MDAnalysis as mda

    universe = mda.Universe(str(job.topology), str(job.trajectory))
    receptor = universe.select_atoms(SELECTION)
    residues = list(receptor.residues)
    if len(residues) != N_RESIDUES:
        raise RuntimeError(f"receptor residue count changed for {job.key}: {len(residues)}")
    sequence = []
    atom_indices = []
    for residue_index, residue in enumerate(residues):
        standard = RESNAME_MAP.get(residue.resname, residue.resname)
        aa = RESTYPE_3TO1.get(standard)
        if aa is None:
            raise RuntimeError(f"unsupported receptor residue {residue.resname} at {residue_index}")
        sequence.append(aa)
        local = {}
        for atom in residue.atoms:
            name = "CD1" if standard == "ILE" and atom.name == "CD" else atom.name
            if name in local:
                raise RuntimeError(f"duplicate atom name {name} in receptor residue {residue_index}")
            local[name] = atom.index
        missing = [name for name in ATOM14[aa] if name not in local]
        if missing:
            raise RuntimeError(f"missing expected atoms at receptor residue {residue_index}: {missing}")
        atom_indices.extend(local[name] for name in ATOM14[aa])
    observed_sequence = "".join(sequence)
    if observed_sequence != job.sequence:
        raise RuntimeError(f"receptor sequence differs from Phase-0 mapping for {job.key}")
    residue_index, backbone_mask, sidechain_mask, records = atom_metadata(observed_sequence)
    atom_z = np.asarray([ATOMIC_NUMBER[name[0]] for aa in observed_sequence for name in ATOM14[aa]], dtype=np.int64)
    if len(atom_indices) != len(atom_z) or len(atom_z) != len(residue_index):
        raise RuntimeError("atom14 heavy-atom mapping length mismatch")
    return universe, np.asarray(atom_indices, dtype=np.int64), atom_z, residue_index, backbone_mask, sidechain_mask, records


def load_frame_batch(universe, atom_indices: np.ndarray, frame_ids: tuple[int, ...] | list[int]) -> np.ndarray:
    coordinates = np.empty((len(frame_ids), len(atom_indices), 3), dtype=np.float32)
    for offset, frame in enumerate(frame_ids):
        universe.trajectory[int(frame)]
        coordinates[offset] = universe.atoms.positions[atom_indices]
    if not np.isfinite(coordinates).all():
        raise RuntimeError("non-finite input coordinates")
    return coordinates


def infer_bs256(model, universe, atom_indices, atom_z, residue_index, backbone_mask, sidechain_mask,
                frame_ids: tuple[int, ...] | list[int], device: str) -> np.ndarray:
    xyz = load_frame_batch(universe, atom_indices, frame_ids)
    readouts = extract_readouts(
        model, xyz, atom_z, residue_index, backbone_mask, sidechain_mask,
        N_RESIDUES, device, BATCH_SIZE,
    )
    bs = np.concatenate((readouts["B"], readouts["S"]), axis=-1).astype(np.float32, copy=False)
    if bs.shape != (len(frame_ids), N_RESIDUES, BS_WIDTH) or not np.isfinite(bs).all():
        raise RuntimeError(f"invalid BS256 output shape/values: {bs.shape}")
    return bs


def replay_metrics(observed: np.ndarray, replayed: np.ndarray) -> dict:
    delta = observed.astype(np.float64) - replayed.astype(np.float64)
    rmse = math.sqrt(float(np.mean(np.square(delta))))
    reference_rms = math.sqrt(float(np.mean(np.square(observed.astype(np.float64)))))
    relative_rmse = rmse / max(reference_rms, np.finfo(np.float64).tiny)
    return {
        "bitwise_identical": bool(np.array_equal(observed, replayed)),
        "max_abs": float(np.max(np.abs(delta))),
        "rmse": rmse,
        "relative_rmse": relative_rmse,
        "limit": REPLAY_RELATIVE_RMSE_LIMIT,
        "passed": bool(relative_rmse < REPLAY_RELATIVE_RMSE_LIMIT),
    }


def load_model(checkpoint: Path, device: str, geom2vec_source: Path):
    if sha256(checkpoint) != EXPECTED_CHECKPOINT_SHA256:
        raise RuntimeError("checkpoint SHA256 differs from frozen candidate")
    source = source_verification(geom2vec_source)
    if not source["source_equivalent"]:
        raise RuntimeError("installed Geom2Vec inference code differs from pinned source")
    if device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA device requested but torch.cuda.is_available() is false")
    model, checkpoint_report = verify_and_strict_load(checkpoint)
    if not checkpoint_report["complete_strict_coverage"]:
        raise RuntimeError("frozen checkpoint strict-load verification failed")
    model = model.to(device).eval()
    return model, {"checkpoint": checkpoint_report, "geom2vec_source": source}


def extract_job(job: Job, model, device: str, source_record: dict, model_record: dict) -> dict:
    target = ensure_output_path(job.cache_path, CACHE_ROOT)
    manifest_path = ensure_output_path(job.manifest_path, REPORT_ROOT)
    target.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.partial.{os.getpid()}")
    ensure_output_path(temporary, CACHE_ROOT)
    started = time.perf_counter()
    universe = None
    output = None
    cache = None
    try:
        universe, atom_indices, atom_z, residue_index, backbone_mask, sidechain_mask, records = topology_mapping(job)
        output = np.lib.format.open_memmap(temporary, mode="w+", dtype=np.float32, shape=job.output_shape)
        total = len(job.frame_ids)
        for start in range(0, total, BATCH_SIZE):
            frame_ids = job.frame_ids[start:start + BATCH_SIZE]
            output[start:start + len(frame_ids)] = infer_bs256(
                model, universe, atom_indices, atom_z, residue_index, backbone_mask, sidechain_mask,
                frame_ids, device,
            )
            output.flush()
            completed = start + len(frame_ids)
            elapsed = time.perf_counter() - started
            print(f"[{job.key}] {completed:4d}/{total} frames ({100.0*completed/total:6.2f}%) elapsed={elapsed:.1f}s", flush=True)
        gly = np.asarray([aa == "G" for aa in job.sequence])
        if np.count_nonzero(output[:, gly, 128:]) != 0:
            raise RuntimeError("Gly side-chain output is not exact zero")
        close_memmap(output)
        output = None
        replay_ids = tuple(frame for frame in REPLAY_FRAME_IDS if frame in job.frame_ids)
        cache = np.load(temporary, mmap_mode="r", allow_pickle=False)
        positions = [job.frame_ids.index(frame) for frame in replay_ids]
        observed = np.asarray(cache[positions])
        replayed = infer_bs256(
            model, universe, atom_indices, atom_z, residue_index, backbone_mask, sidechain_mask,
            replay_ids, device,
        )
        replay = replay_metrics(observed, replayed)
        replay["frame_ids"] = list(replay_ids)
        if not replay["passed"]:
            raise RuntimeError(f"deterministic replay failed for {job.key}: {replay}")
        close_memmap(cache)
        cache = None
        cache_hash = sha256(temporary)
        os.replace(temporary, target)
        identity = expected_manifest_identity(job)
        manifest = {
            **identity,
            "created_at": utc_now(),
            "cache": {"path": str(target), "sha256": cache_hash, "bytes": target.stat().st_size},
            "source": source_record,
            "model_provenance": model_record,
            "atom_mapping": {
                "expected_heavy_atoms": int(len(atom_indices)),
                "backbone_atoms": list(BACKBONE_ATOMS),
                "sidechain_present_residues": int(sum(record["sidechain_present"] for record in records)),
                "glycine_count": job.sequence.count("G"),
            },
            "replay": replay,
            "elapsed_seconds": time.perf_counter() - started,
            "output_scope": "C1-BS256 only; no normalization, RFF, PACER-FKG, contrasts, or biological metrics",
        }
        atomic_json(manifest_path, manifest)
        return {"status": "extracted", "job": job.key, "cache": str(target), "manifest": str(manifest_path),
                "cache_sha256": cache_hash, "shape": list(job.output_shape), "replay": replay,
                "elapsed_seconds": manifest["elapsed_seconds"]}
    except Exception:
        # Only this invocation's incomplete temporary is removable. A final cache
        # is never silently overwritten or deleted.
        close_memmap(cache)
        close_memmap(output)
        if temporary.exists():
            temporary.unlink()
        raise
    finally:
        if universe is not None:
            universe.trajectory.close()


def verify_replay(job: Job, model, device: str) -> dict:
    cache = np.load(job.cache_path, mmap_mode="r", allow_pickle=False)
    manifest = load_json(job.manifest_path)
    replay_ids = tuple(int(x) for x in manifest["replay"]["frame_ids"])
    positions = [job.frame_ids.index(frame) for frame in replay_ids]
    universe, atom_indices, atom_z, residue_index, backbone_mask, sidechain_mask, _ = topology_mapping(job)
    try:
        replayed = infer_bs256(model, universe, atom_indices, atom_z, residue_index, backbone_mask,
                               sidechain_mask, replay_ids, device)
    finally:
        universe.trajectory.close()
    observed = np.asarray(cache[positions])
    close_memmap(cache)
    result = replay_metrics(observed, replayed)
    result["frame_ids"] = list(replay_ids)
    if not result["passed"]:
        raise RuntimeError(f"deterministic replay verification failed for {job.key}")
    return result


def run_extract(mode: str, checkpoint: Path, device: str, geom2vec_source: Path) -> dict:
    frozen = verify_phase0_and_frozen_candidate()
    jobs = phase0_jobs(mode)
    results = []
    model = model_record = None
    for index, job in enumerate(jobs, 1):
        print(f"=== {mode.upper()} {index}/{len(jobs)}: {job.key} ===", flush=True)
        cache_status, source_record = validate_cache(job, verify_source_hashes=True)
        if cache_status["status"] == "valid":
            print(f"[{job.key}] valid cache verified; skipping extraction", flush=True)
            results.append({"job": job.key, "status": "skipped_valid_cache", **cache_status})
            continue
        if model is None:
            model, model_record = load_model(checkpoint, device, geom2vec_source)
        assert source_record is None
        source_record = authenticate_job_inputs(job)
        results.append(extract_job(job, model, device, source_record, model_record))
    return {"version": VERSION, "mode": mode, "created_at": utc_now(), "frozen_provenance": frozen,
            "device": device, "batch_size": BATCH_SIZE, "jobs": results,
            "forbidden_downstream_operations_performed": [], "status": "PASS"}


def discover_verify_jobs() -> list[Job]:
    manifests_root = REPORT_ROOT / "manifests"
    if not manifests_root.is_dir():
        raise RuntimeError("no Phase-1 cache manifests exist to verify")
    available = {(j.mode, j.system, j.replica): j for mode in ("smoke", "full") for j in phase0_jobs(mode)}
    jobs = []
    for path in sorted(manifests_root.rglob("*.json")):
        item = load_json(path)
        key = (item.get("mode"), item.get("system"), int(item.get("replica", -1)))
        job = available.get(key)
        if job is None or path.resolve() != job.manifest_path.resolve():
            raise RuntimeError(f"unexpected Phase-1 manifest: {path}")
        jobs.append(job)
    if not jobs:
        raise RuntimeError("no Phase-1 cache manifests exist to verify")
    cache_files = {path.resolve() for path in CACHE_ROOT.rglob("*") if path.is_file()} if CACHE_ROOT.is_dir() else set()
    expected_cache_files = {job.cache_path.resolve() for job in jobs}
    unexpected = sorted(str(path) for path in cache_files - expected_cache_files)
    if unexpected:
        raise RuntimeError(f"orphan or unexpected files exist in the BS256 cache root: {unexpected}")
    return jobs


def run_verify(checkpoint: Path, device: str, geom2vec_source: Path) -> dict:
    frozen = verify_phase0_and_frozen_candidate()
    jobs = discover_verify_jobs()
    structural = []
    for index, job in enumerate(jobs, 1):
        print(f"=== VERIFY HASH {index}/{len(jobs)}: {job.mode}/{job.key} ===", flush=True)
        status, source = validate_cache(job, verify_source_hashes=True)
        if status["status"] != "valid":
            raise RuntimeError(f"missing cache during verify: {job.key}")
        structural.append({**status, "source": source})
    model, model_record = load_model(checkpoint, device, geom2vec_source)
    replay = []
    for index, job in enumerate(jobs, 1):
        print(f"=== VERIFY REPLAY {index}/{len(jobs)}: {job.mode}/{job.key} ===", flush=True)
        replay.append({"mode": job.mode, "job": job.key, **verify_replay(job, model, device)})
    return {"version": VERSION, "mode": "verify", "created_at": utc_now(), "frozen_provenance": frozen,
            "device": device, "model_provenance": model_record, "structural_verification": structural,
            "deterministic_replay": replay, "forbidden_downstream_operations_performed": [], "status": "PASS"}


def write_run_report(mode: str, report: dict) -> Path:
    stamp = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S.%f%z")
    path = REPORT_ROOT / "runs" / f"phase1_{mode}_{stamp}.json"
    atomic_json(path, report)
    return path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--smoke", action="store_true", help="extract apo/replica_02 frames 0-3 only")
    modes.add_argument("--full", action="store_true", help="extract all 12 authenticated 1000-frame trajectories")
    modes.add_argument("--verify", action="store_true", help="verify all existing caches and replay fixed frames")
    parser.add_argument("--device", default="cuda", help="torch device (default: cuda)")
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--geom2vec-source", type=Path, default=Path(r"C:\projects\geom2vec-source"))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    mode = "smoke" if args.smoke else "full" if args.full else "verify"
    started = time.perf_counter()
    try:
        report = run_verify(args.checkpoint, args.device, args.geom2vec_source) if mode == "verify" else run_extract(
            mode, args.checkpoint, args.device, args.geom2vec_source
        )
        report["elapsed_seconds"] = time.perf_counter() - started
        report["runtime"] = {"python": sys.version, "executable": sys.executable, "platform": platform.platform(),
                             "torch": torch.__version__, "numpy": np.__version__}
        path = write_run_report(mode, report)
        print(json.dumps({"status": "PASS", "mode": mode, "report": str(path)}, indent=2), flush=True)
        return 0
    except Exception as exc:
        failure = {"version": VERSION, "mode": mode, "created_at": utc_now(), "status": "FAIL",
                   "error_type": type(exc).__name__, "error": str(exc), "traceback": traceback.format_exc(),
                   "elapsed_seconds": time.perf_counter() - started,
                   "forbidden_downstream_operations_performed": []}
        try:
            path = write_run_report(mode, failure)
            print(f"failure report: {path}", file=sys.stderr)
        except Exception as report_error:
            print(f"could not write failure report: {report_error}", file=sys.stderr)
        print(f"Phase-1 {mode} failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
