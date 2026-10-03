#!/usr/bin/env python
"""PACER-FKG v02 Phase 2: temporal blocks and frozen R2-only state.

This program deliberately stops before PACER scoring, contrasts, graph
diffusion, or biological analysis.  It never imports or executes ViSNet.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from project.pacer_fkg_v02.path_guard import guard_write_path  # noqa: E402


VERSION = "PACER_FKG_V02_PHASE2_CALIBRATION_v01"
SYSTEMS = (
    "apo",
    "probe_only",
    "compound110__candidate_no_probe",
    "compound110__candidate_probe",
)
REPLICAS = (1, 2, 3)
CALIBRATION_REPLICA = 2
N_FRAMES, N_RESIDUES, BS_WIDTH = 1000, 270, 256
BLOCK_FRAMES, N_BLOCKS = 20, 50
RFF_FEATURES = 512
BASE_SEED = 271828  # historical G2-B/G2-C policy
BRANCHES = {
    "STATE_MOTION": {"width": 512, "seed": BASE_SEED + 512},
    "SIGNED_DRIFT": {"width": 256, "seed": BASE_SEED + 256},
}

CACHE_ROOT = REPO / "project/cache/pacer_fkg_v02_longmd_v01"
BS_ROOT = CACHE_ROOT / "bs256/full"
BLOCK_ROOT = CACHE_ROOT / "phase2_blocks"
RESULT_ROOT = REPO / "project/results/pacer_fkg_v02_longmd_v01/calibration"
MANIFEST_ROOT = REPO / "project/results/pacer_fkg_v02_longmd_v01/encoder/manifests/full"
FREEZE_PATH = RESULT_ROOT / "V02_FREEZE_MANIFEST.json"
RUNNER_PATH = Path(__file__).resolve()
HISTORICAL_G2B = REPO / "project/pacer_dc_training/analyze_pacer_fkg_g2b.py"
HISTORICAL_G2C = REPO / "project/pacer_dc_training/analyze_pacer_fkg_g2c.py"
REGION_DEFINITIONS = REPO / "project/results/pacer_dc_four_context_v01/compound110/G2_REGION_MAP_v02.json"
GRAPH_DEFINITIONS = REPO / "project/results/pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json"
V01_PROTECTION = REPO / "project/pacer_fkg_v02/V01_PROTECTION_MANIFEST_SHA256.txt"


def utc_now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def safe_path(path: Path) -> Path:
    return guard_write_path(path, REPO)


def atomic_json(path: Path, value: Any) -> None:
    target = safe_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.partial.{os.getpid()}")
    safe_path(temporary)
    try:
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_npy(path: Path, value: np.ndarray) -> None:
    target = safe_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.partial.{os.getpid()}.npy")
    safe_path(temporary)
    try:
        np.save(temporary, value, allow_pickle=False)
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def artifact(path: Path) -> dict[str, Any]:
    return {"path": path.relative_to(REPO).as_posix(), "sha256": sha256(path), "bytes": path.stat().st_size}


def cache_path(system: str, replica: int) -> Path:
    return BS_ROOT / system / f"replica_{replica:02d}" / "C1_BS256.npy"


def phase1_manifest_path(system: str, replica: int) -> Path:
    return MANIFEST_ROOT / f"{system}__replica_{replica:02d}.json"


def branch_path(system: str, replica: int, branch: str) -> Path:
    return BLOCK_ROOT / system / f"replica_{replica:02d}" / f"{branch}.npy"


def authenticate_phase1() -> list[dict[str, Any]]:
    """Authenticate exactly the twelve full BS256 cache/manifest pairs."""
    records = []
    for system in SYSTEMS:
        for replica in REPLICAS:
            cache = cache_path(system, replica)
            manifest_path = phase1_manifest_path(system, replica)
            if not cache.is_file() or not manifest_path.is_file():
                raise FileNotFoundError(f"missing Phase-1 pair for {system} replica_{replica:02d}")
            manifest = load_json(manifest_path)
            required = {
                "schema": "pacer_fkg_v02.phase1_bs256_cache_manifest.v1",
                "mode": "full",
                "system": system,
                "replica": replica,
                "shape": [N_FRAMES, N_RESIDUES, BS_WIDTH],
                "dtype": "float32",
                "frame_ids": list(range(N_FRAMES)),
            }
            mismatches = [key for key, value in required.items() if manifest.get(key) != value]
            if mismatches:
                raise RuntimeError(f"Phase-1 manifest mismatch for {system} R{replica}: {mismatches}")
            stated = manifest.get("cache", {})
            actual_hash = sha256(cache)
            if stated.get("sha256") != actual_hash or stated.get("bytes") != cache.stat().st_size:
                raise RuntimeError(f"Phase-1 cache hash/size mismatch for {system} R{replica}")
            if manifest.get("replay", {}).get("passed") is not True:
                raise RuntimeError(f"Phase-1 replay is not authenticated for {system} R{replica}")
            array = np.load(cache, mmap_mode="r", allow_pickle=False)
            if array.shape != (N_FRAMES, N_RESIDUES, BS_WIDTH) or array.dtype != np.float32:
                raise RuntimeError(f"Phase-1 cache array contract failed for {system} R{replica}")
            if not all(np.isfinite(array[start:start + 20]).all() for start in range(0, N_FRAMES, 20)):
                raise RuntimeError(f"non-finite Phase-1 cache for {system} R{replica}")
            records.append({
                "system": system, "replica": replica,
                "cache": artifact(cache), "manifest": artifact(manifest_path),
            })
    if len(records) != 12:
        raise RuntimeError("expected exactly twelve authenticated Phase-1 caches")
    return records


def construct_blocks(frames: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Construct 50 independent contiguous blocks from one trajectory."""
    if frames.shape != (N_FRAMES, N_RESIDUES, BS_WIDTH):
        raise ValueError(f"expected {(N_FRAMES, N_RESIDUES, BS_WIDTH)}, got {frames.shape}")
    blocks = np.asarray(frames, dtype=np.float32).reshape(N_BLOCKS, BLOCK_FRAMES, N_RESIDUES, BS_WIDTH)
    mu_z = blocks.mean(axis=1, dtype=np.float64).astype(np.float32)
    delta = np.diff(blocks, axis=1)
    mu_delta = delta.mean(axis=1, dtype=np.float64).astype(np.float32)
    rms_delta = np.sqrt(np.mean(np.square(delta, dtype=np.float64), axis=1)).astype(np.float32)
    state_motion = np.concatenate((mu_z, rms_delta), axis=-1).astype(np.float32, copy=False)
    return state_motion, mu_delta


def verify_endpoint_identity(frames: np.ndarray, mu_delta: np.ndarray) -> bool:
    blocks = np.asarray(frames).reshape(N_BLOCKS, BLOCK_FRAMES, N_RESIDUES, BS_WIDTH)
    endpoint = (blocks[:, -1] - blocks[:, 0]) / np.float32(BLOCK_FRAMES - 1)
    return bool(np.allclose(mu_delta, endpoint, rtol=2e-5, atol=2e-6))


def fit_channels(samples: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    flat = samples.reshape(-1, samples.shape[-1]).astype(np.float64)
    center = np.median(flat, axis=0)
    scale = np.maximum(1.4826 * np.median(np.abs(flat - center), axis=0), 1e-6)
    return center.astype(np.float64), scale.astype(np.float64)


def shared_bandwidth(data: np.ndarray) -> float:
    x = data.astype(np.float64)
    squared = np.sum(x * x, axis=1)
    d2 = np.maximum(squared[:, None] + squared[None, :] - 2 * x @ x.T, 0)
    upper = d2[np.triu_indices(len(x), 1)]
    positive = upper[upper > 1e-12]
    if not len(positive):
        raise ValueError("degenerate bandwidth calibration")
    return float(np.sqrt(np.median(positive)))


def make_rff(input_dim: int, bandwidth: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    weights = rng.normal(0, 1 / bandwidth, size=(input_dim, RFF_FEATURES)).astype(np.float32)
    bias = rng.uniform(0, 2 * np.pi, size=RFF_FEATURES).astype(np.float32)
    return weights, bias


def verify_v01_protection() -> int:
    lines = [line.strip() for line in V01_PROTECTION.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    checked = 0
    for line in lines:
        expected, raw_path = line.split(None, 1)
        path = Path(raw_path.strip())
        if not path.is_file() or sha256(path).lower() != expected.lower():
            raise RuntimeError(f"v01 protection hash mismatch: {path}")
        checked += 1
    return checked


def definitions() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    method = {
        "schema": "pacer_fkg_v02.phase2_method.v1", "version": VERSION,
        "scope": "block construction and R2-only numerical calibration; no PACER outcomes",
        "block_frames": BLOCK_FRAMES, "blocks_per_trajectory": N_BLOCKS,
        "temporal_branches": {
            "STATE_MOTION": {"definition": "concat(mean(Z), sqrt(mean(diff(Z)^2)))", "width": 512},
            "SIGNED_DRIFT": {"definition": "mean(diff(Z))", "width": 256},
        },
        "normalization": "per-channel median and max(1.4826*MAD, 1e-6)",
        "bandwidth": "historical shared-node RBF median positive pairwise Euclidean distance",
        "bandwidth_population": "five fixed evenly-spaced residues from each of all 200 R2 blocks",
        "rff": {"family": "Gaussian RBF random Fourier features", "features": RFF_FEATURES},
        "historical_sources": [artifact(HISTORICAL_G2B), artifact(HISTORICAL_G2C)],
    }
    block_defs = {
        "schema": "pacer_fkg_v02.block_definitions.v1", "frame_count": N_FRAMES,
        "block_size": BLOCK_FRAMES, "partial_blocks": False,
        "blocks": [{"block_id": i, "start_frame": i * BLOCK_FRAMES,
                    "stop_frame_exclusive": (i + 1) * BLOCK_FRAMES}
                   for i in range(N_BLOCKS)],
    }
    contrasts = {
        "schema": "pacer_fkg_v02.contrast_definitions.v1",
        "status": "DEFINED_BUT_NOT_EVALUATED_IN_PHASE2",
        "contexts": list(SYSTEMS),
        "axes": {
            "synergy_interaction": {SYSTEMS[3]: 1, SYSTEMS[1]: -1, SYSTEMS[2]: -1, SYSTEMS[0]: 1},
            "intrinsic_agonism": {SYSTEMS[2]: 1, SYSTEMS[0]: -1},
            "conditional_pam_effect": {SYSTEMS[3]: 1, SYSTEMS[1]: -1},
        },
    }
    return method, block_defs, contrasts


def run() -> dict[str, Any]:
    if FREEZE_PATH.exists():
        raise FileExistsError("numerical state is already frozen; use --verify")
    v01_before = verify_v01_protection()
    phase1 = authenticate_phase1()
    method, block_defs, contrasts = definitions()
    spec_paths = {
        "method_spec": RESULT_ROOT / "V02_METHOD_SPEC.json",
        "block_definitions": RESULT_ROOT / "V02_BLOCK_DEFINITIONS.json",
        "contrast_definitions": RESULT_ROOT / "V02_CONTRAST_DEFINITIONS.json",
        "seeds": RESULT_ROOT / "V02_SEEDS.json",
    }
    atomic_json(spec_paths["method_spec"], method)
    atomic_json(spec_paths["block_definitions"], block_defs)
    atomic_json(spec_paths["contrast_definitions"], contrasts)
    atomic_json(spec_paths["seeds"], {
        "schema": "pacer_fkg_v02.seeds.v1", "policy": "historical base plus branch input width",
        "historical_base_seed": BASE_SEED,
        "branches": {name: info["seed"] for name, info in BRANCHES.items()},
    })

    derived = []
    calibration_arrays: dict[str, list[np.ndarray]] = {name: [] for name in BRANCHES}
    inventory = []
    for record in phase1:
        system, replica = record["system"], record["replica"]
        frames = np.load(cache_path(system, replica), mmap_mode="r", allow_pickle=False)
        sm, sd = construct_blocks(frames)
        if not verify_endpoint_identity(frames, sd):
            raise RuntimeError(f"muDelta endpoint identity failed for {system} R{replica}")
        for branch, value in (("STATE_MOTION", sm), ("SIGNED_DRIFT", sd)):
            out = branch_path(system, replica, branch)
            atomic_npy(out, value)
            derived.append({"system": system, "replica": replica, "branch": branch,
                            "shape": list(value.shape), "dtype": str(value.dtype), **artifact(out)})
            if replica == CALIBRATION_REPLICA:
                calibration_arrays[branch].append(value)
        inventory.append({"system": system, "replica": replica, "blocks": N_BLOCKS,
                          "frame_ranges": [[i * BLOCK_FRAMES, (i + 1) * BLOCK_FRAMES] for i in range(N_BLOCKS)],
                          "used_for_fitting": replica == CALIBRATION_REPLICA})

    atomic_json(RESULT_ROOT / "V02_CALIBRATION_INVENTORY.json", {
        "schema": "pacer_fkg_v02.calibration_inventory.v1",
        "calibration_replica": CALIBRATION_REPLICA, "calibration_blocks": 200,
        "excluded_replicas": [1, 3], "trajectories": inventory,
    })
    node_ids = np.linspace(0, N_RESIDUES - 1, 5, dtype=int)
    state_artifacts = []
    for branch, branch_info in BRANCHES.items():
        samples = np.concatenate(calibration_arrays[branch], axis=0)
        if samples.shape != (200, N_RESIDUES, branch_info["width"]):
            raise RuntimeError(f"unexpected {branch} calibration population: {samples.shape}")
        center, scale = fit_channels(samples)
        normalized = (samples.astype(np.float64) - center) / scale
        bandwidth_samples = normalized[:, node_ids, :].reshape(-1, branch_info["width"])
        bandwidth = shared_bandwidth(bandwidth_samples)
        weights, bias = make_rff(branch_info["width"], bandwidth, branch_info["seed"])
        paths = {
            "center": RESULT_ROOT / branch / "normalization_center.npy",
            "scale": RESULT_ROOT / branch / "normalization_scale.npy",
            "bandwidth": RESULT_ROOT / branch / "bandwidth.json",
            "weights": RESULT_ROOT / branch / "rff_weights.npy",
            "bias": RESULT_ROOT / branch / "rff_bias.npy",
        }
        atomic_npy(paths["center"], center)
        atomic_npy(paths["scale"], scale)
        atomic_json(paths["bandwidth"], {
            "schema": "pacer_fkg_v02.bandwidth.v1", "branch": branch,
            "family": "RBF", "algorithm": "median positive pairwise Euclidean distance",
            "bandwidth": bandwidth, "calibration_replica": 2,
            "calibration_blocks": 200, "samples": len(bandwidth_samples),
            "residue_indices": node_ids.tolist(), "dtype": "float64",
        })
        atomic_npy(paths["weights"], weights)
        atomic_npy(paths["bias"], bias)
        state_artifacts.append({
            "branch": branch, "input_width": branch_info["width"], "rff_features": RFF_FEATURES,
            "normalization": {"dtype": "float64", "shape": [branch_info["width"]],
                              "center": artifact(paths["center"]), "scale": artifact(paths["scale"])},
            "bandwidth": artifact(paths["bandwidth"]),
            "rff": {"seed": branch_info["seed"], "weights_shape": list(weights.shape),
                    "bias_shape": list(bias.shape), "dtype": "float32",
                    "weights": artifact(paths["weights"]), "bias": artifact(paths["bias"])},
        })

    calibration_inventory_path = RESULT_ROOT / "V02_CALIBRATION_INVENTORY.json"
    freeze = {
        "schema": "pacer_fkg_v02.freeze_manifest.v1", "version": VERSION,
        "created_at": utc_now(), "V02_NUMERICAL_STATE_FROZEN": True,
        "scope_boundary": "No PACER outcome, contrast, graph diffusion, or biological analysis was executed.",
        "phase1_inputs": phase1, "derived_block_arrays": derived,
        "definitions": {
            **{name: artifact(path) for name, path in spec_paths.items()},
            "region_definitions": artifact(REGION_DEFINITIONS),
            "graph_definitions": artifact(GRAPH_DEFINITIONS),
        },
        "calibration_population": artifact(calibration_inventory_path),
        "fitted_state": state_artifacts,
        "source_code": [artifact(RUNNER_PATH), artifact(HISTORICAL_G2B), artifact(HISTORICAL_G2C)],
        "v01_protection_entries_verified_before_and_after": v01_before,
    }
    if verify_v01_protection() != v01_before:
        raise RuntimeError("v01 protection inventory changed during Phase 2")
    atomic_json(FREEZE_PATH, freeze)
    return {"status": "FROZEN", "blocks": 600, "calibration_blocks": 200,
            "freeze_manifest": str(FREEZE_PATH)}


def replay_state(branch: str, inventory: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, float, np.ndarray, np.ndarray]:
    width = BRANCHES[branch]["width"]
    arrays = []
    for item in inventory["trajectories"]:
        if item["used_for_fitting"]:
            arrays.append(np.load(branch_path(item["system"], item["replica"], branch), allow_pickle=False))
    samples = np.concatenate(arrays, axis=0)
    if samples.shape != (200, N_RESIDUES, width):
        raise RuntimeError(f"replay calibration shape failed for {branch}")
    center, scale = fit_channels(samples)
    normalized = (samples.astype(np.float64) - center) / scale
    node_ids = np.linspace(0, N_RESIDUES - 1, 5, dtype=int)
    bandwidth = shared_bandwidth(normalized[:, node_ids, :].reshape(-1, width))
    weights, bias = make_rff(width, bandwidth, BRANCHES[branch]["seed"])
    return center, scale, bandwidth, weights, bias


def verify() -> dict[str, Any]:
    if not FREEZE_PATH.is_file():
        raise FileNotFoundError("freeze manifest does not exist; run --run first")
    freeze = load_json(FREEZE_PATH)
    if freeze.get("V02_NUMERICAL_STATE_FROZEN") is not True:
        raise RuntimeError("numerical-state freeze flag is not true")
    phase1 = authenticate_phase1()
    frozen_phase1 = {(x["system"], x["replica"]): x for x in freeze["phase1_inputs"]}
    for current in phase1:
        old = frozen_phase1.get((current["system"], current["replica"]))
        if old != current:
            raise RuntimeError(f"Phase-1 input changed after freeze: {current['system']} R{current['replica']}")
    for group in (freeze["derived_block_arrays"], freeze["definitions"].values(),
                  [freeze["calibration_population"]], freeze["source_code"]):
        for record in group:
            path = REPO / record["path"]
            if not path.is_file() or sha256(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
                raise RuntimeError(f"frozen artifact mismatch: {path}")
    for state in freeze["fitted_state"]:
        records = [state["normalization"]["center"], state["normalization"]["scale"], state["bandwidth"],
                   state["rff"]["weights"], state["rff"]["bias"]]
        for record in records:
            path = REPO / record["path"]
            if sha256(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
                raise RuntimeError(f"fitted-state hash mismatch: {path}")
    inventory = load_json(REPO / freeze["calibration_population"]["path"])
    fitted = [x for x in inventory["trajectories"] if x["used_for_fitting"]]
    excluded = [x for x in inventory["trajectories"] if not x["used_for_fitting"]]
    if len(fitted) != 4 or any(x["replica"] != 2 for x in fitted):
        raise RuntimeError("calibration inventory is not R2-only")
    if len(excluded) != 8 or {x["replica"] for x in excluded} != {1, 3}:
        raise RuntimeError("R1/R3 exclusion proof failed")
    state_by_branch = {x["branch"]: x for x in freeze["fitted_state"]}
    for branch in BRANCHES:
        center, scale, bandwidth, weights, bias = replay_state(branch, inventory)
        state = state_by_branch[branch]
        if not np.array_equal(center, np.load(REPO / state["normalization"]["center"]["path"])):
            raise RuntimeError(f"deterministic center replay failed for {branch}")
        if not np.array_equal(scale, np.load(REPO / state["normalization"]["scale"]["path"])):
            raise RuntimeError(f"deterministic scale replay failed for {branch}")
        if bandwidth != load_json(REPO / state["bandwidth"]["path"])["bandwidth"]:
            raise RuntimeError(f"deterministic bandwidth replay failed for {branch}")
        if not np.array_equal(weights, np.load(REPO / state["rff"]["weights"]["path"])):
            raise RuntimeError(f"deterministic RFF weight replay failed for {branch}")
        if not np.array_equal(bias, np.load(REPO / state["rff"]["bias"]["path"])):
            raise RuntimeError(f"deterministic RFF bias replay failed for {branch}")
    protected = verify_v01_protection()
    return {"status": "VERIFIED_FROZEN", "phase1_pairs": 12, "trajectory_blocks": 600,
            "calibration_blocks": 200, "fitting_replicas": [2], "excluded_replicas": [1, 3],
            "v01_protected_files": protected}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--verify", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    result = run() if args.run else verify()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
