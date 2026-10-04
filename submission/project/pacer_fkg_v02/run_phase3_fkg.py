#!/usr/bin/env python
"""PACER-FKG v02 Phase 3: frozen long-MD outcome evaluation.

This module is deliberately a consumer of Phase-2 state.  It contains no
fitting, bandwidth-selection, or random-number generation path.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy.sparse import csr_matrix


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from project.pacer_fkg_v02.path_guard import guard_write_path  # noqa: E402
from project.pacer_fkg_v02 import run_phase2_calibration as phase2  # noqa: E402


VERSION = "PACER_FKG_V02_PHASE3_LONGMD_v01"
REQUIRED_FREEZE_SHA256 = "b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd"
FREEZE_PATH = REPO / "project/results/pacer_fkg_v02_longmd_v01/calibration/V02_FREEZE_MANIFEST.json"
RESULT_ROOT = REPO / "project/results/pacer_fkg_v02_longmd_v01"
PHASE3_ROOT = RESULT_ROOT / "phase3"
GRAPH_PATH = REPO / "project/results/pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json"
REGION_PATH = REPO / "project/results/pacer_dc_four_context_v01/compound110/G2_REGION_MAP_v02.json"
CONTRAST_PATH = RESULT_ROOT / "calibration/V02_CONTRAST_DEFINITIONS.json"
PROTECTION_PATH = REPO / "project/pacer_fkg_v02/V01_PROTECTION_MANIFEST_SHA256.txt"

CONTEXTS = {
    "A": "apo",
    "P": "probe_only",
    "C": "compound110__candidate_no_probe",
    "CP": "compound110__candidate_probe",
}
REPLICAS = (1, 2, 3)
REPLICA_LABELS = {
    1: "FROZEN_APPLICATION_REPLICA",
    2: "CALIBRATION_REPLICA",
    3: "FROZEN_APPLICATION_REPLICA",
}
BRANCHES = {
    "STATE_MOTION": {"width": 512, "bandwidth": 32.30188361260893, "seed": 272340},
    "SIGNED_DRIFT": {"width": 256, "bandwidth": 29.29934899925964, "seed": 272084},
}
CONTRASTS = {
    "Delta_PAM": {"CP": 1.0, "P": -1.0},
    "Delta_AGO": {"C": 1.0, "A": -1.0},
    "Delta_INT": {"CP": 1.0, "P": -1.0, "C": -1.0, "A": 1.0},
}
HISTORICAL_AXIS = {
    "Delta_PAM": "conditional_pam_effect",
    "Delta_AGO": "intrinsic_agonism",
    "Delta_INT": "synergy_interaction",
}
ALPHA = 0.65
DIFFUSION_STEPS = 20
N_BLOCKS = 50
N_RESIDUES = 270
RFF_FEATURES = 512


def utc_now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def artifact(path: Path) -> dict[str, Any]:
    return {
        "path": path.relative_to(REPO).as_posix(),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


def staged_artifact(path: Path, stage: Path) -> dict[str, Any]:
    """Describe a staged file at the path it will have after atomic publish."""
    final_path = PHASE3_ROOT / path.relative_to(stage)
    return {
        "path": final_path.relative_to(REPO).as_posix(),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


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


def atomic_npz(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    target = safe_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.partial.{os.getpid()}.npz")
    safe_path(temporary)
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def _verify_record(record: Mapping[str, Any], description: str) -> Path:
    path = REPO / str(record["path"])
    if not path.is_file():
        raise RuntimeError(f"missing {description}: {path}")
    if sha256(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
        raise RuntimeError(f"{description} hash/size mismatch: {path}")
    return path


def verify_frozen_anchor() -> dict[str, Any]:
    """Fail closed unless every Phase-3 prerequisite remains frozen."""
    if not FREEZE_PATH.is_file() or sha256(FREEZE_PATH) != REQUIRED_FREEZE_SHA256:
        raise RuntimeError("required V02_FREEZE_MANIFEST SHA256 mismatch")
    freeze = load_json(FREEZE_PATH)
    if freeze.get("V02_NUMERICAL_STATE_FROZEN") is not True:
        raise RuntimeError("V02_NUMERICAL_STATE_FROZEN is not true")
    # Phase 2's deterministic verifier checks all Phase-1 cache/manifests,
    # block descriptors, normalization, bandwidth, RFF, source provenance,
    # R2-only calibration, R1/R3 exclusion, and the 41 protected v01 files.
    replay = phase2.verify()
    if replay != {
        "status": "VERIFIED_FROZEN", "phase1_pairs": 12,
        "trajectory_blocks": 600, "calibration_blocks": 200,
        "fitting_replicas": [2], "excluded_replicas": [1, 3],
        "v01_protected_files": 41,
    }:
        raise RuntimeError(f"unexpected Phase-2 verification receipt: {replay}")
    state = {item["branch"]: item for item in freeze["fitted_state"]}
    if set(state) != set(BRANCHES):
        raise RuntimeError("frozen branch inventory mismatch")
    for name, expected in BRANCHES.items():
        item = state[name]
        bandwidth = load_json(_verify_record(item["bandwidth"], f"{name} bandwidth"))["bandwidth"]
        if item["input_width"] != expected["width"] or item["rff_features"] != RFF_FEATURES:
            raise RuntimeError(f"{name} dimensions differ from mandatory anchor")
        if item["rff"]["seed"] != expected["seed"] or bandwidth != expected["bandwidth"]:
            raise RuntimeError(f"{name} bandwidth/seed differs from mandatory anchor")
    definitions = freeze["definitions"]
    if definitions["graph_definitions"]["path"] != GRAPH_PATH.relative_to(REPO).as_posix():
        raise RuntimeError("frozen graph provenance path mismatch")
    if definitions["region_definitions"]["path"] != REGION_PATH.relative_to(REPO).as_posix():
        raise RuntimeError("frozen region provenance path mismatch")
    _verify_record(definitions["graph_definitions"], "frozen graph")
    _verify_record(definitions["region_definitions"], "frozen regions")
    return {"freeze": freeze, "phase2_verification": replay, "state": state}


def load_frozen_state(anchor: Mapping[str, Any], branch: str) -> dict[str, Any]:
    item = anchor["state"][branch]
    paths = {
        "center": _verify_record(item["normalization"]["center"], f"{branch} center"),
        "scale": _verify_record(item["normalization"]["scale"], f"{branch} scale"),
        "weights": _verify_record(item["rff"]["weights"], f"{branch} RFF weights"),
        "bias": _verify_record(item["rff"]["bias"], f"{branch} RFF bias"),
    }
    values = {key: np.load(path, allow_pickle=False) for key, path in paths.items()}
    width = BRANCHES[branch]["width"]
    if values["center"].shape != (width,) or values["scale"].shape != (width,):
        raise RuntimeError(f"{branch} frozen normalization shape mismatch")
    if values["weights"].shape != (width, RFF_FEATURES) or values["bias"].shape != (RFF_FEATURES,):
        raise RuntimeError(f"{branch} frozen RFF shape mismatch")
    if not all(np.isfinite(value).all() for value in values.values()) or np.any(values["scale"] <= 0):
        raise RuntimeError(f"{branch} frozen state is non-finite or has invalid scale")
    return {**values, "records": item}


def normalize_frozen(descriptors: np.ndarray, state: Mapping[str, Any]) -> np.ndarray:
    if descriptors.shape[-1] != state["center"].shape[0]:
        raise ValueError("descriptor width does not match frozen normalization")
    return (descriptors.astype(np.float64) - state["center"]) / state["scale"]


def rff_frozen(normalized: np.ndarray, state: Mapping[str, Any]) -> np.ndarray:
    flat = normalized.reshape(-1, normalized.shape[-1]).astype(np.float32)
    mapped = np.cos(flat @ state["weights"] + state["bias"], dtype=np.float32)
    mapped *= np.float32(math.sqrt(2.0 / RFF_FEATURES))
    return mapped.reshape(*normalized.shape[:-1], RFF_FEATURES)


def historical_contrast(context_values: Mapping[str, np.ndarray], name: str) -> np.ndarray:
    """Exact v01 linear signed-representation contrast, using A/P/C/CP aliases."""
    signs = CONTRASTS[name]
    missing = set(signs) - set(context_values)
    if missing:
        raise KeyError(f"missing contrast contexts: {sorted(missing)}")
    result = None
    for context, sign in signs.items():
        term = np.asarray(context_values[context]) * sign
        result = term.copy() if result is None else result + term
    assert result is not None
    return result


def graph_edges(graph: Mapping[str, Any]) -> dict[tuple[int, int], tuple[str, float]]:
    edges: dict[tuple[int, int], tuple[str, float]] = {}
    for edge in graph["edges"]:
        i, j = int(edge["source"]), int(edge["target"])
        if i == j:
            continue
        pair = tuple(sorted((i, j)))
        kind = "backbone" if edge["edge_type"] == "backbone" else "contact"
        weight = 1.0 if kind == "backbone" else float(edge["support_fraction"])
        old = edges.get(pair)
        if old is None or (old[0] != "backbone" and kind == "backbone"):
            edges[pair] = (kind, weight)
        elif old[0] != "backbone":
            edges[pair] = (kind, max(old[1], weight))
    return edges


def graph_transition(graph: Mapping[str, Any]) -> csr_matrix:
    if len(graph["nodes"]) != N_RESIDUES:
        raise RuntimeError("frozen graph node count mismatch")
    row: list[int] = []
    col: list[int] = []
    data: list[float] = []
    for (i, j), (_, weight) in graph_edges(graph).items():
        row.extend((i, j)); col.extend((j, i)); data.extend((weight, weight))
    adjacency = csr_matrix((data, (row, col)), shape=(N_RESIDUES, N_RESIDUES), dtype=np.float64)
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    if np.any(degree <= 0):
        raise RuntimeError("frozen graph has isolated nodes")
    inverse = csr_matrix((1.0 / degree, (np.arange(N_RESIDUES), np.arange(N_RESIDUES))),
                         shape=(N_RESIDUES, N_RESIDUES))
    return inverse @ adjacency


def diffuse(values: np.ndarray, transition: csr_matrix) -> np.ndarray:
    """Historical restart diffusion, independently for each leading sample."""
    if values.shape[-2:] != (N_RESIDUES, RFF_FEATURES):
        raise ValueError("diffusion expects (..., 270, 512)")
    original = np.asarray(values, dtype=np.float32)
    out = original.copy()
    for _ in range(DIFFUSION_STEPS):
        propagated = (transition @ out.transpose(1, 0, 2).reshape(N_RESIDUES, -1)).reshape(
            N_RESIDUES, out.shape[0], RFF_FEATURES).transpose(1, 0, 2)
        out = ((1.0 - ALPHA) * original + ALPHA * propagated).astype(np.float32)
    return out


def region_indices(graph: Mapping[str, Any]) -> dict[str, np.ndarray]:
    regions = {
        name: np.asarray([int(member["embedding_index"]) for member in members], dtype=np.int64)
        for name, members in sorted(graph["regions"].items())
    }
    if not regions or any(len(indices) == 0 for indices in regions.values()):
        raise RuntimeError("empty frozen graph region")
    if any(np.any((indices < 0) | (indices >= N_RESIDUES)) for indices in regions.values()):
        raise RuntimeError("frozen region index out of range")
    return regions


def cosine(a: np.ndarray, b: np.ndarray) -> float | None:
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denominator <= 1e-12:
        return None
    return float(np.clip(np.dot(a, b) / denominator, -1.0, 1.0))


def distribution(values: np.ndarray) -> dict[str, Any]:
    x = np.asarray(values, dtype=np.float64)
    if x.shape != (N_BLOCKS,) or not np.isfinite(x).all():
        raise ValueError("block distribution must contain exactly 50 finite values")
    return {
        "n_correlated_time_blocks": N_BLOCKS,
        "mean": float(np.mean(x)), "median": float(np.median(x)),
        "minimum": float(np.min(x)), "maximum": float(np.max(x)),
        "q025": float(np.percentile(x, 2.5)), "q975": float(np.percentile(x, 97.5)),
        "interpretation": "descriptive only; blocks are correlated time samples, not biological replicates",
    }


def _replica_result(vectors: np.ndarray, vector_key: str) -> dict[str, Any]:
    if vectors.shape != (N_BLOCKS, RFF_FEATURES):
        raise ValueError("regional vectors must have shape (50, 512)")
    mean_vector = vectors.mean(axis=0, dtype=np.float64)
    magnitudes = np.linalg.norm(vectors.astype(np.float64), axis=1)
    return {
        "value": float(np.linalg.norm(mean_vector)),
        "value_definition": "L2 magnitude of the mean signed regional RFF contrast vector",
        "magnitude": float(np.linalg.norm(mean_vector)),
        "sign_direction": "signed vector; no scalar pharmacological sign is defined",
        "signed_direction_vector": mean_vector.tolist(),
        "direction_vector_key": vector_key,
        "block_magnitudes": magnitudes.tolist(),
        "block_distribution": distribution(magnitudes),
    }


def _context_replica_result(vectors: np.ndarray, vector_key: str) -> dict[str, Any]:
    result = _replica_result(vectors, vector_key)
    result["value_definition"] = "L2 magnitude of the mean regional context RFF vector"
    result["sign_direction"] = "context representation vector; contrast sign is not applicable"
    return result


def _verify_phase2_descriptor(anchor: Mapping[str, Any], system: str, replica: int,
                              branch: str, recomputed: np.ndarray) -> dict[str, Any]:
    match = [item for item in anchor["freeze"]["derived_block_arrays"]
             if item["system"] == system and item["replica"] == replica and item["branch"] == branch]
    if len(match) != 1:
        raise RuntimeError("frozen derived descriptor inventory is not unique")
    path = _verify_record(match[0], f"{branch} descriptor {system} R{replica}")
    frozen = np.load(path, mmap_mode="r", allow_pickle=False)
    if frozen.shape != recomputed.shape or not np.array_equal(frozen, recomputed):
        raise RuntimeError(f"Phase-2 branch reconstruction mismatch: {system} R{replica} {branch}")
    return artifact(path)


def evaluate_branch(branch: str, anchor: Mapping[str, Any], graph: Mapping[str, Any],
                    transition: csr_matrix, regions: Mapping[str, np.ndarray], stage: Path
                    ) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    state = load_frozen_state(anchor, branch)
    by_replica: dict[int, dict[str, dict[str, np.ndarray]]] = {}
    block_manifest: list[dict[str, Any]] = []
    for replica in REPLICAS:
        context_nodes: dict[str, np.ndarray] = {}
        source_records = []
        for alias, system in CONTEXTS.items():
            cache = phase2.cache_path(system, replica)
            frames = np.load(cache, mmap_mode="r", allow_pickle=False)
            state_motion, signed_drift = phase2.construct_blocks(frames)
            descriptor = state_motion if branch == "STATE_MOTION" else signed_drift
            frozen_descriptor = _verify_phase2_descriptor(anchor, system, replica, branch, descriptor)
            normalized = normalize_frozen(descriptor, state)
            if not np.isfinite(normalized).all():
                raise RuntimeError(f"non-finite normalized descriptor: {branch} {system} R{replica}")
            mapped = rff_frozen(normalized, state)
            if mapped.shape != (N_BLOCKS, N_RESIDUES, RFF_FEATURES) or not np.isfinite(mapped).all():
                raise RuntimeError(f"RFF contract failed: {branch} {system} R{replica}")
            context_nodes[alias] = mapped
            phase1_record = next(item for item in anchor["freeze"]["phase1_inputs"]
                                 if item["system"] == system and item["replica"] == replica)
            source_records.append({
                "context": alias, "system": system, "phase1_cache": phase1_record["cache"],
                "phase1_manifest": phase1_record["manifest"], "phase2_descriptor": frozen_descriptor,
            })
        context_region_vectors: dict[str, dict[str, np.ndarray]] = {}
        for alias, mapped in context_nodes.items():
            diffused = diffuse(mapped, transition)
            context_region_vectors[alias] = {
                name: diffused[:, indices, :].mean(axis=1, dtype=np.float64).astype(np.float32)
                for name, indices in regions.items()
            }
        contrasts = {
            contrast: {
                region: diffused[:, indices, :].mean(axis=1, dtype=np.float64).astype(np.float32)
                for region, indices in regions.items()
            }
            for contrast in CONTRASTS
            for diffused in [diffuse(historical_contrast(context_nodes, contrast), transition)]
        }
        by_replica[replica] = {"contexts": context_region_vectors, "contrasts": contrasts}
        arrays = {
            **{f"context__{alias}__{region}": vectors
               for alias, region_map in context_region_vectors.items() for region, vectors in region_map.items()},
            **{f"contrast__{name}__{region}": vectors
               for name, region_map in contrasts.items() for region, vectors in region_map.items()},
        }
        replica_dir = stage / branch.lower() / f"replica_{replica:02d}"
        vector_path = replica_dir / "BLOCK_REGIONAL_VECTORS.npz"
        atomic_npz(vector_path, arrays)
        rows = []
        for alias in CONTEXTS:
            for region, vectors in context_region_vectors[alias].items():
                rows.append({
                    "kind": "context", "name": alias, "system": CONTEXTS[alias], "region": region,
                    "block_values": np.linalg.norm(vectors.astype(np.float64), axis=1).tolist(),
                })
        for name, region_map in contrasts.items():
            for region, vectors in region_map.items():
                rows.append({
                    "kind": "contrast", "name": name, "historical_axis": HISTORICAL_AXIS[name],
                    "formula": CONTRASTS[name], "region": region,
                    "block_values": np.linalg.norm(vectors.astype(np.float64), axis=1).tolist(),
                })
        block_path = replica_dir / "BLOCK_LEVEL_VALUES.json"
        atomic_json(block_path, {
            "schema": "pacer_fkg_v02.phase3_block_values.v1",
            "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256, "branch": branch,
            "replica": replica, "replica_label": REPLICA_LABELS[replica],
            "block_count": N_BLOCKS, "blocks_are_correlated_time_samples": True,
            "source_records": source_records, "rows": rows,
        })
        block_manifest.append({
            "branch": branch, "replica": replica, "replica_label": REPLICA_LABELS[replica],
            "block_count": N_BLOCKS, "vectors": staged_artifact(vector_path, stage),
            "values": staged_artifact(block_path, stage),
            "source_records": source_records,
        })

    result: dict[str, Any] = {
        "schema": "pacer_fkg_v02.phase3_branch_results.v1",
        "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256,
        "branch": branch, "branch_input_dimension": BRANCHES[branch]["width"],
        "frozen_bandwidth": BRANCHES[branch]["bandwidth"],
        "frozen_seed": BRANCHES[branch]["seed"], "rff_features": RFF_FEATURES,
        "replica_policy": {"R1": REPLICA_LABELS[1], "R2": REPLICA_LABELS[2], "R3": REPLICA_LABELS[3]},
        "contexts": {}, "contrasts": {},
    }
    for alias, system in CONTEXTS.items():
        result["contexts"][alias] = {"system": system, "regions": {}}
        for region in regions:
            replicas = {}
            for replica in REPLICAS:
                key = f"context__{alias}__{region}"
                replicas[f"R{replica}"] = {
                    "replica_label": REPLICA_LABELS[replica],
                    **_context_replica_result(by_replica[replica]["contexts"][alias][region], key),
                }
            r1 = np.asarray(replicas["R1"]["signed_direction_vector"])
            r3 = np.asarray(replicas["R3"]["signed_direction_vector"])
            agreement = cosine(r1, r3)
            result["contexts"][alias]["regions"][region] = {
                "replicas": replicas, "R1_R3_direction_cosine": agreement,
                "R1_R3_direction_agreement": None if agreement is None else agreement > 0,
                "magnitudes_reported_separately": True,
            }
    for name, signs in CONTRASTS.items():
        result["contrasts"][name] = {
            "historical_axis": HISTORICAL_AXIS[name], "formula": signs, "regions": {},
        }
        for region in regions:
            replicas = {}
            for replica in REPLICAS:
                key = f"contrast__{name}__{region}"
                replicas[f"R{replica}"] = {
                    "replica_label": REPLICA_LABELS[replica],
                    **_replica_result(by_replica[replica]["contrasts"][name][region], key),
                }
            r1 = np.asarray(replicas["R1"]["signed_direction_vector"])
            r3 = np.asarray(replicas["R3"]["signed_direction_vector"])
            agreement = cosine(r1, r3)
            result["contrasts"][name]["regions"][region] = {
                "replicas": replicas, "R1_R3_direction_cosine": agreement,
                "R1_R3_direction_agreement": None if agreement is None else agreement > 0,
                "magnitudes_reported_separately": True,
                "no_magnitude_qualification_gate": True,
            }
    return result, block_manifest


def method_audit(anchor: Mapping[str, Any], graph: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema": "pacer_fkg_v02.phase3_method_audit.v1",
        "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256, "version": VERSION,
        "created_at": utc_now(), "status": "FROZEN_FORMAL_EVALUATION",
        "contexts": CONTEXTS, "replicas": {str(k): v for k, v in REPLICA_LABELS.items()},
        "branches": BRANCHES, "branches_analyzed_separately": True,
        "block_construction": {
            "frames_per_trajectory": 1000, "block_frames": 20, "blocks_per_trajectory": 50,
            "STATE_MOTION": "concat(mean(Z), sqrt(mean(diff(Z)^2)))",
            "SIGNED_DRIFT": "mean(diff(Z))", "lag_frames": 1,
        },
        "frozen_numerical_path": "descriptor -> frozen channel normalization -> frozen Gaussian RBF RFF -> signed context contrast -> frozen graph diffusion -> frozen region mean",
        "software_adaptation": "Phase-2 produces one descriptor per 20-frame block; unlike v01 frame-to-block averaging, the frozen v02 descriptor is mapped once per block. The RBF-RFF, signed contrast, graph diffusion, and regional-mean formulas are unchanged.",
        "contrast_definitions": {name: {"historical_axis": HISTORICAL_AXIS[name], "formula": signs}
                                 for name, signs in CONTRASTS.items()},
        "graph_diffusion": {"alpha": ALPHA, "steps": DIFFUSION_STEPS,
                            "transition": "row-stochastic weighted undirected adjacency with restart"},
        "regions": {name: len(members) for name, members in sorted(graph["regions"].items())},
        "historical_uncertainty": "descriptive distributions of contiguous blocks only; no new bootstrap, p-value, or significance filter",
        "prohibitions_enforced": {
            "no_fitting": True, "no_ViSNet": True, "no_branch_combination": True,
            "no_context_specific_kernel": True, "no_magnitude_gate": True,
            "no_independent_block_claim": True,
        },
        "phase2_verification": anchor["phase2_verification"],
    }


def consistency_report(branch_results: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    rows = []
    for branch, result in branch_results.items():
        for contrast, item in result["contrasts"].items():
            for region, regional in item["regions"].items():
                reps = regional["replicas"]
                rows.append({
                    "branch": branch, "contrast": contrast, "region": region,
                    "R1": {"value": reps["R1"]["value"], "sign_direction": reps["R1"]["sign_direction"]},
                    "R2": {"value": reps["R2"]["value"], "sign_direction": reps["R2"]["sign_direction"],
                           "replica_label": "CALIBRATION_REPLICA"},
                    "R3": {"value": reps["R3"]["value"], "sign_direction": reps["R3"]["sign_direction"]},
                    "R1_R3_direction_cosine": regional["R1_R3_direction_cosine"],
                    "R1_R3_direction_agreement": regional["R1_R3_direction_agreement"],
                    "magnitude_qualification_gate": None,
                })
    return {
        "schema": "pacer_fkg_v02.phase3_replica_consistency.v1",
        "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256,
        "primary_generalization_audit": "R1_vs_R3", "rows": rows,
        "calibration_replica": "R2", "frozen_application_replicas": ["R1", "R3"],
        "branches_combined": False,
    }


def report_markdown(branch_results: Mapping[str, Mapping[str, Any]]) -> str:
    lines = [
        "# PACER-FKG v02 Phase 3 formal frozen long-MD evaluation", "",
        f"freeze_manifest_sha256 = `{REQUIRED_FREEZE_SHA256}`", "",
        "R2 is the **CALIBRATION_REPLICA**. R1 and R3 are the **FROZEN_APPLICATION_REPLICAS**. "
        "All 50 blocks are correlated time samples and are summarized descriptively; they are not treated as biological replicates.", "",
        "STATE_MOTION and SIGNED_DRIFT are reported separately. No branch selection, combination, magnitude gate, p-value, or new bootstrap is used.", "",
    ]
    for branch, result in branch_results.items():
        lines += [f"## {branch}", "", "| Contrast | Region | R1 magnitude | R2 magnitude | R3 magnitude | R1/R3 direction cosine | Agreement |",
                  "|---|---|---:|---:|---:|---:|---|"]
        for contrast, item in result["contrasts"].items():
            for region, regional in item["regions"].items():
                reps = regional["replicas"]
                direction = regional["R1_R3_direction_cosine"]
                direction_text = "undefined" if direction is None else f"{direction:.6g}"
                lines.append(f"| {contrast} | {region} | {reps['R1']['value']:.6g} | {reps['R2']['value']:.6g} | {reps['R3']['value']:.6g} | {direction_text} | {regional['R1_R3_direction_agreement']} |")
        lines += ["", "The full signed direction vectors, context quantities, and all block-level magnitude distributions are serialized in the JSON/NPZ artifacts.", ""]
    lines += ["## Interpretation boundary", "",
              "These are representation-space directional and magnitude results. Weak, null, contradictory, and negative direction agreement is retained unchanged. No pharmacological efficacy or independent-replica significance claim is implied.", ""]
    return "\n".join(lines)


def smoke() -> dict[str, Any]:
    anchor = verify_frozen_anchor()
    graph = load_json(GRAPH_PATH)
    transition = graph_transition(graph)
    regions = region_indices(graph)
    # Deterministic, non-biological 20-frame construction; no output files.
    base = np.arange(20, dtype=np.float32)[:, None, None]
    residue = np.arange(N_RESIDUES, dtype=np.float32)[None, :, None] / 1000.0
    channel = np.arange(256, dtype=np.float32)[None, None, :] / 10000.0
    frames = base + residue + channel
    # Apply the exact frozen definitions to one block without allocating a
    # synthetic 1,000-frame trajectory.
    mu_z = frames.mean(axis=0, dtype=np.float64).astype(np.float32)[None, ...]
    delta = np.diff(frames, axis=0)
    signed_drift = delta.mean(axis=0, dtype=np.float64).astype(np.float32)[None, ...]
    rms_delta = np.sqrt(np.mean(np.square(delta, dtype=np.float64), axis=0)).astype(np.float32)[None, ...]
    state_motion = np.concatenate((mu_z, rms_delta), axis=-1).astype(np.float32, copy=False)
    checks = {}
    for branch, descriptor in (("STATE_MOTION", state_motion), ("SIGNED_DRIFT", signed_drift)):
        state = load_frozen_state(anchor, branch)
        mapped = rff_frozen(normalize_frozen(descriptor[:1], state), state)
        graph_value = diffuse(mapped, transition)
        regional = graph_value[:, next(iter(regions.values())), :].mean(axis=1)
        synthetic = {alias: regional * float(i + 1) for i, alias in enumerate(CONTEXTS)}
        values = {name: historical_contrast(synthetic, name) for name in CONTRASTS}
        if mapped.shape != (1, N_RESIDUES, RFF_FEATURES) or not all(np.isfinite(x).all() for x in values.values()):
            raise RuntimeError(f"{branch} smoke path failed")
        checks[branch] = {"descriptor_shape": list(descriptor.shape), "mapped_shape": list(mapped.shape),
                          "regional_shape": list(regional.shape), "finite": True}
    return {"status": "PHASE3_SMOKE_PASS", "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256,
            "final_biological_interpretation_written": False, "checks": checks}


def run() -> dict[str, Any]:
    anchor = verify_frozen_anchor()
    if PHASE3_ROOT.exists():
        raise FileExistsError(f"refusing silent overwrite of complete or partial Phase-3 root: {PHASE3_ROOT}")
    stage = safe_path(RESULT_ROOT / f".phase3.partial.{os.getpid()}")
    if stage.exists():
        raise FileExistsError(f"staging path already exists: {stage}")
    stage.mkdir(parents=True, exist_ok=False)
    try:
        graph = load_json(GRAPH_PATH)
        transition = graph_transition(graph)
        regions = region_indices(graph)
        branch_results = {}
        manifests = []
        for branch in BRANCHES:
            result, items = evaluate_branch(branch, anchor, graph, transition, regions, stage)
            branch_results[branch] = result
            manifests.extend(items)
        atomic_json(stage / "PHASE3_METHOD_AUDIT.json", method_audit(anchor, graph))
        atomic_json(stage / "STATE_MOTION_RESULTS.json", branch_results["STATE_MOTION"])
        atomic_json(stage / "SIGNED_DRIFT_RESULTS.json", branch_results["SIGNED_DRIFT"])
        atomic_json(stage / "REPLICA_CONSISTENCY.json", consistency_report(branch_results))
        atomic_json(stage / "BLOCK_LEVEL_MANIFEST.json", {
            "schema": "pacer_fkg_v02.phase3_block_manifest.v1",
            "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256, "entries": manifests,
            "expected_entries": 6, "blocks_per_trajectory": N_BLOCKS,
            "trace": "Phase-1 cache -> block index -> frozen branch descriptor -> frozen normalization -> frozen RFF -> graph diffusion -> region -> contrast",
        })
        provenance_dir = stage / "provenance"
        atomic_json(provenance_dir / "FROZEN_INPUT_PROVENANCE.json", {
            "schema": "pacer_fkg_v02.phase3_provenance.v1",
            "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256,
            "freeze_manifest": artifact(FREEZE_PATH), "graph": artifact(GRAPH_PATH),
            "regions": artifact(REGION_PATH), "contrast_definitions": artifact(CONTRAST_PATH),
            "v01_protection_manifest": artifact(PROTECTION_PATH),
            "phase2_verification": anchor["phase2_verification"],
        })
        report_path = stage / "REPORT_PHASE3.md"
        report_path.write_text(report_markdown(branch_results), encoding="utf-8")
        files = sorted(path for path in stage.rglob("*") if path.is_file())
        integrity = {
            "schema": "pacer_fkg_v02.phase3_integrity.v1",
            "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256, "created_at": utc_now(),
            "files": [{"relative_path": path.relative_to(stage).as_posix(),
                       "sha256": sha256(path), "bytes": path.stat().st_size} for path in files],
            "file_count_excluding_self": len(files), "complete": True,
        }
        atomic_json(stage / "PHASE3_INTEGRITY.json", integrity)
        os.replace(stage, safe_path(PHASE3_ROOT))
    except Exception:
        # Preserve failed staging data for diagnosis; never silently replace it.
        raise
    return {"status": "PHASE3_COMPLETE", "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256,
            "output_root": str(PHASE3_ROOT), "branches": list(BRANCHES),
            "replicas": list(REPLICAS), "blocks_per_trajectory": N_BLOCKS}


def verify_results() -> dict[str, Any]:
    anchor = verify_frozen_anchor()
    if not PHASE3_ROOT.is_dir():
        raise FileNotFoundError(f"Phase-3 results do not exist: {PHASE3_ROOT}")
    integrity_path = PHASE3_ROOT / "PHASE3_INTEGRITY.json"
    integrity = load_json(integrity_path)
    if integrity.get("freeze_manifest_sha256") != REQUIRED_FREEZE_SHA256 or integrity.get("complete") is not True:
        raise RuntimeError("invalid Phase-3 integrity header")
    expected = {item["relative_path"]: item for item in integrity["files"]}
    actual = {path.relative_to(PHASE3_ROOT).as_posix(): path for path in PHASE3_ROOT.rglob("*")
              if path.is_file() and path != integrity_path}
    if set(actual) != set(expected):
        raise RuntimeError("Phase-3 file inventory mismatch")
    for relative, path in actual.items():
        record = expected[relative]
        if sha256(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
            raise RuntimeError(f"Phase-3 result hash mismatch: {relative}")
    for name in ("PHASE3_METHOD_AUDIT.json", "STATE_MOTION_RESULTS.json", "SIGNED_DRIFT_RESULTS.json",
                 "REPLICA_CONSISTENCY.json", "BLOCK_LEVEL_MANIFEST.json"):
        if load_json(PHASE3_ROOT / name).get("freeze_manifest_sha256") != REQUIRED_FREEZE_SHA256:
            raise RuntimeError(f"top-level freeze anchor missing or wrong: {name}")
    if REQUIRED_FREEZE_SHA256 not in (PHASE3_ROOT / "REPORT_PHASE3.md").read_text(encoding="utf-8"):
        raise RuntimeError("report freeze anchor missing")
    block_manifest = load_json(PHASE3_ROOT / "BLOCK_LEVEL_MANIFEST.json")
    if len(block_manifest["entries"]) != 6 or any(item["block_count"] != N_BLOCKS for item in block_manifest["entries"]):
        raise RuntimeError("block-level manifest completeness failure")
    branch_results = {
        "STATE_MOTION": load_json(PHASE3_ROOT / "STATE_MOTION_RESULTS.json"),
        "SIGNED_DRIFT": load_json(PHASE3_ROOT / "SIGNED_DRIFT_RESULTS.json"),
    }
    # Replay formulas, aggregate vectors, magnitudes, and descriptive block
    # summaries directly from the immutable serialized regional vectors.
    replayed = 0
    for item in block_manifest["entries"]:
        vector_path = REPO / item["vectors"]["path"]
        values_path = REPO / item["values"]["path"]
        block_values = load_json(values_path)
        if block_values["block_count"] != N_BLOCKS or block_values["replica"] != item["replica"]:
            raise RuntimeError("block-value identity/count mismatch")
        row_index = {(row["kind"], row["name"], row["region"]): row for row in block_values["rows"]}
        result = branch_results[item["branch"]]
        replica_key = f"R{item['replica']}"
        with np.load(vector_path, allow_pickle=False) as arrays:
            regions = sorted({key.split("__", 2)[2] for key in arrays if key.startswith("context__")})
            for region in regions:
                contexts = {alias: arrays[f"context__{alias}__{region}"] for alias in CONTEXTS}
                for alias, vectors in contexts.items():
                    aggregate = result["contexts"][alias]["regions"][region]["replicas"][replica_key]
                    expected = _context_replica_result(vectors, f"context__{alias}__{region}")
                    if aggregate != expected | {"replica_label": REPLICA_LABELS[item["replica"]]}:
                        raise RuntimeError(f"context aggregate replay mismatch: {item['branch']} {replica_key} {alias} {region}")
                    magnitudes = np.linalg.norm(vectors.astype(np.float64), axis=1).tolist()
                    if row_index[("context", alias, region)]["block_values"] != magnitudes:
                        raise RuntimeError(f"context block replay mismatch: {item['branch']} {replica_key} {alias} {region}")
                for contrast in CONTRASTS:
                    stored = arrays[f"contrast__{contrast}__{region}"]
                    if not np.allclose(historical_contrast(contexts, contrast), stored,
                                       rtol=2e-6, atol=2e-7):
                        raise RuntimeError(f"contrast replay mismatch: {item['branch']} R{item['replica']} {contrast} {region}")
                    aggregate = result["contrasts"][contrast]["regions"][region]["replicas"][replica_key]
                    expected = _replica_result(stored, f"contrast__{contrast}__{region}")
                    if aggregate != expected | {"replica_label": REPLICA_LABELS[item["replica"]]}:
                        raise RuntimeError(f"contrast aggregate replay mismatch: {item['branch']} {replica_key} {contrast} {region}")
                    magnitudes = np.linalg.norm(stored.astype(np.float64), axis=1).tolist()
                    if row_index[("contrast", contrast, region)]["block_values"] != magnitudes:
                        raise RuntimeError(f"contrast block replay mismatch: {item['branch']} {replica_key} {contrast} {region}")
                    replayed += 1
    return {"status": "PHASE3_VERIFIED", "freeze_manifest_sha256": REQUIRED_FREEZE_SHA256,
            "phase3_files_verified_excluding_integrity": len(actual), "contrasts_replayed": replayed,
            "phase2_verification": anchor["phase2_verification"]}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--smoke", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--verify", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.smoke:
        result = smoke()
    elif args.run:
        result = run()
    else:
        result = verify_results()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
