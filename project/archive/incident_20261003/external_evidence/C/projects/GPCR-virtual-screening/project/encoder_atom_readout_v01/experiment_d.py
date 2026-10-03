#!/usr/bin/env python
"""Run the fixed PACER-DC Experiment D v01 atom-readout audit once."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch

WORKSPACE = Path(__file__).resolve().parents[2]
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))

from project.encoder_intermediate_v01.core import (  # noqa: E402
    EXPECTED_CHECKPOINT_SHA256,
    PINNED_GEOM2VEC_COMMIT,
    forward_three_states,
    sha256,
    verify_and_strict_load,
)
from project.encoder_intermediate_v01.phase3 import (  # noqa: E402
    dihedral,
    quantiles,
    runtime_report,
    torsion_targets,
    wrapped_angle_error,
)
from project.encoder_intermediate_v01.run_phase_0_2 import (  # noqa: E402
    DEFAULT_AUDIT,
    DEFAULT_CHECKPOINT,
    DEFAULT_INPUT,
    read_sequence,
    source_verification,
)
from project.pacer_dc_training.extract_geom2vec_atom14 import (  # noqa: E402
    ATOM14,
    flatten_atom14,
)

VERSION = "encoder_atom_readout_D_v01"
READOUTS = ("C0-M", "C1-M", "C1-B", "C1-S", "C1-BS", "C1-MC")
CANDIDATE_READOUTS = ("C1-M", "C1-B", "C1-S", "C1-BS", "C1-MC")
STORED_READOUTS = ("M", "B", "S")
BACKBONE_ATOMS = ("N", "CA", "C", "O")
CHI1_GAMMA = {
    "R": "CG",
    "N": "CG",
    "D": "CG",
    "C": "SG",
    "Q": "CG",
    "E": "CG",
    "H": "CG",
    "I": "CG1",
    "L": "CG",
    "K": "CG",
    "M": "CG",
    "F": "CG",
    "P": "CG",
    "S": "OG",
    "T": "OG1",
    "W": "CG",
    "Y": "CG",
    "V": "CG1",
}
CHI1_QUADRUPLETS = {aa: ("N", "CA", "CB", gamma) for aa, gamma in CHI1_GAMMA.items()}

EXPERIMENT_C_ROOT = WORKSPACE / "project/results/encoder_intermediate_C_v01/phase3_run_001"
EXPERIMENT_C_PHASE02 = WORKSPACE / "project/results/encoder_intermediate_C_v01/run_001"
DEFAULT_OUTPUT = WORKSPACE / "project/results/encoder_atom_readout_D_v01/run_001"
DEFAULT_SOURCE = Path(r"C:\projects\geom2vec-source")
SMOKE_RELATIVE = "replica_02/window_000/atom14/apo_w000.npy"
SMOKE_FRAMES = np.asarray([0, 1, 2, 3], dtype=np.int64)
REPLAY_SPECS = (
    ("replica_02/window_000/atom14/apo_w000.npy", (0, 1, 2, 3)),
    ("replica_03/window_000/atom14/apo_w000.npy", (0, 1, 2, 3)),
)
PHASE02_IMPLEMENTATION = WORKSPACE / "project/encoder_intermediate_v01/core.py"

FROZEN_CONFIG = {
    "version": VERSION,
    "candidate_atomic_state": "C1: accumulated x/v entering vis_mp_layers[3]",
    "atomic_descriptor": "concat(x, Euclidean vector-channel norm over xyz), 128D",
    "readouts": {
        "C0-M": "Experiment C production C0 all-heavy mean reference, 128D",
        "C1-M": "all expected heavy atoms mean, 128D",
        "C1-B": "N/CA/C/O mean, 128D",
        "C1-S": "all expected non-N/CA/C/O heavy atoms mean; exact zero for no-sidechain residues, 128D",
        "C1-BS": "concat(C1-B,C1-S), derived during analysis, 256D",
        "C1-MC": "concat(C1-M,exact zeros), derived during analysis, 256D",
    },
    "backbone_atoms": list(BACKBONE_ATOMS),
    "chi1_quadruplets": {aa: list(atoms) for aa, atoms in CHI1_QUADRUPLETS.items()},
    "chi1_excluded_residues": ["A", "G"],
    "extraction_batch_size": 4,
    "stored_cache_arrays": ["M", "B", "S", "sidechain_present", "frame_ids", "residue_index", "sequence"],
    "probe": {
        "family": "linear ridge",
        "alpha": 1.0,
        "regularize_intercept": False,
        "scaling": "per-readout mean/std fitted on R2 backbone-valid union; scale <= 1e-12 replaced by 1",
        "fit": "R2 only",
        "evaluation": "R3 descriptive only",
        "backbone_targets": ["sin(phi)", "cos(phi)", "sin(psi)", "cos(psi)"],
        "sidechain_target": ["sin(chi1)", "cos(chi1)"],
        "aggregate": "residue-balanced mean; phi/psi combined gives phi and psi equal weight",
    },
    "near_constant_rule": "variance <= max(1e-12, median(channel variance)*1e-8)",
    "condition_rule": "positive-spectrum lambda_max/lambda_min above lambda_max*1e-12; also lambda_max*1e-6 regularized surrogate",
    "smoke": {"relative_path": SMOKE_RELATIVE, "frame_ids": SMOKE_FRAMES.tolist()},
    "replay": [{"relative_path": path, "frame_ids": list(frames)} for path, frames in REPLAY_SPECS],
    "success_logic": {
        "accessibility": "C1-BS has lower residue-balanced error than both C1-M and C1-MC on backbone combined angular MAE and chi1 angular MAE, on R2 and R3",
        "not_one_residue": "C1-BS improves per-residue angular MAE for >50% of evaluable residues versus both controls, for both target families and splits",
        "rank": "C1-BS participation rank >= 50% of C1-M on both splits",
        "scaling": "C1-BS has no R2 scaler-collapse channels",
        "numerics": "finite; replay relative RMSE < 1e-4; exact BS blocks and MC zero block",
        "interpretation": "multi-criterion descriptive flag; no scalar winner",
    },
}


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def file_record(path: Path) -> dict:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)}


def parse_metadata(relative: str) -> tuple[int, str, int, str]:
    path = Path(relative)
    replica = int(next(part for part in path.parts if part.startswith("replica_")).split("_")[1])
    split = f"R{replica}"
    window = int(next(part for part in path.parts if part.startswith("window_")).split("_")[1])
    suffix = f"_w{window:03d}"
    condition = path.stem[: -len(suffix)] if path.stem.endswith(suffix) else path.stem
    return replica, split, window, condition


def runtime_identity(runtime: dict) -> dict:
    return {
        key: runtime[key]
        for key in ("python", "python_executable", "platform", "pytorch", "pyg", "numpy", "cuda_build", "cudnn", "device", "device_name", "neighbor_backend")
    }


def verify_provenance(
    checkpoint: Path,
    audit_path: Path,
    input_root: Path,
    source_root: Path,
) -> tuple[torch.nn.Module, list[dict], dict]:
    c_completion = json.loads((EXPERIMENT_C_ROOT / "PHASE3_COMPLETE.json").read_text(encoding="utf-8"))
    c_extraction = json.loads((EXPERIMENT_C_ROOT / "extraction_provenance.json").read_text(encoding="utf-8"))
    c_frozen = json.loads((EXPERIMENT_C_ROOT / "frozen_phase3_config.json").read_text(encoding="utf-8"))
    c_integrity_path = EXPERIMENT_C_ROOT / "integrity_completion_audit.json"
    c_report_path = EXPERIMENT_C_ROOT / "REPORT_PHASE3.md"
    completion_hashes_ok = (
        sha256(c_integrity_path) == c_completion["integrity_audit_sha256"]
        and sha256(c_report_path) == c_completion["report"]["sha256"]
    )
    if not completion_hashes_ok:
        raise RuntimeError("Experiment C completion artifacts changed")

    expected_core = next(
        item for item in c_frozen["phase02_implementation"] if Path(item["path"]).name == "core.py"
    )
    c1_implementation_ok = sha256(PHASE02_IMPLEMENTATION) == expected_core["sha256"]
    if not c1_implementation_ok:
        raise RuntimeError("validated C1 hook implementation changed")

    model, checkpoint_report = verify_and_strict_load(checkpoint)
    if not checkpoint_report["complete_strict_coverage"]:
        raise RuntimeError("checkpoint strict coverage no longer passes")
    source_report = source_verification(source_root)
    if not source_report["source_equivalent"]:
        raise RuntimeError("pinned Geom2Vec source provenance changed")

    current_runtime = runtime_report("cuda")
    runtime_ok = runtime_identity(current_runtime) == runtime_identity(c_extraction["runtime"])
    if not runtime_ok:
        raise RuntimeError("runtime differs from Experiment C Phase 3")
    if c_extraction["checkpoint_sha256"] != EXPECTED_CHECKPOINT_SHA256 or sha256(checkpoint) != EXPECTED_CHECKPOINT_SHA256:
        raise RuntimeError("checkpoint hash differs from Experiment C")
    if c_extraction["source_commit"] != PINNED_GEOM2VEC_COMMIT:
        raise RuntimeError("Experiment C source commit provenance differs")

    archived = json.loads(audit_path.read_text(encoding="utf-8"))
    c_cache_by_relative = {item["relative_path"]: item for item in c_extraction["rows"]}
    rows = []
    for item in archived["rows"]:
        source = input_root / item["input"]
        source_hash = sha256(source)
        if source_hash != item["input_sha256"]:
            raise RuntimeError(f"authenticated atom14 input changed: {source}")
        c_cache = Path(c_cache_by_relative[item["input"]]["cache_path"])
        c_cache_hash = sha256(c_cache)
        if c_cache_hash != c_cache_by_relative[item["input"]]["cache_sha256"]:
            raise RuntimeError(f"Experiment C reference cache changed: {c_cache}")
        replica, split, window, condition = parse_metadata(item["input"])
        rows.append(
            {
                "relative_path": item["input"],
                "source": source,
                "source_sha256": source_hash,
                "sequence_csv": source.with_suffix(".csv"),
                "sequence_sha256": sha256(source.with_suffix(".csv")),
                "frames": int(item["frames"]),
                "replica": replica,
                "split": split,
                "window": window,
                "condition": condition,
                "c_cache": c_cache,
                "c_cache_sha256": c_cache_hash,
            }
        )
    if len(rows) != 40 or {row["split"] for row in rows} != {"R2", "R3"}:
        raise RuntimeError("authenticated inventory is not the expected 40-file R2/R3 set")
    report = {
        "checkpoint": checkpoint_report,
        "c1_implementation": {**file_record(PHASE02_IMPLEMENTATION), "matches_experiment_c": c1_implementation_ok},
        "source": source_report,
        "runtime": current_runtime,
        "runtime_matches_experiment_c": runtime_ok,
        "authenticated_inputs": len(rows),
        "experiment_c_cache_files_verified": len(rows),
        "experiment_c_completion_hashes_verified": completion_hashes_ok,
        "experiment_c_completion": file_record(EXPERIMENT_C_ROOT / "PHASE3_COMPLETE.json"),
        "experiment_c_report": file_record(c_report_path),
        "experiment_c_integrity": file_record(c_integrity_path),
    }
    return model, rows, report


def atom_metadata(sequence: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict]]:
    atom_names = []
    residue_index = []
    records = []
    for residue, aa in enumerate(sequence):
        names = ATOM14[aa]
        backbone = [name for name in names if name in BACKBONE_ATOMS]
        sidechain = [name for name in names if name not in BACKBONE_ATOMS]
        if backbone != list(BACKBONE_ATOMS):
            raise RuntimeError(f"residue {residue} {aa} does not have exact N/CA/C/O backbone ordering")
        atom_names.extend(names)
        residue_index.extend([residue] * len(names))
        records.append(
            {
                "residue_index": residue,
                "residue_type": aa,
                "all_heavy_count": len(names),
                "backbone_count": len(backbone),
                "sidechain_count": len(sidechain),
                "sidechain_present": bool(sidechain),
                "backbone_atoms": backbone,
                "sidechain_atoms": sidechain,
            }
        )
    names = np.asarray(atom_names)
    return (
        np.asarray(residue_index, dtype=np.int64),
        np.isin(names, BACKBONE_ATOMS),
        ~np.isin(names, BACKBONE_ATOMS),
        records,
    )


def prepare_batch(xyz: np.ndarray, atom_z: np.ndarray, device: str):
    pos = torch.as_tensor(xyz, dtype=torch.float32, device=device)
    batch_size, n_atoms, _ = pos.shape
    z0 = torch.as_tensor(atom_z, dtype=torch.long, device=device)
    z = z0.repeat(batch_size)
    batch = torch.arange(batch_size, device=device).repeat_interleave(n_atoms)
    return z, pos.reshape(-1, 3), batch, batch_size, n_atoms


def grouped_mean(
    descriptor: torch.Tensor,
    residue_index: torch.Tensor,
    atom_mask: torch.Tensor,
    n_residues: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    selected = descriptor[atom_mask]
    selected_residue = residue_index[atom_mask]
    output = torch.zeros((n_residues, descriptor.shape[-1]), dtype=descriptor.dtype, device=descriptor.device)
    counts = torch.zeros((n_residues, 1), dtype=descriptor.dtype, device=descriptor.device)
    output.index_add_(0, selected_residue, selected)
    counts.index_add_(0, selected_residue, torch.ones((len(selected_residue), 1), device=descriptor.device))
    return output / counts.clamp_min(1), counts.squeeze(1)


def readouts_from_c1_atomic(
    x: torch.Tensor,
    v: torch.Tensor,
    residue_index: torch.Tensor,
    backbone_mask: torch.Tensor,
    sidechain_mask: torch.Tensor,
    batch_size: int,
    n_atoms: int,
    n_residues: int,
) -> dict[str, torch.Tensor]:
    x = x.reshape(batch_size, n_atoms, 64)
    v = v.reshape(batch_size, n_atoms, 3, 64)
    outputs = {key: [] for key in STORED_READOUTS}
    for frame in range(batch_size):
        descriptor = torch.cat((x[frame], torch.linalg.vector_norm(v[frame], dim=1)), dim=-1)
        mean, _ = grouped_mean(descriptor, residue_index, torch.ones(n_atoms, dtype=torch.bool, device=x.device), n_residues)
        backbone, _ = grouped_mean(descriptor, residue_index, backbone_mask, n_residues)
        sidechain, _ = grouped_mean(descriptor, residue_index, sidechain_mask, n_residues)
        outputs["M"].append(mean)
        outputs["B"].append(backbone)
        outputs["S"].append(sidechain)
    return {key: torch.stack(value) for key, value in outputs.items()}


def extract_readouts(
    model,
    xyz,
    atom_z,
    residue_index,
    backbone_mask,
    sidechain_mask,
    n_residues,
    device,
    batch_size,
) -> dict[str, np.ndarray]:
    pieces = {key: [] for key in STORED_READOUTS}
    ridx = torch.as_tensor(residue_index, dtype=torch.long, device=device)
    bmask = torch.as_tensor(backbone_mask, dtype=torch.bool, device=device)
    smask = torch.as_tensor(sidechain_mask, dtype=torch.bool, device=device)
    with torch.inference_mode():
        for start in range(0, len(xyz), batch_size):
            chunk = xyz[start : start + batch_size]
            z, pos, batch, actual_batch, n_atoms = prepare_batch(chunk, atom_z, device)
            states, _ = forward_three_states(model, z, pos, batch)
            pooled = readouts_from_c1_atomic(
                states["C1"][0], states["C1"][1], ridx, bmask, smask,
                actual_batch, n_atoms, n_residues,
            )
            for key in STORED_READOUTS:
                pieces[key].append(pooled[key].cpu().numpy().astype(np.float32, copy=False))
    return {key: np.concatenate(value, axis=0) for key, value in pieces.items()}


def error_metrics(observed: np.ndarray, reference: np.ndarray) -> dict:
    delta = observed.astype(np.float64) - reference.astype(np.float64)
    reference_rms = math.sqrt(float(np.mean(np.square(reference.astype(np.float64)))))
    rmse = math.sqrt(float(np.mean(np.square(delta))))
    return {
        "max_abs": float(np.abs(delta).max()),
        "rmse": rmse,
        "relative_rmse": rmse / reference_rms,
    }


def derive_readouts(stored: dict[str, np.ndarray], c0: np.ndarray | None = None) -> dict[str, np.ndarray]:
    zero = np.zeros_like(stored["M"])
    result = {
        "C1-M": stored["M"],
        "C1-B": stored["B"],
        "C1-S": stored["S"],
        "C1-BS": np.concatenate((stored["B"], stored["S"]), axis=-1),
        "C1-MC": np.concatenate((stored["M"], zero), axis=-1),
    }
    if c0 is not None:
        result["C0-M"] = c0
    return result


def smoke_test(model, row: dict, device: str) -> dict:
    sequence = read_sequence(row["sequence_csv"])
    coords = np.load(row["source"], mmap_mode="r")
    selected = np.asarray(coords[SMOKE_FRAMES])
    xyz, atom_z, existing_residue_index = flatten_atom14(selected, sequence)
    residue_index, backbone_mask, sidechain_mask, group_records = atom_metadata(sequence)
    if not np.array_equal(existing_residue_index, residue_index):
        raise RuntimeError("atom order differs from validated flatten_atom14 order")
    model = model.to(device).eval()
    first = extract_readouts(model, xyz, atom_z, residue_index, backbone_mask, sidechain_mask, len(sequence), device, len(SMOKE_FRAMES))
    replay = extract_readouts(model, xyz, atom_z, residue_index, backbone_mask, sidechain_mask, len(sequence), device, len(SMOKE_FRAMES))
    with np.load(row["c_cache"]) as archive:
        c_reference = archive["C1"][SMOKE_FRAMES]

    angle = np.deg2rad(37.0)
    axis = np.asarray([1.0, 2.0, -1.0], dtype=np.float64)
    axis /= np.linalg.norm(axis)
    skew = np.asarray([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rotation = np.eye(3) * np.cos(angle) + (1 - np.cos(angle)) * np.outer(axis, axis) + np.sin(angle) * skew
    translation = np.asarray([11.0, -7.0, 3.5])
    rotated = extract_readouts(model, (xyz.astype(np.float64) @ rotation.T).astype(np.float32), atom_z, residue_index, backbone_mask, sidechain_mask, len(sequence), device, len(SMOKE_FRAMES))
    translated = extract_readouts(model, (xyz.astype(np.float64) + translation).astype(np.float32), atom_z, residue_index, backbone_mask, sidechain_mask, len(sequence), device, len(SMOKE_FRAMES))
    derived = derive_readouts(first)
    finite = {name: bool(np.isfinite(value).all()) for name, value in derived.items()}
    shapes = {name: list(value.shape) for name, value in derived.items()}
    exact_blocks = {
        "BS_first_equals_B": bool(np.array_equal(derived["C1-BS"][..., :128], first["B"])),
        "BS_second_equals_S": bool(np.array_equal(derived["C1-BS"][..., 128:], first["S"])),
        "MC_first_equals_M": bool(np.array_equal(derived["C1-MC"][..., :128], first["M"])),
        "MC_second_exact_zero": bool(np.count_nonzero(derived["C1-MC"][..., 128:]) == 0),
    }
    report = {
        "input": str(row["source"].resolve()),
        "input_sha256": row["source_sha256"],
        "frame_ids": SMOKE_FRAMES.tolist(),
        "residue_order": [0, len(sequence) - 1],
        "atom_count": len(atom_z),
        "group_counts": {
            "all_heavy": int(len(atom_z)),
            "backbone": int(backbone_mask.sum()),
            "sidechain": int(sidechain_mask.sum()),
            "sidechain_present_residues": int(sum(record["sidechain_present"] for record in group_records)),
            "sidechain_absent_residues": int(sum(not record["sidechain_present"] for record in group_records)),
            "per_residue": group_records,
        },
        "shapes": shapes,
        "finite": finite,
        "exact_block_checks": exact_blocks,
        "M_vs_experiment_C_C1": error_metrics(first["M"], c_reference),
        "deterministic_replay": {key: error_metrics(replay[key], first[key]) for key in STORED_READOUTS},
        "rotation_invariance": {key: error_metrics(rotated[key], first[key]) for key in STORED_READOUTS},
        "translation_invariance": {key: error_metrics(translated[key], first[key]) for key in STORED_READOUTS},
        "rotation_matrix": rotation.tolist(),
        "translation": translation.tolist(),
    }
    report["passed"] = (
        all(finite.values())
        and shapes == {"C1-M": [4, 270, 128], "C1-B": [4, 270, 128], "C1-S": [4, 270, 128], "C1-BS": [4, 270, 256], "C1-MC": [4, 270, 256]}
        and all(exact_blocks.values())
        and report["M_vs_experiment_C_C1"]["max_abs"] <= 5e-5
        and max(value["relative_rmse"] for value in report["deterministic_replay"].values()) < 1e-4
        and max(value["relative_rmse"] for value in report["rotation_invariance"].values()) < 1e-4
        and max(value["relative_rmse"] for value in report["translation_invariance"].values()) < 1e-4
    )
    return report


def disk_preflight(rows: list[dict], output_root: Path) -> dict:
    frames = sum(row["frames"] for row in rows)
    raw_bytes = frames * 270 * 128 * 4 * len(STORED_READOUTS)
    free = shutil.disk_usage(output_root.parent).free
    required = math.ceil(raw_bytes * 1.10)
    result = {
        "files": len(rows),
        "frames": frames,
        "stored_arrays": list(STORED_READOUTS),
        "expected_raw_cache_bytes": raw_bytes,
        "expected_raw_cache_gib": raw_bytes / 2**30,
        "required_with_10_percent_overhead_bytes": required,
        "free_bytes": free,
        "free_gib": free / 2**30,
        "sufficient": free >= required,
    }
    if not result["sufficient"]:
        raise RuntimeError("insufficient disk space for Experiment D cache")
    return result


def extract_all(model, rows: list[dict], output_root: Path, device: str) -> dict:
    started_all = time.time()
    audit_rows = []
    cache_root = output_root / "cache"
    model = model.to(device).eval()
    for number, row in enumerate(rows, 1):
        coords = np.load(row["source"], mmap_mode="r")
        sequence = read_sequence(row["sequence_csv"])
        xyz, atom_z, existing_residue_index = flatten_atom14(np.asarray(coords), sequence)
        residue_index, backbone_mask, sidechain_mask, group_records = atom_metadata(sequence)
        if not np.array_equal(existing_residue_index, residue_index):
            raise RuntimeError("atom metadata does not match flatten_atom14")
        started = time.time()
        arrays = extract_readouts(
            model, xyz, atom_z, residue_index, backbone_mask, sidechain_mask,
            len(sequence), device, FROZEN_CONFIG["extraction_batch_size"],
        )
        sidechain_present = np.asarray([record["sidechain_present"] for record in group_records], dtype=np.bool_)
        relative_cache = Path(row["relative_path"]).with_suffix(".C1_MBS.npz")
        cache_path = cache_root / relative_cache
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            cache_path,
            M=arrays["M"], B=arrays["B"], S=arrays["S"],
            sidechain_present=sidechain_present,
            frame_ids=np.arange(len(coords), dtype=np.int64),
            residue_index=np.arange(len(sequence), dtype=np.int64),
            sequence=np.asarray(sequence),
        )
        for key in STORED_READOUTS:
            if arrays[key].shape != (row["frames"], 270, 128) or not np.isfinite(arrays[key]).all():
                raise RuntimeError(f"invalid {key} extraction for {row['source']}")
        if np.count_nonzero(arrays["S"][:, ~sidechain_present]) != 0:
            raise RuntimeError("sidechain-absent residues are not exact zero")
        row["cache"] = cache_path
        record = {
            "number": number,
            "relative_path": row["relative_path"],
            "source_path": str(row["source"].resolve()),
            "source_sha256": row["source_sha256"],
            "sequence_sha256": row["sequence_sha256"],
            "replica": row["replica"],
            "split": row["split"],
            "condition": row["condition"],
            "window": row["window"],
            "shape_each_stored_readout": [row["frames"], 270, 128],
            "stored_readouts": list(STORED_READOUTS),
            "sidechain_present_residues": int(sidechain_present.sum()),
            "finite": True,
            "cache_path": str(cache_path.resolve()),
            "cache_bytes": cache_path.stat().st_size,
            "cache_sha256": sha256(cache_path),
            "experiment_c_cache_sha256": row["c_cache_sha256"],
            "elapsed_seconds": round(time.time() - started, 3),
        }
        audit_rows.append(record)
        print(f"[{number:02d}/{len(rows)}] {row['relative_path']} -> {record['cache_bytes']} bytes", flush=True)
    return {
        "method": "one C1 atomic forward per batch; M/B/S derived together; BS/MC not cached",
        "files_expected": len(rows),
        "files_ok": len(audit_rows),
        "frames": sum(row["frames"] for row in rows),
        "elapsed_seconds": round(time.time() - started_all, 3),
        "cache_bytes": sum(row["cache_bytes"] for row in audit_rows),
        "cache_gib": sum(row["cache_bytes"] for row in audit_rows) / 2**30,
        "rows": audit_rows,
    }


def load_representations(row: dict) -> tuple[dict[str, np.ndarray], np.ndarray]:
    with np.load(row["cache"]) as archive:
        stored = {key: archive[key] for key in STORED_READOUTS}
        sidechain_present = archive["sidechain_present"]
    with np.load(row["c_cache"]) as archive:
        c0 = archive["C0"]
    return derive_readouts(stored, c0), sidechain_present


def health_for_chunks(chunks: list[np.ndarray], state: str, split: str, subset: str = "all") -> tuple[dict, list[dict], list[dict]]:
    dimension = chunks[0].shape[-1]
    total_samples = 0
    finite_elements = 0
    total_elements = 0
    sums = np.zeros(dimension, dtype=np.float64)
    sums2 = np.zeros(dimension, dtype=np.float64)
    cross = np.zeros((dimension, dimension), dtype=np.float64)
    norms = []
    absolute = []
    for chunk32 in chunks:
        flat32 = chunk32.reshape(-1, dimension)
        finite_elements += int(np.isfinite(flat32).sum())
        total_elements += flat32.size
        flat = flat32.astype(np.float64)
        total_samples += len(flat)
        sums += flat.sum(axis=0)
        sums2 += np.square(flat).sum(axis=0)
        cross += flat.T @ flat
        norms.append(np.linalg.norm(flat32, axis=1))
        absolute.append(np.abs(flat32).reshape(-1))
    mean = sums / total_samples
    variance = np.maximum(sums2 / total_samples - mean**2, 0)
    covariance = cross / total_samples - np.outer(mean, mean)
    covariance = (covariance + covariance.T) * 0.5
    eigenvalues = np.maximum(np.linalg.eigvalsh(covariance)[::-1], 0)
    eigen_sum = eigenvalues.sum()
    weights = eigenvalues / eigen_sum if eigen_sum else np.zeros_like(eigenvalues)
    participation = eigen_sum**2 / np.square(eigenvalues).sum() if np.square(eigenvalues).sum() else 0.0
    threshold = eigenvalues[0] * 1e-12 if eigenvalues[0] else 0.0
    positive = eigenvalues[eigenvalues > threshold]
    raw_condition = eigenvalues[0] / positive[-1] if len(positive) else math.inf
    epsilon = eigenvalues[0] * 1e-6 if eigenvalues[0] else 1e-12
    regularized_condition = (eigenvalues[0] + epsilon) / (eigenvalues[-1] + epsilon)
    near_threshold = max(1e-12, float(np.median(variance)) * 1e-8)
    near = variance <= near_threshold
    result = {
        "readout": state,
        "split": split,
        "subset": subset,
        "dimension": dimension,
        "observations": total_samples,
        "elements": total_elements,
        "finite_elements": finite_elements,
        "nonfinite_elements": total_elements - finite_elements,
        "feature_rms": math.sqrt(float(sums2.sum() / total_elements)),
        "median_feature_magnitude": float(np.median(np.concatenate(absolute))),
        "feature_norm_distribution": quantiles(np.concatenate(norms)),
        "channel_variance_distribution": quantiles(variance),
        "near_constant_threshold": near_threshold,
        "near_constant_channels": int(near.sum()),
        "covariance_participation_rank": float(participation),
        "covariance_entropy_rank": float(np.exp(-np.sum(weights[weights > 0] * np.log(weights[weights > 0])))),
        "numerical_rank": int(np.linalg.matrix_rank(covariance)),
        "positive_spectrum_eigenvalues": int(len(positive)),
        "positive_spectrum_condition": float(raw_condition),
        "regularized_condition_1e-6": float(regularized_condition),
        "total_256d_rank": int(np.linalg.matrix_rank(covariance)) if dimension == 256 else None,
        "informative_rank": int(np.linalg.matrix_rank(covariance)),
        "zero_block_exact": bool(all(np.count_nonzero(chunk[..., 128:]) == 0 for chunk in chunks)) if state == "C1-MC" else None,
    }
    channel_rows = [
        {
            "readout": state, "split": split, "subset": subset, "channel": channel,
            "mean": float(mean[channel]), "variance": float(variance[channel]),
            "std": float(math.sqrt(variance[channel])), "near_constant": bool(near[channel]),
        }
        for channel in range(dimension)
    ]
    cumulative = np.cumsum(weights)
    spectrum_rows = [
        {
            "readout": state, "split": split, "subset": subset, "eigen_index": index,
            "eigenvalue": float(value), "explained_fraction": float(weights[index]),
            "cumulative_fraction": float(cumulative[index]),
        }
        for index, value in enumerate(eigenvalues)
    ]
    return result, channel_rows, spectrum_rows


def numerical_health(rows: list[dict], output_root: Path) -> dict:
    results = {}
    channel_rows = []
    spectrum_rows = []
    for split in ("R2", "R3"):
        selected = [row for row in rows if row["split"] == split]
        for readout in READOUTS:
            chunks = []
            subset_chunks = []
            for row in selected:
                representations, sidechain_present = load_representations(row)
                chunks.append(representations[readout])
                if readout == "C1-S":
                    subset_chunks.append(representations[readout][:, sidechain_present, :])
            result, channels, spectrum = health_for_chunks(chunks, readout, split)
            results[f"{readout}_{split}"] = result
            channel_rows.extend(channels)
            spectrum_rows.extend(spectrum)
            if readout == "C1-S":
                result, channels, spectrum = health_for_chunks(subset_chunks, readout, split, "sidechain_present")
                results[f"{readout}_{split}_sidechain_present"] = result
                channel_rows.extend(channels)
                spectrum_rows.extend(spectrum)
    for split in ("R2", "R3"):
        results[f"C1-MC_{split}"]["informative_rank"] = results[f"C1-M_{split}"]["numerical_rank"]
        results[f"C1-MC_{split}"]["informative_participation_rank"] = results[f"C1-M_{split}"]["covariance_participation_rank"]
        results[f"C1-BS_{split}"]["informative_participation_rank"] = results[f"C1-BS_{split}"]["covariance_participation_rank"]
    write_json(output_root / "numerical_health.json", results)
    write_csv(output_root / "numerical_health_channels.csv", channel_rows)
    write_csv(output_root / "numerical_health_covariance_spectrum.csv", spectrum_rows)
    return results


def temporal_diagnostics(rows: list[dict], numerical: dict, output_root: Path) -> dict:
    results = {}
    detail = []
    for split in ("R2", "R3"):
        selected = [row for row in rows if row["split"] == split]
        for readout in ("C1-M", "C1-B", "C1-S", "C1-BS"):
            window_centroids = {}
            frame_values = []
            within_values = []
            for row in selected:
                representations, _ = load_representations(row)
                features = representations[readout]
                centroid = features.mean(axis=0, dtype=np.float64)
                trajectory = (row["replica"], row["condition"])
                window_centroids[(trajectory, row["window"])] = centroid
                adjacent = np.sqrt(
                    np.mean(np.square(np.diff(features.astype(np.float64), axis=0)), axis=2)
                )
                frame = float(adjacent.mean())
                within = float(np.sqrt(np.mean(np.square(features.astype(np.float64) - centroid[None]))))
                frame_values.append(frame)
                within_values.append(within)
                detail.append({"readout": readout, "split": split, "condition": row["condition"], "window": row["window"], "frame_displacement": frame, "within_window_rms": within})
            trajectories = sorted({key for key, _window in window_centroids})
            drift_values = []
            for trajectory in trajectories:
                windows = sorted(window for key, window in window_centroids if key == trajectory)
                first = window_centroids[(trajectory, windows[0])]
                last = window_centroids[(trajectory, windows[-1])]
                drift_values.append(float(np.sqrt(np.mean(np.square(last - first)))))
            rms = numerical[f"{readout}_{split}"]["feature_rms"]
            results[f"{readout}_{split}"] = {
                "readout": readout,
                "split": split,
                "windows": len(selected),
                "trajectories": len(trajectories),
                "frame_displacement_rms": float(np.mean(frame_values)),
                "within_window_rms": float(np.mean(within_values)),
                "chronological_first_to_last_rms": float(np.mean(drift_values)),
                "scale_relative_frame_displacement": float(np.mean(frame_values) / rms),
                "scale_relative_within_window": float(np.mean(within_values) / rms),
                "scale_relative_chronological_drift": float(np.mean(drift_values) / rms),
            }
    write_json(output_root / "temporal_diagnostics.json", results)
    write_csv(output_root / "temporal_diagnostics_windows.csv", detail)
    return results


def chi1_targets(coords: np.ndarray, sequence: str) -> dict:
    frames, residues = coords.shape[:2]
    angles = np.full((frames, residues), np.nan, dtype=np.float64)
    defined_residues = []
    for residue, aa in enumerate(sequence):
        quadruplet = CHI1_QUADRUPLETS.get(aa)
        if quadruplet is None:
            continue
        names = ATOM14[aa]
        try:
            slots = [names.index(name) for name in quadruplet]
        except ValueError as exc:
            raise RuntimeError(f"frozen chi1 mapping is unavailable for {aa}: {quadruplet}") from exc
        value, valid = dihedral(*(coords[:, residue, slot, :] for slot in slots))
        angles[:, residue] = np.where(valid, value, np.nan)
        defined_residues.append(residue)
    return {
        "chi1": angles,
        "chi1_mask": np.isfinite(angles),
        "defined_residue_indices": defined_residues,
    }


def r2_scalers(rows: list[dict]) -> dict:
    accum = {}
    for readout in READOUTS:
        dimension = 256 if readout in ("C1-BS", "C1-MC") else 128
        accum[readout] = {"count": 0, "sum": np.zeros(dimension), "sum2": np.zeros(dimension)}
    for row in (row for row in rows if row["split"] == "R2"):
        coords = np.asarray(np.load(row["source"], mmap_mode="r"))
        targets = torsion_targets(coords)
        mask = targets["phi_mask"] | targets["psi_mask"]
        representations, _ = load_representations(row)
        for readout in READOUTS:
            values = representations[readout][mask].astype(np.float64)
            accum[readout]["count"] += len(values)
            accum[readout]["sum"] += values.sum(axis=0)
            accum[readout]["sum2"] += np.square(values).sum(axis=0)
    scalers = {}
    for readout, item in accum.items():
        mean = item["sum"] / item["count"]
        variance = np.maximum(item["sum2"] / item["count"] - mean**2, 0)
        raw_scale = np.sqrt(variance)
        replaced = raw_scale <= 1e-12
        scalers[readout] = {
            "count": item["count"],
            "mean": mean,
            "scale": np.where(replaced, 1.0, raw_scale),
            "replaced_channels": np.flatnonzero(replaced).tolist(),
            "raw_scale_distribution": quantiles(raw_scale),
        }
    return scalers


def fit_models(rows: list[dict], scalers: dict) -> dict:
    torsions = ("phi", "psi", "chi1")
    accum = {}
    for readout in READOUTS:
        dimension = len(scalers[readout]["mean"])
        accum[readout] = {
            torsion: {"xtx": np.zeros((dimension + 1, dimension + 1)), "xty": np.zeros((dimension + 1, 2)), "samples": 0}
            for torsion in torsions
        }
    for row in (row for row in rows if row["split"] == "R2"):
        coords = np.asarray(np.load(row["source"], mmap_mode="r"))
        sequence = read_sequence(row["sequence_csv"])
        targets = {**torsion_targets(coords), **chi1_targets(coords, sequence)}
        representations, _ = load_representations(row)
        for readout in READOUTS:
            standardized = (representations[readout].astype(np.float64) - scalers[readout]["mean"]) / scalers[readout]["scale"]
            for torsion in torsions:
                mask = targets[f"{torsion}_mask"]
                x = standardized[mask]
                angle = targets[torsion][mask]
                y = np.column_stack((np.sin(angle), np.cos(angle)))
                xa = np.column_stack((x, np.ones(len(x))))
                item = accum[readout][torsion]
                item["xtx"] += xa.T @ xa
                item["xty"] += xa.T @ y
                item["samples"] += len(x)
    models = {}
    for readout in READOUTS:
        dimension = len(scalers[readout]["mean"])
        regularizer = np.eye(dimension + 1) * FROZEN_CONFIG["probe"]["alpha"]
        regularizer[-1, -1] = 0
        models[readout] = {}
        for torsion in torsions:
            item = accum[readout][torsion]
            models[readout][torsion] = {
                "samples": item["samples"],
                "coefficient": np.linalg.solve(item["xtx"] + regularizer, item["xty"]),
            }
    return models


def new_metric_accumulator() -> dict:
    return {
        "component_sse": np.zeros(2), "samples": 0,
        "angle_abs": 0.0, "angle_sq": 0.0,
        "per_residue": defaultdict(lambda: {"component_sse": np.zeros(2), "samples": 0, "angle_abs": 0.0, "angle_sq": 0.0}),
    }


def evaluate_models(rows: list[dict], scalers: dict, models: dict, output_root: Path) -> tuple[dict, dict, list[dict], list[dict]]:
    torsions = ("phi", "psi", "chi1")
    accum = {
        (readout, split, torsion, subset): new_metric_accumulator()
        for readout in READOUTS
        for split in ("R2", "R3")
        for torsion in torsions
        for subset in (("all", "sidechain_present") if torsion in ("phi", "psi") else ("valid_chi1",))
    }
    for row in rows:
        coords = np.asarray(np.load(row["source"], mmap_mode="r"))
        sequence = read_sequence(row["sequence_csv"])
        backbone = torsion_targets(coords)
        chi = chi1_targets(coords, sequence)
        targets = {**backbone, **chi}
        representations, sidechain_present = load_representations(row)
        residue_grid = np.broadcast_to(np.arange(270)[None, :], backbone["phi"].shape)
        sidechain_grid = np.broadcast_to(sidechain_present[None, :], backbone["phi"].shape)
        for readout in READOUTS:
            standardized = (representations[readout].astype(np.float64) - scalers[readout]["mean"]) / scalers[readout]["scale"]
            for torsion in torsions:
                base_mask = targets[f"{torsion}_mask"]
                subsets = {"valid_chi1": base_mask} if torsion == "chi1" else {"all": base_mask, "sidechain_present": base_mask & sidechain_grid}
                for subset, mask in subsets.items():
                    x = standardized[mask]
                    observed_angle = targets[torsion][mask]
                    observed = np.column_stack((np.sin(observed_angle), np.cos(observed_angle)))
                    xa = np.column_stack((x, np.ones(len(x))))
                    predicted = xa @ models[readout][torsion]["coefficient"]
                    component_error = predicted - observed
                    angle_error = wrapped_angle_error(np.arctan2(predicted[:, 0], predicted[:, 1]), observed_angle)
                    item = accum[(readout, row["split"], torsion, subset)]
                    item["component_sse"] += np.square(component_error).sum(axis=0)
                    item["samples"] += len(x)
                    item["angle_abs"] += float(angle_error.sum())
                    item["angle_sq"] += float(np.square(angle_error).sum())
                    residues = residue_grid[mask]
                    for residue in np.unique(residues):
                        select = residues == residue
                        ritem = item["per_residue"][int(residue)]
                        ritem["component_sse"] += np.square(component_error[select]).sum(axis=0)
                        ritem["samples"] += int(select.sum())
                        ritem["angle_abs"] += float(angle_error[select].sum())
                        ritem["angle_sq"] += float(np.square(angle_error[select]).sum())
    summaries = {}
    per_residue = []
    for readout in READOUTS:
        for split in ("R2", "R3"):
            summaries[f"{readout}_{split}"] = {"readout": readout, "split": split, "torsions": {}}
            for torsion in torsions:
                subset_names = ("all", "sidechain_present") if torsion in ("phi", "psi") else ("valid_chi1",)
                summaries[f"{readout}_{split}"]["torsions"][torsion] = {}
                for subset in subset_names:
                    item = accum[(readout, split, torsion, subset)]
                    component_mse = item["component_sse"] / item["samples"]
                    residue_component_mse = []
                    residue_angle_mae = []
                    residue_angle_mse = []
                    for residue, ritem in sorted(item["per_residue"].items()):
                        comp = ritem["component_sse"] / ritem["samples"]
                        angle_mae = ritem["angle_abs"] / ritem["samples"]
                        angle_rmse = math.sqrt(ritem["angle_sq"] / ritem["samples"])
                        residue_component_mse.append(comp)
                        residue_angle_mae.append(angle_mae)
                        residue_angle_mse.append(angle_rmse**2)
                        per_residue.append(
                            {
                                "readout": readout, "split": split, "torsion": torsion, "subset": subset,
                                "residue_index": residue, "samples": ritem["samples"],
                                "sin_rmse": math.sqrt(comp[0]), "cos_rmse": math.sqrt(comp[1]),
                                "angular_mae_deg": angle_mae, "angular_rmse_deg": angle_rmse,
                            }
                        )
                    residue_component_mse = np.stack(residue_component_mse)
                    summaries[f"{readout}_{split}"]["torsions"][torsion][subset] = {
                        "samples": item["samples"],
                        "valid_residue_positions": len(item["per_residue"]),
                        "sin_rmse": math.sqrt(component_mse[0]),
                        "cos_rmse": math.sqrt(component_mse[1]),
                        "angular_mae_deg": item["angle_abs"] / item["samples"],
                        "angular_rmse_deg": math.sqrt(item["angle_sq"] / item["samples"]),
                        "residue_balanced_sin_rmse": math.sqrt(float(residue_component_mse[:, 0].mean())),
                        "residue_balanced_cos_rmse": math.sqrt(float(residue_component_mse[:, 1].mean())),
                        "residue_balanced_angular_mae_deg": float(np.mean(residue_angle_mae)),
                        "residue_balanced_angular_rmse_deg": math.sqrt(float(np.mean(residue_angle_mse))),
                    }
            phi = summaries[f"{readout}_{split}"]["torsions"]["phi"]["all"]
            psi = summaries[f"{readout}_{split}"]["torsions"]["psi"]["all"]
            summaries[f"{readout}_{split}"]["backbone_all_combined_angular_mae_deg"] = (phi["residue_balanced_angular_mae_deg"] + psi["residue_balanced_angular_mae_deg"]) / 2
            summaries[f"{readout}_{split}"]["backbone_all_mean_component_rmse"] = np.mean([phi["residue_balanced_sin_rmse"], phi["residue_balanced_cos_rmse"], psi["residue_balanced_sin_rmse"], psi["residue_balanced_cos_rmse"]])
            phi_side = summaries[f"{readout}_{split}"]["torsions"]["phi"]["sidechain_present"]
            psi_side = summaries[f"{readout}_{split}"]["torsions"]["psi"]["sidechain_present"]
            summaries[f"{readout}_{split}"]["backbone_sidechain_present_combined_angular_mae_deg"] = (phi_side["residue_balanced_angular_mae_deg"] + psi_side["residue_balanced_angular_mae_deg"]) / 2

    scaler_json = {
        readout: {
            "count": scaler["count"], "mean": scaler["mean"].tolist(), "scale": scaler["scale"].tolist(),
            "replaced_channels": scaler["replaced_channels"], "raw_scale_distribution": scaler["raw_scale_distribution"],
        }
        for readout, scaler in scalers.items()
    }
    backbone_models = {
        readout: {
            torsion: {"samples": models[readout][torsion]["samples"], "coefficient": models[readout][torsion]["coefficient"].tolist()}
            for torsion in ("phi", "psi")
        }
        for readout in READOUTS
    }
    chi_models = {
        readout: {"samples": models[readout]["chi1"]["samples"], "coefficient": models[readout]["chi1"]["coefficient"].tolist()}
        for readout in READOUTS
    }
    backbone_metrics = {
        key: {
            "readout": value["readout"], "split": value["split"],
            "phi": value["torsions"]["phi"], "psi": value["torsions"]["psi"],
            "all_combined_angular_mae_deg": value["backbone_all_combined_angular_mae_deg"],
            "all_mean_component_rmse": value["backbone_all_mean_component_rmse"],
            "sidechain_present_combined_angular_mae_deg": value["backbone_sidechain_present_combined_angular_mae_deg"],
        }
        for key, value in summaries.items()
    }
    chi_metrics = {
        key: {
            "readout": value["readout"], "split": value["split"],
            **value["torsions"]["chi1"]["valid_chi1"],
        }
        for key, value in summaries.items()
    }
    write_json(output_root / "backbone_torsion_probe.json", {"config": FROZEN_CONFIG["probe"], "scalers": scaler_json, "models": backbone_models, "metrics": backbone_metrics})
    write_json(output_root / "chi1_probe.json", {"mapping": FROZEN_CONFIG["chi1_quadruplets"], "models": chi_models, "metrics": chi_metrics})
    write_csv(output_root / "backbone_torsion_probe_per_residue.csv", [row for row in per_residue if row["torsion"] in ("phi", "psi")])
    write_csv(output_root / "chi1_probe_per_residue.csv", [row for row in per_residue if row["torsion"] == "chi1"])
    return backbone_metrics, chi_metrics, per_residue, scaler_json


def replay_checks(model, rows: list[dict], numerical: dict, output_root: Path, device: str) -> dict:
    lookup = {row["relative_path"]: row for row in rows}
    results = []
    model = model.to(device).eval()
    for relative, frame_tuple in REPLAY_SPECS:
        row = lookup[relative]
        frames = np.asarray(frame_tuple, dtype=np.int64)
        coords = np.load(row["source"], mmap_mode="r")
        sequence = read_sequence(row["sequence_csv"])
        xyz, atom_z, old_ridx = flatten_atom14(np.asarray(coords[frames]), sequence)
        ridx, bmask, smask, _ = atom_metadata(sequence)
        if not np.array_equal(old_ridx, ridx):
            raise RuntimeError("replay atom order mismatch")
        replay = extract_readouts(model, xyz, atom_z, ridx, bmask, smask, len(sequence), device, len(frames))
        with np.load(row["cache"]) as archive:
            cached = {key: archive[key][frames] for key in STORED_READOUTS}
        with np.load(row["c_cache"]) as archive:
            c1_reference = archive["C1"][frames]
        state_errors = {}
        for key in STORED_READOUTS:
            metrics = error_metrics(replay[key], cached[key])
            metrics["relative_to_split_feature_rms"] = metrics["rmse"] / numerical[f"C1-{key}_{row['split']}"]["feature_rms"]
            state_errors[key] = metrics
        results.append(
            {
                "relative_path": relative,
                "split": row["split"],
                "frame_ids": frames.tolist(),
                "source_hash_unchanged": sha256(row["source"]) == row["source_sha256"],
                "cache_replay": state_errors,
                "M_vs_experiment_C_C1": error_metrics(replay["M"], c1_reference),
            }
        )
    report = {"subsets": results}
    write_json(output_root / "replay_audit.json", report)
    return report


def residue_improvement_fraction(per_residue: list[dict], candidate: str, control: str, split: str, family: str) -> dict:
    if family == "backbone":
        def aggregate(readout):
            values = defaultdict(list)
            for row in per_residue:
                if row["readout"] == readout and row["split"] == split and row["torsion"] in ("phi", "psi") and row["subset"] == "all":
                    values[row["residue_index"]].append(row["angular_mae_deg"])
            return {residue: float(np.mean(errors)) for residue, errors in values.items()}
    else:
        def aggregate(readout):
            return {
                row["residue_index"]: row["angular_mae_deg"]
                for row in per_residue
                if row["readout"] == readout and row["split"] == split and row["torsion"] == "chi1"
            }
    candidate_values = aggregate(candidate)
    control_values = aggregate(control)
    common = sorted(set(candidate_values) & set(control_values))
    deltas = np.asarray([candidate_values[index] - control_values[index] for index in common])
    return {
        "common_residue_positions": len(common),
        "fraction_candidate_lower": float(np.mean(deltas < 0)),
        "median_candidate_minus_control_deg": float(np.median(deltas)),
    }


def representation_comparison(
    numerical: dict,
    temporal: dict,
    backbone: dict,
    chi: dict,
    per_residue: list[dict],
    scalers: dict,
    replay: dict,
    continuity: dict,
    output_root: Path,
) -> dict:
    table = []
    for readout in READOUTS:
        for split in ("R2", "R3"):
            n = numerical[f"{readout}_{split}"]
            b = backbone[f"{readout}_{split}"]
            c = chi[f"{readout}_{split}"]
            table.append(
                {
                    "readout": readout, "split": split, "dimension": n["dimension"],
                    "feature_rms": n["feature_rms"], "near_constant_channels": n["near_constant_channels"],
                    "participation_rank": n["covariance_participation_rank"], "numerical_rank": n["numerical_rank"],
                    "positive_condition": n["positive_spectrum_condition"],
                    "backbone_component_rmse": b["all_mean_component_rmse"],
                    "phi_mae_deg": b["phi"]["all"]["residue_balanced_angular_mae_deg"],
                    "psi_mae_deg": b["psi"]["all"]["residue_balanced_angular_mae_deg"],
                    "backbone_combined_mae_deg": b["all_combined_angular_mae_deg"],
                    "chi1_sin_rmse": c["residue_balanced_sin_rmse"], "chi1_cos_rmse": c["residue_balanced_cos_rmse"],
                    "chi1_mae_deg": c["residue_balanced_angular_mae_deg"],
                }
            )
    comparisons = {}
    for split in ("R2", "R3"):
        for control in ("C1-M", "C1-MC"):
            key = f"BS_vs_{control}_{split}"
            comparisons[key] = {
                "backbone_improvement_fraction": 1 - backbone[f"C1-BS_{split}"]["all_combined_angular_mae_deg"] / backbone[f"{control}_{split}"]["all_combined_angular_mae_deg"],
                "chi1_improvement_fraction": 1 - chi[f"C1-BS_{split}"]["residue_balanced_angular_mae_deg"] / chi[f"{control}_{split}"]["residue_balanced_angular_mae_deg"],
                "backbone_per_residue": residue_improvement_fraction(per_residue, "C1-BS", control, split, "backbone"),
                "chi1_per_residue": residue_improvement_fraction(per_residue, "C1-BS", control, split, "chi1"),
            }
    accessibility = all(
        item[metric] > 0
        for item in comparisons.values()
        for metric in ("backbone_improvement_fraction", "chi1_improvement_fraction")
    )
    residue_breadth = all(
        item[family]["fraction_candidate_lower"] > 0.5 and item[family]["median_candidate_minus_control_deg"] < 0
        for item in comparisons.values()
        for family in ("backbone_per_residue", "chi1_per_residue")
    )
    rank_ok = all(
        numerical[f"C1-BS_{split}"]["covariance_participation_rank"] >= 0.5 * numerical[f"C1-M_{split}"]["covariance_participation_rank"]
        for split in ("R2", "R3")
    )
    scaler_ok = len(scalers["C1-BS"]["replaced_channels"]) == 0
    replay_ok = all(
        item["cache_replay"][key]["relative_to_split_feature_rms"] < 1e-4
        for item in replay["subsets"] for key in STORED_READOUTS
    )
    mc_exact = all(numerical[f"C1-MC_{split}"]["zero_block_exact"] for split in ("R2", "R3"))
    criteria = {
        "BS_improves_M_and_MC_both_families_both_splits": accessibility,
        "improvement_broad_across_residues": residue_breadth,
        "rank_acceptable": rank_ok,
        "scaler_acceptable": scaler_ok,
        "replay_acceptable": replay_ok,
        "MC_zero_block_exact": mc_exact,
    }
    result = {
        "continuity": continuity,
        "critical_comparisons": comparisons,
        "success_criteria": criteria,
        "BS_representation_level_promising": all(criteria.values()),
        "no_scalar_winner_selected": True,
        "biological_outcomes_used": False,
        "table": table,
    }
    write_json(output_root / "representation_comparison.json", result)
    write_csv(output_root / "representation_comparison.csv", table)
    return result


def continuity_metrics(backbone: dict, chi: dict, rows: list[dict]) -> dict:
    c_probe = json.loads((EXPERIMENT_C_ROOT / "torsion_probe.json").read_text(encoding="utf-8"))["metrics"]
    result = {"probe_metric_differences": {}, "cache_M_vs_C1": {}}
    for split in ("R2", "R3"):
        result["probe_metric_differences"][split] = {
            "C0_M_combined_mae_minus_experiment_C_C0": backbone[f"C0-M_{split}"]["all_combined_angular_mae_deg"] - c_probe[f"C0_{split}"]["residue_balanced_combined_angular_mae_deg"],
            "C1_M_combined_mae_minus_experiment_C_C1": backbone[f"C1-M_{split}"]["all_combined_angular_mae_deg"] - c_probe[f"C1_{split}"]["residue_balanced_combined_angular_mae_deg"],
        }
        maximum = 0.0
        sumsq = 0.0
        count = 0
        for row in (row for row in rows if row["split"] == split):
            with np.load(row["cache"]) as archive:
                observed = archive["M"].astype(np.float64)
            with np.load(row["c_cache"]) as archive:
                reference = archive["C1"].astype(np.float64)
            delta = observed - reference
            maximum = max(maximum, float(np.abs(delta).max()))
            sumsq += float(np.square(delta).sum())
            count += delta.size
        rmse = math.sqrt(sumsq / count)
        reference_rms = numerical_reference_rms(rows, split, "C1")
        result["cache_M_vs_C1"][split] = {"max_abs": maximum, "rmse": rmse, "relative_rmse": rmse / reference_rms}
    return result


def numerical_reference_rms(rows: list[dict], split: str, state: str) -> float:
    sumsq = 0.0
    count = 0
    for row in (row for row in rows if row["split"] == split):
        with np.load(row["c_cache"]) as archive:
            values = archive[state].astype(np.float64)
        sumsq += float(np.square(values).sum())
        count += values.size
    return math.sqrt(sumsq / count)


def markdown_table(rows: list[dict], columns: list[tuple[str, str]]) -> list[str]:
    lines = ["| " + " | ".join(label for _key, label in columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        values = []
        for key, _label in columns:
            value = row[key]
            values.append(f"{value:.5g}" if isinstance(value, float) else str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def render_report(provenance, disk, smoke, extraction, numerical, temporal, backbone, chi, comparison, replay, integrity) -> str:
    health_rows = []
    probe_rows = []
    temporal_rows = []
    for readout in READOUTS:
        for split in ("R2", "R3"):
            n = numerical[f"{readout}_{split}"]
            b = backbone[f"{readout}_{split}"]
            c = chi[f"{readout}_{split}"]
            health_rows.append({"readout": readout, "split": split, "rms": n["feature_rms"], "near": n["near_constant_channels"], "rank": n["covariance_participation_rank"], "nrank": n["numerical_rank"], "condition": n["positive_spectrum_condition"]})
            probe_rows.append({"readout": readout, "split": split, "component": b["all_mean_component_rmse"], "phi": b["phi"]["all"]["residue_balanced_angular_mae_deg"], "psi": b["psi"]["all"]["residue_balanced_angular_mae_deg"], "combined": b["all_combined_angular_mae_deg"], "chi_sin": c["residue_balanced_sin_rmse"], "chi_cos": c["residue_balanced_cos_rmse"], "chi": c["residue_balanced_angular_mae_deg"]})
    for readout in ("C1-M", "C1-B", "C1-S", "C1-BS"):
        for split in ("R2", "R3"):
            item = temporal[f"{readout}_{split}"]
            temporal_rows.append({"readout": readout, "split": split, "frame": item["scale_relative_frame_displacement"], "within": item["scale_relative_within_window"], "drift": item["scale_relative_chronological_drift"]})
    bs = comparison["BS_representation_level_promising"]
    lines = [
        "# PACER-DC encoder optimization — Experiment D v01",
        "",
        "Experiment D completed the fixed C1 atom-to-residue information-retention audit once. No layer search, learned pooling, biological contrast, long MD, or follow-up attention experiment was performed.",
        "",
        "## Provenance and extraction",
        "",
        f"All provenance checks passed before inference: checkpoint/source/runtime/C1 implementation, 40 atom14 inputs, and all 40 Experiment C reference caches. Disk preflight estimated `{disk['expected_raw_cache_gib']:.3f}` GiB with `{disk['free_gib']:.1f}` GiB free.",
        f"Extraction processed `{extraction['frames']}` frames in `{extraction['elapsed_seconds']:.1f}` seconds and wrote `{extraction['cache_gib']:.3f}` GiB across `{extraction['files_ok']}` caches. Only M/B/S were stored; BS and MC were derived.",
        "",
        "## Smoke test",
        "",
        f"Smoke passed: `{smoke['passed']}`. New M versus Experiment C C1: max absolute `{smoke['M_vs_experiment_C_C1']['max_abs']:.6g}`, relative RMSE `{smoke['M_vs_experiment_C_C1']['relative_rmse']:.6g}`. All exact BS/MC block checks passed.",
        "",
        "## Numerical health",
        "",
        *markdown_table(health_rows, [("readout", "Readout"), ("split", "Split"), ("rms", "Feature RMS"), ("near", "Near-constant"), ("rank", "Participation rank"), ("nrank", "Numerical rank"), ("condition", "Condition")]),
        "",
        "C1-S sidechain-present-only diagnostics and full covariance spectra are stored in the numerical-health artifacts. MC's second 128 channels remain exact zero; its informative rank is reported separately.",
        "",
        "## Structural probes",
        "",
        *markdown_table(probe_rows, [("readout", "Readout"), ("split", "Split"), ("component", "Backbone component RMSE"), ("phi", "Phi MAE"), ("psi", "Psi MAE"), ("combined", "Combined MAE"), ("chi_sin", "Chi1 sin RMSE"), ("chi_cos", "Chi1 cos RMSE"), ("chi", "Chi1 MAE")]),
        "",
        "R2 entries are fits; R3 entries are descriptive evaluations. Metrics are residue-balanced. Sidechain-present backbone subsets and per-residue diagnostics are included in their artifacts.",
        "",
        "## Scale-relative temporal diagnostics",
        "",
        *markdown_table(temporal_rows, [("readout", "Readout"), ("split", "Split"), ("frame", "Frame/RMS"), ("within", "Within/RMS"), ("drift", "Drift/RMS")]),
        "",
        "## Dimension-control interpretation",
        "",
        f"BS representation-level promising under the predeclared rule: `{bs}`. The detailed BS-versus-M and BS-versus-MC improvements and per-residue breadth checks are stored in `representation_comparison.json`.",
        "",
        "## Limitations",
        "",
        "R2 probe errors are in-sample fits. R3 was historically inspected and is descriptive. Frames and windows are correlated. Linear probes measure accessible information rather than all encoded information. Covariance conditioning is descriptive and scale-dependent. Historical runtime bitwise equivalence remains unavailable, and no biological efficacy conclusion follows.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--batch-audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--geom2vec-source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    frozen = dict(FROZEN_CONFIG)
    frozen["checkpoint"] = str(args.checkpoint.resolve())
    frozen["checkpoint_sha256"] = sha256(args.checkpoint)
    frozen["input_root"] = str(args.input_root.resolve())
    frozen["experiment_c_root"] = str(EXPERIMENT_C_ROOT.resolve())
    write_json(args.output_root / "frozen_config.json", frozen)
    write_json(
        args.output_root / "atom_group_definition.json",
        {
            "backbone_atoms": list(BACKBONE_ATOMS),
            "sidechain_definition": "all expected ATOM14 heavy atoms other than N/CA/C/O",
            "sidechain_absent_behavior": "exact zero vector; mask stored but never fed to probes",
            "chi1_quadruplets": FROZEN_CONFIG["chi1_quadruplets"],
            "excluded_chi1_residues": FROZEN_CONFIG["chi1_excluded_residues"],
            "atom14_mapping": ATOM14,
        },
    )

    model, rows, provenance = verify_provenance(args.checkpoint, args.batch_audit, args.input_root, args.geom2vec_source)
    write_json(args.output_root / "provenance.json", provenance)
    disk = disk_preflight(rows, args.output_root)
    write_json(args.output_root / "disk_preflight.json", disk)

    smoke = smoke_test(model, next(row for row in rows if row["relative_path"] == SMOKE_RELATIVE), args.device)
    write_json(args.output_root / "smoke_test.json", smoke)
    if not smoke["passed"]:
        raise SystemExit("Experiment D smoke test failed; stopped before full extraction")

    extraction = extract_all(model, rows, args.output_root, args.device)
    write_json(args.output_root / "extraction_audit.json", extraction)
    numerical = numerical_health(rows, args.output_root)
    temporal = temporal_diagnostics(rows, numerical, args.output_root)
    scalers = r2_scalers(rows)
    models = fit_models(rows, scalers)
    backbone, chi, per_residue, scaler_json = evaluate_models(rows, scalers, models, args.output_root)
    replay = replay_checks(model, rows, numerical, args.output_root, args.device)
    continuity = continuity_metrics(backbone, chi, rows)
    comparison = representation_comparison(numerical, temporal, backbone, chi, per_residue, scaler_json, replay, continuity, args.output_root)

    cache_hashes_ok = all(sha256(Path(item["cache_path"])) == item["cache_sha256"] for item in extraction["rows"])
    input_hashes_ok = all(sha256(row["source"]) == row["source_sha256"] for row in rows)
    c_cache_hashes_ok = all(sha256(row["c_cache"]) == row["c_cache_sha256"] for row in rows)
    core_hash_ok = sha256(PHASE02_IMPLEMENTATION) == provenance["c1_implementation"]["sha256"]
    checkpoint_ok = sha256(args.checkpoint) == EXPECTED_CHECKPOINT_SHA256
    replay_ok = all(item["M_vs_experiment_C_C1"]["max_abs"] <= 5e-5 for item in replay["subsets"])
    integrity = {
        "status": "EXPERIMENT_D_COMPLETE",
        "checkpoint_unchanged": checkpoint_ok,
        "source_inputs_unchanged": input_hashes_ok,
        "c1_hook_implementation_unchanged": core_hash_ok,
        "experiment_c_reference_caches_unchanged": c_cache_hashes_ok,
        "experiment_d_cache_hashes_verified": cache_hashes_ok,
        "M_reproduces_experiment_C_C1": replay_ok,
        "smoke_passed": smoke["passed"],
        "original_archives_overwritten": False,
        "prohibited_outputs_computed": [],
        "biological_labels_used": False,
        "followup_atom_attention_started": False,
        "experiment_d_implementation": file_record(Path(__file__)),
        "extraction": {key: value for key, value in extraction.items() if key != "rows"},
    }
    if not all((checkpoint_ok, input_hashes_ok, core_hash_ok, c_cache_hashes_ok, cache_hashes_ok, replay_ok, smoke["passed"])):
        integrity["status"] = "EXPERIMENT_D_INTEGRITY_FAILURE"
    write_json(args.output_root / "integrity_completion_audit.json", integrity)
    report_path = args.output_root / "REPORT.md"
    report_path.write_text(render_report(provenance, disk, smoke, extraction, numerical, temporal, backbone, chi, comparison, replay, integrity), encoding="utf-8")
    write_json(
        args.output_root / "EXPERIMENT_D_COMPLETE.json",
        {
            "status": integrity["status"],
            "report": file_record(report_path),
            "integrity_audit_sha256": sha256(args.output_root / "integrity_completion_audit.json"),
            "followup_started": False,
        },
    )
    print(json.dumps({"status": integrity["status"], "output": str(args.output_root)}, indent=2))
    if integrity["status"] != "EXPERIMENT_D_COMPLETE":
        raise SystemExit("Experiment D completed with an integrity failure")


if __name__ == "__main__":
    main()
