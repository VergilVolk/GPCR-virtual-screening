#!/usr/bin/env python
"""Experiment C v01 Phase 3: fixed C0/C1/C2 representation audit.

This module intentionally imports the already validated Phase 0--2 core and
does not redefine the representations, hooks, checkpoint loader, or residue
readout.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
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
    pool_batch,
    sha256,
    verify_and_strict_load,
)
from project.encoder_intermediate_v01.run_phase_0_2 import (  # noqa: E402
    DEFAULT_AUDIT,
    DEFAULT_CHECKPOINT,
    DEFAULT_INPUT,
    flatten_atom14,
    package_version,
    production_embed,
    read_sequence,
    runtime_report,
)

VERSION = "encoder_intermediate_C_phase3_v01"
STATES = ("C0", "C1", "C2")
REPLICAS = {2: "R2", 3: "R3"}
DEFAULT_OUTPUT = WORKSPACE / "project/results/encoder_intermediate_C_v01/phase3_run_001"
PHASE02_DIR = WORKSPACE / "project/results/encoder_intermediate_C_v01/run_001"
IMPLEMENTATION_FILES = (
    WORKSPACE / "project/encoder_intermediate_v01/core.py",
    WORKSPACE / "project/encoder_intermediate_v01/run_phase_0_2.py",
)

PHASE3_CONFIG = {
    "version": VERSION,
    "scope": "authenticated historical short R2/R3 atom14 inputs only",
    "states": {
        "C0": "production final normalized x_rep/v_rep",
        "C1": "accumulated x/v entering vis_mp_layers[3]",
        "C2": "final accumulated x/v at out_norm/vec_out_norm pre-hooks",
    },
    "readout": "concat(x, vector_norm(v, dim=xyz)); equal atom mean per residue",
    "output_shape_per_file": [100, 270, 128],
    "extraction_batch_size": 4,
    "cache_dtype": "float32",
    "cache_format": "one uncompressed NPZ containing C0/C1/C2 per authenticated input",
    "probe": {
        "target": ["sin(phi)", "cos(phi)", "sin(psi)", "cos(psi)"],
        "backbone_slots": {"N": 0, "CA": 1, "C": 2},
        "ridge_alpha": 1.0,
        "regularize_intercept": False,
        "feature_scaling": "per-state channel mean/std fit on R2 valid torsion union only; scale<=1e-12 replaced by 1",
        "fit_split": "R2",
        "evaluation_split": "R3 descriptive only",
        "angular_metric": "absolute wrapped atan2 angle error in degrees",
        "aggregate": "residue-balanced: average per-residue metrics, then equal phi/psi average",
    },
    "near_constant_rule": "variance <= max(1e-12, median(channel variance)*1e-8)",
    "effective_rank": {
        "participation_ratio": "sum(lambda)^2/sum(lambda^2)",
        "entropy_rank": "exp(-sum(p*log(p)))",
        "positive_eigen_threshold": "lambda_max*1e-12",
        "regularized_condition_epsilon": "lambda_max*1e-6",
    },
    "temporal": {
        "frame_displacement": "mean RMS feature displacement per residue for adjacent frames; windows equally weighted",
        "within_window": "RMS deviation from per-window/per-residue temporal centroid; windows equally weighted",
        "between_window": "RMS displacement of consecutive window centroids; trajectories equally weighted",
        "within_trajectory": "RMS variation of window centroids around trajectory centroid",
        "between_trajectory": "RMS variation of trajectory centroids around split centroid",
        "chronological_drift": "window-last versus window-first RMS displacement and linear window-index slope RMS",
    },
    "replay_subsets": [
        {"relative_path": "replica_02/window_000/atom14/apo_w000.npy", "frame_ids": [0, 1, 2, 3]},
        {"relative_path": "replica_03/window_000/atom14/apo_w000.npy", "frame_ids": [0, 1, 2, 3]},
    ],
    "promising_candidate_rule": {
        "accessibility": "at least 1% lower residue-balanced combined angular MAE and no higher mean component RMSE than C0 on both R2 and R3",
        "rank": "participation effective rank on both splits >= 50% of C0",
        "numerics": "finite, no severe scaler collapse, and replay/rigid relative RMSE < 1e-4",
        "purpose": "descriptive multi-criterion flag, not model selection",
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


def parse_replica(relative: str) -> tuple[int, str]:
    replica = int(next(part for part in Path(relative).parts if part.startswith("replica_")).split("_")[1])
    return replica, REPLICAS[replica]


def parse_window_condition(relative: str) -> tuple[int, str]:
    path = Path(relative)
    window = int(next(part for part in path.parts if part.startswith("window_")).split("_")[1])
    stem = path.stem
    suffix = f"_w{window:03d}"
    condition = stem[: -len(suffix)] if stem.endswith(suffix) else stem
    return window, condition


def load_authenticated_rows(audit_path: Path, input_root: Path) -> list[dict]:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    rows = []
    for archived in audit["rows"]:
        source = input_root / archived["input"]
        actual = sha256(source)
        if actual != archived["input_sha256"]:
            raise RuntimeError(f"authenticated input changed: {source}")
        replica, split = parse_replica(archived["input"])
        window, condition = parse_window_condition(archived["input"])
        rows.append(
            {
                "relative_path": archived["input"],
                "source": source,
                "source_sha256": actual,
                "sequence_csv": source.with_suffix(".csv"),
                "sequence_sha256": sha256(source.with_suffix(".csv")),
                "frames": int(archived["frames"]),
                "replica": replica,
                "split": split,
                "window": window,
                "condition": condition,
            }
        )
    if len(rows) != 40 or {row["replica"] for row in rows} != {2, 3}:
        raise RuntimeError("expected exactly the authenticated 40-file R2/R3 inventory")
    return rows


def disk_preflight(rows: list[dict], output_root: Path) -> dict:
    frames = sum(row["frames"] for row in rows)
    expected = frames * 270 * 128 * 4 * len(STATES)
    free = shutil.disk_usage(output_root.parent).free
    required = math.ceil(expected * 1.10)
    report = {
        "files": len(rows),
        "frames": frames,
        "expected_raw_cache_bytes": expected,
        "expected_raw_cache_gib": expected / 2**30,
        "required_with_10_percent_overhead_bytes": required,
        "free_bytes_before_extraction": free,
        "free_gib_before_extraction": free / 2**30,
        "sufficient": free >= required,
    }
    if not report["sufficient"]:
        raise RuntimeError("insufficient disk space for predeclared Phase 3 cache")
    return report


def prepare_batch(xyz: np.ndarray, atom_z: np.ndarray, device: str):
    pos = torch.as_tensor(xyz, dtype=torch.float32, device=device)
    batch_size, n_atoms, _ = pos.shape
    z0 = torch.as_tensor(atom_z, dtype=torch.long, device=device)
    z = z0.repeat(batch_size)
    batch = torch.arange(batch_size, device=device).repeat_interleave(n_atoms)
    return z, pos.reshape(-1, 3), batch, batch_size, n_atoms


def extract_array(model, xyz, atom_z, residue_index, n_residues, device, batch_size):
    pieces = {state: [] for state in STATES}
    ridx = torch.as_tensor(residue_index, dtype=torch.long, device=device)
    with torch.inference_mode():
        for start in range(0, len(xyz), batch_size):
            chunk = xyz[start : start + batch_size]
            z, pos, batch, actual_batch, n_atoms = prepare_batch(chunk, atom_z, device)
            states, _ = forward_three_states(model, z, pos, batch)
            for state in STATES:
                pooled = pool_batch(states[state], ridx, n_residues, actual_batch, n_atoms)
                pieces[state].append(pooled.cpu().numpy().astype(np.float32, copy=False))
    return {state: np.concatenate(pieces[state], axis=0) for state in STATES}


def extract_all(model, rows, output_root: Path, device: str) -> dict:
    cache_root = output_root / "cache"
    provenance_rows = []
    started_all = time.time()
    model = model.to(device).eval()
    for number, row in enumerate(rows, 1):
        coords = np.load(row["source"], mmap_mode="r")
        sequence = read_sequence(row["sequence_csv"])
        xyz, atom_z, residue_index = flatten_atom14(np.asarray(coords), sequence)
        started = time.time()
        arrays = extract_array(
            model,
            xyz,
            atom_z,
            residue_index,
            len(sequence),
            device,
            PHASE3_CONFIG["extraction_batch_size"],
        )
        relative_cache = Path(row["relative_path"]).with_suffix(".C012.npz")
        cache_path = cache_root / relative_cache
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            cache_path,
            C0=arrays["C0"],
            C1=arrays["C1"],
            C2=arrays["C2"],
            frame_ids=np.arange(len(coords), dtype=np.int64),
            residue_index=np.arange(len(sequence), dtype=np.int64),
            sequence=np.asarray(sequence),
        )
        for state in STATES:
            if arrays[state].shape != (row["frames"], 270, 128) or not np.isfinite(arrays[state]).all():
                raise RuntimeError(f"invalid {state} extraction for {row['source']}")
        row["cache"] = cache_path
        record = {
            "number": number,
            "relative_path": row["relative_path"],
            "source_path": str(row["source"].resolve()),
            "source_sha256": row["source_sha256"],
            "sequence_path": str(row["sequence_csv"].resolve()),
            "sequence_sha256": row["sequence_sha256"],
            "replica": row["replica"],
            "split": row["split"],
            "condition": row["condition"],
            "window": row["window"],
            "frame_ids": [0, row["frames"] - 1],
            "residue_index": [0, 269],
            "shape": [row["frames"], 270, 128],
            "states": list(STATES),
            "finite": True,
            "cache_path": str(cache_path.resolve()),
            "cache_bytes": cache_path.stat().st_size,
            "cache_sha256": sha256(cache_path),
            "elapsed_seconds": round(time.time() - started, 3),
        }
        provenance_rows.append(record)
        print(f"[{number:02d}/{len(rows)}] {row['relative_path']} -> {record['cache_bytes']} bytes", flush=True)
    return {
        "method": "single-forward C0/C1/C2 extraction using validated fixed hooks",
        "files_expected": len(rows),
        "files_ok": len(provenance_rows),
        "frames": sum(row["frames"] for row in rows),
        "elapsed_seconds": round(time.time() - started_all, 3),
        "cache_bytes": sum(row["cache_bytes"] for row in provenance_rows),
        "cache_gib": sum(row["cache_bytes"] for row in provenance_rows) / 2**30,
        "rows": provenance_rows,
    }


def quantiles(values: np.ndarray) -> dict:
    q = np.quantile(values, [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
    return dict(zip(("min", "p05", "p25", "median", "p75", "p95", "max"), map(float, q)))


def numerical_health(rows: list[dict], output_root: Path) -> tuple[dict, list[dict], list[dict]]:
    result = {}
    channel_rows = []
    spectrum_rows = []
    for split in ("R2", "R3"):
        selected = [row for row in rows if row["split"] == split]
        for state in STATES:
            total_samples = 0
            finite_elements = 0
            total_elements = 0
            sums = np.zeros(128, dtype=np.float64)
            sums2 = np.zeros(128, dtype=np.float64)
            cross = np.zeros((128, 128), dtype=np.float64)
            norms = []
            absolute_chunks = []
            for row in selected:
                with np.load(row["cache"]) as archive:
                    flat32 = archive[state].reshape(-1, 128)
                finite_elements += int(np.isfinite(flat32).sum())
                total_elements += flat32.size
                flat = flat32.astype(np.float64)
                total_samples += len(flat)
                sums += flat.sum(axis=0)
                sums2 += np.square(flat).sum(axis=0)
                cross += flat.T @ flat
                norms.append(np.linalg.norm(flat32, axis=1))
                absolute_chunks.append(np.abs(flat32).reshape(-1))
            mean = sums / total_samples
            variance = np.maximum(sums2 / total_samples - mean**2, 0.0)
            std = np.sqrt(variance)
            covariance = cross / total_samples - np.outer(mean, mean)
            covariance = (covariance + covariance.T) * 0.5
            eigenvalues = np.linalg.eigvalsh(covariance)[::-1]
            eigenvalues = np.maximum(eigenvalues, 0.0)
            eigen_sum = eigenvalues.sum()
            weights = eigenvalues / eigen_sum if eigen_sum > 0 else np.zeros_like(eigenvalues)
            participation = eigen_sum**2 / np.square(eigenvalues).sum() if np.square(eigenvalues).sum() > 0 else 0.0
            positive_threshold = eigenvalues[0] * 1e-12 if eigenvalues[0] > 0 else 0.0
            positive = eigenvalues[eigenvalues > positive_threshold]
            raw_condition = eigenvalues[0] / positive[-1] if len(positive) else math.inf
            epsilon = eigenvalues[0] * 1e-6 if eigenvalues[0] > 0 else 1e-12
            regularized_condition = (eigenvalues[0] + epsilon) / (eigenvalues[-1] + epsilon)
            median_variance = float(np.median(variance))
            near_threshold = max(1e-12, median_variance * 1e-8)
            near_constant = variance <= near_threshold
            all_norms = np.concatenate(norms)
            all_absolute = np.concatenate(absolute_chunks)
            rms = math.sqrt(float(sums2.sum() / total_elements))
            key = f"{state}_{split}"
            result[key] = {
                "state": state,
                "split": split,
                "observations": total_samples,
                "elements": total_elements,
                "finite_elements": finite_elements,
                "nonfinite_elements": total_elements - finite_elements,
                "channel_mean_distribution": quantiles(mean),
                "channel_std_distribution": quantiles(std),
                "channel_variance_distribution": quantiles(variance),
                "near_constant_threshold": near_threshold,
                "near_constant_channels": int(near_constant.sum()),
                "covariance_trace": float(eigen_sum),
                "covariance_rank_numpy": int(np.linalg.matrix_rank(covariance)),
                "effective_rank_participation": float(participation),
                "effective_rank_entropy": float(np.exp(-np.sum(weights[weights > 0] * np.log(weights[weights > 0])))),
                "positive_covariance_eigenvalues": int(len(positive)),
                "condition_positive_spectrum": float(raw_condition),
                "condition_regularized_1e-6": float(regularized_condition),
                "feature_norm_distribution": quantiles(all_norms),
                "median_feature_magnitude": float(np.median(all_absolute)),
                "rms_feature_magnitude": rms,
            }
            cumulative = np.cumsum(weights)
            for index, value in enumerate(eigenvalues):
                spectrum_rows.append(
                    {
                        "state": state,
                        "split": split,
                        "eigen_index": index,
                        "eigenvalue": float(value),
                        "explained_fraction": float(weights[index]),
                        "cumulative_fraction": float(cumulative[index]),
                    }
                )
            for channel in range(128):
                channel_rows.append(
                    {
                        "state": state,
                        "split": split,
                        "channel": channel,
                        "mean": float(mean[channel]),
                        "std": float(std[channel]),
                        "variance": float(variance[channel]),
                        "near_constant": bool(near_constant[channel]),
                    }
                )
    write_json(output_root / "numerical_health.json", result)
    write_csv(output_root / "numerical_health_channels.csv", channel_rows)
    write_csv(output_root / "numerical_health_covariance_spectrum.csv", spectrum_rows)
    return result, channel_rows, spectrum_rows


def rms_delta(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(a.astype(np.float64) - b.astype(np.float64)))))


def temporal_diagnostics(rows: list[dict], output_root: Path) -> tuple[dict, list[dict]]:
    summaries = {}
    detail = []
    for split in ("R2", "R3"):
        selected = [row for row in rows if row["split"] == split]
        for state in STATES:
            window_centroids = {}
            frame_displacements = []
            within_windows = []
            for row in selected:
                with np.load(row["cache"]) as archive:
                    features = archive[state]
                centroid = features.mean(axis=0, dtype=np.float64)
                key = (row["replica"], row["condition"])
                window_centroids[(key, row["window"])] = centroid
                adjacent = np.sqrt(np.mean(np.square(np.diff(features.astype(np.float64), axis=0)), axis=2))
                frame_value = float(adjacent.mean())
                within_value = float(np.sqrt(np.mean(np.square(features.astype(np.float64) - centroid[None, :, :]))))
                frame_displacements.append(frame_value)
                within_windows.append(within_value)
                detail.append(
                    {
                        "state": state,
                        "split": split,
                        "replica": row["replica"],
                        "condition": row["condition"],
                        "window": row["window"],
                        "frame_to_frame_displacement": frame_value,
                        "within_window_rms": within_value,
                    }
                )
            trajectories = sorted({key for key, _window in window_centroids})
            between_consecutive = []
            within_trajectory = []
            chronological_drift = []
            chronological_slope = []
            trajectory_centroids = []
            for trajectory in trajectories:
                windows = sorted(window for key, window in window_centroids if key == trajectory)
                centroids = np.stack([window_centroids[(trajectory, window)] for window in windows])
                trajectory_centroid = centroids.mean(axis=0)
                trajectory_centroids.append(trajectory_centroid)
                between_consecutive.append(float(np.mean([rms_delta(centroids[i], centroids[i - 1]) for i in range(1, len(centroids))])))
                within_trajectory.append(float(np.sqrt(np.mean(np.square(centroids - trajectory_centroid[None])))))
                chronological_drift.append(rms_delta(centroids[-1], centroids[0]))
                index = np.asarray(windows, dtype=np.float64)
                centered_index = index - index.mean()
                slope = np.tensordot(centered_index, centroids, axes=(0, 0)) / np.square(centered_index).sum()
                chronological_slope.append(float(np.sqrt(np.mean(np.square(slope)))))
            trajectory_centroids = np.stack(trajectory_centroids)
            split_centroid = trajectory_centroids.mean(axis=0)
            between_trajectory = float(np.sqrt(np.mean(np.square(trajectory_centroids - split_centroid[None]))))
            within_trajectory_mean = float(np.mean(within_trajectory))
            summaries[f"{state}_{split}"] = {
                "state": state,
                "split": split,
                "windows": len(selected),
                "trajectories": len(trajectories),
                "frame_to_frame_displacement_window_balanced_mean": float(np.mean(frame_displacements)),
                "frame_to_frame_displacement_window_median": float(np.median(frame_displacements)),
                "within_window_rms_window_balanced_mean": float(np.mean(within_windows)),
                "between_consecutive_window_centroid_rms_trajectory_mean": float(np.mean(between_consecutive)),
                "within_trajectory_window_centroid_rms_mean": within_trajectory_mean,
                "between_trajectory_centroid_rms": between_trajectory,
                "between_to_within_trajectory_ratio": between_trajectory / within_trajectory_mean if within_trajectory_mean else math.inf,
                "chronological_first_to_last_rms_trajectory_mean": float(np.mean(chronological_drift)),
                "chronological_window_slope_rms_trajectory_mean": float(np.mean(chronological_slope)),
            }
    write_json(output_root / "temporal_diagnostics.json", summaries)
    write_csv(output_root / "temporal_diagnostics_windows.csv", detail)
    return summaries, detail


def dihedral(p0, p1, p2, p3) -> tuple[np.ndarray, np.ndarray]:
    p0 = p0.astype(np.float64)
    p1 = p1.astype(np.float64)
    p2 = p2.astype(np.float64)
    p3 = p3.astype(np.float64)
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1_norm = np.linalg.norm(b1, axis=-1)
    valid = np.isfinite(np.stack((p0, p1, p2, p3))).all(axis=(0, -1)) & (b1_norm > 1e-8)
    unit = b1 / np.maximum(b1_norm[..., None], 1e-300)
    v = b0 - np.sum(b0 * unit, axis=-1, keepdims=True) * unit
    w = b2 - np.sum(b2 * unit, axis=-1, keepdims=True) * unit
    vnorm = np.linalg.norm(v, axis=-1)
    wnorm = np.linalg.norm(w, axis=-1)
    valid &= (vnorm > 1e-8) & (wnorm > 1e-8)
    x = np.sum(v * w, axis=-1)
    y = np.sum(np.cross(unit, v) * w, axis=-1)
    angle = np.arctan2(y, x)
    valid &= np.isfinite(angle)
    return angle, valid


def torsion_targets(coords: np.ndarray) -> dict:
    n = coords[:, :, 0, :]
    ca = coords[:, :, 1, :]
    c = coords[:, :, 2, :]
    frames, residues = n.shape[:2]
    phi = np.full((frames, residues), np.nan, dtype=np.float64)
    psi = np.full((frames, residues), np.nan, dtype=np.float64)
    phi_values, phi_valid = dihedral(c[:, :-1], n[:, 1:], ca[:, 1:], c[:, 1:])
    psi_values, psi_valid = dihedral(n[:, :-1], ca[:, :-1], c[:, :-1], n[:, 1:])
    phi[:, 1:] = np.where(phi_valid, phi_values, np.nan)
    psi[:, :-1] = np.where(psi_valid, psi_values, np.nan)
    return {
        "phi": phi,
        "psi": psi,
        "phi_mask": np.isfinite(phi),
        "psi_mask": np.isfinite(psi),
    }


def standardizers(rows: list[dict]) -> dict:
    r2 = [row for row in rows if row["split"] == "R2"]
    accum = {state: {"count": 0, "sum": np.zeros(128), "sum2": np.zeros(128)} for state in STATES}
    for row in r2:
        coords = np.asarray(np.load(row["source"], mmap_mode="r"))
        targets = torsion_targets(coords)
        mask = targets["phi_mask"] | targets["psi_mask"]
        with np.load(row["cache"]) as archive:
            for state in STATES:
                x = archive[state][mask].astype(np.float64)
                accum[state]["count"] += len(x)
                accum[state]["sum"] += x.sum(axis=0)
                accum[state]["sum2"] += np.square(x).sum(axis=0)
    result = {}
    for state in STATES:
        count = accum[state]["count"]
        mean = accum[state]["sum"] / count
        variance = np.maximum(accum[state]["sum2"] / count - mean**2, 0)
        raw_scale = np.sqrt(variance)
        replaced = raw_scale <= 1e-12
        scale = np.where(replaced, 1.0, raw_scale)
        result[state] = {
            "count": count,
            "mean": mean,
            "scale": scale,
            "replaced_channels": np.flatnonzero(replaced).tolist(),
            "scale_distribution": quantiles(raw_scale),
        }
    return result


def fit_ridge(rows: list[dict], scalers: dict) -> dict:
    r2 = [row for row in rows if row["split"] == "R2"]
    accum = {
        state: {
            torsion: {"xtx": np.zeros((129, 129)), "xty": np.zeros((129, 2)), "n": 0}
            for torsion in ("phi", "psi")
        }
        for state in STATES
    }
    for row in r2:
        coords = np.asarray(np.load(row["source"], mmap_mode="r"))
        targets = torsion_targets(coords)
        with np.load(row["cache"]) as archive:
            for state in STATES:
                features = archive[state].astype(np.float64)
                standardized = (features - scalers[state]["mean"]) / scalers[state]["scale"]
                for torsion in ("phi", "psi"):
                    mask = targets[f"{torsion}_mask"]
                    x = standardized[mask]
                    angle = targets[torsion][mask]
                    y = np.column_stack((np.sin(angle), np.cos(angle)))
                    xa = np.column_stack((x, np.ones(len(x))))
                    accum[state][torsion]["xtx"] += xa.T @ xa
                    accum[state][torsion]["xty"] += xa.T @ y
                    accum[state][torsion]["n"] += len(x)
    models = {}
    regularizer = np.eye(129) * PHASE3_CONFIG["probe"]["ridge_alpha"]
    regularizer[-1, -1] = 0.0
    for state in STATES:
        models[state] = {}
        for torsion in ("phi", "psi"):
            item = accum[state][torsion]
            coefficient = np.linalg.solve(item["xtx"] + regularizer, item["xty"])
            models[state][torsion] = {"coefficient": coefficient, "training_samples": item["n"]}
    return models


def wrapped_angle_error(predicted: np.ndarray, observed: np.ndarray) -> np.ndarray:
    return np.abs(np.arctan2(np.sin(predicted - observed), np.cos(predicted - observed))) * 180.0 / np.pi


def evaluate_probes(rows: list[dict], scalers: dict, models: dict, output_root: Path) -> tuple[dict, list[dict]]:
    accum = {}
    for state in STATES:
        for split in ("R2", "R3"):
            accum[(state, split)] = {
                torsion: {
                    "component_sse": np.zeros(2),
                    "component_n": 0,
                    "angle_abs": 0.0,
                    "angle_sq": 0.0,
                    "angle_n": 0,
                    "residue": defaultdict(lambda: {"component_sse": np.zeros(2), "n": 0, "angle_abs": 0.0, "angle_sq": 0.0}),
                }
                for torsion in ("phi", "psi")
            }
    for row in rows:
        coords = np.asarray(np.load(row["source"], mmap_mode="r"))
        targets = torsion_targets(coords)
        residue_grid = np.broadcast_to(np.arange(270)[None, :], targets["phi"].shape)
        with np.load(row["cache"]) as archive:
            for state in STATES:
                standardized = (archive[state].astype(np.float64) - scalers[state]["mean"]) / scalers[state]["scale"]
                for torsion in ("phi", "psi"):
                    mask = targets[f"{torsion}_mask"]
                    x = standardized[mask]
                    xa = np.column_stack((x, np.ones(len(x))))
                    observed_angle = targets[torsion][mask]
                    observed = np.column_stack((np.sin(observed_angle), np.cos(observed_angle)))
                    predicted = xa @ models[state][torsion]["coefficient"]
                    component_error = predicted - observed
                    predicted_angle = np.arctan2(predicted[:, 0], predicted[:, 1])
                    angle_error = wrapped_angle_error(predicted_angle, observed_angle)
                    item = accum[(state, row["split"])][torsion]
                    item["component_sse"] += np.square(component_error).sum(axis=0)
                    item["component_n"] += len(predicted)
                    item["angle_abs"] += angle_error.sum()
                    item["angle_sq"] += np.square(angle_error).sum()
                    item["angle_n"] += len(angle_error)
                    residues = residue_grid[mask]
                    for residue in np.unique(residues):
                        select = residues == residue
                        ritem = item["residue"][int(residue)]
                        ritem["component_sse"] += np.square(component_error[select]).sum(axis=0)
                        ritem["n"] += int(select.sum())
                        ritem["angle_abs"] += float(angle_error[select].sum())
                        ritem["angle_sq"] += float(np.square(angle_error[select]).sum())
    summaries = {}
    per_residue = []
    for state in STATES:
        for split in ("R2", "R3"):
            state_summary = {"state": state, "split": split, "torsions": {}}
            component_rmses = []
            angle_maes = []
            for torsion in ("phi", "psi"):
                item = accum[(state, split)][torsion]
                component_mse = item["component_sse"] / item["component_n"]
                residue_component_mse = []
                residue_angle_mae = []
                residue_angle_mse = []
                for residue, ritem in sorted(item["residue"].items()):
                    comp_mse = ritem["component_sse"] / ritem["n"]
                    angle_mae = ritem["angle_abs"] / ritem["n"]
                    angle_rmse = math.sqrt(ritem["angle_sq"] / ritem["n"])
                    residue_component_mse.append(comp_mse)
                    residue_angle_mae.append(angle_mae)
                    residue_angle_mse.append(angle_rmse**2)
                    per_residue.append(
                        {
                            "state": state,
                            "split": split,
                            "torsion": torsion,
                            "residue_index": residue,
                            "samples": ritem["n"],
                            "sin_rmse": math.sqrt(comp_mse[0]),
                            "cos_rmse": math.sqrt(comp_mse[1]),
                            "angular_mae_deg": angle_mae,
                            "angular_rmse_deg": angle_rmse,
                        }
                    )
                residue_component_mse = np.stack(residue_component_mse)
                torsion_summary = {
                    "samples": item["angle_n"],
                    "residues": len(item["residue"]),
                    "sin_mse": float(component_mse[0]),
                    "sin_rmse": float(math.sqrt(component_mse[0])),
                    "cos_mse": float(component_mse[1]),
                    "cos_rmse": float(math.sqrt(component_mse[1])),
                    "angular_mae_deg": item["angle_abs"] / item["angle_n"],
                    "angular_rmse_deg": math.sqrt(item["angle_sq"] / item["angle_n"]),
                    "residue_balanced_sin_rmse": float(math.sqrt(residue_component_mse[:, 0].mean())),
                    "residue_balanced_cos_rmse": float(math.sqrt(residue_component_mse[:, 1].mean())),
                    "residue_balanced_angular_mae_deg": float(np.mean(residue_angle_mae)),
                    "residue_balanced_angular_rmse_deg": float(math.sqrt(np.mean(residue_angle_mse))),
                }
                state_summary["torsions"][torsion] = torsion_summary
                component_rmses.extend([torsion_summary["residue_balanced_sin_rmse"], torsion_summary["residue_balanced_cos_rmse"]])
                angle_maes.append(torsion_summary["residue_balanced_angular_mae_deg"])
            state_summary["residue_balanced_mean_component_rmse"] = float(np.mean(component_rmses))
            state_summary["residue_balanced_combined_angular_mae_deg"] = float(np.mean(angle_maes))
            state_summary["fit_or_evaluation"] = "R2 fit" if split == "R2" else "R3 descriptive evaluation"
            summaries[f"{state}_{split}"] = state_summary
    serializable_scalers = {
        state: {
            "count": scalers[state]["count"],
            "mean": scalers[state]["mean"].tolist(),
            "scale": scalers[state]["scale"].tolist(),
            "replaced_channels": scalers[state]["replaced_channels"],
            "scale_distribution": scalers[state]["scale_distribution"],
        }
        for state in STATES
    }
    serializable_models = {
        state: {
            torsion: {
                "training_samples": models[state][torsion]["training_samples"],
                "coefficient": models[state][torsion]["coefficient"].tolist(),
            }
            for torsion in ("phi", "psi")
        }
        for state in STATES
    }
    write_json(
        output_root / "torsion_probe.json",
        {
            "config": PHASE3_CONFIG["probe"],
            "scalers": serializable_scalers,
            "models": serializable_models,
            "metrics": summaries,
        },
    )
    write_csv(output_root / "torsion_probe_per_residue.csv", per_residue)
    return summaries, per_residue


def replay_checks(model, rows, numerical, output_root: Path, device: str) -> dict:
    results = []
    row_lookup = {row["relative_path"]: row for row in rows}
    model = model.to(device).eval()
    for specification in PHASE3_CONFIG["replay_subsets"]:
        row = row_lookup[specification["relative_path"]]
        frame_ids = np.asarray(specification["frame_ids"], dtype=np.int64)
        coords = np.load(row["source"], mmap_mode="r")
        sequence = read_sequence(row["sequence_csv"])
        xyz, atom_z, residue_index = flatten_atom14(np.asarray(coords[frame_ids]), sequence)
        replay = extract_array(model, xyz, atom_z, residue_index, len(sequence), device, len(frame_ids))
        with np.load(row["cache"]) as archive:
            cached = {state: archive[state][frame_ids] for state in STATES}
        production = production_embed(model, xyz, atom_z, residue_index, len(sequence), len(frame_ids), device)
        state_errors = {}
        for state in STATES:
            delta = replay[state].astype(np.float64) - cached[state].astype(np.float64)
            rms = numerical[f"{state}_{row['split']}"]["rms_feature_magnitude"]
            state_errors[state] = {
                "max_abs": float(np.abs(delta).max()),
                "rmse": float(np.sqrt(np.mean(np.square(delta)))),
                "relative_rmse_to_split_feature_rms": float(np.sqrt(np.mean(np.square(delta))) / rms),
            }
        c0_delta = replay["C0"].astype(np.float64) - production.astype(np.float64)
        results.append(
            {
                "relative_path": row["relative_path"],
                "split": row["split"],
                "frame_ids": frame_ids.tolist(),
                "source_sha256_unchanged": sha256(row["source"]) == row["source_sha256"],
                "state_replay": state_errors,
                "c0_production_consistency": {
                    "max_abs": float(np.abs(c0_delta).max()),
                    "rmse": float(np.sqrt(np.mean(np.square(c0_delta)))),
                    "within_smoke_tolerance_2e-6": float(np.abs(c0_delta).max()) <= 2e-6,
                },
            }
        )
    phase02_smoke = json.loads((PHASE02_DIR / "smoke_test_metrics.json").read_text(encoding="utf-8"))
    rigid_relative = {}
    for state in STATES:
        rms = numerical[f"{state}_R2"]["rms_feature_magnitude"]
        rigid_relative[state] = {
            "phase02_rotation_rmse": phase02_smoke["rigid_transform"]["rotation_invariance"][state]["rmse"],
            "phase02_rotation_relative_to_phase3_R2_feature_rms": phase02_smoke["rigid_transform"]["rotation_invariance"][state]["rmse"] / rms,
            "phase02_translation_rmse": phase02_smoke["rigid_transform"]["translation_invariance"][state]["rmse"],
            "phase02_translation_relative_to_phase3_R2_feature_rms": phase02_smoke["rigid_transform"]["translation_invariance"][state]["rmse"] / rms,
        }
    report = {"subsets": results, "phase02_rigid_errors_relative_to_phase3_feature_rms": rigid_relative}
    write_json(output_root / "replay_integrity.json", report)
    return report


def comparison(numerical, temporal, probe, replay, output_root: Path) -> tuple[dict, list[dict]]:
    table = []
    for state in STATES:
        for split in ("R2", "R3"):
            health = numerical[f"{state}_{split}"]
            dynamics = temporal[f"{state}_{split}"]
            accessibility = probe[f"{state}_{split}"]
            table.append(
                {
                    "state": state,
                    "split": split,
                    "rms_feature_magnitude": health["rms_feature_magnitude"],
                    "near_constant_channels": health["near_constant_channels"],
                    "effective_rank_participation": health["effective_rank_participation"],
                    "effective_rank_entropy": health["effective_rank_entropy"],
                    "regularized_condition": health["condition_regularized_1e-6"],
                    "frame_displacement": dynamics["frame_to_frame_displacement_window_balanced_mean"],
                    "within_window_rms": dynamics["within_window_rms_window_balanced_mean"],
                    "between_window_centroid_rms": dynamics["between_consecutive_window_centroid_rms_trajectory_mean"],
                    "chronological_drift_rms": dynamics["chronological_first_to_last_rms_trajectory_mean"],
                    "torsion_component_rmse": accessibility["residue_balanced_mean_component_rmse"],
                    "torsion_angular_mae_deg": accessibility["residue_balanced_combined_angular_mae_deg"],
                    "phi_angular_mae_deg": accessibility["torsions"]["phi"]["residue_balanced_angular_mae_deg"],
                    "psi_angular_mae_deg": accessibility["torsions"]["psi"]["residue_balanced_angular_mae_deg"],
                }
            )
    flags = {}
    for candidate in ("C1", "C2"):
        access = {}
        rank = {}
        for split in ("R2", "R3"):
            candidate_probe = probe[f"{candidate}_{split}"]
            baseline_probe = probe[f"C0_{split}"]
            access[split] = {
                "angular_improvement_fraction": 1 - candidate_probe["residue_balanced_combined_angular_mae_deg"] / baseline_probe["residue_balanced_combined_angular_mae_deg"],
                "component_rmse_no_higher": candidate_probe["residue_balanced_mean_component_rmse"] <= baseline_probe["residue_balanced_mean_component_rmse"],
            }
            rank[split] = numerical[f"{candidate}_{split}"]["effective_rank_participation"] / numerical[f"C0_{split}"]["effective_rank_participation"]
        replay_relative = max(
            subset["state_replay"][candidate]["relative_rmse_to_split_feature_rms"] for subset in replay["subsets"]
        )
        rigid_relative = max(
            replay["phase02_rigid_errors_relative_to_phase3_feature_rms"][candidate]["phase02_rotation_relative_to_phase3_R2_feature_rms"],
            replay["phase02_rigid_errors_relative_to_phase3_feature_rms"][candidate]["phase02_translation_relative_to_phase3_R2_feature_rms"],
        )
        scaler_collapse = len(json.loads((output_root / "torsion_probe.json").read_text())["scalers"][candidate]["replaced_channels"]) > 0
        criteria = {
            "R2_accessibility": access["R2"]["angular_improvement_fraction"] >= 0.01 and access["R2"]["component_rmse_no_higher"],
            "R3_accessibility": access["R3"]["angular_improvement_fraction"] >= 0.01 and access["R3"]["component_rmse_no_higher"],
            "rank_retention": rank["R2"] >= 0.5 and rank["R3"] >= 0.5,
            "no_scaler_collapse": not scaler_collapse,
            "numerical_invariance": replay_relative < 1e-4 and rigid_relative < 1e-4,
        }
        flags[candidate] = {
            "accessibility": access,
            "effective_rank_fraction_of_C0": rank,
            "max_replay_relative_rmse": replay_relative,
            "max_rigid_relative_rmse": rigid_relative,
            "criteria": criteria,
            "promising_structural_accessibility_advantage": all(criteria.values()),
        }
    result = {
        "interpretation_rule": PHASE3_CONFIG["promising_candidate_rule"],
        "candidate_assessments": flags,
        "no_single_scalar_winner_selected": True,
        "biological_outcomes_used": False,
    }
    write_json(output_root / "representation_comparison.json", {"summary": result, "table": table})
    write_csv(output_root / "representation_comparison.csv", table)
    return result, table


def markdown_table(rows: list[dict], columns: list[tuple[str, str]], digits: int = 5) -> list[str]:
    headers = [label for _key, label in columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        values = []
        for key, _label in columns:
            value = row[key]
            values.append(f"{value:.{digits}g}" if isinstance(value, float) else str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def render_report(integrity, numerical, temporal, probe, comparison_result, comparison_table, replay) -> str:
    health_rows = []
    temporal_rows = []
    probe_rows = []
    for state in STATES:
        for split in ("R2", "R3"):
            h = numerical[f"{state}_{split}"]
            t = temporal[f"{state}_{split}"]
            p = probe[f"{state}_{split}"]
            health_rows.append({"state": state, "split": split, "rms": h["rms_feature_magnitude"], "median_abs": h["median_feature_magnitude"], "near": h["near_constant_channels"], "rank": h["effective_rank_participation"], "cond": h["condition_regularized_1e-6"]})
            temporal_rows.append({"state": state, "split": split, "frame": t["frame_to_frame_displacement_window_balanced_mean"], "within": t["within_window_rms_window_balanced_mean"], "windows": t["between_consecutive_window_centroid_rms_trajectory_mean"], "between": t["between_trajectory_centroid_rms"], "drift": t["chronological_first_to_last_rms_trajectory_mean"]})
            probe_rows.append({"state": state, "split": split, "component": p["residue_balanced_mean_component_rmse"], "phi": p["torsions"]["phi"]["residue_balanced_angular_mae_deg"], "psi": p["torsions"]["psi"]["residue_balanced_angular_mae_deg"], "combined": p["residue_balanced_combined_angular_mae_deg"]})
    c1 = comparison_result["candidate_assessments"]["C1"]
    c2 = comparison_result["candidate_assessments"]["C2"]
    lines = [
        "# PACER-DC encoder optimization — Experiment C v01 Phase 3",
        "",
        "Phase 3 completed the fixed C0/C1/C2 representation audit once. No biological contrasts, PACER outcomes, layer sweep, or follow-up representation were evaluated.",
        "",
        "## Extraction and provenance",
        "",
        f"All `{integrity['extraction']['files_ok']}` authenticated short R2/R3 files (`{integrity['extraction']['frames']}` frames) were extracted successfully. C0/C1/C2 were captured in the same forward pass. The cache occupies `{integrity['extraction']['cache_gib']:.3f}` GiB. Checkpoint, input, source, cache, and validated implementation hashes are recorded in the integrity files.",
        "",
        "Historical runtime bitwise equivalence remains unestablished; the Phase 3 runtime itself is recorded exactly.",
        "",
        "## Numerical health",
        "",
        *markdown_table(health_rows, [("state", "State"), ("split", "Split"), ("rms", "Feature RMS"), ("median_abs", "Median |feature|"), ("near", "Near-constant"), ("rank", "Effective rank"), ("cond", "Reg. condition")]),
        "",
        "No pre-diagnostic feature normalization was applied. Full channel statistics and covariance spectra are stored in the numerical-health artifacts.",
        "",
        "## Temporal diagnostics",
        "",
        *markdown_table(temporal_rows, [("state", "State"), ("split", "Split"), ("frame", "Frame Δ"), ("within", "Within-window"), ("windows", "Between-window"), ("between", "Between-trajectory"), ("drift", "Chronological drift")]),
        "",
        "Frames and windows were retained as correlated observations; window and trajectory summaries were balanced at their stated grouping level.",
        "",
        "## R2-fit/R3-descriptive torsion probe",
        "",
        "The fixed probe used ridge alpha `1.0`. Each state used the identical R2-only scaling procedure and identical coordinate-derived masks.",
        "",
        *markdown_table(probe_rows, [("state", "State"), ("split", "Split"), ("component", "Mean component RMSE"), ("phi", "Phi MAE (deg)"), ("psi", "Psi MAE (deg)"), ("combined", "Combined MAE (deg)")]),
        "",
        "## Relative replay and rigid-transform diagnostics",
        "",
    ]
    for state in STATES:
        replay_max = max(x["state_replay"][state]["relative_rmse_to_split_feature_rms"] for x in replay["subsets"])
        rigid = replay["phase02_rigid_errors_relative_to_phase3_feature_rms"][state]
        lines.append(f"- {state}: max replay relative RMSE `{replay_max:.6g}`; rotation relative RMSE `{rigid['phase02_rotation_relative_to_phase3_R2_feature_rms']:.6g}`; translation relative RMSE `{rigid['phase02_translation_relative_to_phase3_R2_feature_rms']:.6g}`.")
    lines += [
        "",
        "## Representation comparison",
        "",
        f"C1 promising under the predeclared multi-criterion rule: `{c1['promising_structural_accessibility_advantage']}`. C2 promising: `{c2['promising_structural_accessibility_advantage']}`.",
        "",
        "The assessment considers R2 and R3 torsion accessibility, effective-rank retention, scaler behavior, replay error, and rigid-transform error together. It does not select a representation from a single metric.",
        "",
        "## Limitations",
        "",
        "R2 is the fit set, so its probe error is optimistic. R3 has already been historically inspected and is only a descriptive engineering evaluation. Frames and windows are correlated, the linear probe measures accessibility rather than total encoded information, runtime bitwise equivalence to the historical extraction is unknown, and these results do not support any biological efficacy claim.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--batch-audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    # This is written before extraction or representation metrics are viewed.
    frozen_config = dict(PHASE3_CONFIG)
    frozen_config["checkpoint"] = str(args.checkpoint.resolve())
    frozen_config["checkpoint_sha256"] = sha256(args.checkpoint)
    frozen_config["input_root"] = str(args.input_root.resolve())
    frozen_config["batch_audit"] = file_record(args.batch_audit)
    frozen_config["phase02_implementation"] = [file_record(path) for path in IMPLEMENTATION_FILES]
    write_json(args.output_root / "frozen_phase3_config.json", frozen_config)

    model, checkpoint_report = verify_and_strict_load(args.checkpoint)
    if not checkpoint_report["complete_strict_coverage"]:
        raise RuntimeError("Phase 0 checkpoint gate no longer passes")
    phase02_verification = json.loads((PHASE02_DIR / "verification_report.json").read_text(encoding="utf-8"))
    if not phase02_verification["source"]["source_equivalent"]:
        raise RuntimeError("Phase 0 source gate no longer passes")
    rows = load_authenticated_rows(args.batch_audit, args.input_root)
    disk = disk_preflight(rows, args.output_root)
    write_json(args.output_root / "disk_preflight.json", disk)

    extraction = extract_all(model, rows, args.output_root, args.device)
    write_json(
        args.output_root / "extraction_provenance.json",
        {
            "version": VERSION,
            "checkpoint_sha256": sha256(args.checkpoint),
            "source_commit": PINNED_GEOM2VEC_COMMIT,
            "runtime": runtime_report(args.device),
            **extraction,
        },
    )
    numerical, _, _ = numerical_health(rows, args.output_root)
    temporal, _ = temporal_diagnostics(rows, args.output_root)
    scalers = standardizers(rows)
    models = fit_ridge(rows, scalers)
    probe, _ = evaluate_probes(rows, scalers, models, args.output_root)
    replay = replay_checks(model, rows, numerical, args.output_root, args.device)
    comparison_result, comparison_table = comparison(numerical, temporal, probe, replay, args.output_root)

    cache_hash_valid = all(sha256(Path(item["cache_path"])) == item["cache_sha256"] for item in extraction["rows"])
    source_hash_valid = all(sha256(row["source"]) == row["source_sha256"] for row in rows)
    implementation_after = [file_record(path) for path in IMPLEMENTATION_FILES]
    implementation_unchanged = implementation_after == frozen_config["phase02_implementation"]
    checkpoint_unchanged = sha256(args.checkpoint) == EXPECTED_CHECKPOINT_SHA256
    production_consistent = all(item["c0_production_consistency"]["within_smoke_tolerance_2e-6"] for item in replay["subsets"])
    integrity = {
        "status": "PHASE3_COMPLETE",
        "phase3_followup_started": False,
        "checkpoint_unchanged": checkpoint_unchanged,
        "source_inputs_unchanged": source_hash_valid,
        "cache_hashes_verified": cache_hash_valid,
        "phase02_implementation_unchanged": implementation_unchanged,
        "phase02_implementation_after": implementation_after,
        "c0_production_consistency_within_smoke_tolerance": production_consistent,
        "extraction": {key: value for key, value in extraction.items() if key != "rows"},
        "prohibited_outputs_computed": [],
        "biological_labels_used": False,
    }
    complete = all((checkpoint_unchanged, source_hash_valid, cache_hash_valid, implementation_unchanged, production_consistent))
    if not complete:
        integrity["status"] = "PHASE3_INTEGRITY_FAILURE"
    write_json(args.output_root / "integrity_completion_audit.json", integrity)
    report_path = args.output_root / "REPORT_PHASE3.md"
    report_path.write_text(
        render_report(integrity, numerical, temporal, probe, comparison_result, comparison_table, replay),
        encoding="utf-8",
    )
    completion = {
        "status": integrity["status"],
        "report": file_record(report_path),
        "integrity_audit_sha256": sha256(args.output_root / "integrity_completion_audit.json"),
        "phase3_followup_started": False,
    }
    write_json(args.output_root / "PHASE3_COMPLETE.json", completion)
    print(json.dumps({"status": integrity["status"], "output": str(args.output_root)}, indent=2))
    if not complete:
        raise SystemExit("Phase 3 completed with an integrity failure")


if __name__ == "__main__":
    main()

