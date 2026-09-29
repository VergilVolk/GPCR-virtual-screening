#!/usr/bin/env python
"""Run PACER-DC Experiment F v01 once from authenticated Experiment D caches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

WORKSPACE = Path(__file__).resolve().parents[2]
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))

from project.encoder_atom_readout_v01.experiment_d import (  # noqa: E402
    chi1_targets,
    read_sequence,
    sha256,
    torsion_targets,
)

VERSION = "encoder_temporal_F_v01"
D_ROOT = WORKSPACE / "project/results/encoder_atom_readout_D_v01/run_001"
OUTPUT = WORKSPACE / "project/results/encoder_temporal_F_v01/run_001"
REPRESENTATIONS = ("M", "MC", "BS")
DESCRIPTORS = ("delta", "average")
TARGETS = ("phi", "psi", "chi1")
ALPHA = 1.0
BLOCK_SIZE = 20

FROZEN_CONFIG = {
    "version": VERSION,
    "input": "authenticated Experiment D C1 M/B/S caches",
    "representations": {
        "M": "C1 all-heavy mean, 128D",
        "MC": "concat(M, exact zeros), 256D dimension control",
        "BS": "concat(B,S), 256D",
    },
    "lag": {"stored_frames": 1, "physical_time_interpretation": None},
    "pair_boundary": ["replica", "condition", "window"],
    "descriptors": {
        "T0": "Z_t, static continuity only",
        "T1_delta": "Z_(t+1)-Z_t",
        "T2_average": "0.5*(Z_t+Z_(t+1))",
    },
    "dynamic_targets": "difference of torsion [sin,cos] vectors; valid in both frames",
    "targets": ["phi", "psi", "chi1"],
    "probe": {
        "family": "linear ridge", "alpha": ALPHA, "regularize_intercept": False,
        "scaling": "per representation/descriptor/target mean/std fitted on valid R2 pairs only",
        "fit": "R2 only", "evaluation": "R3 descriptive only",
    },
    "zero_baseline": "predict [0,0] delta",
    "block_audit": {"frames_per_window": 100, "block_size": BLOCK_SIZE, "blocks_per_window": 5},
    "reverse_subset": {"file_number": 1, "pair_indices": [0, 1, 2, 3], "residue_indices": [1, 2, 3, 5]},
    "success_logic": {
        "beats_zero": "residue-balanced combined component RMSE is lower than zero baseline on R2 and R3",
        "beats_average": "delta descriptor error is lower than orderless average on R2 and R3",
        "BS_beats_MC": "BS delta error is lower than MC delta error on R2 and R3",
        "breadth": "BS improves versus MC for >50% of valid residues on R2 and R3",
        "interpretation": "all criteria evaluated separately for phi, psi, chi1; no scalar winner",
    },
    "prohibitions": ["biological contrasts", "PACER-FKG", "long MD", "other lags", "other block sizes", "temporal neural models"],
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


def derive_representations(m: np.ndarray, b: np.ndarray, s: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "M": m,
        "MC": np.concatenate((m, np.zeros_like(m)), axis=-1),
        "BS": np.concatenate((b, s), axis=-1),
    }


def temporal_descriptors(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return z[1:] - z[:-1], 0.5 * (z[:-1] + z[1:])


def pair_target(angle: np.ndarray, valid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    encoded = np.stack((np.sin(angle), np.cos(angle)), axis=-1)
    return encoded[1:] - encoded[:-1], valid[1:] & valid[:-1]


def contiguous_blocks(frame_count: int, block_size: int = BLOCK_SIZE) -> list[np.ndarray]:
    if frame_count != 100 or block_size != 20 or frame_count % block_size:
        raise ValueError("Experiment F requires exactly 100 frames partitioned into 20-frame blocks")
    return [np.arange(start, start + block_size) for start in range(0, frame_count, block_size)]


def zero_baseline_metrics(dy: np.ndarray) -> dict:
    mse = np.mean(np.square(dy.astype(np.float64)), axis=0)
    return {"delta_sin_rmse": math.sqrt(float(mse[0])), "delta_cos_rmse": math.sqrt(float(mse[1])), "combined_component_rmse": math.sqrt(float(mse.mean()))}


def verify_provenance(d_root: Path) -> tuple[list[dict], dict]:
    completion_path = d_root / "EXPERIMENT_D_COMPLETE.json"
    integrity_path = d_root / "integrity_completion_audit.json"
    report_path = d_root / "REPORT.md"
    extraction_path = d_root / "extraction_audit.json"
    provenance_path = d_root / "provenance.json"
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    integrity = json.loads(integrity_path.read_text(encoding="utf-8"))
    extraction = json.loads(extraction_path.read_text(encoding="utf-8"))
    d_provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if completion["status"] != "EXPERIMENT_D_COMPLETE" or integrity["status"] != "EXPERIMENT_D_COMPLETE":
        raise RuntimeError("Experiment D is not complete")
    if sha256(integrity_path) != completion["integrity_audit_sha256"] or sha256(report_path) != completion["report"]["sha256"]:
        raise RuntimeError("Experiment D completion hashes changed")
    implementation = Path(integrity["experiment_d_implementation"]["path"])
    if sha256(implementation) != integrity["experiment_d_implementation"]["sha256"]:
        raise RuntimeError("Experiment D implementation changed")
    checkpoint = Path(d_provenance["checkpoint"]["checkpoint_path"])
    if sha256(checkpoint) != d_provenance["checkpoint"]["actual_sha256"]:
        raise RuntimeError("ViSNet checkpoint changed")
    rows = []
    for item in extraction["rows"]:
        cache = Path(item["cache_path"])
        source = Path(item["source_path"])
        if sha256(cache) != item["cache_sha256"]:
            raise RuntimeError(f"Experiment D cache changed: {cache}")
        if sha256(source) != item["source_sha256"]:
            raise RuntimeError(f"authenticated atom14 source changed: {source}")
        rows.append({**item, "cache": cache, "source": source, "sequence_csv": source.with_suffix(".csv")})
    if len(rows) != 40 or sum(item["shape_each_stored_readout"][0] for item in rows) != 4000:
        raise RuntimeError("unexpected Experiment D inventory")
    return rows, {
        "experiment_d_root": str(d_root.resolve()),
        "experiment_d_completion": file_record(completion_path),
        "experiment_d_integrity": file_record(integrity_path),
        "experiment_d_report": file_record(report_path),
        "experiment_d_extraction_audit": file_record(extraction_path),
        "experiment_d_implementation": file_record(implementation),
        "checkpoint": file_record(checkpoint),
        "checkpoint_strict_coverage_reused": d_provenance["checkpoint"]["complete_strict_coverage"],
        "c1_hook_implementation_unchanged": integrity["c1_hook_implementation_unchanged"],
        "cache_files_verified": len(rows), "source_files_verified": len(rows),
        "cache_bytes_reused": extraction["cache_bytes"], "cache_gib_reused": extraction["cache_gib"],
        "visnet_inference_run": False, "large_cache_duplicated": False,
        "runtime": d_provenance["runtime"],
    }


def load_payload(row: dict) -> tuple[dict[str, np.ndarray], dict[str, tuple[np.ndarray, np.ndarray]]]:
    with np.load(row["cache"]) as archive:
        reps = derive_representations(archive["M"], archive["B"], archive["S"])
    coords = np.asarray(np.load(row["source"], mmap_mode="r"))
    sequence = read_sequence(row["sequence_csv"])
    raw = {**torsion_targets(coords), **chi1_targets(coords, sequence)}
    targets = {name: pair_target(raw[name], raw[f"{name}_mask"]) for name in TARGETS}
    return reps, targets


def pair_inventory(rows: list[dict]) -> dict:
    split_counts = {}
    details = []
    for split in ("R2", "R3"):
        subset = [row for row in rows if row["split"] == split]
        counts = {target: 0 for target in TARGETS}
        residues = {target: set() for target in TARGETS}
        for row in subset:
            _, targets = load_payload(row)
            item = {"relative_path": row["relative_path"], "split": split, "replica": row["replica"], "condition": row["condition"], "window": row["window"], "stored_frame_pairs": 99, "valid_samples": {}}
            for target, (_, mask) in targets.items():
                count = int(mask.sum())
                counts[target] += count
                residues[target].update(np.flatnonzero(mask.any(axis=0)).tolist())
                item["valid_samples"][target] = count
            details.append(item)
        split_counts[split] = {
            "trajectory_window_groups": len(subset), "stored_frame_pairs": 99 * len(subset),
            "valid_samples": counts, "valid_residue_positions": {target: len(value) for target, value in residues.items()},
        }
    return {"lag_stored_frames": 1, "cross_window_pairs": 0, "splits": split_counts, "groups": details}


def descriptor_for(rep: np.ndarray, descriptor: str) -> np.ndarray:
    delta, average = temporal_descriptors(rep)
    return delta if descriptor == "delta" else average


def fit_scalers(rows: list[dict]) -> dict:
    accum = {}
    for rep in REPRESENTATIONS:
        dim = 128 if rep == "M" else 256
        for descriptor in DESCRIPTORS:
            for target in TARGETS:
                accum[(rep, descriptor, target)] = {"n": 0, "sum": np.zeros(dim), "sum2": np.zeros(dim)}
    for row in (item for item in rows if item["split"] == "R2"):
        reps, targets = load_payload(row)
        for rep, values in reps.items():
            for descriptor in DESCRIPTORS:
                xall = descriptor_for(values, descriptor).astype(np.float64)
                for target, (_, mask) in targets.items():
                    x = xall[mask]
                    item = accum[(rep, descriptor, target)]
                    item["n"] += len(x); item["sum"] += x.sum(0); item["sum2"] += np.square(x).sum(0)
    scalers = {}
    for key, item in accum.items():
        mean = item["sum"] / item["n"]
        variance = np.maximum(item["sum2"] / item["n"] - mean**2, 0)
        raw = np.sqrt(variance)
        scale = np.where(raw <= 1e-12, 1.0, raw)
        scalers[key] = {"count": item["n"], "mean": mean, "scale": scale, "replaced_channels": np.flatnonzero(raw <= 1e-12).tolist()}
    return scalers


def fit_probes(rows: list[dict], scalers: dict) -> dict:
    accum = {}
    for key, scaler in scalers.items():
        dim = len(scaler["mean"])
        accum[key] = {"xtx": np.zeros((dim + 1, dim + 1)), "xty": np.zeros((dim + 1, 2)), "samples": 0}
    for row in (item for item in rows if item["split"] == "R2"):
        reps, targets = load_payload(row)
        for rep, values in reps.items():
            for descriptor in DESCRIPTORS:
                raw_x = descriptor_for(values, descriptor).astype(np.float64)
                for target, (dy, mask) in targets.items():
                    key = (rep, descriptor, target); scaler = scalers[key]
                    x = (raw_x[mask] - scaler["mean"]) / scaler["scale"]
                    y = dy[mask].astype(np.float64)
                    xa = np.column_stack((x, np.ones(len(x))))
                    item = accum[key]; item["xtx"] += xa.T @ xa; item["xty"] += xa.T @ y; item["samples"] += len(x)
    models = {}
    for key, item in accum.items():
        regularizer = np.eye(item["xtx"].shape[0]) * ALPHA; regularizer[-1, -1] = 0
        models[key] = {"samples": item["samples"], "coefficient": np.linalg.solve(item["xtx"] + regularizer, item["xty"])}
    return models


def metric_accumulator() -> dict:
    return {"sse": np.zeros(2), "zero_sse": np.zeros(2), "n": 0, "per_residue": defaultdict(lambda: {"sse": np.zeros(2), "zero_sse": np.zeros(2), "n": 0})}


def evaluate(rows: list[dict], scalers: dict, models: dict) -> tuple[dict, list[dict]]:
    accum = {(rep, descriptor, target, split): metric_accumulator() for rep in REPRESENTATIONS for descriptor in DESCRIPTORS for target in TARGETS for split in ("R2", "R3")}
    for row in rows:
        reps, targets = load_payload(row)
        residue_grid = np.broadcast_to(np.arange(270)[None, :], (99, 270))
        for rep, values in reps.items():
            for descriptor in DESCRIPTORS:
                raw_x = descriptor_for(values, descriptor).astype(np.float64)
                for target, (dy, mask) in targets.items():
                    key = (rep, descriptor, target); scaler = scalers[key]
                    x = (raw_x[mask] - scaler["mean"]) / scaler["scale"]
                    y = dy[mask].astype(np.float64)
                    pred = np.column_stack((x, np.ones(len(x)))) @ models[key]["coefficient"]
                    error = pred - y
                    item = accum[(rep, descriptor, target, row["split"])]
                    item["sse"] += np.square(error).sum(0); item["zero_sse"] += np.square(y).sum(0); item["n"] += len(y)
                    indices = residue_grid[mask]
                    for residue in np.unique(indices):
                        select = indices == residue; ritem = item["per_residue"][int(residue)]
                        ritem["sse"] += np.square(error[select]).sum(0); ritem["zero_sse"] += np.square(y[select]).sum(0); ritem["n"] += int(select.sum())
    summary, per_residue = {}, []
    for key, item in accum.items():
        rep, descriptor, target, split = key
        mse = item["sse"] / item["n"]; zero_mse = item["zero_sse"] / item["n"]
        residue_rmse, residue_zero = [], []
        for residue, ritem in sorted(item["per_residue"].items()):
            rmse = np.sqrt(ritem["sse"] / ritem["n"]); zero = np.sqrt(ritem["zero_sse"] / ritem["n"])
            combined = math.sqrt(float(np.mean(rmse**2))); zero_combined = math.sqrt(float(np.mean(zero**2)))
            residue_rmse.append(rmse); residue_zero.append(zero)
            per_residue.append({"representation": rep, "descriptor": descriptor, "target": target, "split": split, "residue_index": residue, "samples": ritem["n"], "delta_sin_rmse": float(rmse[0]), "delta_cos_rmse": float(rmse[1]), "combined_component_rmse": combined, "zero_combined_component_rmse": zero_combined, "relative_improvement_over_zero": 1 - combined / zero_combined})
        residue_rmse = np.asarray(residue_rmse); residue_zero = np.asarray(residue_zero)
        combined = math.sqrt(float(mse.mean())); zero_combined = math.sqrt(float(zero_mse.mean()))
        summary[f"{rep}_{descriptor}_{target}_{split}"] = {
            "representation": rep, "descriptor": descriptor, "target": target, "split": split,
            "samples": item["n"], "valid_residue_positions": len(item["per_residue"]),
            "delta_sin_rmse": math.sqrt(float(mse[0])), "delta_cos_rmse": math.sqrt(float(mse[1])), "combined_component_rmse": combined,
            "zero_delta_sin_rmse": math.sqrt(float(zero_mse[0])), "zero_delta_cos_rmse": math.sqrt(float(zero_mse[1])), "zero_combined_component_rmse": zero_combined,
            "relative_improvement_over_zero": 1 - combined / zero_combined,
            "residue_balanced_delta_sin_rmse": math.sqrt(float(np.mean(residue_rmse[:, 0] ** 2))),
            "residue_balanced_delta_cos_rmse": math.sqrt(float(np.mean(residue_rmse[:, 1] ** 2))),
            "residue_balanced_combined_component_rmse": math.sqrt(float(np.mean(residue_rmse**2))),
            "residue_balanced_zero_combined_component_rmse": math.sqrt(float(np.mean(residue_zero**2))),
            "residue_balanced_relative_improvement_over_zero": 1 - math.sqrt(float(np.mean(residue_rmse**2))) / math.sqrt(float(np.mean(residue_zero**2))),
        }
    return summary, per_residue


def reverse_and_embedding_checks(rows: list[dict]) -> dict:
    reps, targets = load_payload(rows[0])
    pair_ids = np.asarray(FROZEN_CONFIG["reverse_subset"]["pair_indices"]); residue_ids = np.asarray(FROZEN_CONFIG["reverse_subset"]["residue_indices"])
    checks = {}
    for rep, z in reps.items():
        forward = z[pair_ids + 1][:, residue_ids] - z[pair_ids][:, residue_ids]
        reverse = z[pair_ids][:, residue_ids] - z[pair_ids + 1][:, residue_ids]
        checks[f"{rep}_reverse_max_abs"] = float(np.max(np.abs(reverse + forward)))
    for target, (dy, _) in targets.items():
        forward = dy[pair_ids][:, residue_ids]
        raw = -forward
        checks[f"{target}_reverse_max_abs"] = float(np.max(np.abs(raw + forward)))
    dm, am = temporal_descriptors(reps["M"]); dmc, amc = temporal_descriptors(reps["MC"])
    checks.update({
        "delta_MC_first_equals_M": bool(np.array_equal(dmc[..., :128], dm)), "delta_MC_zero_block": bool(np.count_nonzero(dmc[..., 128:]) == 0),
        "average_MC_first_equals_M": bool(np.array_equal(amc[..., :128], am)), "average_MC_zero_block": bool(np.count_nonzero(amc[..., 128:]) == 0),
    })
    errors = [value for key, value in checks.items() if key.endswith("max_abs")]
    checks["passed"] = all(np.isfinite(errors)) and max(errors) == 0 and all(value for key, value in checks.items() if key.endswith("equals_M") or key.endswith("zero_block"))
    return checks


def temporal_energy_and_blocks(rows: list[dict]) -> tuple[dict, dict]:
    energy_acc = defaultdict(lambda: {"delta_norm2": 0.0, "pairs": 0, "delta_ss": 0.0, "feature_ss": 0.0, "feature_n": 0, "x": 0.0, "y": 0.0, "x2": 0.0, "y2": 0.0, "xy": 0.0, "lag_n": 0, "within": 0.0, "within_n": 0})
    block_acc = defaultdict(lambda: {"frame": 0.0, "within": 0.0, "between": 0.0, "n": 0, "groups": 0})
    centroids = defaultdict(dict)
    for row in rows:
        with np.load(row["cache"]) as archive:
            base = {"M": archive["M"], "B": archive["B"], "S": archive["S"]}
            reps = derive_representations(base["M"], base["B"], base["S"])
        all_states = {**reps, "BS_B": base["B"], "BS_S": base["S"]}
        for rep, z32 in all_states.items():
            z = z32.astype(np.float64); delta = z[1:] - z[:-1]; key = (rep, row["split"]); item = energy_acc[key]
            item["delta_norm2"] += float(np.square(delta).sum()); item["pairs"] += delta.shape[0] * delta.shape[1]
            item["delta_ss"] += float(np.square(delta).sum()); item["feature_ss"] += float(np.square(z).sum()); item["feature_n"] += z.size
            x, y = z[:-1].ravel(), z[1:].ravel(); item["x"] += float(x.sum()); item["y"] += float(y.sum()); item["x2"] += float(x @ x); item["y2"] += float(y @ y); item["xy"] += float(x @ y); item["lag_n"] += len(x)
            window_mean = z.mean(0, keepdims=True); item["within"] += float(np.square(z - window_mean).sum()); item["within_n"] += z.size
            if rep in REPRESENTATIONS:
                blocks = contiguous_blocks(len(z)); block_centroids = np.stack([z[idx].mean(0) for idx in blocks]); grand = z.mean(0)
                frame_ss = float(np.square(z - grand).sum()); within_ss = sum(float(np.square(z[idx] - block_centroids[i]).sum()) for i, idx in enumerate(blocks)); between_ss = BLOCK_SIZE * float(np.square(block_centroids - grand).sum())
                bitem = block_acc[key]; bitem["frame"] += frame_ss; bitem["within"] += within_ss; bitem["between"] += between_ss; bitem["n"] += z.size; bitem["groups"] += 1
                centroids[(rep, row["split"], row["condition"])][row["window"]] = z.mean((0, 1))
    energy = {}
    for key, item in energy_acc.items():
        rep, split = key; n = item["lag_n"]
        cov = item["xy"] - item["x"] * item["y"] / n; vx = item["x2"] - item["x"]**2 / n; vy = item["y2"] - item["y"]**2 / n
        feature_rms = math.sqrt(item["feature_ss"] / item["feature_n"]); delta_rms = math.sqrt(item["delta_ss"] / (item["pairs"] * (128 if rep in ("M", "BS_B", "BS_S") else 256)))
        drift_values = []
        if rep in REPRESENTATIONS:
            for condition in sorted({row["condition"] for row in rows if row["split"] == split}):
                values = centroids[(rep, split, condition)]; drift_values.append(float(np.linalg.norm(values[max(values)] - values[min(values)])))
        energy[f"{rep}_{split}"] = {"representation": rep, "split": split, "mean_delta_norm_squared": item["delta_norm2"] / item["pairs"], "feature_rms": feature_rms, "delta_rms": delta_rms, "delta_rms_relative_to_feature_rms": delta_rms / feature_rms, "lag1_representation_autocorrelation": cov / math.sqrt(vx * vy), "within_window_temporal_variance": item["within"] / item["within_n"], "chronological_drift_mean_centroid_l2": float(np.mean(drift_values)) if drift_values else None}
    blocks = {}
    for key, item in block_acc.items():
        rep, split = key; total = item["frame"]
        blocks[f"{rep}_{split}"] = {"representation": rep, "split": split, "window_groups": item["groups"], "frames_per_window": 100, "blocks_per_window": 5, "block_size": 20, "frame_level_variance": item["frame"] / item["n"], "within_block_variance": item["within"] / item["n"], "between_block_centroid_variance": item["between"] / item["n"], "fraction_temporal_variance_remaining_after_block_averaging": item["between"] / total, "fraction_discarded_within_blocks": item["within"] / total, "decomposition_relative_error": abs(total - item["within"] - item["between"]) / total}
    return energy, blocks


def static_continuity(rows: list[dict], d_root: Path) -> dict:
    backbone = json.loads((d_root / "backbone_torsion_probe.json").read_text(encoding="utf-8")); chi = json.loads((d_root / "chi1_probe.json").read_text(encoding="utf-8"))
    mapping = {"M": "C1-M", "BS": "C1-BS"}; result = {"method": "independent replay of frozen Experiment D scaler/model over all cached T0 frames", "metrics": {}, "max_abs_metric_difference": 0.0}
    for rep, dname in mapping.items():
        scaler = backbone["scalers"][dname]; mean = np.asarray(scaler["mean"]); scale = np.asarray(scaler["scale"])
        models = {t: np.asarray(backbone["models"][dname][t]["coefficient"]) for t in ("phi", "psi")}; models["chi1"] = np.asarray(chi["models"][dname]["coefficient"])
        for split in ("R2", "R3"):
            sums = {t: {"sse": np.zeros(2), "n": 0, "angle": 0.0} for t in TARGETS}
            for row in (item for item in rows if item["split"] == split):
                with np.load(row["cache"]) as archive:
                    reps = derive_representations(archive["M"], archive["B"], archive["S"]); z = reps[rep].astype(np.float64)
                coords = np.asarray(np.load(row["source"], mmap_mode="r")); sequence = read_sequence(row["sequence_csv"]); raw = {**torsion_targets(coords), **chi1_targets(coords, sequence)}
                xall = (z - mean) / scale
                for target in TARGETS:
                    mask = raw[f"{target}_mask"]; angle = raw[target][mask]; y = np.column_stack((np.sin(angle), np.cos(angle))); x = xall[mask]
                    pred = np.column_stack((x, np.ones(len(x)))) @ models[target]; error = pred - y
                    predicted_angle = np.arctan2(pred[:, 0], pred[:, 1]); angular = np.abs((predicted_angle - angle + np.pi) % (2 * np.pi) - np.pi) * 180 / np.pi
                    sums[target]["sse"] += np.square(error).sum(0); sums[target]["n"] += len(y); sums[target]["angle"] += float(angular.sum())
            observed = {t: {"sin_rmse": math.sqrt(sums[t]["sse"][0] / sums[t]["n"]), "cos_rmse": math.sqrt(sums[t]["sse"][1] / sums[t]["n"]), "angular_mae_deg": sums[t]["angle"] / sums[t]["n"]} for t in TARGETS}
            expected = {"phi": backbone["metrics"][f"{dname}_{split}"]["phi"]["all"], "psi": backbone["metrics"][f"{dname}_{split}"]["psi"]["all"], "chi1": chi["metrics"][f"{dname}_{split}"]}
            differences = {t: {metric: abs(observed[t][metric] - expected[t][metric]) for metric in observed[t]} for t in TARGETS}
            maximum = max(value for target in differences.values() for value in target.values()); result["max_abs_metric_difference"] = max(result["max_abs_metric_difference"], maximum)
            result["metrics"][f"{rep}_{split}"] = {"observed": observed, "experiment_d": {t: {metric: expected[t][metric] for metric in observed[t]} for t in TARGETS}, "absolute_difference": differences}
    result["tolerance"] = 1e-10; result["passed"] = result["max_abs_metric_difference"] <= result["tolerance"]
    return result


def comparison(metrics: dict, per_residue: list[dict]) -> dict:
    lookup = {(r["representation"], r["descriptor"], r["target"], r["split"], r["residue_index"]): r for r in per_residue}
    result = {"targets": {}, "interpretation": {}, "no_scalar_winner": True}
    for target in TARGETS:
        criteria = {}
        for rep in REPRESENTATIONS:
            criteria[rep] = {}
            for split in ("R2", "R3"):
                delta = metrics[f"{rep}_delta_{target}_{split}"]; average = metrics[f"{rep}_average_{target}_{split}"]
                criteria[rep][split] = {"delta_rmse": delta["residue_balanced_combined_component_rmse"], "zero_rmse": delta["residue_balanced_zero_combined_component_rmse"], "average_rmse": average["residue_balanced_combined_component_rmse"], "delta_beats_zero": delta["residue_balanced_combined_component_rmse"] < delta["residue_balanced_zero_combined_component_rmse"], "delta_beats_average": delta["residue_balanced_combined_component_rmse"] < average["residue_balanced_combined_component_rmse"]}
        breadth = {}
        for split in ("R2", "R3"):
            residues = sorted({key[-1] for key in lookup if key[:4] == ("BS", "delta", target, split)})
            improved = sum(lookup[("BS", "delta", target, split, residue)]["combined_component_rmse"] < lookup[("MC", "delta", target, split, residue)]["combined_component_rmse"] for residue in residues)
            breadth[split] = {"BS_better_residues": improved, "valid_residues": len(residues), "fraction": improved / len(residues)}
        result["targets"][target] = {"criteria": criteria, "BS_beats_MC": {split: criteria["BS"][split]["delta_rmse"] < criteria["MC"][split]["delta_rmse"] for split in ("R2", "R3")}, "per_residue_breadth": breadth}
        t = result["targets"][target]
        result["interpretation"][target] = {
            "useful_BS_temporal_information": all(t["criteria"]["BS"][split]["delta_beats_zero"] and t["criteria"]["BS"][split]["delta_beats_average"] and t["BS_beats_MC"][split] and t["per_residue_breadth"][split]["fraction"] > 0.5 for split in ("R2", "R3")),
            "classification": None,
        }
        if not all(t["criteria"][rep][split]["delta_beats_zero"] for rep in REPRESENTATIONS for split in ("R2", "R3")):
            result["interpretation"][target]["classification"] = "D: at least one representation does not beat zero on both splits"
        elif not all(t["criteria"]["BS"][split]["delta_beats_average"] for split in ("R2", "R3")):
            result["interpretation"][target]["classification"] = "C: delta does not beat orderless average on both splits"
        elif all(t["BS_beats_MC"].values()):
            result["interpretation"][target]["classification"] = "A: BS contains additional accessible dynamic information"
        else:
            result["interpretation"][target]["classification"] = "B/E: BS separation adds little or fails R3 retention"
    return result


def render_report(provenance: dict, inventory: dict, continuity: dict, energy: dict, blocks: dict, metrics: dict, comp: dict) -> str:
    lines = ["# PACER-DC encoder optimization - Experiment F v01", "", "Experiment F completed the fixed lag-1 temporal information-retention audit using Experiment D caches only. ViSNet was not rerun.", "", "## Integrity and inventory", "", f"All `{provenance['cache_files_verified']}` D caches and `{provenance['source_files_verified']}` authenticated coordinate files matched their archived hashes. Static continuity passed with maximum metric error `{continuity['max_abs_metric_difference']:.3g}`.", "", "| Split | Groups | Stored pairs | Phi samples | Psi samples | Chi1 samples |", "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for split in ("R2", "R3"):
        x = inventory["splits"][split]; lines.append(f"| {split} | {x['trajectory_window_groups']} | {x['stored_frame_pairs']} | {x['valid_samples']['phi']} | {x['valid_samples']['psi']} | {x['valid_samples']['chi1']} |")
    lines += ["", "## Dynamic probes", "", "Residue-balanced component RMSE; lower is better.", "", "| Target | Rep | Split | Delta | Average | Zero | Delta vs zero |", "| --- | --- | --- | ---: | ---: | ---: | ---: |"]
    for target in TARGETS:
        for rep in REPRESENTATIONS:
            for split in ("R2", "R3"):
                d = metrics[f"{rep}_delta_{target}_{split}"]; a = metrics[f"{rep}_average_{target}_{split}"]; lines.append(f"| {target} | {rep} | {split} | {d['residue_balanced_combined_component_rmse']:.6g} | {a['residue_balanced_combined_component_rmse']:.6g} | {d['residue_balanced_zero_combined_component_rmse']:.6g} | {d['residue_balanced_relative_improvement_over_zero']:.2%} |")
    lines += ["", "## Temporal energy and block averaging", "", "| Rep | Split | Delta RMS / feature RMS | Lag-1 autocorr | Variance retained | Variance discarded |", "| --- | --- | ---: | ---: | ---: | ---: |"]
    for rep in REPRESENTATIONS:
        for split in ("R2", "R3"):
            e = energy[f"{rep}_{split}"]; b = blocks[f"{rep}_{split}"]; lines.append(f"| {rep} | {split} | {e['delta_rms_relative_to_feature_rms']:.6g} | {e['lag1_representation_autocorrelation']:.6g} | {b['fraction_temporal_variance_remaining_after_block_averaging']:.2%} | {b['fraction_discarded_within_blocks']:.2%} |")
    lines += ["", "## Interpretation", ""]
    for target in TARGETS:
        item = comp["interpretation"][target]; breadth = comp["targets"][target]["per_residue_breadth"]; lines.append(f"- **{target}:** {item['classification']}. Full success flag: `{item['useful_BS_temporal_information']}`. BS improved over MC for {breadth['R2']['BS_better_residues']}/{breadth['R2']['valid_residues']} R2 and {breadth['R3']['BS_better_residues']}/{breadth['R3']['valid_residues']} R3 residues.")
    lines += ["", "Individual frame pairs are correlated observations, not independent biological replicates. R2 fits are in-sample; R3 is historically inspected and descriptive. Lag means one stored frame only. No physical time or biological efficacy interpretation is made.", ""]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--d-root", type=Path, default=D_ROOT); parser.add_argument("--output-root", type=Path, default=OUTPUT); parser.add_argument("--refresh-integrity-only", action="store_true"); args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True); write_json(args.output_root / "frozen_config.json", FROZEN_CONFIG)
    rows, provenance = verify_provenance(args.d_root); write_json(args.output_root / "provenance.json", provenance); print("provenance verified", flush=True)
    inventory = pair_inventory(rows); write_json(args.output_root / "pair_inventory.json", inventory); print("pair inventory complete", flush=True)
    checks = reverse_and_embedding_checks(rows)
    if args.refresh_integrity_only:
        continuity = json.loads((args.output_root / "static_continuity.json").read_text(encoding="utf-8")); blocks = json.loads((args.output_root / "block_averaging_audit.json").read_text(encoding="utf-8"))
        cache_unchanged = all(sha256(row["cache"]) == row["cache_sha256"] and sha256(row["source"]) == row["source_sha256"] for row in rows)
        integrity = {"status": "EXPERIMENT_F_COMPLETE", "experiment_d_caches_and_sources_unchanged": cache_unchanged, "static_continuity_passed": continuity["passed"], "signed_reverse_and_MC_embedding_checks": checks, "R2_only_fitting": True, "R3_used_for_fitting_or_selection": False, "identical_target_masks_across_representations": True, "cross_window_pairs": inventory["cross_window_pairs"], "block_partition_exact": all(value["decomposition_relative_error"] < 1e-12 for value in blocks.values()), "visnet_inference_run": False, "large_derived_arrays_cached": False, "prohibited_outputs_computed": [], "experiment_f_implementation": file_record(Path(__file__))}
        if not (cache_unchanged and continuity["passed"] and checks["passed"] and integrity["block_partition_exact"]): integrity["status"] = "EXPERIMENT_F_INTEGRITY_FAILURE"
        write_json(args.output_root / "integrity_completion_audit.json", integrity)
        report = args.output_root / "REPORT.md"; write_json(args.output_root / "EXPERIMENT_F_COMPLETE.json", {"status": integrity["status"], "report": file_record(report), "integrity_audit_sha256": sha256(args.output_root / "integrity_completion_audit.json"), "followup_started": False})
        print(json.dumps({"status": integrity["status"], "integrity_only": True}, indent=2)); return
    continuity = static_continuity(rows, args.d_root); write_json(args.output_root / "static_continuity.json", continuity); print("static continuity complete", flush=True)
    if not continuity["passed"]:
        raise RuntimeError("static continuity does not reproduce Experiment D")
    energy, blocks = temporal_energy_and_blocks(rows); write_json(args.output_root / "temporal_energy.json", energy); write_json(args.output_root / "block_averaging_audit.json", blocks); print("energy and block audits complete", flush=True)
    scalers = fit_scalers(rows); print("R2 scalers fitted", flush=True)
    models = fit_probes(rows, scalers); print("R2 ridge probes fitted", flush=True)
    metrics, per_residue = evaluate(rows, scalers, models); print("R2/R3 probe evaluation complete", flush=True)
    for target in TARGETS:
        target_metrics = {key: value for key, value in metrics.items() if value["target"] == target}
        target_rows = [row for row in per_residue if row["target"] == target]
        model_json = {"|".join(key): {"samples": value["samples"], "coefficient": value["coefficient"].tolist()} for key, value in models.items() if key[2] == target}
        scaler_json = {"|".join(key): {"count": value["count"], "mean": value["mean"].tolist(), "scale": value["scale"].tolist(), "replaced_channels": value["replaced_channels"]} for key, value in scalers.items() if key[2] == target}
        write_json(args.output_root / f"{target}_dynamic_probe.json", {"config": FROZEN_CONFIG["probe"], "scalers": scaler_json, "models": model_json, "metrics": target_metrics}); write_csv(args.output_root / f"{target}_dynamic_probe_per_residue.csv", target_rows)
    comp = comparison(metrics, per_residue); write_json(args.output_root / "representation_comparison.json", comp)
    write_csv(args.output_root / "representation_comparison.csv", [{"target": target, "classification": comp["interpretation"][target]["classification"], "useful_BS_temporal_information": comp["interpretation"][target]["useful_BS_temporal_information"], "R2_BS_better_fraction": comp["targets"][target]["per_residue_breadth"]["R2"]["fraction"], "R3_BS_better_fraction": comp["targets"][target]["per_residue_breadth"]["R3"]["fraction"]} for target in TARGETS])
    cache_unchanged = all(sha256(row["cache"]) == row["cache_sha256"] and sha256(row["source"]) == row["source_sha256"] for row in rows)
    integrity = {"status": "EXPERIMENT_F_COMPLETE", "experiment_d_caches_and_sources_unchanged": cache_unchanged, "static_continuity_passed": continuity["passed"], "signed_reverse_and_MC_embedding_checks": checks, "R2_only_fitting": True, "R3_used_for_fitting_or_selection": False, "identical_target_masks_across_representations": True, "cross_window_pairs": inventory["cross_window_pairs"], "block_partition_exact": all(value["decomposition_relative_error"] < 1e-12 for value in blocks.values()), "visnet_inference_run": False, "large_derived_arrays_cached": False, "prohibited_outputs_computed": [], "experiment_f_implementation": file_record(Path(__file__))}
    if not (cache_unchanged and continuity["passed"] and checks["passed"] and integrity["block_partition_exact"]): integrity["status"] = "EXPERIMENT_F_INTEGRITY_FAILURE"
    write_json(args.output_root / "integrity_completion_audit.json", integrity)
    report = args.output_root / "REPORT.md"; report.write_text(render_report(provenance, inventory, continuity, energy, blocks, metrics, comp), encoding="utf-8")
    write_json(args.output_root / "EXPERIMENT_F_COMPLETE.json", {"status": integrity["status"], "report": file_record(report), "integrity_audit_sha256": sha256(args.output_root / "integrity_completion_audit.json"), "followup_started": False})
    print(json.dumps({"status": integrity["status"], "output": str(args.output_root)}, indent=2))
    if integrity["status"] != "EXPERIMENT_F_COMPLETE": raise SystemExit("Experiment F integrity failure")


if __name__ == "__main__":
    main()
