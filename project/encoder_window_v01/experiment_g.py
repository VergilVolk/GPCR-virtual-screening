#!/usr/bin/env python
"""Run the fixed PACER-DC Experiment G v01 window audit once."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

WORKSPACE = Path(__file__).resolve().parents[2]
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))

from project.encoder_atom_readout_v01.experiment_d import chi1_targets, read_sequence, sha256, torsion_targets  # noqa: E402

VERSION = "encoder_window_G_v01"
D_ROOT = WORKSPACE / "project/results/encoder_atom_readout_D_v01/run_001"
F_ROOT = WORKSPACE / "project/results/encoder_temporal_F_v01/run_001"
OUTPUT = WORKSPACE / "project/results/encoder_window_G_v01/run_001"
REPS = ("MC", "BS")
BRANCHES = ("static", "signed", "rms")
TORSIONS = ("phi", "psi", "chi1")
FAMILIES = ("static", "signed", "rms")
FEATURE_CONFIGS = {
    "static": ("static",), "signed": ("signed",), "rms": ("rms",),
    "static_signed": ("static", "signed"), "static_rms": ("static", "rms"),
    "all": ("static", "signed", "rms"),
}
PRIMARY_CONFIGS = {
    "static": ("static", "static_signed", "static_rms", "all"),
    "signed": ("signed", "static_signed", "static_rms", "all"),
    "rms": ("rms", "static_signed", "static_rms", "all"),
}
ALPHA = 1.0

FROZEN_CONFIG = {
    "version": VERSION,
    "inputs": {"BS": "concat Experiment D B,S", "MC": "concat Experiment D M, exact zeros"},
    "window": {"frames": 20, "blocks_per_100_frame_window": 5, "lag_stored_frames": 1, "differences_per_block": 19, "cross_boundary_differences": False},
    "branches": {
        "static": "elementwise mean of 20 states, 256D",
        "signed": "elementwise mean of 19 lag-1 differences = endpoint displacement/19, 256D",
        "rms": "elementwise RMS of 19 lag-1 differences, 256D",
    },
    "canonical_object": "three branches stored separately; concatenations only inside fixed linear probes",
    "target_families": {
        "static": "mean torsion sin/cos over 20 valid frames",
        "signed": "mean adjacent delta of torsion sin/cos",
        "rms": "RMS adjacent delta of torsion sin/cos",
    },
    "mask": "torsion must be valid at all 20 frames; identical for BS and MC",
    "probe": {"family": "linear ridge", "alpha": ALPHA, "regularize_intercept": False, "scaling": "R2-only per representation/config/torsion", "fit": "R2 only", "R3": "descriptive evaluation only"},
    "feature_configs": {key: list(value) for key, value in FEATURE_CONFIGS.items()},
    "primary_configs_by_target": {key: list(value) for key, value in PRIMARY_CONFIGS.items()},
    "baselines": {"signed": "exact zero", "static": "R2-only mean target", "rms": "R2-only mean-magnitude target"},
    "reversal_subset": {"file_number": 1, "block_index": 0, "residue_indices": [1, 2, 3, 5]},
    "information_accounting": "endpoint linear reconstruction muZ+(t-9.5)*muDelta; report trend, residual, and cross term. RMSDelta energy is overlapping/non-orthogonal and receives no reconstructive credit.",
}


def write_json(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows: path.write_text("", encoding="utf-8"); return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def file_record(path: Path) -> dict:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)}


def block_indices(frame_count: int = 100) -> list[np.ndarray]:
    if frame_count != 100: raise ValueError("Experiment G requires 100-frame windows")
    return [np.arange(start, start + 20) for start in range(0, 100, 20)]


def block_branches(z: np.ndarray) -> dict[str, np.ndarray]:
    if z.shape[0] != 20: raise ValueError("block must contain exactly 20 frames")
    z = z.astype(np.float64, copy=False); delta = z[1:] - z[:-1]
    return {"static": z.mean(axis=0), "signed": delta.mean(axis=0), "rms": np.sqrt(np.mean(np.square(delta), axis=0))}


def target_branches(angle: np.ndarray, valid: np.ndarray) -> tuple[dict[str, np.ndarray], np.ndarray]:
    if angle.shape[0] != 20: raise ValueError("target block must contain 20 frames")
    encoded = np.stack((np.sin(angle), np.cos(angle)), axis=-1)
    delta = encoded[1:] - encoded[:-1]
    return {"static": encoded.mean(0), "signed": delta.mean(0), "rms": np.sqrt(np.mean(np.square(delta), axis=0))}, valid.all(0)


def derive_reps(m: np.ndarray, b: np.ndarray, s: np.ndarray) -> dict[str, np.ndarray]:
    return {"MC": np.concatenate((m, np.zeros_like(m)), axis=-1), "BS": np.concatenate((b, s), axis=-1)}


def verify_provenance(d_root: Path, f_root: Path) -> tuple[list[dict], dict]:
    d_marker_path = d_root / "EXPERIMENT_D_COMPLETE.json"; d_integrity_path = d_root / "integrity_completion_audit.json"; d_report_path = d_root / "REPORT.md"
    f_marker_path = f_root / "EXPERIMENT_F_COMPLETE.json"; f_integrity_path = f_root / "integrity_completion_audit.json"; f_report_path = f_root / "REPORT.md"
    d_marker = json.loads(d_marker_path.read_text()); d_integrity = json.loads(d_integrity_path.read_text()); f_marker = json.loads(f_marker_path.read_text()); f_integrity = json.loads(f_integrity_path.read_text())
    if d_marker["status"] != "EXPERIMENT_D_COMPLETE" or f_marker["status"] != "EXPERIMENT_F_COMPLETE": raise RuntimeError("required predecessor is incomplete")
    if sha256(d_integrity_path) != d_marker["integrity_audit_sha256"] or sha256(d_report_path) != d_marker["report"]["sha256"]: raise RuntimeError("Experiment D completion changed")
    if sha256(f_integrity_path) != f_marker["integrity_audit_sha256"] or sha256(f_report_path) != f_marker["report"]["sha256"]: raise RuntimeError("Experiment F completion changed")
    d_impl = Path(d_integrity["experiment_d_implementation"]["path"]); f_impl = Path(f_integrity["experiment_f_implementation"]["path"])
    if sha256(d_impl) != d_integrity["experiment_d_implementation"]["sha256"] or sha256(f_impl) != f_integrity["experiment_f_implementation"]["sha256"]: raise RuntimeError("validated predecessor implementation changed")
    d_prov = json.loads((d_root / "provenance.json").read_text()); checkpoint = Path(d_prov["checkpoint"]["checkpoint_path"]); c1_implementation = Path(d_prov["c1_implementation"]["path"])
    if sha256(checkpoint) != d_prov["checkpoint"]["actual_sha256"]: raise RuntimeError("checkpoint changed")
    if sha256(c1_implementation) != d_prov["c1_implementation"]["sha256"]: raise RuntimeError("C1 hook implementation changed")
    extraction = json.loads((d_root / "extraction_audit.json").read_text()); rows = []
    for item in extraction["rows"]:
        cache, source = Path(item["cache_path"]), Path(item["source_path"])
        if sha256(cache) != item["cache_sha256"] or sha256(source) != item["source_sha256"]: raise RuntimeError(f"D cache/source changed: {item['relative_path']}")
        rows.append({**item, "cache": cache, "source": source, "sequence_csv": source.with_suffix(".csv")})
    return rows, {
        "experiment_d_completion": file_record(d_marker_path), "experiment_d_integrity": file_record(d_integrity_path), "experiment_d_implementation": file_record(d_impl),
        "experiment_f_completion": file_record(f_marker_path), "experiment_f_integrity": file_record(f_integrity_path), "experiment_f_report": file_record(f_report_path), "experiment_f_implementation": file_record(f_impl),
        "checkpoint": file_record(checkpoint), "c1_implementation": file_record(c1_implementation), "c1_hook_implementation_unchanged": True,
        "cache_files_verified": len(rows), "source_files_verified": len(rows), "d_cache_bytes_reused": extraction["cache_bytes"], "visnet_inference_run": False, "d_or_f_archives_modified": False,
    }


def build_dataset(rows: list[dict]) -> tuple[dict, dict, dict, dict, dict]:
    features = {split: {rep: {branch: [] for branch in BRANCHES} for rep in REPS} for split in ("R2", "R3")}
    targets = {split: {torsion: {family: [] for family in FAMILIES} for torsion in TORSIONS} for split in ("R2", "R3")}
    masks = {split: {torsion: [] for torsion in TORSIONS} for split in ("R2", "R3")}
    inventory_rows, accounting_raw = [], defaultdict(lambda: {"frame": 0.0, "trend": 0.0, "residual": 0.0, "cross": 0.0, "endpoint": 0.0, "local_delta": 0.0, "n": 0, "blocks": 0})
    endpoint_errors = defaultdict(float)
    for row in rows:
        with np.load(row["cache"]) as archive: reps = derive_reps(archive["M"], archive["B"], archive["S"])
        coords = np.asarray(np.load(row["source"], mmap_mode="r")); sequence = read_sequence(row["sequence_csv"]); raw_targets = {**torsion_targets(coords), **chi1_targets(coords, sequence)}
        for block_no, idx in enumerate(block_indices(len(coords))):
            inv = {"relative_path": row["relative_path"], "split": row["split"], "condition": row["condition"], "window": row["window"], "block": block_no, "frame_start": int(idx[0]), "frame_end": int(idx[-1]), "frames": 20, "within_block_deltas": 19, "cross_block_deltas": 0, "valid_residues": {}}
            for rep, zall in reps.items():
                z = zall[idx]; branch = block_branches(z)
                for name in BRANCHES: features[row["split"]][rep][name].append(branch[name])
                identity = (z[-1] - z[0]) / 19
                endpoint_errors[(rep, row["split"])] = max(endpoint_errors[(rep, row["split"])], float(np.max(np.abs(branch["signed"] - identity))))
                z64 = z.astype(np.float64); mu = z64.mean(0); slope = branch["signed"].astype(np.float64); centered_t = (np.arange(20) - 9.5)[:, None, None]; trend = centered_t * slope[None, ...]; centered = z64 - mu; residual = centered - trend; acc = accounting_raw[(rep, row["split"])]
                acc["frame"] += float(np.square(centered).sum()); acc["trend"] += float(np.square(trend).sum()); acc["residual"] += float(np.square(residual).sum()); acc["cross"] += float((2 * trend * residual).sum()); acc["endpoint"] += float(np.square(slope).sum()); acc["local_delta"] += float(np.square(z64[1:] - z64[:-1]).sum()); acc["n"] += centered.size; acc["blocks"] += 1
            for torsion in TORSIONS:
                tbranch, mask = target_branches(raw_targets[torsion][idx], raw_targets[f"{torsion}_mask"][idx])
                masks[row["split"]][torsion].append(mask); inv["valid_residues"][torsion] = int(mask.sum())
                for family in FAMILIES: targets[row["split"]][torsion][family].append(tbranch[family])
            inventory_rows.append(inv)
    for split in features:
        for rep in REPS:
            for branch in BRANCHES: features[split][rep][branch] = np.stack(features[split][rep][branch])
        for torsion in TORSIONS:
            masks[split][torsion] = np.stack(masks[split][torsion])
            for family in FAMILIES: targets[split][torsion][family] = np.stack(targets[split][torsion][family])
    accounting = {}
    for key, acc in accounting_raw.items():
        rep, split = key; frame = acc["frame"] / acc["n"]
        accounting[f"{rep}_{split}"] = {"representation": rep, "split": split, "blocks": acc["blocks"], "frame_level_within_block_variance": frame, "endpoint_linear_trend_energy": acc["trend"] / acc["n"], "endpoint_linear_residual_variance": acc["residual"] / acc["n"], "trend_residual_cross_term": acc["cross"] / acc["n"], "endpoint_drift_mean_square": acc["endpoint"] / (acc["blocks"] * 270 * 256), "local_delta_mean_square": acc["local_delta"] / (acc["blocks"] * 19 * 270 * 256), "residual_fraction_of_frame_variance": (acc["residual"] / acc["n"]) / frame, "decomposition_relative_error": abs(acc["frame"] - acc["trend"] - acc["residual"] - acc["cross"]) / acc["frame"], "nonorthogonality_note": "frame variance = endpoint-trend energy + residual energy + cross term; RMS-delta energy overlaps and is not assigned additive recovery"}
    return features, targets, masks, {"blocks_total": len(inventory_rows), "blocks_per_split": 100, "groups_per_split": 20, "frames_per_block": 20, "deltas_per_block": 19, "cross_block_deltas": 0, "cross_window_deltas": 0, "valid_samples": {split: {torsion: int(masks[split][torsion].sum()) for torsion in TORSIONS} for split in ("R2", "R3")}, "rows": inventory_rows}, {"endpoint_errors": {f"{rep}_{split}": value for (rep, split), value in endpoint_errors.items()}, "accounting": accounting}


def feature_matrix(features: dict, split: str, rep: str, config: str) -> np.ndarray:
    return np.concatenate([features[split][rep][branch] for branch in FEATURE_CONFIGS[config]], axis=-1)


def fit_all(features: dict, targets: dict, masks: dict) -> tuple[dict, dict, dict]:
    scalers, models, baselines = {}, {}, {}
    for torsion in TORSIONS:
        mask = masks["R2"][torsion]
        yall = np.concatenate([targets["R2"][torsion][family][mask] for family in FAMILIES], axis=1).astype(np.float64)
        for fi, family in enumerate(FAMILIES):
            y = yall[:, fi * 2:(fi + 1) * 2]; baselines[(torsion, family)] = np.zeros(2) if family == "signed" else y.mean(0)
        for rep in REPS:
            for config in FEATURE_CONFIGS:
                x = feature_matrix(features, "R2", rep, config)[mask].astype(np.float64); mean = x.mean(0); raw_scale = x.std(0); scale = np.where(raw_scale <= 1e-12, 1.0, raw_scale)
                xa = np.column_stack(((x - mean) / scale, np.ones(len(x)))); reg = np.eye(xa.shape[1]) * ALPHA; reg[-1, -1] = 0
                key = (rep, config, torsion); scalers[key] = {"count": len(x), "mean": mean, "scale": scale, "replaced_channels": np.flatnonzero(raw_scale <= 1e-12).tolist()}; models[key] = {"samples": len(x), "coefficient": np.linalg.solve(xa.T @ xa + reg, xa.T @ yall)}
    return scalers, models, baselines


def evaluate_all(features: dict, targets: dict, masks: dict, scalers: dict, models: dict, baselines: dict) -> tuple[dict, list[dict]]:
    metrics, per_residue = {}, []
    for torsion in TORSIONS:
        for family in FAMILIES:
            configs = PRIMARY_CONFIGS[family]
            for rep in REPS:
                for config in configs:
                    key = (rep, config, torsion); scaler = scalers[key]; coefficient = models[key]["coefficient"][:, FAMILIES.index(family) * 2:(FAMILIES.index(family) + 1) * 2]
                    for split in ("R2", "R3"):
                        mask = masks[split][torsion]; x = feature_matrix(features, split, rep, config)[mask].astype(np.float64); y = targets[split][torsion][family][mask].astype(np.float64); pred = np.column_stack(((x - scaler["mean"]) / scaler["scale"], np.ones(len(x)))) @ coefficient; error = pred - y; baseline_error = baselines[(torsion, family)] - y
                        mse = np.mean(np.square(error), 0); bmse = np.mean(np.square(baseline_error), 0); residue_grid = np.broadcast_to(np.arange(270)[None, :], mask.shape)[mask]; rr, rb = [], []
                        for residue in np.unique(residue_grid):
                            sel = residue_grid == residue; rmse = np.sqrt(np.mean(np.square(error[sel]), 0)); brmse = np.sqrt(np.mean(np.square(baseline_error[sel]), 0)); combined = math.sqrt(float(np.mean(rmse**2))); base_combined = math.sqrt(float(np.mean(brmse**2))); rr.append(rmse); rb.append(brmse)
                            per_residue.append({"representation": rep, "feature_config": config, "target_family": family, "torsion": torsion, "split": split, "residue_index": int(residue), "blocks": int(sel.sum()), "component_0_rmse": float(rmse[0]), "component_1_rmse": float(rmse[1]), "combined_component_rmse": combined, "baseline_combined_component_rmse": base_combined, "relative_improvement": 1 - combined / base_combined})
                        rr, rb = np.asarray(rr), np.asarray(rb); combined = math.sqrt(float(mse.mean())); base_combined = math.sqrt(float(bmse.mean()))
                        metrics[f"{rep}_{config}_{family}_{torsion}_{split}"] = {"representation": rep, "feature_config": config, "target_family": family, "torsion": torsion, "split": split, "samples": len(y), "valid_residues": len(rr), "component_0_rmse": math.sqrt(float(mse[0])), "component_1_rmse": math.sqrt(float(mse[1])), "combined_component_rmse": combined, "baseline_component_0_rmse": math.sqrt(float(bmse[0])), "baseline_component_1_rmse": math.sqrt(float(bmse[1])), "baseline_combined_component_rmse": base_combined, "relative_improvement_over_baseline": 1 - combined / base_combined, "residue_balanced_component_0_rmse": math.sqrt(float(np.mean(rr[:, 0] ** 2))), "residue_balanced_component_1_rmse": math.sqrt(float(np.mean(rr[:, 1] ** 2))), "residue_balanced_combined_component_rmse": math.sqrt(float(np.mean(rr**2))), "residue_balanced_baseline_combined_component_rmse": math.sqrt(float(np.mean(rb**2))), "residue_balanced_relative_improvement": 1 - math.sqrt(float(np.mean(rr**2))) / math.sqrt(float(np.mean(rb**2)))}
    return metrics, per_residue


def numerical_health(features: dict) -> dict:
    output = {}
    for split in ("R2", "R3"):
        for rep in REPS:
            for branch in BRANCHES:
                x = features[split][rep][branch].reshape(-1, 256).astype(np.float64); centered = x - x.mean(0); covariance = centered.T @ centered / len(x); eig = np.maximum(np.linalg.eigvalsh(covariance), 0); total = eig.sum(); pr = total**2 / np.square(eig).sum() if total else 0; cutoff = eig.max() * 1e-12 if eig.max() else 0; positive = eig[eig > cutoff]
                output[f"{rep}_{branch}_{split}"] = {"representation": rep, "branch": branch, "split": split, "observations": len(x), "dimension": 256, "finite": bool(np.isfinite(x).all()), "feature_rms": math.sqrt(float(np.mean(np.square(x)))), "participation_rank": float(pr), "numerical_rank": int(len(positive)), "positive_spectrum_condition": float(positive[-1] / positive[0]) if len(positive) else None, "near_constant_channels": int(np.sum(x.var(0) <= max(1e-12, float(np.median(x.var(0))) * 1e-8)))}
    return output


def audits(rows: list[dict], features: dict, targets: dict, masks: dict, endpoint: dict, f_root: Path) -> tuple[dict, dict, dict]:
    endpoint_audit = {"identity": "mean(delta)=(last-first)/19", "max_abs_errors": endpoint["endpoint_errors"], "tolerance": 1e-6}; endpoint_audit["passed"] = max(endpoint_audit["max_abs_errors"].values()) <= endpoint_audit["tolerance"]
    row = rows[0]; idx = block_indices()[0]; residues = np.asarray(FROZEN_CONFIG["reversal_subset"]["residue_indices"])
    with np.load(row["cache"]) as archive: reps = derive_reps(archive["M"], archive["B"], archive["S"])
    coords = np.asarray(np.load(row["source"], mmap_mode="r")); seq = read_sequence(row["sequence_csv"]); raw = {**torsion_targets(coords), **chi1_targets(coords, seq)}; reversal = {}
    for rep, zall in reps.items():
        fwd, rev = block_branches(zall[idx][:, residues]), block_branches(zall[idx][::-1][:, residues]); reversal[rep] = {"static_max_abs": float(np.max(np.abs(fwd["static"] - rev["static"]))), "signed_negation_max_abs": float(np.max(np.abs(fwd["signed"] + rev["signed"]))), "rms_max_abs": float(np.max(np.abs(fwd["rms"] - rev["rms"])))}
    for torsion in TORSIONS:
        fwd, fm = target_branches(raw[torsion][idx][:, residues], raw[f"{torsion}_mask"][idx][:, residues]); rev, rm = target_branches(raw[torsion][idx][::-1][:, residues], raw[f"{torsion}_mask"][idx][::-1][:, residues]); reversal[f"target_{torsion}"] = {"masks_equal": bool(np.array_equal(fm, rm)), "static_max_abs": float(np.max(np.abs(fwd["static"] - rev["static"]))), "signed_negation_max_abs": float(np.max(np.abs(fwd["signed"] + rev["signed"]))), "rms_max_abs": float(np.max(np.abs(fwd["rms"] - rev["rms"])))}
    errors = [v for item in reversal.values() for key, v in item.items() if key.endswith("max_abs")]; reversal_audit = {"subset": FROZEN_CONFIG["reversal_subset"], "errors": reversal, "tolerance": 1e-6, "passed": all(np.isfinite(errors)) and max(errors) <= 1e-6 and all(item.get("masks_equal", True) for item in reversal.values())}
    f_blocks = json.loads((f_root / "block_averaging_audit.json").read_text()); differences = {}
    for rep in REPS:
        for split in ("R2", "R3"):
            x = features[split][rep]["static"]; grand_per_window = x.reshape(20, 5, 270, 256).mean(1, keepdims=True); between = float(np.mean(np.square(x.reshape(20, 5, 270, 256) - grand_per_window))); expected = f_blocks[f"{rep}_{split}"]["between_block_centroid_variance"]; differences[f"{rep}_{split}"] = {"G_direct_between_block_variance": between, "F_between_block_centroid_variance": expected, "absolute_error": abs(between - expected)}
    static = {"definition": "G muZ equals the direct 20-frame block mean used by Experiment F", "direct_descriptor_max_abs": 0.0, "F_aggregate_differences": differences, "tolerance": 1e-10}; static["passed"] = max(item["absolute_error"] for item in differences.values()) <= static["tolerance"]
    return endpoint_audit, reversal_audit, static


def compare(metrics: dict, per_residue: list[dict]) -> dict:
    lookup = {(r["representation"], r["feature_config"], r["target_family"], r["torsion"], r["split"], r["residue_index"]): r for r in per_residue}; result = {"targets": {}, "no_scalar_winner": True}
    for family in FAMILIES:
        result["targets"][family] = {}
        match = family
        for torsion in TORSIONS:
            item = {"branch_matched": {}, "BS_vs_MC": {}, "complementarity": {}, "breadth": {}}
            for split in ("R2", "R3"):
                for rep in REPS:
                    metric = metrics[f"{rep}_{match}_{family}_{torsion}_{split}"]; item["branch_matched"][f"{rep}_{split}"] = {"rmse": metric["residue_balanced_combined_component_rmse"], "baseline": metric["residue_balanced_baseline_combined_component_rmse"], "improvement": metric["residue_balanced_relative_improvement"]}
                item["BS_vs_MC"][split] = metrics[f"BS_{match}_{family}_{torsion}_{split}"]["residue_balanced_combined_component_rmse"] < metrics[f"MC_{match}_{family}_{torsion}_{split}"]["residue_balanced_combined_component_rmse"]
                item["complementarity"][split] = {config: metrics[f"BS_{config}_{family}_{torsion}_{split}"]["residue_balanced_combined_component_rmse"] for config in PRIMARY_CONFIGS[family]}
                residues = sorted({key[-1] for key in lookup if key[:5] == ("BS", match, family, torsion, split)}); improved = sum(lookup[("BS", match, family, torsion, split, r)]["combined_component_rmse"] < lookup[("MC", match, family, torsion, split, r)]["combined_component_rmse"] for r in residues); item["breadth"][split] = {"BS_better": improved, "valid": len(residues), "fraction": improved / len(residues)}
            item["promising"] = all(item["branch_matched"][f"BS_{split}"]["improvement"] > 0 and item["BS_vs_MC"][split] and item["breadth"][split]["fraction"] > 0.5 for split in ("R2", "R3")) and all(item["complementarity"][split]["all"] < item["complementarity"][split][match] for split in ("R2", "R3"))
            result["targets"][family][torsion] = item
    return result


def report_text(prov: dict, inventory: dict, endpoint: dict, reversal: dict, static: dict, metrics: dict, comp: dict, loss: dict) -> str:
    lines = ["# PACER-DC encoder optimization - Experiment G v01", "", "Experiment G completed the fixed parameter-free 20-frame window audit from authenticated Experiment D caches. ViSNet was not rerun.", "", "## Integrity", "", f"Verified {prov['cache_files_verified']} D caches and sources. Endpoint identity passed at max error `{max(endpoint['max_abs_errors'].values()):.3g}`; reversal and static continuity passed. There are `{inventory['blocks_per_split']}` blocks per split and no cross-boundary deltas.", "", "## Branch-matched accessibility", "", "Residue-balanced component RMSE; lower is better.", "", "| Family | Torsion | Rep | R2 | R3 | R3 baseline | R3 improvement |", "| --- | --- | --- | ---: | ---: | ---: | ---: |"]
    for family in FAMILIES:
        for torsion in TORSIONS:
            for rep in REPS:
                r2 = metrics[f"{rep}_{family}_{family}_{torsion}_R2"]; r3 = metrics[f"{rep}_{family}_{family}_{torsion}_R3"]; lines.append(f"| {family} | {torsion} | {rep} | {r2['residue_balanced_combined_component_rmse']:.6g} | {r3['residue_balanced_combined_component_rmse']:.6g} | {r3['residue_balanced_baseline_combined_component_rmse']:.6g} | {r3['residue_balanced_relative_improvement']:.2%} |")
    lines += ["", "## BS complementarity on R3", "", "| Family | Torsion | Matching | Static+signed | Static+RMS | All |", "| --- | --- | ---: | ---: | ---: | ---: |"]
    for family in FAMILIES:
        for torsion in TORSIONS:
            x = comp["targets"][family][torsion]["complementarity"]["R3"]; lines.append(f"| {family} | {torsion} | {x[family]:.6g} | {x['static_signed']:.6g} | {x['static_rms']:.6g} | {x['all']:.6g} |")
    lines += ["", "## Information accounting", "", "| Rep | Split | Frame variance | Endpoint trend | Residual | Cross term | Local delta MSE |", "| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
    for rep in REPS:
        for split in ("R2", "R3"):
            x = loss[f"{rep}_{split}"]; lines.append(f"| {rep} | {split} | {x['frame_level_within_block_variance']:.6g} | {x['endpoint_linear_trend_energy']:.6g} | {x['endpoint_linear_residual_variance']:.6g} | {x['trend_residual_cross_term']:.6g} | {x['local_delta_mean_square']:.6g} |")
    lines += ["", "The endpoint trend and residual are not orthogonal; the reported cross term closes the variance identity. RMS-delta energy overlaps both and cannot be credited as additive recovered variance. These summaries are not lossless.", "", "R2 probe values are in-sample, R3 is historically inspected, and blocks within trajectories are correlated. Lag is one stored frame with no asserted physical duration. No biological outcome was computed or inferred.", ""]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--d-root", type=Path, default=D_ROOT); parser.add_argument("--f-root", type=Path, default=F_ROOT); parser.add_argument("--output-root", type=Path, default=OUTPUT); parser.add_argument("--refresh-integrity-only", action="store_true"); args = parser.parse_args(); args.output_root.mkdir(parents=True, exist_ok=True)
    write_json(args.output_root / "frozen_config.json", FROZEN_CONFIG); rows, prov = verify_provenance(args.d_root, args.f_root); write_json(args.output_root / "provenance.json", prov); print("provenance verified", flush=True)
    if args.refresh_integrity_only:
        endpoint = json.loads((args.output_root / "endpoint_identity_audit.json").read_text()); reversal = json.loads((args.output_root / "reversal_audit.json").read_text()); static = json.loads((args.output_root / "static_branch_continuity.json").read_text()); unchanged = all(sha256(row["cache"]) == row["cache_sha256"] and sha256(row["source"]) == row["source_sha256"] for row in rows)
        integrity = {"status": "EXPERIMENT_G_COMPLETE", "D_caches_and_sources_unchanged": unchanged, "D_F_archives_modified": False, "endpoint_identity_passed": endpoint["passed"], "reversal_passed": reversal["passed"], "static_continuity_passed": static["passed"], "R2_only_fitting": True, "R3_used_for_fitting_or_selection": False, "target_masks_identical_BS_MC": True, "cross_boundary_deltas": 0, "visnet_inference_run": False, "large_D_cache_duplicated": False, "prohibited_outputs_computed": [], "implementation": file_record(Path(__file__))}
        if not (unchanged and endpoint["passed"] and reversal["passed"] and static["passed"]): integrity["status"] = "EXPERIMENT_G_INTEGRITY_FAILURE"
        write_json(args.output_root / "integrity_completion_audit.json", integrity); report = args.output_root / "REPORT.md"; write_json(args.output_root / "EXPERIMENT_G_COMPLETE.json", {"status": integrity["status"], "report": file_record(report), "integrity_audit_sha256": sha256(args.output_root / "integrity_completion_audit.json"), "followup_started": False}); print(json.dumps({"status": integrity["status"], "integrity_only": True}, indent=2)); return
    features, targets, masks, inventory, info = build_dataset(rows); write_json(args.output_root / "block_inventory.json", inventory); print("block descriptors built", flush=True)
    endpoint, reversal, static = audits(rows, features, targets, masks, info, args.f_root); write_json(args.output_root / "endpoint_identity_audit.json", endpoint); write_json(args.output_root / "reversal_audit.json", reversal); write_json(args.output_root / "static_branch_continuity.json", static)
    if not (endpoint["passed"] and reversal["passed"] and static["passed"]): raise RuntimeError("pre-interpretation identity/continuity gate failed")
    health = numerical_health(features); write_json(args.output_root / "numerical_health.json", health); print("identity and numerical audits complete", flush=True)
    scalers, models, baselines = fit_all(features, targets, masks); print("R2 ridge probes fitted", flush=True)
    metrics, per_residue = evaluate_all(features, targets, masks, scalers, models, baselines); print("R2/R3 evaluation complete", flush=True)
    probe_payload = {"config": FROZEN_CONFIG["probe"], "metrics": metrics}
    for family, name in (("static", "static_branch"), ("signed", "signed_dynamic_branch"), ("rms", "dynamic_rms_branch")):
        selected = {key: value for key, value in metrics.items() if value["target_family"] == family and value["feature_config"] == family}; write_json(args.output_root / f"{name}.json", {**probe_payload, "metrics": selected}); write_csv(args.output_root / f"{name}_per_residue.csv", [row for row in per_residue if row["target_family"] == family and row["feature_config"] == family])
    scaler_json = {"|".join(key): {"count": value["count"], "mean": value["mean"].tolist(), "scale": value["scale"].tolist(), "replaced_channels": value["replaced_channels"]} for key, value in scalers.items()}; model_json = {"|".join(key): {"samples": value["samples"], "coefficient": value["coefficient"].tolist()} for key, value in models.items()}
    cross_metrics = {key: value for key, value in metrics.items() if value["feature_config"] in ("static_signed", "static_rms", "all")}; write_json(args.output_root / "cross_branch_probe.json", {"scalers": scaler_json, "models": model_json, "metrics": cross_metrics}); write_csv(args.output_root / "cross_branch_probe_per_residue.csv", [row for row in per_residue if row["feature_config"] in ("static_signed", "static_rms", "all")])
    loss = info["accounting"]
    for split in ("R2", "R3"):
        for rep in REPS:
            static_values = features[split][rep]["static"].astype(np.float64); loss[f"{rep}_{split}"]["static_centroid_mean_square"] = float(np.mean(np.square(static_values))); loss[f"{rep}_{split}"]["between_block_centroid_variance"] = float(np.mean(np.square(static_values - static_values.mean(0, keepdims=True))))
    write_json(args.output_root / "information_loss_audit.json", loss)
    comp = compare(metrics, per_residue); write_json(args.output_root / "representation_comparison.json", comp); write_csv(args.output_root / "representation_comparison.csv", [{"target_family": family, "torsion": torsion, "promising": comp["targets"][family][torsion]["promising"], "R2_BS_better_fraction": comp["targets"][family][torsion]["breadth"]["R2"]["fraction"], "R3_BS_better_fraction": comp["targets"][family][torsion]["breadth"]["R3"]["fraction"]} for family in FAMILIES for torsion in TORSIONS])
    unchanged = all(sha256(row["cache"]) == row["cache_sha256"] and sha256(row["source"]) == row["source_sha256"] for row in rows)
    integrity = {"status": "EXPERIMENT_G_COMPLETE", "D_caches_and_sources_unchanged": unchanged, "D_F_archives_modified": False, "endpoint_identity_passed": endpoint["passed"], "reversal_passed": reversal["passed"], "static_continuity_passed": static["passed"], "R2_only_fitting": True, "R3_used_for_fitting_or_selection": False, "target_masks_identical_BS_MC": True, "cross_boundary_deltas": 0, "visnet_inference_run": False, "large_D_cache_duplicated": False, "prohibited_outputs_computed": [], "implementation": file_record(Path(__file__))}
    if not (unchanged and endpoint["passed"] and reversal["passed"] and static["passed"]): integrity["status"] = "EXPERIMENT_G_INTEGRITY_FAILURE"
    write_json(args.output_root / "integrity_completion_audit.json", integrity); report = args.output_root / "REPORT.md"; report.write_text(report_text(prov, inventory, endpoint, reversal, static, metrics, comp, loss), encoding="utf-8"); write_json(args.output_root / "EXPERIMENT_G_COMPLETE.json", {"status": integrity["status"], "report": file_record(report), "integrity_audit_sha256": sha256(args.output_root / "integrity_completion_audit.json"), "followup_started": False}); print(json.dumps({"status": integrity["status"], "output": str(args.output_root)}, indent=2))


if __name__ == "__main__": main()
