#!/usr/bin/env python
"""Run the fixed Experiment E v01 compact B/S fusion audit exactly once."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
from torch import nn

WORKSPACE = Path(__file__).resolve().parents[2]
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))

from project.encoder_atom_readout_v01.experiment_d import (  # noqa: E402
    CHI1_QUADRUPLETS,
    chi1_targets,
    file_record,
    sha256,
    write_csv,
    write_json,
)
from project.encoder_intermediate_v01.phase3 import (  # noqa: E402
    quantiles,
    runtime_report,
    torsion_targets,
    wrapped_angle_error,
)

VERSION = "encoder_bs_fusion_E_v01"
SEEDS = (17, 43, 101)
MODEL_KINDS = ("M", "AVG", "STATIC", "DYNAMIC")
LEARNED_KINDS = ("STATIC", "DYNAMIC")
DEFAULT_D_ROOT = WORKSPACE / "project/results/encoder_atom_readout_D_v01/run_001"
DEFAULT_OUTPUT = WORKSPACE / "project/results/encoder_bs_fusion_E_v01/run_001"
D_IMPLEMENTATION = WORKSPACE / "project/encoder_atom_readout_v01/experiment_d.py"

FROZEN_CONFIG = {
    "version": VERSION,
    "input": "authenticated Experiment D C1 M/B/S caches; no ViSNet inference",
    "references": {"M": "C1 all-heavy mean, 128D", "BS": "concat(B,S), 256D information reference"},
    "candidates": {
        "M": "C1-M, 128D",
        "AVG": "(B+S)/2, 128D",
        "STATIC": "channelwise sigmoid gate constant across samples, 128D",
        "DYNAMIC": "sigmoid(W2(tanh(W1 concat(B,S)+b1)+b2)); hidden=16, 128D",
    },
    "dynamic_hidden_dimension": 16,
    "decoder": "single affine Linear(128,256) with bias; identical family and initialization rule for M/AVG/STATIC/DYNAMIC",
    "target": "raw frozen concat(B,S), no target scaling",
    "loss": "mean squared error over all 256 reconstruction channels",
    "optimizer": "Adam",
    "learning_rate": 0.001,
    "weight_decay": 0.0,
    "batch_size": 8192,
    "epochs": 20,
    "seeds": list(SEEDS),
    "stopping_rule": "exactly 20 epochs; final epoch retained; no early stopping",
    "shuffle": "full R2 permutation each epoch from seed+epoch; identical order across model kinds for a seed",
    "decoder_initialization": "Xavier uniform from seed+1000; identical decoder start across model kinds for a seed; zero bias",
    "dynamic_initialization": "Xavier uniform from seed+2000/+2001; zero biases",
    "static_initialization": "gate logits exactly zero",
    "training_split": "R2 only",
    "evaluation_split": "R3 descriptive only; never used for training/stopping/selection",
    "structural_probe": {
        "method": "exact Experiment D linear ridge protocol",
        "alpha": 1.0,
        "scaling": "per-representation R2-only mean/std; scale<=1e-12 replaced by 1",
        "targets": ["sin(phi)", "cos(phi)", "sin(psi)", "cos(psi)", "sin(chi1)", "cos(chi1)"],
        "chi1_mapping": {aa: list(atoms) for aa, atoms in CHI1_QUADRUPLETS.items()},
    },
    "success_logic": {
        "beats_M": "every seed has lower backbone-combined and chi1 residue-balanced MAE than M on R2 and R3",
        "substantial_BS_recovery": "every seed recovers >=50% of the M-to-BS advantage for both target families on R2 and R3",
        "R3_reconstruction": "every seed has R3 total reconstruction MSE <=1.25 times its R2 MSE",
        "seed_stability": "cross-seed CV <=5% for R3 reconstruction, backbone MAE, and chi1 MAE",
        "static_beats_AVG": "median STATIC is lower than AVG for both structural families on both splits and median R3 reconstruction is lower than AVG median",
        "dynamic_justified": "for every paired seed DYNAMIC is lower than STATIC on R3 reconstruction, backbone MAE, and chi1 MAE",
        "no_single_score": True,
    },
}


class FusionDecoder(nn.Module):
    def __init__(self, kind: str, seed: int) -> None:
        super().__init__()
        if kind not in MODEL_KINDS:
            raise ValueError(kind)
        self.kind = kind
        if kind == "STATIC":
            self.gate_logits = nn.Parameter(torch.zeros(128))
        elif kind == "DYNAMIC":
            self.gate_in = nn.Linear(256, 16)
            self.gate_out = nn.Linear(16, 128)
            gate_generator = torch.Generator(device="cpu").manual_seed(seed + 2000)
            nn.init.xavier_uniform_(self.gate_in.weight, generator=gate_generator)
            nn.init.zeros_(self.gate_in.bias)
            gate_generator = torch.Generator(device="cpu").manual_seed(seed + 2001)
            nn.init.xavier_uniform_(self.gate_out.weight, generator=gate_generator)
            nn.init.zeros_(self.gate_out.bias)
        self.decoder = nn.Linear(128, 256)
        decoder_generator = torch.Generator(device="cpu").manual_seed(seed + 1000)
        nn.init.xavier_uniform_(self.decoder.weight, generator=decoder_generator)
        nn.init.zeros_(self.decoder.bias)

    def fuse(self, m: torch.Tensor, b: torch.Tensor, s: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor | None]:
        if self.kind == "M":
            return m, None
        if self.kind == "AVG":
            return (b + s) * 0.5, None
        if self.kind == "STATIC":
            gate = torch.sigmoid(self.gate_logits).expand_as(b)
            return gate * b + (1 - gate) * s, gate
        gate = torch.sigmoid(self.gate_out(torch.tanh(self.gate_in(torch.cat((b, s), dim=-1)))))
        return gate * b + (1 - gate) * s, gate

    def forward(self, m: torch.Tensor, b: torch.Tensor, s: torch.Tensor):
        z, gate = self.fuse(m, b, s)
        return self.decoder(z), z, gate


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False


def runtime_identity(runtime: dict) -> dict:
    return {
        key: runtime[key]
        for key in ("python", "python_executable", "platform", "pytorch", "pyg", "numpy", "cuda_build", "cudnn", "device", "device_name", "neighbor_backend")
    }


def verify_provenance(d_root: Path, output_root: Path) -> tuple[list[dict], dict]:
    completion_path = d_root / "EXPERIMENT_D_COMPLETE.json"
    integrity_path = d_root / "integrity_completion_audit.json"
    report_path = d_root / "REPORT.md"
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    integrity = json.loads(integrity_path.read_text(encoding="utf-8"))
    if completion["status"] != "EXPERIMENT_D_COMPLETE" or integrity["status"] != "EXPERIMENT_D_COMPLETE":
        raise RuntimeError("Experiment D is not complete")
    if sha256(report_path) != completion["report"]["sha256"] or sha256(integrity_path) != completion["integrity_audit_sha256"]:
        raise RuntimeError("Experiment D completion artifacts changed")
    if sha256(D_IMPLEMENTATION) != integrity["experiment_d_implementation"]["sha256"]:
        raise RuntimeError("Experiment D atom grouping/readout implementation changed")
    d_frozen = json.loads((d_root / "frozen_config.json").read_text(encoding="utf-8"))
    d_groups = json.loads((d_root / "atom_group_definition.json").read_text(encoding="utf-8"))
    if d_frozen["chi1_quadruplets"] != FROZEN_CONFIG["structural_probe"]["chi1_mapping"]:
        raise RuntimeError("Experiment D chi1 mapping differs")
    if d_groups["backbone_atoms"] != ["N", "CA", "C", "O"]:
        raise RuntimeError("Experiment D backbone atom definition differs")
    d_provenance = json.loads((d_root / "provenance.json").read_text(encoding="utf-8"))
    current_runtime = runtime_report("cuda")
    if runtime_identity(current_runtime) != runtime_identity(d_provenance["runtime"]):
        raise RuntimeError("runtime differs from Experiment D")
    extraction = json.loads((d_root / "extraction_audit.json").read_text(encoding="utf-8"))
    rows = []
    for item in extraction["rows"]:
        cache = Path(item["cache_path"])
        if sha256(cache) != item["cache_sha256"]:
            raise RuntimeError(f"Experiment D cache changed: {cache}")
        source = Path(item["source_path"])
        if sha256(source) != item["source_sha256"]:
            raise RuntimeError(f"authenticated source changed: {source}")
        rows.append(
            {
                "relative_path": item["relative_path"],
                "source": source,
                "source_sha256": item["source_sha256"],
                "sequence_csv": source.with_suffix(".csv"),
                "cache": cache,
                "cache_sha256": item["cache_sha256"],
                "split": item["split"],
                "replica": item["replica"],
                "condition": item["condition"],
                "window": item["window"],
            }
        )
    if len(rows) != 40:
        raise RuntimeError("expected 40 Experiment D caches")
    provenance = {
        "experiment_d_root": str(d_root.resolve()),
        "experiment_d_completion": file_record(completion_path),
        "experiment_d_integrity": file_record(integrity_path),
        "experiment_d_report": file_record(report_path),
        "experiment_d_frozen_config": file_record(d_root / "frozen_config.json"),
        "experiment_d_atom_groups": file_record(d_root / "atom_group_definition.json"),
        "experiment_d_implementation": file_record(D_IMPLEMENTATION),
        "runtime": current_runtime,
        "cache_files_verified": len(rows),
        "source_files_verified": len(rows),
        "visnet_inference_run": False,
        "large_cache_duplicated": False,
    }
    return rows, provenance


def load_data(rows: list[dict]) -> tuple[dict[str, dict[str, torch.Tensor]], dict[str, list[dict]], dict[str, dict[str, np.ndarray]]]:
    data = {}
    file_layout = {"R2": [], "R3": []}
    targets = {}
    for split in ("R2", "R3"):
        selected = [row for row in rows if row["split"] == split]
        m_parts, b_parts, s_parts = [], [], []
        phi_parts, psi_parts, chi_parts, residue_parts = [], [], [], []
        offset = 0
        for row in selected:
            with np.load(row["cache"]) as archive:
                m = archive["M"].reshape(-1, 128)
                b = archive["B"].reshape(-1, 128)
                s = archive["S"].reshape(-1, 128)
                sequence = str(archive["sequence"])
                frames = len(archive["frame_ids"])
            coords = np.asarray(np.load(row["source"], mmap_mode="r"))
            backbone = torsion_targets(coords)
            chi = chi1_targets(coords, sequence)
            count = len(m)
            m_parts.append(m)
            b_parts.append(b)
            s_parts.append(s)
            phi_parts.append(backbone["phi"].reshape(-1))
            psi_parts.append(backbone["psi"].reshape(-1))
            chi_parts.append(chi["chi1"].reshape(-1))
            residue_parts.append(np.tile(np.arange(270, dtype=np.int64), frames))
            file_layout[split].append({**row, "start": offset, "stop": offset + count, "frames": frames})
            offset += count
        data[split] = {
            "M": torch.from_numpy(np.concatenate(m_parts).astype(np.float32, copy=False)),
            "B": torch.from_numpy(np.concatenate(b_parts).astype(np.float32, copy=False)),
            "S": torch.from_numpy(np.concatenate(s_parts).astype(np.float32, copy=False)),
        }
        targets[split] = {
            "phi": np.concatenate(phi_parts),
            "psi": np.concatenate(psi_parts),
            "chi1": np.concatenate(chi_parts),
            "residue_index": np.concatenate(residue_parts),
        }
    return data, file_layout, targets


def batches(total: int, batch_size: int, order: torch.Tensor | None = None):
    if order is None:
        order = torch.arange(total)
    for start in range(0, total, batch_size):
        yield order[start : start + batch_size]


def train_one(kind: str, seed: int, train: dict[str, torch.Tensor], device: str, output_root: Path) -> tuple[dict, Path]:
    set_seed(seed)
    model = FusionDecoder(kind, seed).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=FROZEN_CONFIG["learning_rate"], weight_decay=0.0)
    history = []
    total = len(train["M"])
    for epoch in range(FROZEN_CONFIG["epochs"]):
        generator = torch.Generator(device="cpu").manual_seed(seed + epoch)
        order = torch.randperm(total, generator=generator)
        sum_squared = 0.0
        elements = 0
        model.train()
        for index in batches(total, FROZEN_CONFIG["batch_size"], order):
            m = train["M"][index].to(device)
            b = train["B"][index].to(device)
            s = train["S"][index].to(device)
            target = torch.cat((b, s), dim=-1)
            prediction, _, _ = model(m, b, s)
            loss = torch.mean(torch.square(prediction - target))
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            sum_squared += float(loss.detach()) * target.numel()
            elements += target.numel()
        history.append({"epoch": epoch + 1, "mse": sum_squared / elements})
    checkpoint = output_root / "checkpoints" / f"{kind.lower()}_seed{seed}.pt"
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"kind": kind, "seed": seed, "state_dict": model.state_dict()}, checkpoint)
    return {"kind": kind, "seed": seed, "epochs": history}, checkpoint


def load_model(checkpoint: Path, device: str) -> FusionDecoder:
    payload = torch.load(checkpoint, map_location="cpu", weights_only=True)
    model = FusionDecoder(payload["kind"], int(payload["seed"]))
    model.load_state_dict(payload["state_dict"], strict=True)
    return model.to(device).eval()


def reconstruction_metrics(model: FusionDecoder, data: dict[str, torch.Tensor], device: str) -> dict:
    sse_b = 0.0
    sse_s = 0.0
    count = 0
    with torch.inference_mode():
        for index in batches(len(data["M"]), FROZEN_CONFIG["batch_size"]):
            m = data["M"][index].to(device)
            b = data["B"][index].to(device)
            s = data["S"][index].to(device)
            prediction, _, _ = model(m, b, s)
            b_hat, s_hat = prediction.chunk(2, dim=-1)
            sse_b += float(torch.square(b_hat - b).sum())
            sse_s += float(torch.square(s_hat - s).sum())
            count += b.numel()
    b_mse = sse_b / count
    s_mse = sse_s / count
    return {
        "samples": len(data["M"]),
        "B_mse": b_mse,
        "S_mse": s_mse,
        "total_mse": (b_mse + s_mse) / 2,
        "B_rmse": math.sqrt(b_mse),
        "S_rmse": math.sqrt(s_mse),
        "total_rmse": math.sqrt((b_mse + s_mse) / 2),
    }


def derive_representation(
    kind: str,
    data: dict[str, torch.Tensor],
    device: str,
    checkpoint: Path | None = None,
    return_gates: bool = False,
) -> tuple[np.ndarray, np.ndarray | None]:
    total = len(data["M"])
    output = np.empty((total, 128), dtype=np.float32)
    gates = np.empty((total, 128), dtype=np.float32) if return_gates else None
    model = load_model(checkpoint, device) if checkpoint else None
    with torch.inference_mode():
        for index in batches(total, FROZEN_CONFIG["batch_size"]):
            m = data["M"][index].to(device)
            b = data["B"][index].to(device)
            s = data["S"][index].to(device)
            if kind == "M":
                z, gate = m, None
            elif kind == "AVG":
                z, gate = (b + s) * 0.5, None
            else:
                z, gate = model.fuse(m, b, s)
            output[index.numpy()] = z.cpu().numpy()
            if gates is not None and gate is not None:
                gates[index.numpy()] = gate.cpu().numpy()
    return output, gates


def health(array: np.ndarray, name: str, split: str) -> tuple[dict, list[dict], list[dict]]:
    dimension = array.shape[1]
    sums = np.zeros(dimension, dtype=np.float64)
    sums2 = np.zeros(dimension, dtype=np.float64)
    cross = np.zeros((dimension, dimension), dtype=np.float64)
    norms = []
    finite = 0
    for start in range(0, len(array), 50000):
        chunk32 = array[start : start + 50000]
        finite += int(np.isfinite(chunk32).sum())
        chunk = chunk32.astype(np.float64)
        sums += chunk.sum(axis=0)
        sums2 += np.square(chunk).sum(axis=0)
        cross += chunk.T @ chunk
        norms.append(np.linalg.norm(chunk32, axis=1))
    mean = sums / len(array)
    variance = np.maximum(sums2 / len(array) - mean**2, 0)
    covariance = cross / len(array) - np.outer(mean, mean)
    covariance = (covariance + covariance.T) * 0.5
    eigenvalues = np.maximum(np.linalg.eigvalsh(covariance)[::-1], 0)
    eigen_sum = eigenvalues.sum()
    weights = eigenvalues / eigen_sum
    participation = eigen_sum**2 / np.square(eigenvalues).sum()
    threshold = eigenvalues[0] * 1e-12
    positive = eigenvalues[eigenvalues > threshold]
    near_threshold = max(1e-12, float(np.median(variance)) * 1e-8)
    result = {
        "representation": name, "split": split, "dimension": dimension,
        "observations": len(array), "elements": array.size,
        "finite_elements": finite, "nonfinite_elements": array.size - finite,
        "feature_rms": math.sqrt(float(sums2.sum() / array.size)),
        "feature_norm_distribution": quantiles(np.concatenate(norms)),
        "channel_variance_distribution": quantiles(variance),
        "near_constant_channels": int((variance <= near_threshold).sum()),
        "participation_rank": float(participation),
        "entropy_rank": float(np.exp(-np.sum(weights[weights > 0] * np.log(weights[weights > 0])))),
        "numerical_rank": int(np.linalg.matrix_rank(covariance)),
        "positive_spectrum_condition": float(eigenvalues[0] / positive[-1]),
        "regularized_condition_1e-6": float((eigenvalues[0] + eigenvalues[0] * 1e-6) / (eigenvalues[-1] + eigenvalues[0] * 1e-6)),
    }
    channels = [
        {"representation": name, "split": split, "channel": i, "mean": float(mean[i]), "variance": float(variance[i]), "std": math.sqrt(float(variance[i]))}
        for i in range(dimension)
    ]
    cumulative = np.cumsum(weights)
    spectrum = [
        {"representation": name, "split": split, "eigen_index": i, "eigenvalue": float(eigenvalues[i]), "explained_fraction": float(weights[i]), "cumulative_fraction": float(cumulative[i])}
        for i in range(dimension)
    ]
    return result, channels, spectrum


def temporal(array: np.ndarray, layout: list[dict], rms: float, name: str, split: str) -> tuple[dict, list[dict]]:
    frame_values, within_values, detail = [], [], []
    window_centroids = {}
    for item in layout:
        features = array[item["start"] : item["stop"]].reshape(item["frames"], 270, -1)
        centroid = features.mean(axis=0, dtype=np.float64)
        adjacent = np.sqrt(np.mean(np.square(np.diff(features.astype(np.float64), axis=0)), axis=2))
        frame = float(adjacent.mean())
        within = float(np.sqrt(np.mean(np.square(features.astype(np.float64) - centroid[None]))))
        trajectory = (item["replica"], item["condition"])
        window_centroids[(trajectory, item["window"])] = centroid
        frame_values.append(frame)
        within_values.append(within)
        detail.append({"representation": name, "split": split, "condition": item["condition"], "window": item["window"], "frame_displacement": frame, "within_window_rms": within})
    trajectories = sorted({key for key, _window in window_centroids})
    drifts = []
    for trajectory in trajectories:
        windows = sorted(window for key, window in window_centroids if key == trajectory)
        delta = window_centroids[(trajectory, windows[-1])] - window_centroids[(trajectory, windows[0])]
        drifts.append(float(np.sqrt(np.mean(np.square(delta)))))
    result = {
        "representation": name, "split": split,
        "frame_displacement": float(np.mean(frame_values)),
        "within_window_rms": float(np.mean(within_values)),
        "chronological_drift": float(np.mean(drifts)),
        "scale_relative_frame_displacement": float(np.mean(frame_values) / rms),
        "scale_relative_within_window": float(np.mean(within_values) / rms),
        "scale_relative_chronological_drift": float(np.mean(drifts) / rms),
    }
    return result, detail


def fit_probe(array: np.ndarray, targets: dict[str, np.ndarray]) -> tuple[dict, dict]:
    union = np.isfinite(targets["phi"]) | np.isfinite(targets["psi"])
    values = array[union].astype(np.float64)
    mean = values.mean(axis=0)
    raw_scale = values.std(axis=0)
    replaced = raw_scale <= 1e-12
    scale = np.where(replaced, 1.0, raw_scale)
    standardized = (array.astype(np.float64) - mean) / scale
    models = {}
    for torsion in ("phi", "psi", "chi1"):
        mask = np.isfinite(targets[torsion])
        x = standardized[mask]
        angle = targets[torsion][mask]
        y = np.column_stack((np.sin(angle), np.cos(angle)))
        xa = np.column_stack((x, np.ones(len(x))))
        regularizer = np.eye(xa.shape[1])
        regularizer[-1, -1] = 0
        coefficient = np.linalg.solve(xa.T @ xa + regularizer, xa.T @ y)
        models[torsion] = {"coefficient": coefficient, "samples": len(x)}
    scaler = {"mean": mean, "scale": scale, "replaced_channels": np.flatnonzero(replaced).tolist(), "count": len(values)}
    return scaler, models


def evaluate_probe(
    array: np.ndarray,
    targets: dict[str, np.ndarray],
    scaler: dict,
    models: dict,
    name: str,
    split: str,
) -> tuple[dict, list[dict]]:
    standardized = (array.astype(np.float64) - scaler["mean"]) / scaler["scale"]
    residue_index = targets["residue_index"]
    result = {"representation": name, "split": split, "torsions": {}}
    per_residue = []
    for torsion in ("phi", "psi", "chi1"):
        mask = np.isfinite(targets[torsion])
        x = standardized[mask]
        angle = targets[torsion][mask]
        y = np.column_stack((np.sin(angle), np.cos(angle)))
        xa = np.column_stack((x, np.ones(len(x))))
        prediction = xa @ models[torsion]["coefficient"]
        component = prediction - y
        angular = wrapped_angle_error(np.arctan2(prediction[:, 0], prediction[:, 1]), angle)
        residues = residue_index[mask]
        residue_component, residue_angle_mae, residue_angle_mse = [], [], []
        for residue in np.unique(residues):
            select = residues == residue
            comp_mse = np.square(component[select]).mean(axis=0)
            angle_mae = float(angular[select].mean())
            angle_rmse = float(np.sqrt(np.square(angular[select]).mean()))
            residue_component.append(comp_mse)
            residue_angle_mae.append(angle_mae)
            residue_angle_mse.append(angle_rmse**2)
            per_residue.append({"representation": name, "split": split, "torsion": torsion, "residue_index": int(residue), "samples": int(select.sum()), "sin_rmse": math.sqrt(float(comp_mse[0])), "cos_rmse": math.sqrt(float(comp_mse[1])), "angular_mae_deg": angle_mae, "angular_rmse_deg": angle_rmse})
        residue_component = np.asarray(residue_component)
        result["torsions"][torsion] = {
            "samples": len(angle), "valid_residue_positions": len(residue_angle_mae),
            "sin_rmse": math.sqrt(float(np.square(component[:, 0]).mean())),
            "cos_rmse": math.sqrt(float(np.square(component[:, 1]).mean())),
            "residue_balanced_sin_rmse": math.sqrt(float(residue_component[:, 0].mean())),
            "residue_balanced_cos_rmse": math.sqrt(float(residue_component[:, 1].mean())),
            "angular_mae_deg": float(angular.mean()),
            "residue_balanced_angular_mae_deg": float(np.mean(residue_angle_mae)),
            "residue_balanced_angular_rmse_deg": math.sqrt(float(np.mean(residue_angle_mse))),
        }
    phi, psi = result["torsions"]["phi"], result["torsions"]["psi"]
    result["backbone_combined_mae_deg"] = (phi["residue_balanced_angular_mae_deg"] + psi["residue_balanced_angular_mae_deg"]) / 2
    result["backbone_mean_component_rmse"] = float(np.mean([phi["residue_balanced_sin_rmse"], phi["residue_balanced_cos_rmse"], psi["residue_balanced_sin_rmse"], psi["residue_balanced_cos_rmse"]]))
    result["chi1_mae_deg"] = result["torsions"]["chi1"]["residue_balanced_angular_mae_deg"]
    return result, per_residue


def gate_diagnostics(
    kind: str,
    seed: int,
    gate_by_split: dict[str, np.ndarray] | None,
    checkpoint: Path,
    layout: dict[str, list[dict]],
) -> tuple[dict, list[dict]]:
    model = load_model(checkpoint, "cpu")
    rows = []
    if kind == "STATIC":
        gate = torch.sigmoid(model.gate_logits).detach().numpy()
        entropy = -(gate * np.log(gate) + (1 - gate) * np.log(1 - gate))
        result = {
            "kind": kind, "seed": seed,
            "gate_distribution": quantiles(gate),
            "gate_entropy_mean": float(entropy.mean()),
            "variation_across_frames": 0.0,
            "variation_across_residues": 0.0,
            "R2_R3_gate_mean_rmse": 0.0,
        }
        for channel, value in enumerate(gate):
            rows.append({"kind": kind, "seed": seed, "split": "constant", "channel": channel, "mean": float(value), "std": 0.0})
        return result, rows
    split_stats = {}
    for split in ("R2", "R3"):
        gates = gate_by_split[split]
        channel_mean = gates.mean(axis=0, dtype=np.float64)
        channel_std = gates.std(axis=0, dtype=np.float64)
        entropy = -(gates.astype(np.float64) * np.log(np.clip(gates, 1e-12, 1)) + (1 - gates.astype(np.float64)) * np.log(np.clip(1 - gates, 1e-12, 1)))
        frame_variation, residue_variation = [], []
        for item in layout[split]:
            window = gates[item["start"] : item["stop"]].reshape(item["frames"], 270, 128)
            frame_variation.append(float(window.std(axis=0, dtype=np.float64).mean()))
            residue_variation.append(float(window.std(axis=1, dtype=np.float64).mean()))
        split_stats[split] = {
            "gate_mean_distribution": quantiles(channel_mean),
            "gate_std_distribution": quantiles(channel_std),
            "global_gate_mean": float(gates.mean()),
            "global_gate_std": float(gates.std()),
            "gate_entropy_mean": float(entropy.mean()),
            "variation_across_frames": float(np.mean(frame_variation)),
            "variation_across_residues": float(np.mean(residue_variation)),
            "channel_mean": channel_mean,
        }
        for channel in range(128):
            rows.append({"kind": kind, "seed": seed, "split": split, "channel": channel, "mean": float(channel_mean[channel]), "std": float(channel_std[channel])})
    delta = split_stats["R3"]["channel_mean"] - split_stats["R2"]["channel_mean"]
    result = {
        "kind": kind, "seed": seed,
        "R2": {key: value for key, value in split_stats["R2"].items() if key != "channel_mean"},
        "R3": {key: value for key, value in split_stats["R3"].items() if key != "channel_mean"},
        "R2_R3_channel_mean_rmse": float(np.sqrt(np.mean(np.square(delta)))),
        "R2_R3_channel_mean_max_abs": float(np.abs(delta).max()),
        "causal_interpretation": False,
    }
    return result, rows


def cv(values: list[float]) -> float:
    mean = float(np.mean(values))
    return float(np.std(values, ddof=1) / abs(mean)) if mean else 0.0


def compare(
    reconstruction: list[dict],
    probes: dict[str, dict],
    output_root: Path,
) -> dict:
    m = {split: probes[f"M_{split}"] for split in ("R2", "R3")}
    bs = {split: probes[f"BS_{split}"] for split in ("R2", "R3")}
    assessment = {}
    for kind in ("AVG", "STATIC", "DYNAMIC"):
        names = ["AVG"] if kind == "AVG" else [f"{kind}_seed{seed}" for seed in SEEDS]
        entries = []
        for name in names:
            item = {"representation": name, "splits": {}}
            for split in ("R2", "R3"):
                probe = probes[f"{name}_{split}"]
                item["splits"][split] = {
                    "backbone_mae": probe["backbone_combined_mae_deg"],
                    "chi1_mae": probe["chi1_mae_deg"],
                    "backbone_beats_M": probe["backbone_combined_mae_deg"] < m[split]["backbone_combined_mae_deg"],
                    "chi1_beats_M": probe["chi1_mae_deg"] < m[split]["chi1_mae_deg"],
                    "backbone_BS_advantage_recovered": (m[split]["backbone_combined_mae_deg"] - probe["backbone_combined_mae_deg"]) / (m[split]["backbone_combined_mae_deg"] - bs[split]["backbone_combined_mae_deg"]),
                    "chi1_BS_advantage_recovered": (m[split]["chi1_mae_deg"] - probe["chi1_mae_deg"]) / (m[split]["chi1_mae_deg"] - bs[split]["chi1_mae_deg"]),
                }
            entries.append(item)
        assessment[kind] = entries
    recon_lookup = {(row["kind"], row["seed"], row["split"]): row for row in reconstruction if row["kind"] != "BS"}
    family = {}
    for kind in ("STATIC", "DYNAMIC"):
        stable = {
            "R3_reconstruction_cv": cv([recon_lookup[(kind, seed, "R3")]["total_mse"] for seed in SEEDS]),
            "R3_backbone_cv": cv([probes[f"{kind}_seed{seed}_R3"]["backbone_combined_mae_deg"] for seed in SEEDS]),
            "R3_chi1_cv": cv([probes[f"{kind}_seed{seed}_R3"]["chi1_mae_deg"] for seed in SEEDS]),
        }
        all_entries = assessment[kind]
        criteria = {
            "all_seeds_beat_M_both_targets_both_splits": all(
                split_data[metric]
                for entry in all_entries for split_data in entry["splits"].values()
                for metric in ("backbone_beats_M", "chi1_beats_M")
            ),
            "all_seeds_recover_at_least_half_BS_advantage": all(
                split_data[metric] >= 0.5
                for entry in all_entries for split_data in entry["splits"].values()
                for metric in ("backbone_BS_advantage_recovered", "chi1_BS_advantage_recovered")
            ),
            "R3_reconstruction_retained_all_seeds": all(
                recon_lookup[(kind, seed, "R3")]["total_mse"] <= 1.25 * recon_lookup[(kind, seed, "R2")]["total_mse"]
                for seed in SEEDS
            ),
            "seed_CV_within_5_percent": all(value <= 0.05 for value in stable.values()),
        }
        family[kind] = {"seed_stability": stable, "criteria": criteria, "promising": all(criteria.values())}
    avg_probe = {split: probes[f"AVG_{split}"] for split in ("R2", "R3")}
    static_medians = {
        split: {
            "backbone": float(np.median([probes[f"STATIC_seed{seed}_{split}"]["backbone_combined_mae_deg"] for seed in SEEDS])),
            "chi1": float(np.median([probes[f"STATIC_seed{seed}_{split}"]["chi1_mae_deg"] for seed in SEEDS])),
        }
        for split in ("R2", "R3")
    }
    avg_recon_median = float(np.median([recon_lookup[("AVG", seed, "R3")]["total_mse"] for seed in SEEDS]))
    static_recon_median = float(np.median([recon_lookup[("STATIC", seed, "R3")]["total_mse"] for seed in SEEDS]))
    static_beats_avg = (
        all(static_medians[split][metric] < avg_probe[split]["backbone_combined_mae_deg" if metric == "backbone" else "chi1_mae_deg"] for split in ("R2", "R3") for metric in ("backbone", "chi1"))
        and static_recon_median < avg_recon_median
    )
    paired_dynamic = {
        str(seed): {
            "R3_reconstruction_lower": recon_lookup[("DYNAMIC", seed, "R3")]["total_mse"] < recon_lookup[("STATIC", seed, "R3")]["total_mse"],
            "R3_backbone_lower": probes[f"DYNAMIC_seed{seed}_R3"]["backbone_combined_mae_deg"] < probes[f"STATIC_seed{seed}_R3"]["backbone_combined_mae_deg"],
            "R3_chi1_lower": probes[f"DYNAMIC_seed{seed}_R3"]["chi1_mae_deg"] < probes[f"STATIC_seed{seed}_R3"]["chi1_mae_deg"],
        }
        for seed in SEEDS
    }
    dynamic_justified = all(all(values.values()) for values in paired_dynamic.values())
    bs_reference_retained = all(
        bs[split][metric] <= probes[f"{name}_{split}"][metric]
        for split in ("R2", "R3")
        for metric in ("backbone_combined_mae_deg", "chi1_mae_deg")
        for name in ["M", "AVG", *[f"STATIC_seed{s}" for s in SEEDS], *[f"DYNAMIC_seed{s}" for s in SEEDS]]
    )
    result = {
        "candidate_assessments": assessment,
        "family_assessments": family,
        "static_beats_AVG": static_beats_avg,
        "static_medians": static_medians,
        "paired_dynamic_vs_static": paired_dynamic,
        "dynamic_gating_justified": dynamic_justified,
        "BS_remains_information_rich_reference": bs_reference_retained,
        "stop_learned_fusion_development": not family["STATIC"]["promising"] and not family["DYNAMIC"]["promising"],
        "no_seed_selected": True,
        "no_biological_outputs_used": True,
    }
    write_json(output_root / "representation_comparison.json", result)
    return result


def render_report(
    training: list[dict],
    reconstruction: list[dict],
    numerical: dict,
    temporal_metrics: dict,
    probes: dict,
    gate_metrics: dict,
    comparison: dict,
    integrity: dict,
) -> str:
    recon_rows = [row for row in reconstruction if row["kind"] != "BS"]
    probe_rows = []
    for name in ["M", "AVG", *[f"STATIC_seed{s}" for s in SEEDS], *[f"DYNAMIC_seed{s}" for s in SEEDS], "BS"]:
        for split in ("R2", "R3"):
            item = probes[f"{name}_{split}"]
            probe_rows.append({"name": name, "split": split, "phi": item["torsions"]["phi"]["residue_balanced_angular_mae_deg"], "psi": item["torsions"]["psi"]["residue_balanced_angular_mae_deg"], "backbone": item["backbone_combined_mae_deg"], "chi1": item["chi1_mae_deg"]})
    compact_names = ["M", "AVG", *[f"STATIC_seed{s}" for s in SEEDS], *[f"DYNAMIC_seed{s}" for s in SEEDS]]
    health_rows = []
    for name in compact_names:
        for split in ("R2", "R3"):
            n = numerical[f"{name}_{split}"]
            t = temporal_metrics[f"{name}_{split}"]
            health_rows.append({"name": name, "split": split, "rms": n["feature_rms"], "rank": n["participation_rank"], "nrank": n["numerical_rank"], "condition": n["positive_spectrum_condition"], "frame": t["scale_relative_frame_displacement"], "within": t["scale_relative_within_window"], "drift": t["scale_relative_chronological_drift"]})

    def table(rows, columns):
        output = ["| " + " | ".join(label for _key, label in columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
        for row in rows:
            values = []
            for key, _label in columns:
                value = row[key]
                values.append(f"{value:.5g}" if isinstance(value, float) else str(value))
            output.append("| " + " | ".join(values) + " |")
        return output

    lines = [
        "# PACER-DC encoder optimization — Experiment E v01",
        "",
        "Experiment E completed the fixed compact B/S fusion comparison once using Experiment D caches only. ViSNet was not rerun. No structural target or biological outcome was used to train fusion.",
        "",
        "## Training and reconstruction",
        "",
        f"All M/AVG/STATIC/DYNAMIC models used the same linear decoder, Adam `1e-3`, batch size `8192`, 20 fixed epochs, and seeds `{list(SEEDS)}`. No early stopping or seed selection was used.",
        "",
        *table(recon_rows, [("kind", "Kind"), ("seed", "Seed"), ("split", "Split"), ("total_mse", "Total MSE"), ("B_mse", "B MSE"), ("S_mse", "S MSE")]),
        "",
        "BS identity reconstruction is exactly zero and is retained only as the information reference.",
        "",
        "## Structural accessibility",
        "",
        *table(probe_rows, [("name", "Representation"), ("split", "Split"), ("phi", "Phi MAE"), ("psi", "Psi MAE"), ("backbone", "Combined MAE"), ("chi1", "Chi1 MAE")]),
        "",
        "## Numerical and temporal health",
        "",
        *table(health_rows, [("name", "Representation"), ("split", "Split"), ("rms", "RMS"), ("rank", "Participation rank"), ("nrank", "Numerical rank"), ("condition", "Condition"), ("frame", "Frame/RMS"), ("within", "Within/RMS"), ("drift", "Drift/RMS")]),
        "",
        "## Interpretation",
        "",
        f"BS remains the information-rich reference: `{comparison['BS_remains_information_rich_reference']}`.",
        f"STATIC promising: `{comparison['family_assessments']['STATIC']['promising']}`. DYNAMIC promising: `{comparison['family_assessments']['DYNAMIC']['promising']}`.",
        f"STATIC beats AVG under the frozen rule: `{comparison['static_beats_AVG']}`. Dynamic gating reproducibly justified over STATIC: `{comparison['dynamic_gating_justified']}`.",
        f"Stop learned fusion development: `{comparison['stop_learned_fusion_development']}`.",
        "",
        "Gate statistics are descriptive engineering diagnostics and are not causal atom-importance estimates.",
        "",
        "## Limitations",
        "",
        "Fusion reconstruction and structural probes share the same frozen representations but structural targets never enter fusion training. R2 probe values are in-sample ridge fits; R3 is historically inspected and descriptive. Frames and windows are correlated. The three seeds characterize initialization stability but are not biological replicates. No biological efficacy inference follows.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--d-root", type=Path, default=DEFAULT_D_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--resume-analysis", action="store_true", help="Reuse the complete frozen checkpoints/history after an artifact-writing interruption.")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    frozen = dict(FROZEN_CONFIG)
    frozen["experiment_d_root"] = str(args.d_root.resolve())
    write_json(args.output_root / "frozen_config.json", frozen)

    rows, provenance = verify_provenance(args.d_root, args.output_root)
    write_json(args.output_root / "provenance.json", provenance)
    data, layout, targets = load_data(rows)

    checkpoints = {(kind, seed): args.output_root / "checkpoints" / f"{kind.lower()}_seed{seed}.pt" for kind in MODEL_KINDS for seed in SEEDS}
    if args.resume_analysis:
        missing = [str(path) for path in checkpoints.values() if not path.is_file()]
        if missing:
            raise RuntimeError(f"Cannot resume: missing frozen checkpoints: {missing}")
        training_history = json.loads((args.output_root / "training_history.json").read_text(encoding="utf-8"))["runs"]
        reconstruction = json.loads((args.output_root / "reconstruction_metrics.json").read_text(encoding="utf-8"))
    else:
        training_history = []
        reconstruction = []
        started = time.time()
        for kind in MODEL_KINDS:
            for seed in SEEDS:
                history, checkpoint = train_one(kind, seed, data["R2"], args.device, args.output_root)
                checkpoints[(kind, seed)] = checkpoint
                training_history.append(history)
                model = load_model(checkpoint, args.device)
                for split in ("R2", "R3"):
                    metrics = reconstruction_metrics(model, data[split], args.device)
                    reconstruction.append({"kind": kind, "seed": seed, "split": split, **metrics})
                print(f"trained {kind} seed {seed}", flush=True)
        reconstruction.append({"kind": "BS", "seed": None, "split": "R2", "samples": len(data["R2"]["M"]), "B_mse": 0.0, "S_mse": 0.0, "total_mse": 0.0, "B_rmse": 0.0, "S_rmse": 0.0, "total_rmse": 0.0})
        reconstruction.append({"kind": "BS", "seed": None, "split": "R3", "samples": len(data["R3"]["M"]), "B_mse": 0.0, "S_mse": 0.0, "total_mse": 0.0, "B_rmse": 0.0, "S_rmse": 0.0, "total_rmse": 0.0})
        write_json(args.output_root / "training_history.json", {"config": FROZEN_CONFIG, "elapsed_seconds": round(time.time() - started, 3), "runs": training_history})
        write_json(args.output_root / "reconstruction_metrics.json", reconstruction)
        write_csv(args.output_root / "reconstruction_metrics.csv", reconstruction)

    d_numerical = json.loads((args.d_root / "numerical_health.json").read_text(encoding="utf-8"))
    d_temporal = json.loads((args.d_root / "temporal_diagnostics.json").read_text(encoding="utf-8"))
    d_backbone = json.loads((args.d_root / "backbone_torsion_probe.json").read_text(encoding="utf-8"))["metrics"]
    d_chi = json.loads((args.d_root / "chi1_probe.json").read_text(encoding="utf-8"))["metrics"]
    d_per_residue = []
    for filename in ("backbone_torsion_probe_per_residue.csv", "chi1_probe_per_residue.csv"):
        with (args.d_root / filename).open(encoding="utf-8", newline="") as handle:
            d_per_residue.extend(list(csv.DictReader(handle)))

    numerical = {}
    temporal_metrics = {}
    channel_rows, spectrum_rows, temporal_rows = [], [], []
    probes = {}
    per_residue_rows = []
    scaler_records = {}
    model_records = {}
    gate_reports = {}
    gate_channel_rows = []
    representations = [("M", "M", None), ("AVG", "AVG", None)]
    representations += [(f"STATIC_seed{seed}", "STATIC", checkpoints[("STATIC", seed)]) for seed in SEEDS]
    representations += [(f"DYNAMIC_seed{seed}", "DYNAMIC", checkpoints[("DYNAMIC", seed)]) for seed in SEEDS]

    for name, kind, checkpoint in representations:
        arrays = {}
        gates = {}
        for split in ("R2", "R3"):
            arrays[split], gate = derive_representation(kind, data[split], args.device, checkpoint, return_gates=kind == "DYNAMIC")
            if gate is not None:
                gates[split] = gate
            h, channels, spectrum = health(arrays[split], name, split)
            numerical[f"{name}_{split}"] = h
            channel_rows.extend(channels)
            spectrum_rows.extend(spectrum)
            t, detail = temporal(arrays[split], layout[split], h["feature_rms"], name, split)
            temporal_metrics[f"{name}_{split}"] = t
            temporal_rows.extend(detail)
        scaler, ridge_models = fit_probe(arrays["R2"], targets["R2"])
        scaler_records[name] = {"count": scaler["count"], "mean": scaler["mean"].tolist(), "scale": scaler["scale"].tolist(), "replaced_channels": scaler["replaced_channels"]}
        model_records[name] = {torsion: {"samples": item["samples"], "coefficient": item["coefficient"].tolist()} for torsion, item in ridge_models.items()}
        for split in ("R2", "R3"):
            metric, residue = evaluate_probe(arrays[split], targets[split], scaler, ridge_models, name, split)
            probes[f"{name}_{split}"] = metric
            per_residue_rows.extend(residue)
        if kind in LEARNED_KINDS:
            seed = int(name.split("seed")[1])
            diagnostic, diagnostic_rows = gate_diagnostics(kind, seed, gates if kind == "DYNAMIC" else None, checkpoint, layout)
            gate_reports[name] = diagnostic
            gate_channel_rows.extend(diagnostic_rows)
        del arrays, gates
        print(f"analyzed {name}", flush=True)

    for reference, d_name in (("M", "C1-M"), ("BS", "C1-BS")):
        for split in ("R2", "R3"):
            dn = d_numerical[f"{d_name}_{split}"]
            numerical[f"{reference}_{split}"] = {
                "representation": reference, "split": split, "dimension": dn["dimension"],
                "observations": dn["observations"], "elements": dn["elements"], "finite_elements": dn["finite_elements"], "nonfinite_elements": dn["nonfinite_elements"],
                "feature_rms": dn["feature_rms"], "feature_norm_distribution": dn["feature_norm_distribution"], "channel_variance_distribution": dn["channel_variance_distribution"],
                "near_constant_channels": dn["near_constant_channels"], "participation_rank": dn["covariance_participation_rank"], "entropy_rank": dn["covariance_entropy_rank"],
                "numerical_rank": dn["numerical_rank"], "positive_spectrum_condition": dn["positive_spectrum_condition"], "regularized_condition_1e-6": dn["regularized_condition_1e-6"],
                "reused_from_experiment_d": True,
            }
            if d_name in ("C1-M", "C1-BS"):
                dt = d_temporal[f"{d_name}_{split}"]
                temporal_metrics[f"{reference}_{split}"] = {**dt, "representation": reference, "reused_from_experiment_d": True}
            db = d_backbone[f"{d_name}_{split}"]
            dc = d_chi[f"{d_name}_{split}"]
            probes[f"{reference}_{split}"] = {
                "representation": reference, "split": split,
                "torsions": {"phi": db["phi"]["all"], "psi": db["psi"]["all"], "chi1": dc},
                "backbone_combined_mae_deg": db["all_combined_angular_mae_deg"],
                "backbone_mean_component_rmse": db["all_mean_component_rmse"],
                "chi1_mae_deg": dc["residue_balanced_angular_mae_deg"],
                "reused_from_experiment_d": True,
            }
    for row in d_per_residue:
        if row["readout"] in ("C1-M", "C1-BS") and row.get("subset", "valid_chi1") in ("all", "valid_chi1"):
            per_residue_rows.append({
                "representation": "M" if row["readout"] == "C1-M" else "BS",
                "split": row["split"], "torsion": row["torsion"],
                "residue_index": int(row["residue_index"]), "samples": int(row["samples"]),
                "sin_rmse": float(row["sin_rmse"]), "cos_rmse": float(row["cos_rmse"]),
                "angular_mae_deg": float(row["angular_mae_deg"]), "angular_rmse_deg": float(row["angular_rmse_deg"]),
            })

    write_json(args.output_root / "gate_diagnostics.json", gate_reports)
    write_csv(args.output_root / "gate_diagnostics_channels.csv", gate_channel_rows)
    write_json(args.output_root / "numerical_health.json", numerical)
    write_csv(args.output_root / "numerical_health_channels.csv", channel_rows)
    write_csv(args.output_root / "numerical_health_covariance_spectrum.csv", spectrum_rows)
    write_json(args.output_root / "temporal_diagnostics.json", temporal_metrics)
    write_csv(args.output_root / "temporal_diagnostics_windows.csv", temporal_rows)
    backbone_output = {key: {"representation": value["representation"], "split": value["split"], "phi": value["torsions"]["phi"], "psi": value["torsions"]["psi"], "combined_mae_deg": value["backbone_combined_mae_deg"], "mean_component_rmse": value["backbone_mean_component_rmse"]} for key, value in probes.items()}
    chi_output = {key: {"representation": value["representation"], "split": value["split"], **value["torsions"]["chi1"]} for key, value in probes.items()}
    write_json(args.output_root / "backbone_probe.json", {"config": FROZEN_CONFIG["structural_probe"], "scalers": scaler_records, "models": model_records, "metrics": backbone_output})
    write_csv(args.output_root / "backbone_probe_per_residue.csv", [row for row in per_residue_rows if row["torsion"] in ("phi", "psi")])
    write_json(args.output_root / "chi1_probe.json", {"mapping": FROZEN_CONFIG["structural_probe"]["chi1_mapping"], "metrics": chi_output})
    write_csv(args.output_root / "chi1_probe_per_residue.csv", [row for row in per_residue_rows if row["torsion"] == "chi1"])

    comparison = compare(reconstruction, probes, args.output_root)
    checkpoint_records = [file_record(path) for path in checkpoints.values()]
    d_hashes_unchanged = all(sha256(row["cache"]) == row["cache_sha256"] for row in rows)
    integrity = {
        "status": "EXPERIMENT_E_COMPLETE",
        "experiment_d_caches_unchanged": d_hashes_unchanged,
        "experiment_d_archives_modified": False,
        "visnet_inference_run": False,
        "all_seeds_reported": all(any(row["seed"] == seed for row in reconstruction if row["kind"] == kind) for kind in MODEL_KINDS for seed in SEEDS),
        "structural_targets_used_for_fusion_training": False,
        "R3_used_for_training_or_selection": False,
        "favorable_seed_selected": False,
        "followup_attention_started": False,
        "prohibited_outputs_computed": [],
        "checkpoints": checkpoint_records,
        "experiment_e_implementation": file_record(Path(__file__)),
    }
    if not d_hashes_unchanged:
        integrity["status"] = "EXPERIMENT_E_INTEGRITY_FAILURE"
    write_json(args.output_root / "integrity_completion_audit.json", integrity)
    report_path = args.output_root / "REPORT.md"
    report_path.write_text(render_report(training_history, reconstruction, numerical, temporal_metrics, probes, gate_reports, comparison, integrity), encoding="utf-8")
    write_json(args.output_root / "EXPERIMENT_E_COMPLETE.json", {"status": integrity["status"], "report": file_record(report_path), "integrity_audit_sha256": sha256(args.output_root / "integrity_completion_audit.json"), "followup_started": False})
    print(json.dumps({"status": integrity["status"], "output": str(args.output_root)}, indent=2))
    if integrity["status"] != "EXPERIMENT_E_COMPLETE":
        raise SystemExit("Experiment E integrity failure")


if __name__ == "__main__":
    main()
