#!/usr/bin/env python3
"""Encoder-free physical collective variables for PACER four-context MD.

The module deliberately avoids learned embeddings.  It turns atom14 receptor
coordinates into distances, side-chain distances and region compactness, then
forms factorial contrasts only after each trajectory/window is summarized.
Frames are never treated as independent compounds or replicas.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping

import numpy as np


BACKBONE_SLOTS = (0, 1, 2, 3)
CONTEXTS = ("apo", "probe_only", "candidate_no_probe", "candidate_probe")
ATOM14_COUNTS = {
    "A": 5, "R": 11, "N": 8, "D": 8, "C": 5, "Q": 9, "E": 9, "G": 4,
    "H": 10, "I": 8, "L": 8, "K": 9, "M": 8, "F": 11, "P": 7, "S": 6,
    "T": 7, "W": 14, "Y": 12, "V": 7,
}


def load_mapping(path: Path) -> dict[int, int]:
    rows = list(csv.DictReader(path.open(encoding="utf-8-sig", newline="")))
    mapping = {int(row["structure_resid"]): int(row["embedding_index"]) for row in rows}
    if len(mapping) != len(rows):
        raise ValueError("duplicate structure residue in mapping")
    return mapping


def load_sequence(path: Path) -> str:
    rows = list(csv.DictReader(path.open(encoding="utf-8-sig", newline="")))
    if len(rows) != 1 or not rows[0].get("seqres"):
        raise ValueError(f"missing or ambiguous seqres: {path}")
    return rows[0]["seqres"].strip()


def _valid_atoms(residue: np.ndarray) -> np.ndarray:
    finite = np.isfinite(residue).all(axis=-1)
    nonzero = np.linalg.norm(residue, axis=-1) > 1e-6
    return finite & nonzero


def sidechain_centroid(atom14: np.ndarray, sequence: str) -> np.ndarray:
    """Return (frames,residues,3), falling back to CA for Gly/missing side chain."""
    out = np.empty(atom14.shape[:2] + (3,), dtype=np.float64)
    if len(sequence) != atom14.shape[1]:
        raise ValueError("sequence length does not match atom14 residues")
    for resid, aa in enumerate(sequence):
        atoms = np.asarray(atom14[:, resid], dtype=np.float64)
        count = ATOM14_COUNTS.get(aa)
        if count is None:
            raise ValueError(f"unsupported residue {aa!r}")
        side_slots = np.arange(4, count)
        for frame in range(atom14.shape[0]):
            out[frame, resid] = atoms[frame, side_slots].mean(axis=0) if side_slots.size else atoms[frame, 1]
    return out


def _distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.linalg.norm(a - b, axis=-1)


def extract_frame_features(atom14: np.ndarray, sequence: str, mapping: Mapping[int, int], config: Mapping) -> tuple[list[str], np.ndarray]:
    if atom14.ndim != 4 or atom14.shape[2:] != (14, 3):
        raise ValueError(f"expected (frames,residues,14,3), got {atom14.shape}")
    if not np.isfinite(atom14).all():
        raise ValueError("non-finite atom14 coordinates")
    ca = np.asarray(atom14[:, :, 1, :], dtype=np.float64)
    side = sidechain_centroid(atom14, sequence)
    names: list[str] = []
    columns: list[np.ndarray] = []

    def idx(resid: int) -> int:
        if int(resid) not in mapping:
            raise KeyError(f"M4 residue {resid} absent from mapping")
        return mapping[int(resid)]

    for i, j in config["ca_distance_pairs"]:
        names.append(f"ca_dist_{i}_{j}")
        columns.append(_distance(ca[:, idx(i)], ca[:, idx(j)]))
    for i, j in config["sidechain_distance_pairs"]:
        names.append(f"sc_dist_{i}_{j}")
        columns.append(_distance(side[:, idx(i)], side[:, idx(j)]))

    region_centroids = {}
    for region, residues in config["regions"].items():
        ri = [idx(r) for r in residues]
        xyz = ca[:, ri, :]
        centroid = xyz.mean(axis=1)
        region_centroids[region] = centroid
        rg = np.sqrt(np.mean(np.sum((xyz - centroid[:, None, :]) ** 2, axis=-1), axis=1))
        names.append(f"ca_rg_{region}")
        columns.append(rg)
    for left, right in config["region_centroid_pairs"]:
        names.append(f"centroid_dist_{left}__{right}")
        columns.append(_distance(region_centroids[left], region_centroids[right]))

    values = np.stack(columns, axis=1)
    if not np.isfinite(values).all():
        raise ValueError("non-finite physical features")
    return names, values


def summarize_trace(values: np.ndarray) -> dict[str, np.ndarray]:
    if values.ndim != 2 or values.shape[0] < 3:
        raise ValueError("need at least three frames")
    centered = values - values.mean(axis=0, keepdims=True)
    denom = np.sum(centered[:-1] ** 2, axis=0)
    lag1 = np.divide(
        np.sum(centered[:-1] * centered[1:], axis=0), denom,
        out=np.zeros(values.shape[1]), where=denom > 1e-12,
    )
    return {
        "mean": values.mean(axis=0),
        "std": values.std(axis=0, ddof=1),
        "q10": np.quantile(values, 0.10, axis=0),
        "q90": np.quantile(values, 0.90, axis=0),
        "lag1": np.clip(lag1, -1.0, 1.0),
    }


def factorial_contrast(context_vectors: Mapping[str, np.ndarray], signs: Mapping[str, float]) -> np.ndarray:
    missing = set(signs) - set(context_vectors)
    if missing:
        raise KeyError(f"missing contexts: {sorted(missing)}")
    shape = np.asarray(next(iter(context_vectors.values()))).shape
    out = np.zeros(shape, dtype=np.float64)
    for context, sign in signs.items():
        value = np.asarray(context_vectors[context], dtype=np.float64)
        if value.shape != shape:
            raise ValueError("context vectors have inconsistent shapes")
        out += float(sign) * value
    return out


def robust_location_scale(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    center = np.median(matrix, axis=0)
    mad = np.median(np.abs(matrix - center), axis=0)
    scale = 1.4826 * mad
    fallback = np.std(matrix, axis=0, ddof=1)
    scale = np.where(scale > 1e-8, scale, np.where(fallback > 1e-8, fallback, 1.0))
    return center, scale


def cosine(a: np.ndarray, b: np.ndarray) -> float | None:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.clip(np.dot(a, b) / denom, -1.0, 1.0)) if denom > 1e-12 else None


def load_config(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))
