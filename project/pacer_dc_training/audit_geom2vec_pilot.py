#!/usr/bin/env python
"""Audit Geom2Vec feasibility before PACER-DC adopts it as an encoder."""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np

from extract_geom2vec_atom14 import (
    create_model_from_checkpoint,
    embed,
    flatten_atom14,
    install_cpu_radius_graph_fallback,
    read_sequence,
)


CONTEXTS = ("apo", "probe_only", "candidate_no_probe", "candidate_probe")


def cosine(a, b):
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a.ravel(), b.ravel()) / denom) if denom else float("nan")


def load_regions(path: Path, n_residues: int):
    data = json.loads(path.read_text(encoding="utf-8"))["regions"]
    regions = {
        name: np.asarray([entry["embedding_index"] for entry in entries], dtype=int)
        for name, entries in data.items()
        if entries
    }
    regions["global"] = np.arange(n_residues, dtype=int)
    return regions


def region_summary(features, indices):
    # frames x residues x channels -> frames x channels
    frames = features[:, indices, :].mean(axis=1)
    centroid = frames.mean(axis=0)
    temporal_rms = float(np.sqrt(np.mean(np.sum((frames - centroid) ** 2, axis=1))))
    return frames, centroid, temporal_rms


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-dir", type=Path, required=True)
    p.add_argument("--atom14-dir", type=Path, required=True)
    p.add_argument("--region-map", type=Path, required=True)
    p.add_argument("--checkpoint", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    arrays = {}
    for context in CONTEXTS:
        path = args.embedding_dir / f"{context}_stride20.npz"
        with np.load(path) as z:
            arrays[context] = z["residue_features"].astype(np.float64)
    shapes = {tuple(v.shape) for v in arrays.values()}
    if len(shapes) != 1:
        raise ValueError(f"context embedding shapes differ: {shapes}")
    n_frames, n_residues, channels = next(iter(shapes))
    regions = load_regions(args.region_map, n_residues)

    region_rows = []
    for name, indices in regions.items():
        summaries = {c: region_summary(x, indices) for c, x in arrays.items()}
        centroids = {c: summaries[c][1] for c in CONTEXTS}
        within = float(np.mean([summaries[c][2] for c in CONTEXTS]))
        between = float(np.median([
            np.linalg.norm(centroids[a] - centroids[b]) for a, b in combinations(CONTEXTS, 2)
        ]))
        d_ago = centroids["candidate_no_probe"] - centroids["apo"]
        d_pam = centroids["candidate_probe"] - centroids["probe_only"]
        d_int = d_pam - d_ago
        region_rows.append({
            "region": name,
            "n_residues": int(len(indices)),
            "within_context_temporal_rms": within,
            "median_between_context_distance": between,
            "between_within_ratio": between / within if within else None,
            "dAGO_norm": float(np.linalg.norm(d_ago)),
            "dPAM_norm": float(np.linalg.norm(d_pam)),
            "dINT_norm": float(np.linalg.norm(d_int)),
            "dPAM_dAGO_cosine": cosine(d_pam, d_ago),
        })

    # Determinism and SE(3) invariance on one real apo frame.
    atom_path = args.atom14_dir / "apo_w000.npy"
    sequence = read_sequence(args.atom14_dir / "apo_w000.csv")
    coords = np.load(atom_path, mmap_mode="r")[:1]
    xyz, z, ridx = flatten_atom14(np.asarray(coords), sequence)
    fallback = install_cpu_radius_graph_fallback()
    model = create_model_from_checkpoint(str(args.checkpoint), device="cpu")
    base = embed(model, xyz, z, ridx, len(sequence), 1, "cpu")
    repeat = embed(model, xyz, z, ridx, len(sequence), 1, "cpu")
    translated = embed(model, xyz + np.asarray([11.25, -7.5, 3.0], np.float32), z, ridx, len(sequence), 1, "cpu")
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], np.float32)
    rotated = embed(model, xyz @ rotation.T, z, ridx, len(sequence), 1, "cpu")

    def comparison(other):
        return {
            "max_abs_delta": float(np.max(np.abs(base - other))),
            "relative_l2_delta": float(np.linalg.norm(base - other) / np.linalg.norm(base)),
            "cosine": cosine(base, other),
        }

    report = {
        "method": "PACER-DC frozen Geom2Vec ViSNet feasibility audit",
        "evidence_level": "encoder_feasibility_pilot",
        "input": {
            "candidate": "compound110",
            "replica": 1,
            "window": 0,
            "frames_per_context": n_frames,
            "residues": n_residues,
            "channels": channels,
            "note": "Five strided frames are descriptive, not independent replicates.",
        },
        "se3_and_determinism": {
            "repeat": comparison(repeat),
            "translation": comparison(translated),
            "rotation": comparison(rotated),
            "cpu_radius_graph_fallback": fallback,
        },
        "four_context_descriptive": region_rows,
        "decision": {
            "engineering_usable": bool(
                np.isfinite(base).all()
                and comparison(repeat)["relative_l2_delta"] < 1e-7
                and comparison(translated)["relative_l2_delta"] < 1e-5
                and comparison(rotated)["relative_l2_delta"] < 1e-5
            ),
            "cross_replica_qualified": False,
            "reason": "Only replica 1 atom14 tensors are committed; R2/R3 must run this exact frozen extractor before replacement is accepted.",
        },
        "claim_boundary": "Shows executable invariant encoding and non-collapsed four-context signal only; does not establish replica stability, PAM efficacy, or superiority to OneProt.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
