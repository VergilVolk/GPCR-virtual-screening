#!/usr/bin/env python
"""PACER-FKG: distributional four-context contrasts on a residue graph.

This is a training-free phase-1 estimator.  It does not pair frames from
independent trajectories.  Instead, each context is represented by an RKHS
kernel mean, and the 2x2 candidate-by-probe interaction is measured as
mu(C+A) - mu(A) - mu(C) + mu(0).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


CONTEXTS = ("apo", "probe_only", "candidate_no_probe", "candidate_probe")


def load_contexts(root: Path, pattern: str):
    arrays = {}
    for context in CONTEXTS:
        path = root / pattern.format(context=context)
        with np.load(path) as z:
            arrays[context] = z["residue_features"].astype(np.float64)
    shapes = {a.shape[1:] for a in arrays.values()}
    if len(shapes) != 1:
        raise ValueError(f"incompatible residue/channel shapes: {shapes}")
    return arrays


def preprocess(arrays, mode):
    values = arrays
    if mode == "global_residual_channel_robust":
        # Remove frame-wide encoder drift before asking which residue network
        # changed.  Global shift remains available as an explicit baseline,
        # but cannot masquerade as a localized allosteric path.
        values = {k: v - v.mean(axis=1, keepdims=True) for k, v in arrays.items()}
    elif mode != "channel_robust":
        raise ValueError(f"unknown preprocessing mode: {mode}")
    pooled = np.concatenate(list(values.values()), axis=0)
    flat = pooled.reshape(-1, pooled.shape[-1])
    median = np.median(flat, axis=0).reshape(1, 1, -1)
    mad = np.median(np.abs(flat - median.reshape(1, -1)), axis=0).reshape(1, 1, -1)
    scale = np.maximum(1.4826 * mad, 1e-6)
    return {k: (v - median) / scale for k, v in values.items()}


def median_bandwidth(x):
    d2 = np.sum((x[:, None, :] - x[None, :, :]) ** 2, axis=-1)
    values = d2[np.triu_indices(len(x), 1)]
    values = values[values > 0]
    return float(np.sqrt(np.median(values))) if len(values) else 1.0


def signed_kernel_stat(groups, signs, bandwidth=None):
    xs, weights = [], []
    for name, sign in signs.items():
        x = groups[name]
        xs.append(x)
        weights.extend([sign / len(x)] * len(x))
    x = np.concatenate(xs, axis=0)
    w = np.asarray(weights)
    bandwidth = bandwidth or median_bandwidth(x)
    d2 = np.sum((x[:, None, :] - x[None, :, :]) ** 2, axis=-1)
    kernel = np.exp(-d2 / (2.0 * bandwidth * bandwidth))
    biased_squared = float(w @ kernel @ w)

    # U-statistic: exclude the always-one diagonal within each context.  The
    # ordinary empirical RKHS norm is badly upward biased when n is small.
    unbiased_squared = 0.0
    offset_i = 0
    names = list(signs)
    for i, name_i in enumerate(names):
        ni = len(groups[name_i])
        block_i = slice(offset_i, offset_i + ni)
        offset_j = 0
        for j, name_j in enumerate(names):
            nj = len(groups[name_j])
            block_j = slice(offset_j, offset_j + nj)
            block = kernel[block_i, block_j]
            if i == j:
                if ni < 2:
                    raise ValueError("unbiased kernel statistic requires >=2 frames per context")
                estimate = (block.sum() - np.trace(block)) / (ni * (ni - 1))
            else:
                estimate = block.mean()
            unbiased_squared += signs[name_i] * signs[name_j] * float(estimate)
            offset_j += nj
        offset_i += ni
    return {
        "unbiased_squared": unbiased_squared,
        "unbiased_positive_root": max(unbiased_squared, 0.0) ** 0.5,
        "biased_squared": biased_squared,
        "bandwidth": bandwidth,
    }


def linear_norm(groups, signs):
    delta = sum(sign * groups[name].mean(axis=0) for name, sign in signs.items())
    return float(np.linalg.norm(delta))


def graph_diffuse(scores, graph, alpha=0.65, steps=20):
    n = len(graph["nodes"])
    adjacency = np.zeros((n, n), dtype=np.float64)
    for edge in graph["edges"]:
        i, j = edge["source"], edge["target"]
        weight = 1.0 if edge["edge_type"] == "backbone" else edge["support_fraction"]
        adjacency[i, j] = max(adjacency[i, j], weight)
        adjacency[j, i] = max(adjacency[j, i], weight)
    degree = adjacency.sum(axis=1, keepdims=True)
    transition = adjacency / np.maximum(degree, 1e-12)
    out = scores.copy()
    for _ in range(steps):
        out = (1.0 - alpha) * scores + alpha * transition @ out
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-root", type=Path, required=True)
    p.add_argument("--pattern", default="{context}_stride20.npz")
    p.add_argument("--graph", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--preprocess", choices=("channel_robust", "global_residual_channel_robust"),
                   default="global_residual_channel_robust")
    args = p.parse_args()

    graph = json.loads(args.graph.read_text(encoding="utf-8"))
    raw = load_contexts(args.embedding_root, args.pattern)
    arrays = preprocess(raw, args.preprocess)
    signs = {
        "synergy_interaction": {"candidate_probe": 1, "probe_only": -1, "candidate_no_probe": -1, "apo": 1},
        "intrinsic_agonism": {"candidate_no_probe": 1, "apo": -1},
        "conditional_pam_effect": {"candidate_probe": 1, "probe_only": -1},
    }

    node_rows = []
    node_scores = {axis: np.zeros(raw["apo"].shape[1]) for axis in signs}
    for index in range(raw["apo"].shape[1]):
        groups = {name: value[:, index, :] for name, value in arrays.items()}
        row = {"embedding_index": index, "site": graph["nodes"][index]["site"]}
        for axis, contrast in signs.items():
            estimate = signed_kernel_stat(groups, contrast)
            score = estimate["unbiased_positive_root"]
            node_scores[axis][index] = score
            row[f"{axis}_kernel_u2"] = estimate["unbiased_squared"]
            row[f"{axis}_kernel"] = score
            row[f"{axis}_kernel_biased2"] = estimate["biased_squared"]
            row[f"{axis}_linear"] = linear_norm(groups, contrast)
            row[f"{axis}_bandwidth"] = estimate["bandwidth"]
        node_rows.append(row)

    diffused = {axis: graph_diffuse(values, graph) for axis, values in node_scores.items()}
    region_rows = []
    for region, members in graph["regions"].items():
        idx = np.asarray([m["embedding_index"] for m in members], dtype=int)
        if len(idx) == 0:
            continue
        groups = {name: value[:, idx, :].reshape(len(value), -1) for name, value in arrays.items()}
        row = {"region": region, "n_residues": int(len(idx))}
        for axis, contrast in signs.items():
            estimate = signed_kernel_stat(groups, contrast)
            row[f"{axis}_kernel_u2"] = estimate["unbiased_squared"]
            row[f"{axis}_kernel"] = estimate["unbiased_positive_root"]
            row[f"{axis}_kernel_biased2"] = estimate["biased_squared"]
            row[f"{axis}_linear"] = linear_norm(groups, contrast)
            row[f"{axis}_bandwidth"] = estimate["bandwidth"]
            row[f"{axis}_graph_diffused_mean"] = float(diffused[axis][idx].mean())
        region_rows.append(row)

    top_nodes = {}
    for axis, values in diffused.items():
        order = np.argsort(values)[::-1][:15]
        top_nodes[axis] = [
            {"site": graph["nodes"][int(i)]["site"], "score": float(values[i])} for i in order
        ]

    report = {
        "method": "PACER-FKG_v01 (frozen Geom2Vec + factorial kernel mean + consensus residue graph)",
        "evidence_level": "single_replica_method_smoke",
        "preprocessing": args.preprocess,
        "frames_per_context": {k: int(len(v)) for k, v in raw.items()},
        "axes": {
            "synergy_interaction": "mu(C+A)-mu(A)-mu(C)+mu(0)",
            "intrinsic_agonism": "mu(C)-mu(0)",
            "conditional_pam_effect": "mu(C+A)-mu(A)",
        },
        "region_scores": region_rows,
        "top_graph_diffused_nodes": top_nodes,
        "node_scores": node_rows,
        "validation_required": [
            "repeat on >=3 independent replicas",
            "estimate uncertainty by replica/block bootstrap, never by treating frames as independent",
            "compare with linear mean, global pooling, fixed-region mean and non-graph kernel baselines",
            "require PAM-vs-compound110/inactive directional separation and distal-control specificity",
        ],
        "claim_boundary": "Smoke test only. No p-value, PAM classification, convergence or performance claim from one replica/five frames.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(args.output), "regions": len(region_rows), "nodes": len(node_rows)}, indent=2))


if __name__ == "__main__":
    main()
