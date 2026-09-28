from __future__ import annotations

import numpy as np


AXIS_SIGNS = {
    "dAGO": {
        "candidate": 1.0,
        "apo": -1.0,
    },
    "dPAM": {
        "candidate_probe": 1.0,
        "probe": -1.0,
    },
    "dINT": {
        "candidate_probe": 1.0,
        "probe": -1.0,
        "candidate": -1.0,
        "apo": 1.0,
    },
}


def preprocess_features(features, median, scale):
    """
    Frozen G2-C preprocessing:
      1. per-frame mean-center over the residue axis;
      2. frozen G2-B channel median/MAD scaling;
      3. return float32 for the RFF map.
    """
    a = np.asarray(features, dtype=np.float32)
    med = np.asarray(median, dtype=np.float64)
    scl = np.asarray(scale, dtype=np.float64)

    if a.ndim != 3:
        raise ValueError("features must be frames x residues x channels")
    if med.shape != (a.shape[-1],) or scl.shape != (a.shape[-1],):
        raise ValueError("median/scale shape mismatch")
    if not np.isfinite(a).all():
        raise ValueError("nonfinite features")
    if not np.isfinite(med).all() or not np.isfinite(scl).all():
        raise ValueError("nonfinite preprocessing parameters")
    if (scl <= 0).any():
        raise ValueError("scale must be positive")

    a = a - a.mean(axis=1, keepdims=True)
    return ((a - med) / scl).astype(np.float32)


def rff_map(features, weights, bias):
    """
    Exact frozen G2-C RFF form:
      sqrt(2/D) * cos(x W + b)
    """
    a = np.asarray(features, dtype=np.float32)
    w = np.asarray(weights, dtype=np.float32)
    b = np.asarray(bias, dtype=np.float32)

    if a.ndim != 3:
        raise ValueError("features must be frames x residues x channels")
    if w.ndim != 2 or w.shape[0] != a.shape[-1]:
        raise ValueError("RFF weight shape mismatch")
    if b.shape != (w.shape[1],):
        raise ValueError("RFF bias shape mismatch")

    x = a.reshape(-1, a.shape[-1])
    z = np.sqrt(2.0 / len(b)) * np.cos(x @ w + b)

    return z.reshape(
        a.shape[0],
        a.shape[1],
        len(b),
    ).astype(np.float32)


def temporal_block_mean(mapped, block_frames):
    a = np.asarray(mapped)

    if a.ndim != 3:
        raise ValueError("mapped tensor must be frames x residues x RFF")
    if block_frames < 1 or a.shape[0] % block_frames:
        raise ValueError("invalid block partition")

    n_blocks = a.shape[0] // block_frames

    return a.reshape(
        n_blocks,
        block_frames,
        a.shape[1],
        a.shape[2],
    ).mean(axis=1)


def signed_block_contrast(context_blocks, axis):
    """
    Context arrays must already be aligned:
      block index b in every context refers to the same paired time block.
    """
    if axis not in AXIS_SIGNS:
        raise ValueError(f"unknown axis: {axis}")

    signs = AXIS_SIGNS[axis]
    missing = set(signs) - set(context_blocks)
    if missing:
        raise ValueError(f"missing contexts: {sorted(missing)}")

    arrays = {
        c: np.asarray(context_blocks[c])
        for c in signs
    }

    shapes = {a.shape for a in arrays.values()}
    if len(shapes) != 1:
        raise ValueError("paired context block shapes differ")

    out = np.zeros(
        next(iter(arrays.values())).shape,
        dtype=np.float32,
    )

    for context, sign in signs.items():
        out += np.float32(sign) * arrays[context].astype(
            np.float32,
            copy=False,
        )

    return out


def region_block_vectors(nodeblocks, residue_indices):
    """
    Reference-compatible regional vectors:
      mean over selected residues for each temporal block.
    """
    a = np.asarray(nodeblocks)
    ix = np.asarray(residue_indices, dtype=np.int64)

    if a.ndim != 3:
        raise ValueError("nodeblocks must be blocks x residues x RFF")
    if ix.ndim != 1 or len(ix) == 0:
        raise ValueError("empty/invalid residue index")

    return a[:, ix, :].mean(axis=1, dtype=np.float64)


def pooled_replica_vector(block_vectors):
    """
    One independent replica vector:
      mean over its temporal blocks.
    Blocks are not independent biological replicates.
    """
    a = np.asarray(block_vectors, dtype=np.float64)

    if a.ndim != 2:
        raise ValueError("block_vectors must be blocks x RFF")

    return a.mean(axis=0, dtype=np.float64)


def cosine(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    den = np.linalg.norm(a) * np.linalg.norm(b)
    if den <= 0:
        return None

    return float(np.dot(a, b) / den)


def summarize_three_replicas(replica_vectors):
    expected = {1, 2, 3}
    if set(replica_vectors) != expected:
        raise ValueError("expected exactly replicas R1/R2/R3")

    v = {
        r: np.asarray(replica_vectors[r], dtype=np.float64)
        for r in expected
    }

    if len({x.shape for x in v.values()}) != 1:
        raise ValueError("replica vector shape mismatch")

    mean_vector = np.mean(
        np.stack([v[1], v[2], v[3]]),
        axis=0,
        dtype=np.float64,
    )

    return {
        "replica_norms": {
            r: float(np.linalg.norm(v[r]))
            for r in expected
        },
        "mean_replica_vector_norm": float(
            np.mean([np.linalg.norm(v[r]) for r in expected])
        ),
        "norm_of_mean_vector": float(np.linalg.norm(mean_vector)),
        "pairwise_direction_cosines": {
            "R1_vs_R2": cosine(v[1], v[2]),
            "R1_vs_R3": cosine(v[1], v[3]),
            "R2_vs_R3": cosine(v[2], v[3]),
        },
        "mean_vector": mean_vector,
        "qualification_gate": "NOT_DEFINED",
    }
