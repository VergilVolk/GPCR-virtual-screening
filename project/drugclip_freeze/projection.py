"""Frozen projection-head apply (Model B).

The deployed weights are materialised two-layer MLPs:

    linear1: (512, 512)   linear2: (128, 512)
    forward(x) = L2_normalize(linear2(relu(linear1(x))))

This is byte-for-byte the forward pass used when the weights were produced.
It is re-expressed here functionally so that no *training* module
(`finetune_drugclip_*`, `train_drugclip_*`) is ever imported into the
production path.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

PROJECTION_KEYS = ("linear1.weight", "linear1.bias", "linear2.weight", "linear2.bias")
EXPECTED_INPUT_DIM = 512
EXPECTED_OUTPUT_DIM = 128


class ProjectionError(ValueError):
    """Raised when a projection state is missing, malformed or inconsistent."""


def _load_state_arrays(path: Path) -> dict[str, np.ndarray]:
    """Load a projection `.pt` without importing torch.

    The deployed files are plain zip archives containing a pickled dict of
    tensors; `torch.load` is the canonical reader.  Torch is used when
    available and the call is strictly read-only (`weights_only=False` is
    required because the payload is a dict, not a bare state_dict).
    """
    import torch  # local import: this module stays importable without torch

    payload = torch.load(str(path), map_location="cpu", weights_only=False)
    if not isinstance(payload, dict):
        raise ProjectionError(f"projection payload is not a dict: {path}")
    return payload


def load_head(payload: dict[str, Any], key: str) -> dict[str, np.ndarray]:
    state = payload.get(key)
    if state is None:
        raise ProjectionError(f"projection payload lacks {key!r}")
    missing = [k for k in PROJECTION_KEYS if k not in state]
    if missing:
        raise ProjectionError(f"{key} is missing tensor keys: {missing}")
    unexpected = sorted(set(state) - set(PROJECTION_KEYS))
    if unexpected:
        raise ProjectionError(f"{key} has unexpected tensor keys: {unexpected}")
    out: dict[str, np.ndarray] = {}
    for name in PROJECTION_KEYS:
        tensor = state[name]
        array = tensor.detach().cpu().numpy().astype(np.float32, copy=False)
        out[name] = array
    w1, b1 = out["linear1.weight"], out["linear1.bias"]
    w2, b2 = out["linear2.weight"], out["linear2.bias"]
    if w1.shape != (EXPECTED_INPUT_DIM, EXPECTED_INPUT_DIM):
        raise ProjectionError(f"{key}.linear1.weight shape {w1.shape} != {(EXPECTED_INPUT_DIM, EXPECTED_INPUT_DIM)}")
    if b1.shape != (EXPECTED_INPUT_DIM,):
        raise ProjectionError(f"{key}.linear1.bias shape {b1.shape} != {(EXPECTED_INPUT_DIM,)}")
    if w2.shape != (EXPECTED_OUTPUT_DIM, EXPECTED_INPUT_DIM):
        raise ProjectionError(f"{key}.linear2.weight shape {w2.shape} != {(EXPECTED_OUTPUT_DIM, EXPECTED_INPUT_DIM)}")
    if b2.shape != (EXPECTED_OUTPUT_DIM,):
        raise ProjectionError(f"{key}.linear2.bias shape {b2.shape} != {(EXPECTED_OUTPUT_DIM,)}")
    if not all(np.isfinite(v).all() for v in out.values()):
        raise ProjectionError(f"{key} contains non-finite weights")
    return out


def load_projection_file(path: Path) -> dict[str, dict[str, np.ndarray]]:
    payload = _load_state_arrays(Path(path))
    return {"mol_project": load_head(payload, "mol_project"),
            "pocket_project": load_head(payload, "pocket_project"),
            "metadata": {k: v for k, v in payload.items() if k not in ("mol_project", "pocket_project")}}


def forward(representations: np.ndarray, head: dict[str, np.ndarray]) -> np.ndarray:
    """Frozen forward pass: L2_normalize(linear2(relu(linear1(x))))."""
    x = np.asarray(representations, dtype=np.float32)
    if x.ndim != 2 or x.shape[1] != EXPECTED_INPUT_DIM:
        raise ProjectionError(f"expected representations of shape (n, {EXPECTED_INPUT_DIM}), got {x.shape}")
    hidden = np.maximum(x @ head["linear1.weight"].T + head["linear1.bias"], 0.0).astype(np.float32)
    out = (hidden @ head["linear2.weight"].T + head["linear2.bias"]).astype(np.float32)
    norm = np.linalg.norm(out, axis=1, keepdims=True)
    norm[norm == 0.0] = 1.0
    return (out / norm).astype(np.float32)
