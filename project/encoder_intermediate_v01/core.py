"""Strict model loading and fixed C0/C1/C2 extraction for Experiment C v01."""

from __future__ import annotations

import hashlib
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import torch

EXPECTED_CHECKPOINT_SHA256 = (
    "b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417"
)
PINNED_GEOM2VEC_COMMIT = "371d642ec1061664f16e49fcac702d07fc8d0b51"

FROZEN_CONFIG = {
    "hidden_channels": 64,
    "num_layers": 6,
    "num_heads": 8,
    "num_rbf": 64,
    "cutoff": 7.5,
    "max_num_neighbors": 32,
    "lmax": 1,
    "max_z": 100,
    "vecnorm_type": "max_min",
    "trainable_vecnorm": True,
    "trainable_rbf": False,
    "vertex": False,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def construct_visnet() -> torch.nn.Module:
    """Construct ViSNet explicitly; the permissive installed factory is not used."""
    from geom2vec.models.representation.visnet import ViSNet

    return ViSNet(**FROZEN_CONFIG)


def verify_and_strict_load(checkpoint_path: Path) -> tuple[torch.nn.Module, dict]:
    actual_hash = sha256(checkpoint_path)
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    if not isinstance(checkpoint, dict) or not checkpoint:
        raise TypeError("checkpoint is not a nonempty state-dict mapping")
    if not all(isinstance(key, str) and torch.is_tensor(value) for key, value in checkpoint.items()):
        raise TypeError("checkpoint is not a flat string-to-tensor state dict")

    model = construct_visnet().cpu()
    model_state = model.state_dict()
    parameter_keys = {name for name, _ in model.named_parameters()}
    buffer_keys = {name for name, _ in model.named_buffers()}
    rows = []
    for key in sorted(set(model_state) | set(checkpoint)):
        expected = model_state.get(key)
        observed = checkpoint.get(key)
        row = {
            "key": key,
            "kind": "parameter" if key in parameter_keys else "persistent_buffer" if key in buffer_keys else "unexpected",
            "model_shape": list(expected.shape) if expected is not None else None,
            "checkpoint_shape": list(observed.shape) if observed is not None else None,
            "model_dtype": str(expected.dtype) if expected is not None else None,
            "checkpoint_dtype": str(observed.dtype) if observed is not None else None,
            "checkpoint_finite": bool(torch.isfinite(observed).all()) if observed is not None else None,
        }
        row["shape_match"] = expected is not None and observed is not None and expected.shape == observed.shape
        row["dtype_match"] = expected is not None and observed is not None and expected.dtype == observed.dtype
        rows.append(row)

    missing = sorted(set(model_state) - set(checkpoint))
    unexpected = sorted(set(checkpoint) - set(model_state))
    shape_mismatches = sorted(
        key for key in set(model_state) & set(checkpoint) if model_state[key].shape != checkpoint[key].shape
    )
    dtype_mismatches = sorted(
        key for key in set(model_state) & set(checkpoint) if model_state[key].dtype != checkpoint[key].dtype
    )
    nonfinite = sorted(key for key, value in checkpoint.items() if not bool(torch.isfinite(value).all()))
    matched = set(model_state) & set(checkpoint) - set(shape_mismatches) - set(dtype_mismatches) - set(nonfinite)
    matched_parameters = parameter_keys & matched
    matched_state = set(model_state) & matched

    def numel(keys: set[str], mapping: dict[str, torch.Tensor]) -> int:
        return sum(mapping[key].numel() for key in keys)

    strict_error = None
    strict_ok = False
    try:
        model.load_state_dict(checkpoint, strict=True)
        strict_ok = True
    except RuntimeError as exc:
        strict_error = str(exc)

    complete = (
        actual_hash == EXPECTED_CHECKPOINT_SHA256
        and not missing
        and not unexpected
        and not shape_mismatches
        and not dtype_mismatches
        and not nonfinite
        and strict_ok
    )
    report = {
        "checkpoint_path": str(checkpoint_path.resolve()),
        "expected_sha256": EXPECTED_CHECKPOINT_SHA256,
        "actual_sha256": actual_hash,
        "sha256_match": actual_hash == EXPECTED_CHECKPOINT_SHA256,
        "state_dict_structure": "flat string-to-tensor mapping",
        "checkpoint_key_count": len(checkpoint),
        "model_state_key_count": len(model_state),
        "missing_keys": missing,
        "unexpected_keys": unexpected,
        "shape_mismatches": shape_mismatches,
        "dtype_mismatches": dtype_mismatches,
        "nonfinite_keys": nonfinite,
        "parameter_count_coverage": {
            "matched": len(matched_parameters),
            "total": len(parameter_keys),
            "fraction": len(matched_parameters) / len(parameter_keys),
        },
        "parameter_element_coverage": {
            "matched": numel(matched_parameters, model_state),
            "total": numel(parameter_keys, model_state),
            "fraction": numel(matched_parameters, model_state) / numel(parameter_keys, model_state),
        },
        "state_element_coverage": {
            "matched": numel(matched_state, model_state),
            "total": numel(set(model_state), model_state),
            "fraction": numel(matched_state, model_state) / numel(set(model_state), model_state),
        },
        "strict_load_attempted": True,
        "strict_load_ok": strict_ok,
        "strict_load_error": strict_error,
        "complete_strict_coverage": complete,
        "keys": rows,
    }
    return model, report


def invariant_residue_features(
    x: torch.Tensor,
    v: torch.Tensor,
    residue_index: torch.Tensor,
    n_residues: int,
) -> torch.Tensor:
    """The unchanged production [scalar, vector norm] equal-atom-mean readout."""
    invariant = torch.cat((x, torch.linalg.vector_norm(v, dim=1)), dim=-1)
    output = torch.zeros((n_residues, invariant.shape[-1]), dtype=invariant.dtype, device=invariant.device)
    counts = torch.zeros((n_residues, 1), dtype=invariant.dtype, device=invariant.device)
    output.index_add_(0, residue_index, invariant)
    counts.index_add_(0, residue_index, torch.ones((len(residue_index), 1), device=invariant.device))
    return output / counts.clamp_min(1.0)


def pool_batch(
    state: tuple[torch.Tensor, torch.Tensor],
    residue_index: torch.Tensor,
    n_residues: int,
    batch_size: int,
    n_atoms: int,
) -> torch.Tensor:
    x, v = state
    x = x.reshape(batch_size, n_atoms, 64)
    v = v.reshape(batch_size, n_atoms, 3, 64)
    return torch.stack(
        [invariant_residue_features(x[i], v[i], residue_index, n_residues) for i in range(batch_size)]
    )


@contextmanager
def capture_fixed_states(model: torch.nn.Module) -> Iterator[dict[str, tuple[torch.Tensor, torch.Tensor]]]:
    """Capture exactly C1 and C2 with pre-hooks, and always remove hooks."""
    representation = model.representation_model
    captured: dict[str, torch.Tensor | tuple[torch.Tensor, torch.Tensor]] = {}

    def midpoint_pre_hook(_module, inputs):
        captured["C1"] = (inputs[0], inputs[1])

    def scalar_pre_hook(_module, inputs):
        captured["C2_x"] = inputs[0]

    def vector_pre_hook(_module, inputs):
        captured["C2_v"] = inputs[0]

    handles = [
        representation.vis_mp_layers[3].register_forward_pre_hook(midpoint_pre_hook),
        representation.out_norm.register_forward_pre_hook(scalar_pre_hook),
        representation.vec_out_norm.register_forward_pre_hook(vector_pre_hook),
    ]
    try:
        yield captured  # tensors remain valid for pooling inside inference mode
    finally:
        for handle in handles:
            handle.remove()


def forward_three_states(
    model: torch.nn.Module, z: torch.Tensor, pos: torch.Tensor, batch: torch.Tensor
) -> tuple[dict[str, tuple[torch.Tensor, torch.Tensor]], tuple[torch.Tensor, ...]]:
    """Extract C0/C1/C2 during one model forward."""
    with capture_fixed_states(model) as captured:
        returned = model(z=z, pos=pos, batch=batch)
        states = {
            "C0": (returned[0], returned[1]),
            "C1": captured["C1"],
            "C2": (captured["C2_x"], captured["C2_v"]),
        }
    return states, returned


def hook_count(model: torch.nn.Module) -> int:
    representation = model.representation_model
    targets = [representation.vis_mp_layers[3], representation.out_norm, representation.vec_out_norm]
    return sum(len(module._forward_pre_hooks) for module in targets)

