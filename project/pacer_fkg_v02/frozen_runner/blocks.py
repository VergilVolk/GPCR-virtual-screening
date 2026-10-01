"""Parameterised temporal-block construction for arbitrary observation windows.

The frozen Phase-2 constructor (`run_phase2_calibration.construct_blocks`) is
hard-wired to the 1000-frame / 50-block long-MD window; the frozen Stage-B driver
carries an equivalent 400-frame / 20-block implementation.  Prospective runs need
the same arithmetic at a configurable window length.

This module therefore re-expresses that arithmetic with explicit geometry
parameters and PROVES bit-identity against both frozen implementations
(`prove_equivalence`) before it may be used.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from project.pacer_fkg_v02.frozen_runner import engines

ENDPOINT_RTOL = 2e-5
ENDPOINT_ATOL = 2e-6


def construct_blocks(
    frames: np.ndarray,
    *,
    block_frames: int,
    n_blocks: int,
    n_residues: int | None = None,
    bs_width: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return `(STATE_MOTION, SIGNED_DRIFT)` descriptors for one trajectory.

    Arithmetic is element-for-element the frozen Phase-2 / Stage-B definition:
    float32 storage, float64 reductions, one-frame internal differences only.
    """
    arr = np.asarray(frames, dtype=np.float32)
    residues = int(n_residues if n_residues is not None else engines.N_RESIDUES)
    width = int(bs_width if bs_width is not None else engines.P2_BS_WIDTH)
    expected = (n_blocks * block_frames, residues, width)
    if arr.shape != expected:
        raise ValueError(f"expected BS256 frames of shape {expected}, got {arr.shape}")
    blocks = arr.reshape(n_blocks, block_frames, residues, width)
    mu_z = blocks.mean(axis=1, dtype=np.float64).astype(np.float32)
    delta = np.diff(blocks, axis=1)
    mu_delta = delta.mean(axis=1, dtype=np.float64).astype(np.float32)
    rms_delta = np.sqrt(np.mean(np.square(delta, dtype=np.float64), axis=1)).astype(np.float32)
    state_motion = np.concatenate((mu_z, rms_delta), axis=-1).astype(np.float32, copy=False)
    return state_motion, mu_delta


def verify_endpoint_identity(frames: np.ndarray, mu_delta: np.ndarray, *, block_frames: int, n_blocks: int,
                             n_residues: int | None = None, bs_width: int | None = None) -> bool:
    """Frozen identity: mean(DeltaZ) == (Z_last - Z_first) / (block_frames - 1)."""
    residues = int(n_residues if n_residues is not None else engines.N_RESIDUES)
    width = int(bs_width if bs_width is not None else engines.P2_BS_WIDTH)
    blocks = np.asarray(frames).reshape(n_blocks, block_frames, residues, width)
    endpoint = (blocks[:, -1] - blocks[:, 0]) / np.float32(block_frames - 1)
    return bool(np.allclose(mu_delta, endpoint, rtol=ENDPOINT_RTOL, atol=ENDPOINT_ATOL))


def prove_equivalence(*, n_residues: int | None = None, bs_width: int | None = None,
                      seed: int = 20261002) -> dict[str, Any]:
    """Bit-identity proof against the two frozen block implementations."""
    residues = int(n_residues if n_residues is not None else engines.N_RESIDUES)
    width = int(bs_width if bs_width is not None else engines.P2_BS_WIDTH)
    rng = np.random.default_rng(seed)
    report: dict[str, Any] = {"n_residues": residues, "bs_width": width, "checks": {}}

    # 1) long-MD geometry: 1000 frames / 50 blocks vs frozen Phase-2
    long_frames = rng.standard_normal((engines.P2_N_FRAMES, residues, width), dtype=np.float32)
    sm, sd = construct_blocks(long_frames, block_frames=engines.P2_BLOCK_FRAMES, n_blocks=engines.P2_N_BLOCKS,
                              n_residues=residues, bs_width=width)
    ref_sm, ref_sd = engines.construct_blocks(long_frames)
    report["checks"]["phase2_1000x50_STATE_MOTION_bit_identical"] = bool(np.array_equal(sm, ref_sm))
    report["checks"]["phase2_1000x50_SIGNED_DRIFT_bit_identical"] = bool(np.array_equal(sd, ref_sd))
    report["checks"]["phase2_endpoint_identity"] = bool(
        engines.verify_endpoint_identity(long_frames, sd))

    # 2) 20 ns geometry: 400 frames / 20 blocks vs the frozen Stage-B driver
    short_frames = rng.standard_normal((400, residues, width), dtype=np.float32)
    sm2, sd2 = construct_blocks(short_frames, block_frames=20, n_blocks=20,
                                n_residues=residues, bs_width=width)
    try:
        from project.pacer_fkg_v02.run_cm00734_stage_b_phase2a_frozen_apply_v01 import construct as stage_b_construct

        ref = stage_b_construct(short_frames)
        report["checks"]["stage_b_400x20_STATE_MOTION_bit_identical"] = bool(
            np.array_equal(sm2, ref["STATE_MOTION"]))
        report["checks"]["stage_b_400x20_SIGNED_DRIFT_bit_identical"] = bool(
            np.array_equal(sd2, ref["SIGNED_DRIFT"]))
        report["stage_b_driver_available"] = True
    except Exception as exc:  # pragma: no cover - reported, never silently skipped
        report["stage_b_driver_available"] = False
        report["stage_b_driver_error"] = f"{type(exc).__name__}: {exc}"

    report["checks"]["stage_b_endpoint_identity"] = bool(
        verify_endpoint_identity(short_frames, sd2, block_frames=20, n_blocks=20,
                                 n_residues=residues, bs_width=width))
    report["all_identical"] = all(
        v for k, v in report["checks"].items() if k.endswith("bit_identical")
    ) if report.get("stage_b_driver_available") else False
    return report
