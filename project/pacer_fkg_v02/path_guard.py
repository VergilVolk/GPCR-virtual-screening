"""Explicit write-path isolation for PACER-FKG v02.

Call :func:`guard_write_path` immediately before every v02 write.  The check is
based on resolved paths and ``os.path.commonpath``; it deliberately does not use
string-prefix matching.
"""

from __future__ import annotations

import os
from pathlib import Path


ALLOWED_RELATIVE_ROOTS = (
    "project/results/pacer_fkg_v02_longmd_v01",
    "project/cache/pacer_fkg_v02_longmd_v01",
    "project/pacer_fkg_v02",
)

# Includes the seven initially manifested result roots and all additional v01
# implementation/input/calibration roots discovered during the Phase-0 audit.
PROTECTED_RELATIVE_ROOTS = (
    "project/results/pacer_dc_fkg_G1_diagnostics_v01",
    "project/results/pacer_dc_fkg_G2A_PBC_v01",
    "project/results/pacer_dc_fkg_G2A_v01",
    "project/results/pacer_dc_fkg_G2A_v02",
    "project/results/pacer_dc_fkg_G2B_robustness_v01",
    "project/results/pacer_dc_fkg_G2B_v01",
    "project/results/pacer_dc_fkg_G2C_v01",
    "project/results/pacer_dc_fkg_R2R3_v01",
    "project/results/pacer_dc_fkg_U2_erratum_v02",
    "project/results/pacer_dc_four_context_v01",
    "project/results/pacer_dc_geom2vec_pilot_v01",
    "project/results/pacer_dc_geom2vec_R2R3_full_v01",
    "project/results/pacer_dc_G2A_source_v01",
    "project/results/pacer_dc_production_v01",
    "project/results/pacer_dc_membrane_reference_v01",
    "project/results/pacer_dc_restraint_release_v01",
    "project/results/pacer_dc_common_protein_v01",
)


class UnsafeWritePath(ValueError):
    """Raised when a requested output is protected or outside v02 roots."""


def _resolved(path: os.PathLike[str] | str) -> Path:
    return Path(path).expanduser().resolve(strict=False)


def _within(path: Path, root: Path) -> bool:
    """Case-normalized, component-aware containment test."""
    try:
        return os.path.normcase(os.path.commonpath((str(path), str(root)))) == os.path.normcase(str(root))
    except ValueError:  # different Windows drives
        return False


def guard_write_path(path: os.PathLike[str] | str, repo_root: os.PathLike[str] | str) -> Path:
    """Return the resolved safe path or reject it before any write occurs."""
    candidate = _resolved(path)
    root = _resolved(repo_root)
    protected = tuple(_resolved(root / item) for item in PROTECTED_RELATIVE_ROOTS)
    allowed = tuple(_resolved(root / item) for item in ALLOWED_RELATIVE_ROOTS)
    if any(_within(candidate, item) for item in protected):
        raise UnsafeWritePath(f"protected historical PACER-FKG path: {candidate}")
    if not any(_within(candidate, item) for item in allowed):
        raise UnsafeWritePath(f"path is outside the approved v02 output roots: {candidate}")
    return candidate
