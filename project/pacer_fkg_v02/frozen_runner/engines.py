"""Import-only re-export layer over the frozen PACER-FKG v02 engines.

NOTHING in this file redefines a frozen numeric.  Every name below is an alias
to the object created by the frozen module.  `assert_frozen_numerics()` proves
at runtime that the imported values still equal the values recorded in the
frozen v02 freeze manifest / method spec; any drift aborts the run.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from project.pacer_fkg_v02 import run_phase2_calibration as _p2
from project.pacer_fkg_v02 import run_phase3_fkg as _p3
from project.pacer_fkg_v02.path_guard import PROTECTED_RELATIVE_ROOTS, UnsafeWritePath

# --- frozen Phase-2 (calibration / block construction) --------------------
construct_blocks = _p2.construct_blocks
verify_endpoint_identity = _p2.verify_endpoint_identity
P2_N_FRAMES = _p2.N_FRAMES
P2_N_RESIDUES = _p2.N_RESIDUES
P2_BS_WIDTH = _p2.BS_WIDTH
P2_BLOCK_FRAMES = _p2.BLOCK_FRAMES
P2_N_BLOCKS = _p2.N_BLOCKS
P2_RFF_FEATURES = _p2.RFF_FEATURES
P2_BASE_SEED = _p2.BASE_SEED
P2_BRANCHES = _p2.BRANCHES

# --- frozen Phase-3 (apply / graph / regions) -----------------------------
REQUIRED_FREEZE_SHA256 = _p3.REQUIRED_FREEZE_SHA256
FREEZE_PATH = _p3.FREEZE_PATH
GRAPH_PATH = _p3.GRAPH_PATH
REGION_PATH = _p3.REGION_PATH
CONTRAST_PATH = _p3.CONTRAST_PATH
CONTEXTS = _p3.CONTEXTS
REPLICAS = _p3.REPLICAS
REPLICA_LABELS = _p3.REPLICA_LABELS
BRANCHES = _p3.BRANCHES
CONTRASTS = _p3.CONTRASTS
HISTORICAL_AXIS = _p3.HISTORICAL_AXIS
ALPHA = _p3.ALPHA
DIFFUSION_STEPS = _p3.DIFFUSION_STEPS
N_BLOCKS = _p3.N_BLOCKS
N_RESIDUES = _p3.N_RESIDUES
RFF_FEATURES = _p3.RFF_FEATURES

verify_frozen_anchor = _p3.verify_frozen_anchor

_ANCHOR_CACHE: dict[str, Any] | None = None


def cached_anchor(force: bool = False) -> dict[str, Any]:
    """Verify the frozen anchor once per process.

    `verify_frozen_anchor` re-hashes roughly 3.8 GB of Phase-1/Phase-2 caches
    and the 41 protected v01 files.  The anchor is immutable for the lifetime of
    a process, so repeated calls are pure cost; the first call is still strict
    and fail-closed.
    """
    global _ANCHOR_CACHE
    if _ANCHOR_CACHE is None or force:
        assert_frozen_numerics()
        _ANCHOR_CACHE = verify_frozen_anchor()
    return _ANCHOR_CACHE
load_frozen_state = _p3.load_frozen_state
normalize_frozen = _p3.normalize_frozen
rff_frozen = _p3.rff_frozen
historical_contrast = _p3.historical_contrast
graph_edges = _p3.graph_edges
graph_transition = _p3.graph_transition
diffuse = _p3.diffuse
region_indices = _p3.region_indices
cosine = _p3.cosine
distribution = _p3.distribution
sha256 = _p3.sha256
load_json = _p3.load_json

# --- frozen Phase-1 (encoder) --------------------------------------------
phase1 = None  # populated lazily: importing torch/ViSNet is expensive.


def load_phase1_module():
    """Import the frozen Phase-1 encoder module on demand."""
    global phase1
    if phase1 is None:
        from project.pacer_fkg_v02 import run_phase1_bs256 as _p1

        phase1 = _p1
    return phase1


# --- frozen value snapshot ------------------------------------------------
FROZEN_EXPECTED = {
    "REQUIRED_FREEZE_SHA256": "b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd",
    "REPLICAS": (1, 2, 3),
    "CONTEXTS": {"A": "apo", "P": "probe_only", "C": "compound110__candidate_no_probe", "CP": "compound110__candidate_probe"},
    "CONTRASTS": {
        "Delta_PAM": {"CP": 1.0, "P": -1.0},
        "Delta_AGO": {"C": 1.0, "A": -1.0},
        "Delta_INT": {"CP": 1.0, "P": -1.0, "C": -1.0, "A": 1.0},
    },
    "HISTORICAL_AXIS": {
        "Delta_PAM": "conditional_pam_effect",
        "Delta_AGO": "intrinsic_agonism",
        "Delta_INT": "synergy_interaction",
    },
    "ALPHA": 0.65,
    "DIFFUSION_STEPS": 20,
    "N_RESIDUES": 270,
    "RFF_FEATURES": 512,
    "BRANCHES": {
        "STATE_MOTION": {"width": 512, "bandwidth": 32.30188361260893, "seed": 272340},
        "SIGNED_DRIFT": {"width": 256, "bandwidth": 29.29934899925964, "seed": 272084},
    },
    "P2_N_FRAMES": 1000,
    "P2_BLOCK_FRAMES": 20,
    "P2_N_BLOCKS": 50,
    "P2_RFF_FEATURES": 512,
    "P2_BASE_SEED": 271828,
    "P2_BRANCHES": {"STATE_MOTION": {"width": 512, "seed": 272340}, "SIGNED_DRIFT": {"width": 256, "seed": 272084}},
}


def frozen_numerics_report() -> dict[str, Any]:
    """Return the live imported values alongside the expected frozen values."""
    live = {
        "REQUIRED_FREEZE_SHA256": REQUIRED_FREEZE_SHA256,
        "REPLICAS": tuple(REPLICAS),
        "CONTEXTS": dict(CONTEXTS),
        "CONTRASTS": {k: dict(v) for k, v in CONTRASTS.items()},
        "HISTORICAL_AXIS": dict(HISTORICAL_AXIS),
        "ALPHA": float(ALPHA),
        "DIFFUSION_STEPS": int(DIFFUSION_STEPS),
        "N_RESIDUES": int(N_RESIDUES),
        "RFF_FEATURES": int(RFF_FEATURES),
        "BRANCHES": {k: dict(v) for k, v in BRANCHES.items()},
        "P2_N_FRAMES": int(P2_N_FRAMES),
        "P2_BLOCK_FRAMES": int(P2_BLOCK_FRAMES),
        "P2_N_BLOCKS": int(P2_N_BLOCKS),
        "P2_RFF_FEATURES": int(P2_RFF_FEATURES),
        "P2_BASE_SEED": int(P2_BASE_SEED),
        "P2_BRANCHES": {k: dict(v) for k, v in P2_BRANCHES.items()},
    }
    mismatches = [key for key, expected in FROZEN_EXPECTED.items() if live[key] != expected]
    return {"live": live, "expected": FROZEN_EXPECTED, "mismatches": mismatches, "ok": not mismatches}


def assert_frozen_numerics() -> None:
    """Fail closed when any imported frozen numeric has drifted."""
    report = frozen_numerics_report()
    if not report["ok"]:
        raise RuntimeError(f"frozen numeric drift detected: {report['mismatches']}")


FROZEN_MODULE_FILES = (
    "project/pacer_fkg_v02/run_phase1_bs256.py",
    "project/pacer_fkg_v02/run_phase2_calibration.py",
    "project/pacer_fkg_v02/run_phase3_fkg.py",
    "project/pacer_fkg_v02/path_guard.py",
)


def frozen_module_ledger(repo_root: Path) -> dict[str, Any]:
    """LF-normalised SHA256 of the frozen modules (EOL-insensitive identity)."""

    def lf_sha256(path: Path) -> str:
        raw = path.read_bytes().replace(b"\r\n", b"\n")
        return hashlib.sha256(raw).hexdigest()

    entries = {}
    for rel in FROZEN_MODULE_FILES:
        p = repo_root / rel
        entries[rel] = {
            "lf_sha256": lf_sha256(p),
            "bytes_on_disk": p.stat().st_size,
            "lf_bytes": len(p.read_bytes().replace(b"\r\n", b"\n")),
        }
    return entries
