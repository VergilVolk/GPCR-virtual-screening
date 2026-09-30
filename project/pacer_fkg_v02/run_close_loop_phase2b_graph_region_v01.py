#!/usr/bin/env python
"""PACER-DC closure Track-B matched-20ns graph/region analysis.

Two matched four-context panels:
1. compound110 retrospective reference
2. LY2119620 validation

Both use the same historical apo/probe controls, first 20 ns only,
and the frozen PACER-FKG-v02 numerical state, graph and regions.

No fitting, calibration, threshold selection, or biological gate is performed.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from project.pacer_fkg_v02 import run_phase2_calibration as phase2
from project.pacer_fkg_v02.run_phase3_fkg import (
    REQUIRED_FREEZE_SHA256,
    GRAPH_PATH,
    CONTRASTS,
    HISTORICAL_AXIS,
    RFF_FEATURES,
    diffuse,
    graph_transition,
    historical_contrast,
    load_frozen_state,
    normalize_frozen,
    region_indices,
    rff_frozen,
    verify_frozen_anchor,
)

VERSION = "PACER_DC_CLOSE_LOOP_PHASE2B_GRAPH_REGION_v01"

CLOSURE = REPO / "project/results/pacer_dc_close_loop_20ns_v01"

PHASE2A_RECEIPT = (
    CLOSURE
    / "track_b_fkg_v02"
    / "phase2a_frozen_apply"
    / "PHASE2A_FREEZE_RECEIPT_v01.json"
)

EXPECTED_PHASE2A_SHA = (
    "e42d3445b3374a3837a6c12d6b1c5e2"
    "e97b54824fd2747d08902394faad6744f"
)

LY_RFF_ROOT = (
    REPO
    / "project/cache/pacer_dc_close_loop_20ns_v01"
    / "track_b_fkg_v02"
    / "phase2a_frozen_apply"
    / "full"
)

OUTPUT_ROOT = (
    CLOSURE
    / "track_b_fkg_v02"
    / "phase2b_graph_region"
)

RESULT_PATH = OUTPUT_ROOT / "TRACK_B_MATCHED_20NS_RESULTS_v01.json"

HISTORICAL_CONTEXTS = {
    "A": "apo",
    "P": "probe_only",
    "C": "compound110__candidate_no_probe",
    "CP": "compound110__candidate_probe",
}

LY_CONTEXTS = {
    "A": "apo",
    "P": "probe_only",
    "C": "LY2119620__candidate_no_probe",
    "CP": "LY2119620__candidate_probe",
}

REPLICAS = (1, 2, 3)
BRANCHES = ("STATE_MOTION", "SIGNED_DRIFT")

N_FRAMES = 400
N_BLOCKS = 20
BLOCK_FRAMES = 20
N_RESIDUES = 270
BS_WIDTH = 256


def now():
    return datetime.now(timezone.utc).astimezone().isoformat()


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def atomic_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.partial.{os.getpid()}")
    tmp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def atomic_npz(path: Path, arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.partial.{os.getpid()}.npz")
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, path)


def construct_20ns_blocks(frames):
    frames = np.asarray(frames, dtype=np.float32)

    if frames.shape != (400, 270, 256):
        raise ValueError(f"unexpected matched-20ns shape: {frames.shape}")

    blocks = frames.reshape(
        N_BLOCKS, BLOCK_FRAMES, N_RESIDUES, BS_WIDTH
    )

    mu_z = blocks.mean(axis=1, dtype=np.float64).astype(np.float32)

    delta = np.diff(blocks, axis=1)

    signed_drift = delta.mean(
        axis=1, dtype=np.float64
    ).astype(np.float32)

    rms_delta = np.sqrt(
        np.mean(
            np.square(delta, dtype=np.float64),
            axis=1,
        )
    ).astype(np.float32)

    state_motion = np.concatenate(
        (mu_z, rms_delta),
        axis=-1,
    ).astype(np.float32, copy=False)

    endpoint = (
        blocks[:, -1] - blocks[:, 0]
    ) / np.float32(BLOCK_FRAMES - 1)

    if not np.allclose(
        signed_drift,
        endpoint,
        rtol=2e-5,
        atol=2e-6,
    ):
        raise RuntimeError("SIGNED_DRIFT endpoint identity failed")

    return {
        "STATE_MOTION": state_motion,
        "SIGNED_DRIFT": signed_drift,
    }


def historical_rff(anchor, system, replica, branch):
    path = phase2.cache_path(system, replica)

    frames = np.load(
        path,
        mmap_mode="r",
        allow_pickle=False,
    )

    if frames.shape != (1000, 270, 256):
        raise RuntimeError(
            f"historical cache shape mismatch: {system} R{replica}"
        )

    descriptors = construct_20ns_blocks(
        np.asarray(frames[:400])
    )

    state = load_frozen_state(anchor, branch)

    mapped = rff_frozen(
        normalize_frozen(descriptors[branch], state),
        state,
    )

    if mapped.shape != (20, 270, 512):
        raise RuntimeError("historical RFF shape mismatch")

    if not np.isfinite(mapped).all():
        raise RuntimeError("historical RFF non-finite")

    return mapped


def ly_rff(system, replica, branch):
    filename = f"{branch}_RFF.npy"

    path = (
        LY_RFF_ROOT
        / system
        / f"replica_{replica:02d}"
        / filename
    )

    arr = np.load(
        path,
        mmap_mode="r",
        allow_pickle=False,
    )

    if arr.shape != (20, 270, 512):
        raise RuntimeError(
            f"LY RFF shape mismatch: {system} R{replica} {branch}"
        )

    if not np.isfinite(arr).all():
        raise RuntimeError("LY RFF non-finite")

    return np.asarray(arr)


def cosine(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    denom = float(np.linalg.norm(a) * np.linalg.norm(b))

    if denom <= 1e-12:
        return None

    return float(
        np.clip(np.dot(a, b) / denom, -1.0, 1.0)
    )


def summarize(vectors):
    vectors = np.asarray(vectors, dtype=np.float64)

    mean_vector = vectors.mean(axis=0)
    block_magnitude = np.linalg.norm(vectors, axis=1)

    q10, q50, q90 = np.quantile(
        block_magnitude,
        [0.10, 0.50, 0.90],
    )

    return {
        "mean_vector_norm": float(np.linalg.norm(mean_vector)),
        "block_magnitude_mean": float(block_magnitude.mean()),
        "block_magnitude_q10": float(q10),
        "block_magnitude_median": float(q50),
        "block_magnitude_q90": float(q90),
    }


def evaluate_panel(
    panel_name,
    branch,
    replica,
    context_nodes,
    transition,
    regions,
):
    context_region = {}

    for alias, values in context_nodes.items():
        diffused = diffuse(values, transition)

        context_region[alias] = {
            region: diffused[:, indices, :].mean(
                axis=1,
                dtype=np.float64,
            ).astype(np.float32)
            for region, indices in regions.items()
        }

    contrast_region = {}

    for contrast in CONTRASTS:
        contrast_nodes = historical_contrast(
            context_nodes,
            contrast,
        )

        diffused = diffuse(
            contrast_nodes,
            transition,
        )

        contrast_region[contrast] = {
            region: diffused[:, indices, :].mean(
                axis=1,
                dtype=np.float64,
            ).astype(np.float32)
            for region, indices in regions.items()
        }

    arrays = {}

    for alias, region_map in context_region.items():
        for region, vectors in region_map.items():
            arrays[f"context__{alias}__{region}"] = vectors

    for contrast, region_map in contrast_region.items():
        for region, vectors in region_map.items():
            arrays[f"contrast__{contrast}__{region}"] = vectors

    out = (
        OUTPUT_ROOT
        / branch.lower()
        / panel_name
        / f"replica_{replica:02d}"
        / "BLOCK_REGIONAL_VECTORS.npz"
    )

    atomic_npz(out, arrays)

    return contrast_region, {
        "path": out.relative_to(REPO).as_posix(),
        "sha256": sha256(out),
        "bytes": out.stat().st_size,
    }


def main():
    if RESULT_PATH.exists():
        raise FileExistsError(
            f"refusing overwrite: {RESULT_PATH}"
        )

    if sha256(PHASE2A_RECEIPT) != EXPECTED_PHASE2A_SHA:
        raise RuntimeError("Phase-2A receipt SHA mismatch")

    receipt = load_json(PHASE2A_RECEIPT)

    if receipt.get("status") != "PHASE2A_FROZEN":
        raise RuntimeError("Phase-2A is not frozen")

    # Single mandatory historical integrity gate.
    anchor = verify_frozen_anchor()

    if REQUIRED_FREEZE_SHA256 != (
        "b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd"
    ):
        raise RuntimeError("historical freeze SHA changed")

    graph = load_json(GRAPH_PATH)
    transition = graph_transition(graph)
    regions = region_indices(graph)

    all_results = {}
    vectors_for_comparison = {}
    vector_artifacts = []

    for branch in BRANCHES:
        print(f"=== {branch} ===", flush=True)

        all_results[branch] = {
            "compound110_reference": {},
            "LY2119620": {},
        }

        vectors_for_comparison[branch] = {
            "compound110_reference": {},
            "LY2119620": {},
        }

        for replica in REPLICAS:
            print(f"R{replica}", flush=True)

            hist = {
                alias: historical_rff(
                    anchor,
                    system,
                    replica,
                    branch,
                )
                for alias, system in HISTORICAL_CONTEXTS.items()
            }

            ly = {
                "A": hist["A"],
                "P": hist["P"],
                "C": ly_rff(
                    LY_CONTEXTS["C"],
                    replica,
                    branch,
                ),
                "CP": ly_rff(
                    LY_CONTEXTS["CP"],
                    replica,
                    branch,
                ),
            }

            for panel_name, contexts in (
                ("compound110_reference", hist),
                ("LY2119620", ly),
            ):
                contrasts, artifact = evaluate_panel(
                    panel_name,
                    branch,
                    replica,
                    contexts,
                    transition,
                    regions,
                )

                vector_artifacts.append(
                    {
                        "branch": branch,
                        "panel": panel_name,
                        "replica": replica,
                        **artifact,
                    }
                )

                vectors_for_comparison[branch][panel_name][replica] = contrasts

        for panel_name in (
            "compound110_reference",
            "LY2119620",
        ):
            panel_result = {}

            for contrast in CONTRASTS:
                panel_result[contrast] = {
                    "historical_axis":
                        HISTORICAL_AXIS[contrast],
                    "regions": {},
                }

                for region in regions:
                    reps = {}

                    for replica in REPLICAS:
                        v = vectors_for_comparison[
                            branch
                        ][panel_name][replica][contrast][region]

                        reps[f"R{replica}"] = summarize(v)

                    r1 = vectors_for_comparison[
                        branch
                    ][panel_name][1][contrast][region].mean(axis=0)

                    r3 = vectors_for_comparison[
                        branch
                    ][panel_name][3][contrast][region].mean(axis=0)

                    panel_result[contrast]["regions"][region] = {
                        "replicas": reps,
                        "R1_R3_direction_cosine": cosine(r1, r3),
                    }

            all_results[branch][panel_name] = panel_result

    comparison = {}

    for branch in BRANCHES:
        comparison[branch] = {}

        for contrast in CONTRASTS:
            comparison[branch][contrast] = {}

            for region in regions:
                replicas = {}

                for replica in REPLICAS:
                    ref = vectors_for_comparison[
                        branch
                    ]["compound110_reference"][replica][contrast][region].mean(axis=0)

                    ly = vectors_for_comparison[
                        branch
                    ]["LY2119620"][replica][contrast][region].mean(axis=0)

                    ref_norm = float(np.linalg.norm(ref))
                    ly_norm = float(np.linalg.norm(ly))

                    replicas[f"R{replica}"] = {
                        "direction_cosine_LY_vs_compound110":
                            cosine(ly, ref),

                        "LY_mean_vector_norm":
                            ly_norm,

                        "compound110_mean_vector_norm":
                            ref_norm,

                        "magnitude_ratio_LY_over_compound110":
                            (
                                ly_norm / ref_norm
                                if ref_norm > 1e-12
                                else None
                            ),
                    }

                comparison[branch][contrast][region] = {
                    "replicas": replicas
                }

    report = {
        "schema":
            "pacer_dc.close_loop_20ns.track_b."
            "matched_graph_region.v1",

        "version": VERSION,
        "created_at": now(),
        "status": "TRACK_B_MATCHED_20NS_COMPLETE",

        "upstream": {
            "phase2a_freeze_receipt_sha256":
                EXPECTED_PHASE2A_SHA,

            "historical_fkg_v02_freeze_sha256":
                REQUIRED_FREEZE_SHA256,
        },

        "matched_design": {
            "frames_per_trajectory": 400,
            "frame_spacing_ps": 50,
            "duration_ns": 20,
            "block_frames": 20,
            "blocks_per_trajectory": 20,

            "shared_controls": {
                "A": "apo first 20 ns",
                "P": "probe_only first 20 ns",
            },

            "compound110_reference": {
                "C":
                    "compound110__candidate_no_probe first 20 ns",
                "CP":
                    "compound110__candidate_probe first 20 ns",
            },

            "LY2119620": {
                "C":
                    "LY2119620__candidate_no_probe 20 ns",
                "CP":
                    "LY2119620__candidate_probe 20 ns",
            },
        },

        "contrasts": {
            name: {
                "historical_axis": HISTORICAL_AXIS[name],
                "formula": formula,
            }
            for name, formula in CONTRASTS.items()
        },

        "graph": {
            "path": GRAPH_PATH.relative_to(REPO).as_posix(),
            "sha256": sha256(GRAPH_PATH),
            "diffusion_alpha": 0.65,
            "diffusion_steps": 20,
            "regions": {
                name: int(len(indices))
                for name, indices in regions.items()
            },
        },

        "results": all_results,

        "LY_vs_compound110_descriptive_comparison":
            comparison,

        "vector_artifacts": vector_artifacts,

        "claim_boundary": (
            "Retrospective matched-20ns descriptive mechanism comparison. "
            "No threshold tuning, p-values, classifier fitting, efficacy "
            "prediction, or independent-block inference."
        ),

        "forbidden_operations_performed": {
            "calibration": False,
            "normalization_refit": False,
            "bandwidth_refit": False,
            "rff_refit": False,
            "graph_refit": False,
            "region_refit": False,
            "threshold_selection": False,
            "outcome_driven_tuning": False,
        },
    }

    atomic_json(RESULT_PATH, report)

    print(
        json.dumps(
            {
                "status": report["status"],
                "branches": list(BRANCHES),
                "panels": [
                    "compound110_reference",
                    "LY2119620",
                ],
                "replicas": list(REPLICAS),
                "blocks_per_trajectory": 20,
                "vector_artifacts":
                    len(vector_artifacts),
                "output": str(RESULT_PATH),
                "output_sha256":
                    sha256(RESULT_PATH),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()