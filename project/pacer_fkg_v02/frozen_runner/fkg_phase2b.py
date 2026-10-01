"""Generic Phase-2B driver: frozen graph diffusion + frozen region pooling +
frozen signed contrasts.

Both temporal branches are evaluated separately and are never combined.  The
primary cross-replica statistic is the R1/R3 direction cosine; the calibration
replica (R2) is still mapped and reported as a descriptive replica result so the
Stage-A/Stage-B semantics are preserved exactly.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from project.pacer_fkg_v02.frozen_runner import engines, paths as pathlib_policy
from project.pacer_fkg_v02.frozen_runner.blocks import construct_blocks
from project.pacer_fkg_v02.frozen_runner.io_utils import check_npy, load_npy_copy
from project.pacer_fkg_v02.frozen_runner.spec import RunSpec, sha256_file
from project.pacer_fkg_v02 import run_phase2_calibration as phase2

VERSION = "PACER_FKG_V02_FROZEN_RUNNER_PHASE2B_v01"


def _now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def _atomic_json(path: Path, payload: Any) -> None:
    tmp = path.with_name(f".{path.name}.partial.{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _atomic_npz(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    tmp = path.with_name(f".{path.name}.partial.{os.getpid()}.npz")
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, path)


def summarize(values: np.ndarray) -> dict[str, float]:
    """Frozen Stage-A/Stage-B block-magnitude summary."""
    x = np.asarray(values, dtype=np.float64)
    mean_vector = x.mean(axis=0)
    magnitude = np.linalg.norm(x, axis=1)
    q10, q50, q90 = np.quantile(magnitude, [0.10, 0.50, 0.90])
    return {"mean_vector_norm": float(np.linalg.norm(mean_vector)),
            "block_magnitude_mean": float(magnitude.mean()),
            "block_magnitude_q10": float(q10),
            "block_magnitude_median": float(q50),
            "block_magnitude_q90": float(q90)}


def candidate_rff(spec: RunSpec, system: str, replica: int, branch: str) -> np.ndarray:
    path = spec.paths.phase2a_cache_root / system / f"replica_{replica:02d}" / f"{branch}_RFF.npy"
    expected = (spec.blocks_per_trajectory, engines.N_RESIDUES, engines.RFF_FEATURES)
    try:
        return check_npy(path, expected)
    except ValueError as exc:
        raise RuntimeError(f"candidate RFF invalid: {path} ({exc})") from exc


def reference_branches(spec: RunSpec, system: str, replica: int) -> dict[str, np.ndarray]:
    """Frozen historical reference descriptors, read and blocked once per replica."""
    cache = phase2.cache_path(system, replica)
    if not cache.is_file():
        raise FileNotFoundError(f"frozen long-MD reference cache missing: {cache}")
    frames = load_npy_copy(cache)  # owned copy: never close a mapping still referenced
    if frames.shape != (phase2.N_FRAMES, phase2.N_RESIDUES, phase2.BS_WIDTH):
        raise RuntimeError(f"frozen reference cache shape mismatch: {cache} {frames.shape}")
    window = np.array(frames[: spec.n_frames], dtype=np.float32, copy=True)
    del frames
    state_motion, signed_drift = construct_blocks(window, block_frames=spec.block_frames,
                                                  n_blocks=spec.blocks_per_trajectory)
    return {"STATE_MOTION": state_motion, "SIGNED_DRIFT": signed_drift}


def reference_rff(spec: RunSpec, system: str, replica: int, branch: str, *,
                  anchor: Mapping[str, Any] | None = None,
                  descriptors: Mapping[str, np.ndarray] | None = None) -> np.ndarray:
    """Apply the frozen normalization / RFF state to the historical reference."""
    values = (descriptors if descriptors is not None else reference_branches(spec, system, replica))[branch]
    state = engines.load_frozen_state(anchor if anchor is not None else engines.cached_anchor(), branch)
    mapped = engines.rff_frozen(engines.normalize_frozen(values, state), state)
    expected = (spec.blocks_per_trajectory, engines.N_RESIDUES, engines.RFF_FEATURES)
    if mapped.shape != expected or not np.isfinite(mapped).all():
        raise RuntimeError("historical reference RFF contract failed")
    return mapped


def _panel(name: str, branch: str, replica: int, contexts: Mapping[str, np.ndarray],
           transition, regions: Mapping[str, np.ndarray], output_root: Path, root: Path,
           allowed: tuple[Path, ...]) -> tuple[dict[str, Any], dict[str, Any]]:
    context_regions: dict[str, Any] = {}
    for alias, value in contexts.items():
        diffused = engines.diffuse(value, transition)
        context_regions[alias] = {region: diffused[:, idx, :].mean(axis=1, dtype=np.float64).astype(np.float32)
                                  for region, idx in regions.items()}
    contrast_regions: dict[str, Any] = {}
    for contrast in engines.CONTRASTS:
        diffused = engines.diffuse(engines.historical_contrast(contexts, contrast), transition)
        contrast_regions[contrast] = {region: diffused[:, idx, :].mean(axis=1, dtype=np.float64).astype(np.float32)
                                      for region, idx in regions.items()}
    arrays = {**{f"context__{a}__{r}": v for a, m in context_regions.items() for r, v in m.items()},
              **{f"contrast__{c}__{r}": v for c, m in contrast_regions.items() for r, v in m.items()}}
    path = pathlib_policy.guard_new_output_noclobber(
        output_root / branch.lower() / name / f"replica_{replica:02d}" / "BLOCK_REGIONAL_VECTORS.npz",
        root, allowed)
    path.parent.mkdir(parents=True, exist_ok=True)
    _atomic_npz(path, arrays)
    return contrast_regions, {"path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}


def panel_contexts(spec: RunSpec, replica: int, branch: str, anchor, candidate: Mapping[str, str],
                   reference: Mapping[str, str]) -> tuple[dict[str, dict[str, np.ndarray]],
                                                         dict[str, dict[str, np.ndarray]]]:
    """Build the reference and candidate context RFFs for one (replica, branch)."""
    alias_systems = {**reference, "A": candidate["A"], "P": candidate["P"]}
    ref_contexts: dict[str, np.ndarray] = {}
    for alias, system in alias_systems.items():
        descriptors = reference_branches(spec, system, replica)
        ref_contexts[alias] = reference_rff(spec, system, replica, branch, anchor=anchor,
                                            descriptors=descriptors)
        del descriptors
    cand_contexts = {alias: candidate_rff(spec, system, replica, branch)
                     for alias, system in candidate.items()}
    return {"reference_panel": ref_contexts, "candidate": cand_contexts}


def run(spec: RunSpec, *, repo_root: Path) -> dict[str, Any]:
    root = repo_root.resolve()
    result_path = pathlib_policy.guard_new_output_noclobber(
        spec.paths.report_root / "phase2b" / "FROZEN_MATCHED_RESULTS_V01.json", root,
        spec.paths.allowed_output_roots)

    anchor = engines.cached_anchor()
    graph = engines.load_json(engines.GRAPH_PATH)
    transition = engines.graph_transition(graph)
    regions = engines.region_indices(graph)

    candidate = spec.contexts
    reference = spec.reference_systems
    results: dict[str, Any] = {}
    vectors: dict[str, Any] = {}
    artifacts: list[dict[str, Any]] = []

    for branch in ("STATE_MOTION", "SIGNED_DRIFT"):
        results[branch] = {"reference_panel": {}, "candidate": {}}
        vectors[branch] = {"reference_panel": {}, "candidate": {}}
        for replica in spec.replicas:
            contexts = panel_contexts(spec, replica, branch, anchor, candidate, reference)
            for name, panel_context in contexts.items():
                contrast_regions, artifact = _panel(name, branch, replica, panel_context, transition,
                                                    regions, spec.paths.report_root / "phase2b", root,
                                                    spec.paths.allowed_output_roots)
                vectors[branch][name][replica] = contrast_regions
                artifacts.append({"branch": branch, "panel": name, "replica": replica, **artifact})
            del contexts
        for name in ("reference_panel", "candidate"):
            panel: dict[str, Any] = {}
            for contrast in engines.CONTRASTS:
                panel[contrast] = {"historical_axis": engines.HISTORICAL_AXIS[contrast], "regions": {}}
                for region in regions:
                    per_replica = {f"R{r}": summarize(vectors[branch][name][r][contrast][region])
                                   for r in spec.replicas}
                    r1 = vectors[branch][name][1][contrast][region].mean(axis=0)
                    r3 = vectors[branch][name][3][contrast][region].mean(axis=0)
                    panel[contrast]["regions"][region] = {
                        "replicas": per_replica,
                        "R1_R3_direction_cosine": engines.cosine(r1, r3),
                    }
            results[branch][name] = panel

    report = {
        "schema": "pacer.final.frozen_runner.matched_results.v1", "version": VERSION,
        "created_at": _now(), "status": "FROZEN_MATCHED_COMPLETE", "run_id": spec.run_id,
        "candidate": {"candidate_id": spec.molecule.candidate_id,
                      "canonical_smiles": spec.molecule.canonical_smiles,
                      "role": spec.molecule.role, "chemotype": spec.molecule.chemotype},
        "upstream": {"historical_fkg_v02_freeze_sha256": engines.REQUIRED_FREEZE_SHA256,
                     "phase2a_receipt_sha256": sha256_file(
                         spec.paths.report_root / "phase2a" / "PHASE2A_FREEZE_RECEIPT_v01.json"),
                     "run_spec_sha256": sha256_file(spec.source_path)},
        "matched_design": {"n_frames": spec.n_frames, "frame_spacing_ps": 50,
                           "duration_ns": spec.n_frames / 20.0, "block_frames": spec.block_frames,
                           "blocks_per_trajectory": spec.blocks_per_trajectory,
                           "replicas": list(spec.replicas),
                           "seeds": {str(k): v for k, v in spec.seeds.items()},
                           "calibration_replica": 2, "application_replicas": [1, 3],
                           "shared_controls": {"A": "frozen long-MD apo", "P": "frozen long-MD probe_only"},
                           "reference_panel": {k: v for k, v in reference.items()}},
        "contrasts": {n: {"historical_axis": engines.HISTORICAL_AXIS[n], "formula": f}
                      for n, f in engines.CONTRASTS.items()},
        "graph": {"path": str(engines.GRAPH_PATH), "sha256": sha256_file(engines.GRAPH_PATH),
                  "diffusion_alpha": engines.ALPHA, "diffusion_steps": engines.DIFFUSION_STEPS,
                  "regions": {n: int(len(i)) for n, i in regions.items()}},
        "results": results,
        "vector_artifacts": artifacts,
        "claim_boundary": ("Prospective representation-space mechanistic evidence only. No threshold, "
                           "no PAM probability, no classifier and no efficacy claim. Temporal blocks are "
                           "correlated time samples and are never treated as biological replicates."),
        "forbidden_operations_performed": {"calibration": False, "normalization_refit": False,
                                           "bandwidth_refit": False, "rff_refit": False, "graph_refit": False,
                                           "region_refit": False, "threshold_selection": False,
                                           "outcome_driven_tuning": False, "branch_combination": False},
    }
    _atomic_json(result_path, report)
    return {"status": report["status"], "output": str(result_path), "output_sha256": sha256_file(result_path),
            "branches": ["STATE_MOTION", "SIGNED_DRIFT"], "vector_artifacts": len(artifacts)}


def verify(spec: RunSpec, *, repo_root: Path) -> dict[str, Any]:
    root = repo_root.resolve()
    engines.cached_anchor()
    result_path = spec.paths.report_root / "phase2b" / "FROZEN_MATCHED_RESULTS_V01.json"
    report = json.loads(result_path.read_text(encoding="utf-8-sig"))
    if report.get("status") != "FROZEN_MATCHED_COMPLETE":
        raise RuntimeError("invalid Phase-2B result status")
    for artifact in report["vector_artifacts"]:
        path = Path(artifact["path"])
        if not path.is_file() or sha256_file(path) != artifact["sha256"]:
            raise RuntimeError(f"vector artifact verification failed: {path}")
    return {"status": "FROZEN_MATCHED_VERIFIED", "output_sha256": sha256_file(result_path),
            "vector_artifacts": len(report["vector_artifacts"])}
