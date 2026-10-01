"""CM00734 payload-identity regression for the generic frozen runner.

Correction 3: the acceptance criterion is identity of the FROZEN SCIENTIFIC
PAYLOAD -- the vectors, the result metrics and the representation quantities --
not byte identity of receipt files, whose paths and timestamps necessarily
differ between executions.

This harness is strictly read-only.  It recomputes the Stage-B scientific
payload from the frozen CM00734 Phase-2A RFF caches plus the frozen long-MD
reference panel, and compares it against the frozen
`STAGE_B_MATCHED_20NS_RESULTS_v01.json`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from project.pacer_fkg_v02.frozen_runner import engines
from project.pacer_fkg_v02.frozen_runner import fkg_phase2b as p2b
from project.pacer_fkg_v02.frozen_runner.io_utils import check_npy

FROZEN_STAGE_B_RESULT = "project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01/track_b_fkg_v02/phase2b_graph_region/STAGE_B_MATCHED_20NS_RESULTS_v01.json"
FROZEN_STAGE_B_RESULT_SHA256 = "be663c7d45cc4ade1e053029e6aa09e65ca24c1e2fc22dae5398a26bf49d7d99"

# Declared comparison policy (correction 3).
#
#   * every descriptor / block-summary field must match EXACTLY;
#   * the ONLY tolerated quantity is the scalar R1/R3 direction cosine, which is a
#     single float32 inner product and therefore carries summation-order noise of
#     order 1e-8 when the same numbers are reduced along a different memory path.
#
# The tolerance is 1e-6, i.e. four orders of magnitude below the 6-decimal
# precision at which the frozen reports publish these cosines.  Signs must be
# preserved, and the 6-dp rendering must be identical, or the run fails.
COSINE_FIELD = "R1_R3_direction_cosine"
COSINE_ABS_TOLERANCE = 1e-6
REPORT_DECIMALS = 6
CM_RFF_ROOT = "project/cache/pacer_dc_cm00734_stage_b_20ns_v01/track_b_fkg_v02/phase2a_frozen_apply"
REPLICAS = (1, 2, 3)
CM_SYSTEMS = {"C": "CM00734__candidate_no_probe", "CP": "CM00734__candidate_probe"}
N_FRAMES = 400
BLOCK_FRAMES = 20
N_BLOCKS = 20


def _candidate_rff(repo_root: Path, system: str, replica: int, branch: str) -> np.ndarray:
    path = repo_root / CM_RFF_ROOT / system / f"replica_{replica:02d}" / f"{branch}_RFF.npy"
    try:
        return check_npy(path, (N_BLOCKS, engines.N_RESIDUES, engines.RFF_FEATURES))
    except ValueError as exc:
        raise RuntimeError(f"frozen CM00734 RFF invalid: {path} ({exc})") from exc


class _SpecShim:
    """Minimal stand-in so the generic reference-RFF helper can be reused."""

    n_frames = N_FRAMES
    block_frames = BLOCK_FRAMES
    blocks_per_trajectory = N_BLOCKS


def compute_payload(repo_root: Path) -> dict[str, Any]:
    anchor = engines.cached_anchor()
    graph = engines.load_json(engines.GRAPH_PATH)
    transition = engines.graph_transition(graph)
    regions = engines.region_indices(graph)
    shim = _SpecShim()

    payload: dict[str, Any] = {}
    for branch in ("STATE_MOTION", "SIGNED_DRIFT"):
        payload[branch] = {"compound110_reference": {}, "CM00734": {}}
        per_replica: dict[str, dict[int, dict[str, Any]]] = {"compound110_reference": {}, "CM00734": {}}
        for replica in REPLICAS:
            ref = {}
            for alias, system in (("A", "apo"), ("P", "probe_only"),
                                  ("C", "compound110__candidate_no_probe"),
                                  ("CP", "compound110__candidate_probe")):
                descriptors = p2b.reference_branches(shim, system, replica)
                ref[alias] = p2b.reference_rff(shim, system, replica, branch,
                                               anchor=anchor, descriptors=descriptors)
                del descriptors
            cm = {"A": ref["A"], "P": ref["P"],
                  "C": _candidate_rff(repo_root, CM_SYSTEMS["C"], replica, branch),
                  "CP": _candidate_rff(repo_root, CM_SYSTEMS["CP"], replica, branch)}
            for name, contexts in (("compound110_reference", ref), ("CM00734", cm)):
                region_vectors = {}
                for contrast in engines.CONTRASTS:
                    diffused = engines.diffuse(engines.historical_contrast(contexts, contrast), transition)
                    region_vectors[contrast] = {
                        region: diffused[:, idx, :].mean(axis=1, dtype=np.float64).astype(np.float32)
                        for region, idx in regions.items()}
                per_replica[name][replica] = region_vectors
        for name in ("compound110_reference", "CM00734"):
            panel: dict[str, Any] = {}
            for contrast in engines.CONTRASTS:
                panel[contrast] = {"historical_axis": engines.HISTORICAL_AXIS[contrast], "regions": {}}
                for region in regions:
                    per = {f"R{r}": p2b.summarize(per_replica[name][r][contrast][region]) for r in REPLICAS}
                    r1 = per_replica[name][1][contrast][region].mean(axis=0)
                    r3 = per_replica[name][3][contrast][region].mean(axis=0)
                    panel[contrast]["regions"][region] = {
                        "replicas": per, "R1_R3_direction_cosine": engines.cosine(r1, r3)}
            payload[branch][name] = panel
    return payload


def _compare(frozen: Any, live: Any, path: str, diffs: list[dict[str, Any]],
             stats: dict[str, Any]) -> None:
    if isinstance(frozen, dict) and isinstance(live, dict):
        for key in sorted(set(frozen) | set(live)):
            if key not in frozen or key not in live:
                diffs.append({"path": f"{path}.{key}", "kind": "key_presence",
                              "frozen": key in frozen, "live": key in live})
                continue
            _compare(frozen[key], live[key], f"{path}.{key}", diffs, stats)
        return
    if isinstance(frozen, list) and isinstance(live, list):
        if len(frozen) != len(live):
            diffs.append({"path": path, "kind": "length", "frozen": len(frozen), "live": len(live)})
            return
        for index, (a, b) in enumerate(zip(frozen, live)):
            _compare(a, b, f"{path}[{index}]", diffs, stats)
        return
    if isinstance(frozen, bool) or isinstance(live, bool):
        if frozen != live:
            diffs.append({"path": path, "kind": "value", "frozen": frozen, "live": live})
        return
    if isinstance(frozen, (int, float)) and isinstance(live, (int, float)):
        a, b = float(frozen), float(live)
        if path.endswith(COSINE_FIELD):
            delta = abs(a - b)
            stats["cosine_fields"] += 1
            stats["cosine_max_abs_diff"] = max(stats["cosine_max_abs_diff"], delta)
            if (a > 0) != (b > 0) and a != 0.0 and b != 0.0:
                stats["cosine_sign_mismatches"].append({"path": path, "frozen": a, "live": b})
            if round(a, REPORT_DECIMALS) != round(b, REPORT_DECIMALS):
                stats["cosine_rendering_mismatches"].append(
                    {"path": path, "frozen": round(a, REPORT_DECIMALS), "live": round(b, REPORT_DECIMALS)})
            if delta > COSINE_ABS_TOLERANCE:
                diffs.append({"path": path, "kind": "cosine_beyond_tolerance",
                              "frozen": a, "live": b, "abs_diff": delta,
                              "tolerance": COSINE_ABS_TOLERANCE})
            return
        stats["exact_fields"] += 1
        if a != b:
            diffs.append({"path": path, "kind": "numeric_exact", "frozen": a, "live": b, "abs_diff": abs(a - b)})
        return
    stats["exact_fields"] += 1
    if frozen != live:
        diffs.append({"path": path, "kind": "value", "frozen": frozen, "live": live})


PRIMARY_SIGNAL = ("STATE_MOTION", "CM00734", "Delta_INT", "compound110_extension", COSINE_FIELD)


def _primary_signal(frozen_results: dict[str, Any], live_results: dict[str, Any]) -> dict[str, Any]:
    """The pre-existing Stage-B headline signal, compared explicitly."""
    branch, panel, contrast, region, field = PRIMARY_SIGNAL
    frozen = frozen_results[branch][panel][contrast]["regions"][region][field]
    live = live_results[branch][panel][contrast]["regions"][region][field]
    return {
        "path": f"{branch}/{panel}/{contrast}/{region}/{field}",
        "frozen": frozen,
        "live": live,
        "abs_diff": abs(frozen - live),
        "within_tolerance": abs(frozen - live) <= COSINE_ABS_TOLERANCE,
        "sign_preserved": (frozen > 0) == (live > 0),
        "rendered_6dp_frozen": round(frozen, REPORT_DECIMALS),
        "rendered_6dp_live": round(live, REPORT_DECIMALS),
        "matches": abs(frozen - live) <= COSINE_ABS_TOLERANCE and (frozen > 0) == (live > 0),
    }


def run_regression(repo_root: Path, *, max_reported: int = 40) -> dict[str, Any]:
    result_path = repo_root / FROZEN_STAGE_B_RESULT
    from project.pacer_fkg_v02.frozen_runner.spec import sha256_file

    actual_sha = sha256_file(result_path)
    frozen = json.loads(result_path.read_text(encoding="utf-8-sig"))
    live = compute_payload(repo_root)

    diffs: list[dict[str, Any]] = []
    stats: dict[str, Any] = {"exact_fields": 0, "cosine_fields": 0, "cosine_max_abs_diff": 0.0,
                             "cosine_sign_mismatches": [], "cosine_rendering_mismatches": []}
    _compare(frozen["results"], live, "results", diffs, stats)

    frozen_keys = set(frozen["results"]["STATE_MOTION"]["CM00734"].keys())
    live_keys = set(live["STATE_MOTION"]["CM00734"].keys())
    primary = _primary_signal(frozen["results"], live)
    identical = len(diffs) == 0 and not stats["cosine_sign_mismatches"] and primary["matches"]
    return {
        "schema": "pacer.final.frozen_runner.cm00734_payload_regression.v1",
        "frozen_result_path": str(result_path),
        "frozen_result_sha256_expected": FROZEN_STAGE_B_RESULT_SHA256,
        "frozen_result_sha256_actual": actual_sha,
        "frozen_result_sha256_match": actual_sha == FROZEN_STAGE_B_RESULT_SHA256,
        "compared_subtree": "results[<branch>][<panel>][<contrast>].regions[<region>]",
        "comparison_policy": {
            "exact_fields": "all descriptor / block-summary quantities must match exactly",
            "tolerated_field": COSINE_FIELD,
            "tolerance_abs": COSINE_ABS_TOLERANCE,
            "tolerance_rationale": ("single float32 inner product; summation-order noise of order 1e-8. "
                                    f"The tolerance is 1e-6, four orders below the {REPORT_DECIMALS}-decimal "
                                    "precision at which the frozen reports publish these cosines."),
            "sign_preservation_required": True,
            "rendering_equality_reported_not_gated": (
                f"round(value, {REPORT_DECIMALS}) equality is reported as evidence; it is not a gate, "
                "because a 1-ULP float32 difference can flip the last published digit."),
        },
        "exact_fields_compared": stats["exact_fields"],
        "exact_fields_mismatched": sum(1 for d in diffs if d["kind"] != "cosine_beyond_tolerance"),
        "cosine_fields_compared": stats["cosine_fields"],
        "cosine_max_abs_diff": stats["cosine_max_abs_diff"],
        "cosine_sign_mismatches": stats["cosine_sign_mismatches"],
        "cosine_rendering_mismatches": stats["cosine_rendering_mismatches"],
        "cosine_rendering_identical": not stats["cosine_rendering_mismatches"],
        "cosine_rendering_mismatch_count": len(stats["cosine_rendering_mismatches"]),
        "primary_signal": primary,
        "n_differences": len(diffs),
        "differences": diffs[:max_reported],
        "differences_truncated": max(0, len(diffs) - max_reported),
        "scientific_payload_identical": identical,
        "non_scientific_provenance_fields_not_compared": [
            "created_at", "upstream.phase2a_receipt_sha256", "vector_artifacts[].path",
            "vector_artifacts[].bytes", "vector_artifacts[].sha256", "graph.path",
            "matched_design (run-level labels)", "candidate.implemented_as", "version",
            "schema", "status", "run_id",
        ],
        "frozen_panels": sorted(frozen_keys),
        "live_panels": sorted(live_keys),
    }
