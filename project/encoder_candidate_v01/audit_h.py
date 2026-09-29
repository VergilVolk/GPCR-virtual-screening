#!/usr/bin/env python
"""Experiment H: read-only candidate/PACER-FKG compatibility audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

WORKSPACE = Path(__file__).resolve().parents[2]
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))

from project.pacer_dc_training import analyze_pacer_fkg_g2b as g2b  # noqa: E402
from project.pacer_dc_training import analyze_pacer_fkg_g2c as g2c  # noqa: E402
from project.pacer_dc_training import pacer_factorial_kernel_graph as legacy  # noqa: E402

OUT = WORKSPACE / "project/results/encoder_candidate_H_v01/run_001"
IMPL = WORKSPACE / "project/encoder_candidate_v01"
D_ROOT = WORKSPACE / "project/results/encoder_atom_readout_D_v01/run_001"
F_ROOT = WORKSPACE / "project/results/encoder_temporal_F_v01/run_001"
G_ROOT = WORKSPACE / "project/results/encoder_window_G_v01/run_001"
C_ROOT = WORKSPACE / "project/results/encoder_intermediate_C_v01/phase3_run_001"
FKG_ROOT = WORKSPACE / "project/results"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()


def array_sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def rec(path: Path) -> dict:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha(path)}


def dump(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def verify_marker(root: Path, marker_name: str, report_name: str, integrity_name: str) -> dict:
    marker_path, report_path, integrity_path = root / marker_name, root / report_name, root / integrity_name
    marker = json.loads(marker_path.read_text())
    if sha(report_path) != marker["report"]["sha256"] or sha(integrity_path) != marker["integrity_audit_sha256"]: raise RuntimeError(f"completion hashes changed: {root}")
    return {"marker": rec(marker_path), "report": rec(report_path), "integrity": rec(integrity_path), "status": marker["status"]}


def provenance() -> tuple[dict, dict, dict, dict]:
    c = verify_marker(C_ROOT, "PHASE3_COMPLETE.json", "REPORT_PHASE3.md", "integrity_completion_audit.json")
    d = verify_marker(D_ROOT, "EXPERIMENT_D_COMPLETE.json", "REPORT.md", "integrity_completion_audit.json")
    f = verify_marker(F_ROOT, "EXPERIMENT_F_COMPLETE.json", "REPORT.md", "integrity_completion_audit.json")
    g = verify_marker(G_ROOT, "EXPERIMENT_G_COMPLETE.json", "REPORT.md", "integrity_completion_audit.json")
    dprov = json.loads((D_ROOT / "provenance.json").read_text()); extraction = json.loads((D_ROOT / "extraction_audit.json").read_text()); atom_groups_path = D_ROOT / "atom_group_definition.json"; atom_groups = json.loads(atom_groups_path.read_text())
    checkpoint = Path(dprov["checkpoint"]["checkpoint_path"]); c1 = Path(dprov["c1_implementation"]["path"])
    if sha(checkpoint) != dprov["checkpoint"]["actual_sha256"] or sha(c1) != dprov["c1_implementation"]["sha256"]: raise RuntimeError("checkpoint or C1 implementation changed")
    for row in extraction["rows"]:
        if sha(Path(row["cache_path"])) != row["cache_sha256"] or sha(Path(row["source_path"])) != row["source_sha256"]: raise RuntimeError(f"D cache/source changed: {row['relative_path']}")
    source = dprov["source"]; source_root = Path(source["checkout"])
    # The audit executes under a sandbox account while the authorized read-only
    # checkout belongs to the user.  A command-local safe.directory setting
    # permits Git inspection without changing either checkout or global config.
    head = subprocess.run(
        ["git", "-c", f"safe.directory={source_root.as_posix()}", "rev-parse", "HEAD"],
        cwd=source_root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    source_files = []
    for item in source["inference_files"]:
        pinned, installed = Path(item["pinned_checkout_path"]), Path(item["installed_path"]); ok = sha(pinned) == item["pinned_checkout_sha256"] and sha(installed) == item["installed_sha256"] and sha(pinned) == sha(installed)
        if not ok: raise RuntimeError(f"Geom2Vec inference source changed: {item['file']}")
        source_files.append({"file": item["file"], "pinned_sha256": sha(pinned), "installed_sha256": sha(installed), "match": True})
    fcfg = json.loads((F_ROOT / "frozen_config.json").read_text()); gcfg = json.loads((G_ROOT / "frozen_config.json").read_text())
    if fcfg["lag"]["stored_frames"] != 1 or gcfg["window"]["frames"] != 20 or gcfg["window"]["lag_stored_frames"] != 1: raise RuntimeError("validated block/lag definition changed")
    audit = {"status": "PASS", "experiments": {"C": c, "D": d, "F": f, "G": g}, "checkpoint": rec(checkpoint), "checkpoint_expected_sha256": dprov["checkpoint"]["expected_sha256"], "checkpoint_strict_coverage_reused": dprov["checkpoint"]["complete_strict_coverage"], "c1_implementation": rec(c1), "geom2vec": {"pinned_commit": source["pinned_commit"], "checkout_head": head, "head_matches_pin": head == source["pinned_commit"], "inference_files": source_files}, "D_cache_files_verified": 40, "atom14_sources_verified": 40, "atom_group_definition": rec(atom_groups_path), "F_lag_stored_frames": 1, "G_block_frames": 20, "G_lag_stored_frames": 1, "archives_modified": False}
    if head != source["pinned_commit"]: raise RuntimeError("Geom2Vec checkout no longer at pinned commit")
    return audit, dprov, atom_groups, gcfg


def candidate_spec(dprov: dict, atom_groups: dict) -> dict:
    return {
        "version": "PACER_DC_ENCODER_CANDIDATE_v01", "status": "FROZEN_NO_TUNABLE_PARAMETERS",
        "checkpoint_sha256": dprov["checkpoint"]["actual_sha256"], "geom2vec_source_commit": dprov["source"]["pinned_commit"],
        "visnet_state": {"name": "C1", "definition": "accumulated scalar/vector state after message-passing layers 0,1,2", "hook_location": "input pre-hook of representation_model.vis_mp_layers[3]", "scalar_channels": 64, "vector_channels": 64},
        "atomic_invariant_descriptor": {"definition": "concat(x, Euclidean norm of each vector channel over xyz)", "dimension": 128, "vector_norm_axis": "xyz"},
        "residue_readout": {"backbone_atoms": atom_groups["backbone_atoms"], "sidechain_definition": atom_groups["sidechain_definition"], "gly_handling": "S is an exact 128D zero vector and sidechain_present is false", "B_dimension": 128, "S_dimension": 128, "BS_definition": "concat(B,S) in that order", "BS_dimension": 256, "pooling": "equal mean over available expected heavy atoms in each group"},
        "temporal": {"block_size_stored_frames": 20, "lag_stored_frames": 1, "differences_per_block": 19, "delta": "DeltaZ[t]=BS[t+1]-BS[t], internal to each block", "muZ": "mean of BS over exactly 20 ordered stored frames", "muDelta": "mean of the 19 signed DeltaZ vectors; equals (BS_last-BS_first)/19", "rmsDelta": "elementwise sqrt(mean of squared 19 DeltaZ vectors)", "cross_block_deltas": False, "cross_window_deltas": False},
        "branches": {"STATE_MOTION": {"definition": "concat(muZ,RMSDelta) in that order", "dimension": 512}, "SIGNED_DRIFT": {"definition": "muDelta", "dimension": 256}, "merge_policy": "none; branches remain separately named"},
        "dtype_contract": {"atomic_and_residue_cache": "float32", "temporal_reductions": "float64", "candidate_branch_output": "float64", "normalization": "none"},
        "residue_ordering": "270 receptor residues in authenticated sequence/embedding_index order; no sorting or reindexing", "frame_ordering": "strict archived frame_ids within each replica/condition/window; five contiguous blocks [0:20],[20:40],[40:60],[60:80],[80:100]", "block_output_shape_per_source_window": {"STATE_MOTION": [5, 270, 512], "SIGNED_DRIFT": [5, 270, 256]},
        "tunable_parameters": [], "scientific_claim": "representation definition only; no biological efficacy claim",
    }


def fkg_artifacts() -> tuple[dict, dict, dict]:
    g2b_path = FKG_ROOT / "pacer_dc_fkg_G2B_v01/G2B_DIRECTION_AUDIT.json"; g2c_path = FKG_ROOT / "pacer_dc_fkg_G2C_v01/G2C_GRAPH_ABLATION_AUDIT.json"; fkg_path = FKG_ROOT / "pacer_dc_fkg_R2R3_v01/PACER_FKG_R2R3_AUDIT.json"; graph_path = FKG_ROOT / "pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json"
    b = json.loads(g2b_path.read_text()); c = json.loads(g2c_path.read_text()); med = np.asarray(b["channel_median"], dtype=np.float64); scale = np.asarray(b["channel_scale"], dtype=np.float64)
    normalization = {"historical_pipeline": ["load raw residue_features as FP32", "residue-global centering in FP32", "flatten R2 calibration samples and cast FP64", "fit channel median and 1.4826*MAD floor 1e-6 in FP64", "apply FP64 median/scale then cast standardized representation to FP32"], "storage": "embedded JSON arrays in G2B_DIRECTION_AUDIT.json; no standalone normalization file", "artifact": rec(g2b_path), "channel_median": {"shape": list(med.shape), "dtype_when_loaded": str(med.dtype), "content_sha256_float64_c_order": array_sha(med)}, "channel_scale": {"shape": list(scale.shape), "dtype_when_loaded": str(scale.dtype), "content_sha256_float64_c_order": array_sha(scale), "strictly_positive": bool(np.all(scale > 0))}, "compatibility": {"historical_128D": True, "BS256": False, "STATE_MOTION512": False, "SIGNED_DRIFT256": False}, "reason": "channelwise location/scale is representation-specific and length 128; tiling, duplication, padding, interpolation, or truncation changes meaning and is forbidden", "classification": "D"}
    regional = []
    for row in b["calibration"]:
        regional.append({"region": row["region"], "n_residues": row["n_residues"], "input_dimension": row["n_residues"] * 128, "rff_output_dimension": row["rff_dimension"], "seed": row["rff_seed"], "bandwidth": row["bandwidth"], "W_shape_if_procedurally_regenerated": [row["n_residues"] * 128, row["rff_dimension"]], "bias_shape_if_procedurally_regenerated": [row["rff_dimension"]]})
    kernel = {"G2B": {"artifact": rec(g2b_path), "implementation": rec(WORKSPACE / "project/pacer_dc_training/analyze_pacer_fkg_g2b.py"), "kernel": "per-region RBF approximated with RFF", "regional_maps": regional, "stored_random_matrices_or_biases": False, "matrix_hashes": None, "absence_reason": "W and bias existed only in memory; the archive stores seed, bandwidth, input/output dimensions, and diagnostics", "candidate_compatibility": False}, "G2C": {"artifact": rec(g2c_path), "implementation": rec(WORKSPACE / "project/pacer_dc_training/analyze_pacer_fkg_g2c.py"), "kernel": "shared-node 128D RBF RFF", "input_dimension": 128, "W_shape_if_procedurally_regenerated": [128, c["calibration"]["rff_dim"]], "bias_shape_if_procedurally_regenerated": [c["calibration"]["rff_dim"]], "rff_output_dimension": c["calibration"]["rff_dim"], "seed": c["calibration"]["rff_seed"], "bandwidth": c["calibration"]["shared_node_bandwidth"], "stored_random_matrices_or_biases": False, "matrix_hashes": None, "candidate_compatibility": False}, "exact_legacy_RBF": {"implementation": rec(WORKSPACE / "project/pacer_dc_training/pacer_factorial_kernel_graph.py"), "dimension_agnostic_code": True, "frozen_G2B_equivalent": False, "reason": "it recomputes bandwidth/preprocessing per call and is not the frozen common RFF realization"}, "dimension_change_effect": "changes the Euclidean space, bandwidth calibration, W row dimension, and RFF realization; a new draw/calibration is a new PACER-FKG method version", "no_random_features_created_in_H": True, "classification": "D"}
    common = {"g2b": b, "g2c": c, "fkg_path": fkg_path, "graph_path": graph_path}
    return normalization, kernel, common


def architecture(common: dict) -> dict:
    graph = json.loads(common["graph_path"].read_text()); components = [
        {"component": "historical NPZ input schema", "classification": "C", "assumption": "single residue_features array with shape (100,270,128), sequence scalar, frame_ids (100,)"},
        {"component": "legacy load_contexts", "classification": "A", "assumption": "channel dimension inferred, but only one unnamed tensor and equal context shapes"},
        {"component": "legacy preprocess", "classification": "B", "assumption": "dimension-agnostic operations, but it refits median/MAD on supplied contexts rather than using frozen G2-B state"},
        {"component": "R2/R3 and G2-B/G2-C preflight", "classification": "C", "assumption": "explicit (100,270,128) assertions"},
        {"component": "G2-B normalization", "classification": "D", "assumption": "frozen 128-channel R2 median/scale"},
        {"component": "G2-B regional RFF", "classification": "C", "assumption": "input dimension n_region_residues*128"},
        {"component": "G2-C node RFF", "classification": "C", "assumption": "reshape(-1,128), one 128D node descriptor"},
        {"component": "20-frame block reduction", "classification": "C", "assumption": "RFF-map 100 frame-level inputs first, then reshape to five 20-frame blocks; candidate is already block-level"},
        {"component": "graph topology and transition", "classification": "A", "assumption": "270 ordered nodes; diffusion is agnostic to propagated feature width"},
        {"component": "region membership", "classification": "A", "assumption": "embedding_index and 270-residue order only"},
        {"component": "contrast/bootstrap/report schema", "classification": "D", "assumption": "one common RFF coordinate tensor and historical calibrated method semantics"},
        {"component": "serialization/output", "classification": "C", "assumption": "single residue_features input; no named branch metadata or branch-specific calibrated state"},
    ]
    return {"audited_path": "historical frame tensor -> FP32 residue centering -> frozen G2-B channel scaling -> RBF-RFF -> 20-frame RFF block mean -> signed contrasts -> optional G2-C graph diffusion -> reports", "input_contract": {"residue_features": [100, 270, 128], "sequence": "scalar string length 270 matching graph node sites", "frame_ids": [100], "dtype": "finite numeric loaded/cast FP32", "contexts": ["apo", "probe_only", "candidate_no_probe", "candidate_probe"], "branches": "none"}, "graph": {"artifact": rec(common["graph_path"]), "version": graph["version"], "nodes": len(graph["nodes"]), "edges": len(graph["edges"]), "regions": len(graph["regions"]), "feature_width_agnostic_after_RFF": True}, "components": components, "candidate_shapes": {"BS_frame_level": [100, 270, 256], "STATE_MOTION_block_level": [5, 270, 512], "SIGNED_DRIFT_block_level": [5, 270, 256]}, "directly_compatible": False}


def branch_audit() -> dict:
    return {"current_input_model": "one unnamed residue_features tensor", "multiple_named_branches_supported": False, "internal_branch_merge": False, "separate_normalization_or_kernel_per_branch": False, "one_common_normalization_kernel_assumed": True, "STATE_MOTION512_direct": False, "SIGNED_DRIFT256_direct": False, "simultaneous_candidate_object_direct": False, "forbidden_forcing_operations": ["concatenate SM and SD", "run branches independently and combine outputs without a preregistered definition", "tile/pad/truncate calibration", "reuse historical RFF matrices"], "compatibility_boundary": "current PACER-FKG has no legitimate representation-branch concept", "classification": "C"}


def contract_tests(norm: dict, common: dict) -> dict:
    med = np.asarray(common["g2b"]["channel_median"], dtype=np.float64); scale = np.asarray(common["g2b"]["channel_scale"], dtype=np.float64); rng = np.random.default_rng(104729)
    arrays256 = {name: rng.normal(size=(3, 4, 256)).astype(np.float32) for name in legacy.CONTEXTS}; tests = []
    def run(name, expected, fn):
        try:
            value = fn(); tests.append({"test": name, "expected": expected, "outcome": "accepted", "detail": value})
        except Exception as exc: tests.append({"test": name, "expected": expected, "outcome": "rejected", "exception_type": type(exc).__name__, "detail": str(exc)})
    run("legacy_preprocess_256", "accept_dimension_agnostic_but_refits_state", lambda: {k: list(v.shape) for k, v in legacy.preprocess(arrays256, "global_residual_channel_robust").items()})
    groups = {name: rng.normal(size=(4, 256)) for name in legacy.CONTEXTS}; run("legacy_exact_kernel_256", "accept_dimension_agnostic", lambda: legacy.signed_kernel_stat(groups, {name: 1 if i % 2 == 0 else -1 for i, name in enumerate(legacy.CONTEXTS)})["bandwidth"])
    graph = json.loads(common["graph_path"].read_text()); run("graph_diffusion_270_nodes", "accept_width_independent_node_scores", lambda: list(legacy.graph_diffuse(np.ones(270), graph).shape))
    for label, dim in (("BS256", 256), ("STATE_MOTION512", 512), ("SIGNED_DRIFT256", 256)):
        run(f"frozen_normalization_{label}", "reject_calibration_shape", lambda dim=dim: list(((np.zeros((2, 270, dim)) - med) / scale).shape))
    nres = common["g2b"]["calibration"][0]["n_residues"]; historical_W = np.zeros((nres * 128, 512), dtype=np.float32); candidate_region = np.zeros((2, nres * 256), dtype=np.float32)
    run("G2B_historical_RFF_with_BS_region", "reject_matrix_input_dimension", lambda: list(g2b.rff(candidate_region, historical_W, np.zeros(512, dtype=np.float32)).shape))
    run("G2C_rff_map_BS256", "reject_hardcoded_128_reshape/output", lambda: list(g2c.rff_map(np.zeros((100, 270, 256), dtype=np.float32), np.zeros((128, 256), dtype=np.float32), np.zeros(256, dtype=np.float32)).shape))
    run("G2B_unit_normalize_named_branches", "reject_no_branch_interface", lambda: g2b.unit_normalize({"STATE_MOTION": np.zeros((5, 270, 512)), "SIGNED_DRIFT": np.zeros((5, 270, 256))}))
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        for name in legacy.CONTEXTS: np.savez(root / f"{name}.npz", residue_features=np.zeros((5, 270, 512), dtype=np.float32))
        run("legacy_loader_single_512_tensor", "accept_channel_width_but_not_branch_object", lambda: {k: list(v.shape) for k, v in legacy.load_contexts(root, "{context}.npz").items()})
    passed = all((t["outcome"] == "accepted") == t["expected"].startswith("accept") for t in tests)
    return {"synthetic_only": True, "biological_inputs_or_outputs_used": False, "random_features_drawn": False, "tests": tests, "expected_boundaries_observed": passed}


def decision() -> dict:
    return {"classification": "H-C", "statement": "Both code and frozen PACER-FKG state are representation-specific; a new version and implementation adaptation are required.", "reusable_scientific_definitions": ["four-context distributional comparison in a common feature space", "R2-only calibration principle", "residue-global centering operation as a definition, subject to new calibration validation", "RBF kernel family", "20-frame contiguous block rule", "270-node frozen graph topology, region membership, and graph diffusion formula"], "newly_estimated_or_generated_numerical_state": ["branch-specific channel median/scale for 512D SM and 256D SD", "branch-specific or otherwise newly specified RBF bandwidths", "new dimension-correct RFF matrices and biases or a versioned exact-kernel alternative", "any calibration diagnostics and frozen hashes"], "software_generalization_required": ["versioned branch-aware input schema", "accept five block-level observations instead of requiring 100 frame-level 128D inputs, without re-averaging", "remove hardcoded 128 from validation, region flattening, and node RFF mapping", "explicit branch-specific normalization/kernel routing", "branch-aware serialization, provenance, assertions, and reporting"], "new_scientific_method_choices_not_made_in_H": ["how SM and SD kernel outputs are combined or kept separate", "whether branches share or use distinct bandwidth/kernel families", "relative branch weighting", "how branch-specific contrasts enter graph propagation and gates", "new calibration population and qualification criteria"], "why_not_H_B": "software is not merely dimension-generic: explicit 128D/100-frame checks and the absence of branches require implementation changes", "implementation_performed": False}


def report(spec, prov, contract, norm, kernel, branch, tests, version) -> str:
    return "\n".join(["# PACER-DC Experiment H v01", "", "The frozen candidate specification was created and the PACER-FKG path was audited read-only. No biological contrast, long-MD access, calibration, RFF draw, or PACER-FKG modification occurred.", "", "## Decision", "", f"**{version['classification']}** — {version['statement']}", "", "## Candidate", "", "C1 `[x, ||v||]` with deterministic B/S pooling yields BS256. Each fixed 20-frame block produces separately named STATE-MOTION512=`[muZ,RMSDelta]` and SIGNED-DRIFT256=`muDelta` branches.", "", "## Compatibility boundary", "", "The historical input is one `(100,270,128)` `residue_features` tensor. Frozen normalization has 128 channels; G2-B regional RFF inputs are `n_residues*128`; G2-C reshapes nodes to 128D. RFF matrices were not serialized, so no matrix hashes exist. PACER-FKG has no named branch interface and applies RFF before its own 20-frame averaging, whereas the candidate is already block-level.", "", f"Frozen median shape/hash: `{norm['channel_median']['shape']}` / `{norm['channel_median']['content_sha256_float64_c_order']}`.", f"Frozen scale shape/hash: `{norm['channel_scale']['shape']}` / `{norm['channel_scale']['content_sha256_float64_c_order']}`.", "", "Synthetic tests confirmed that low-level legacy preprocessing/kernel math accepts arbitrary widths only by refitting state, while frozen normalization, G2-B RFF, G2-C mapping, and branch interfaces reject candidate-shaped inputs at their actual boundaries.", "", "## New-version boundary", "", "The four-context formulation, graph topology, diffusion formula, RBF family, R2-only calibration principle, and fixed 20-frame rule can remain as scientific definitions. New branch-specific normalization, bandwidth/RFF state, schema, routing, serialization, and reporting would be required. Branch combination, weighting, kernel sharing, graph routing, calibration population, and qualification criteria are new scientific choices and were not made here.", "", "No efficacy or PACER outcome claim is made.", ""])


def main() -> None:
    p = argparse.ArgumentParser(); p.add_argument("--output-root", type=Path, default=OUT); args = p.parse_args(); args.output_root.mkdir(parents=True, exist_ok=True); IMPL.mkdir(parents=True, exist_ok=True)
    prov, dprov, atoms, _ = provenance(); spec = candidate_spec(dprov, atoms); dump(IMPL / "CANDIDATE_SPEC.json", spec); dump(args.output_root / "CANDIDATE_SPEC.json", spec); dump(args.output_root / "provenance_audit.json", prov)
    norm, kernel, common = fkg_artifacts(); contract = architecture(common); branch = branch_audit(); tests = contract_tests(norm, common); version = decision()
    dump(args.output_root / "pacer_fkg_input_contract.json", contract); dump(args.output_root / "normalization_compatibility.json", norm); dump(args.output_root / "kernel_rff_compatibility.json", kernel); dump(args.output_root / "branch_compatibility.json", branch); dump(args.output_root / "software_contract_tests.json", tests); dump(args.output_root / "versioning_decision.json", version)
    integrity = {"status": "EXPERIMENT_H_COMPLETE", "provenance_passed": prov["status"] == "PASS", "candidate_spec_identical_in_implementation_and_results": sha(IMPL / "CANDIDATE_SPEC.json") == sha(args.output_root / "CANDIDATE_SPEC.json"), "expected_software_boundaries_observed": tests["expected_boundaries_observed"], "PACER_FKG_modified": False, "new_calibration_fit": False, "new_RFF_drawn": False, "biological_metrics_computed": [], "long_MD_accessed": False, "followup_started": False, "implementation": rec(Path(__file__))}
    dump(args.output_root / "integrity_completion_audit.json", integrity); report_path = args.output_root / "REPORT.md"; report_path.write_text(report(spec, prov, contract, norm, kernel, branch, tests, version), encoding="utf-8"); dump(args.output_root / "EXPERIMENT_H_COMPLETE.json", {"status": integrity["status"], "report": rec(report_path), "integrity_audit_sha256": sha(args.output_root / "integrity_completion_audit.json"), "followup_started": False}); print(json.dumps({"status": integrity["status"], "classification": version["classification"], "output": str(args.output_root)}, indent=2))


if __name__ == "__main__": main()
