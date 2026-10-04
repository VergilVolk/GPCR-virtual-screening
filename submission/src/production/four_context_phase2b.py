"""Stage4 Phase2b: candidate-specific matched contrasts on frozen graph/regions.

Replicas stay separate. All Stage4 replicas are evaluation-only, including R2.
The output is a descriptive evidence package, without pharmacological labels.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from project.pacer_fkg_v02 import stage4_prospective_common_v01 as c
from project.pacer_fkg_v02 import run_stage4_prospective_phase2a_frozen_apply_v01 as p2

OUTPUT = c.RESULT_ROOT / "phase2b_graph_region"
RESULT = OUTPUT / "STAGE4_PROSPECTIVE_RESULTS_v01.json"
RECEIPT = OUTPUT / "PHASE2B_FREEZE_RECEIPT_v01.json"
# The three formal historical names/formulas are inherited exactly. These two
# additional names describe requested contrasts, with no biological threshold.
ADDITIONAL_CONTRASTS = {"Probe_effect": {"P": 1.0, "A": -1.0},
                        "Probe_background_effect": {"CP": 1.0, "C": -1.0}}


def all_contrasts(engine):
    return {**engine.CONTRASTS, **ADDITIONAL_CONTRASTS}


def contrast(contexts, name, engine):
    if set(contexts) != set(c.ALIASES):
        raise ValueError("requires matched A/P/C/CP")
    if name in engine.CONTRASTS:
        return engine.historical_contrast(contexts, name)
    signs = ADDITIONAL_CONTRASTS[name]
    out = None
    for alias, sign in signs.items():
        term = contexts[alias] * sign
        out = term.copy() if out is None else out + term
    return out


def summarize(v):
    values = np.asarray(v, dtype=np.float64)
    if values.shape != (10, 512) or not np.isfinite(values).all():
        raise ValueError("requires 10 correlated finite time blocks, each 512D")
    vector = values.mean(axis=0)
    mag = np.linalg.norm(values, axis=1)
    return {"mean_vector_norm": float(np.linalg.norm(vector)), "signed_direction_vector": vector.tolist(),
            "block_magnitude_mean": float(mag.mean()), "block_magnitude_q10": float(np.quantile(mag, .1)),
            "block_magnitude_median": float(np.median(mag)), "block_magnitude_q90": float(np.quantile(mag, .9)),
            "n_correlated_time_blocks": 10, "replica_role": "FROZEN_STATE_EVALUATION_ONLY"}


def vector_path(candidate, branch, replica, mode="full"):
    return OUTPUT / mode / candidate / branch / f"replica_{replica:02d}" / "BLOCK_REGIONAL_VECTORS.npz"


def evaluate_candidate(candidate, engine, transition, regions, mode):
    mapping = c.context_map(candidate)
    reports, artifacts = {}, []
    for branch in c.BRANCHES:
        contrasts_by_rep = {}
        for rep in c.REPLICAS:
            contexts = {}
            upstream = []
            for job in c.expected_jobs():
                if job.candidate != candidate or job.replica != rep:
                    continue
                m = p2.validate_pair(job, branch, "full")
                if m is None:
                    raise RuntimeError("missing authenticated Phase2a job")
                contexts[job.context] = np.load(c.verify_artifact(m["rff"]), allow_pickle=False)
                upstream.append(c.artifact(p2.paths(job, branch, "full")[2]))
            if set(contexts) != set(c.ALIASES):
                raise RuntimeError("candidate-specific four contexts missing")
            cr = {}
            for alias, val in contexts.items():
                dv = engine.diffuse(val, transition)
                cr[alias] = {r: dv[:, idx, :].mean(axis=1, dtype=np.float64).astype(np.float32) for r, idx in regions.items()}
            ctr = {}
            for name in all_contrasts(engine):
                dv = engine.diffuse(contrast(contexts, name, engine), transition)
                ctr[name] = {r: dv[:, idx, :].mean(axis=1, dtype=np.float64).astype(np.float32) for r, idx in regions.items()}
            arrays = {**{f"context__{a}__{r}": v for a, region in cr.items() for r, v in region.items()},
                      **{f"contrast__{n}__{r}": v for n, region in ctr.items() for r, v in region.items()}}
            vp = vector_path(candidate, branch, rep, mode)
            # Resume only after comparing the actual numeric arrays, not a
            # newly serialized zip's non-scientific metadata.
            if vp.exists():
                with np.load(vp, allow_pickle=False) as old:
                    if set(old.files) != set(arrays) or any(not np.array_equal(old[k], v) for k, v in arrays.items()):
                        raise RuntimeError(f"existing graph vectors differ: {vp}")
            else:
                c.write_array(vp, arrays, compressed=True)
            artifacts.append({"candidate": candidate, "cluster": c.CANDIDATE_CLUSTERS[candidate],
                              "branch": branch, "replica": rep, "seed": c.SEEDS[rep],
                              "contexts": mapping, "artifact": c.artifact(vp), "upstream": upstream})
            contrasts_by_rep[rep] = ctr
        report = {}
        for name in all_contrasts(engine):
            region_reports = {}
            for region in regions:
                reps = {f"R{rep}": summarize(contrasts_by_rep[rep][name][region]) for rep in c.REPLICAS}
                vectors = {rep: np.asarray(reps[f"R{rep}"]["signed_direction_vector"]) for rep in c.REPLICAS}
                region_reports[region] = {
                    "replicas": reps,
                    "R1_R3_direction_cosine": engine.cosine(vectors[1], vectors[3]),
                    "all_pair_direction_cosines": {f"R{a}_R{b}": engine.cosine(vectors[a], vectors[b])
                                                   for a, b in ((1, 2), (1, 3), (2, 3))}}
            report[name] = {"formula": all_contrasts(engine)[name], "regions": region_reports}
        reports[branch] = report
    return reports, artifacts


def run(source, mode="full"):
    engine, anchor = c.verify_frozen_anchor()
    p2.verify(source)
    graph_path = c.FROZEN_ROOT / anchor["freeze"]["definitions"]["graph_definitions"]["path"]
    graph = c.load_json(graph_path)
    transition, regions = engine.graph_transition(graph), engine.region_indices(graph)
    candidates = list(c.CANDIDATE_CLUSTERS) if mode == "full" else ["PACER0073"]
    result, artifacts = {}, []
    for candidate in candidates:
        report, arts = evaluate_candidate(candidate, engine, transition, regions, mode)
        result[candidate] = {"cluster": c.CANDIDATE_CLUSTERS[candidate], "contexts": c.context_map(candidate),
                             "branches": report}
        artifacts += arts
        print(f"matched four-context evaluation complete: {candidate}", flush=True)
    payload = {"schema": "pacer.stage4.prospective.frozen.graph_region.v1",
               "status": f"STAGE4_{mode.upper()}_FROZEN_EVALUATION_COMPLETE", "results": result,
               "historical_anchor": anchor["receipt"], "phase2a_receipt": c.artifact(p2.RECEIPT),
               "temporal_contract": c.TEMPORAL_CONTRACT, "vector_artifacts": artifacts,
               "graph": c.artifact(graph_path), "graph_alpha": engine.ALPHA,
               "graph_steps": engine.DIFFUSION_STEPS, "regions": {r: list(map(int, idx)) for r, idx in regions.items()},
               "prospective_replica_policy": {f"R{r}": {"seed": c.SEEDS[r], "role": "EVALUATION_ONLY"} for r in c.REPLICAS},
               "historical_calibration_source": "frozen historical long-MD R2; Stage4 R2 never fitted",
               "implementation": c.code_identity(), "forbidden_operations_performed": c.FORBIDDEN,
               "claim_boundary": c.CLAIM}
    path = RESULT if mode == "full" else OUTPUT / "STAGE4_SMOKE_RESULTS_v01.json"
    c.write_json(path, payload)
    c.write_json(RECEIPT if mode == "full" else OUTPUT / "SMOKE_RECEIPT_v01.json",
                 {"status": payload["status"], "result": c.artifact(path), "phase2a_receipt": c.artifact(p2.RECEIPT),
                  "historical_freeze_sha256": c.REQUIRED_FREEZE_SHA256,
                  "implementation": c.code_identity(), "forbidden_operations_performed": c.FORBIDDEN})
    return {"status": payload["status"], "candidates": len(candidates), "vector_artifacts": len(artifacts), "result": str(path)}


def verify(source=c.SOURCE_ROOT):
    c.verify_frozen_anchor()
    p2.verify(source)
    receipt = c.load_json(RECEIPT)
    if receipt.get("historical_freeze_sha256") != c.REQUIRED_FREEZE_SHA256 or receipt.get("implementation") != c.code_identity():
        raise RuntimeError("Phase2b freeze/implementation drift")
    c.verify_artifact(receipt["phase2a_receipt"])
    if Path(receipt["result"]["path"]).resolve() != RESULT.resolve():
        raise RuntimeError("unexpected Phase2b result path")
    c.verify_artifact(receipt["result"])
    result = c.load_json(RESULT)
    if result.get("status") != "STAGE4_FULL_FROZEN_EVALUATION_COMPLETE" or set(result["results"]) != set(c.CANDIDATE_CLUSTERS):
        raise RuntimeError("Phase2b result panel mismatch")
    if result.get("temporal_contract") != c.TEMPORAL_CONTRACT or len(result["vector_artifacts"]) != 18:
        raise RuntimeError("Phase2b temporal/vector inventory mismatch")
    c.verify_artifact(result["graph"])
    for candidate, data in result["results"].items():
        if data["contexts"] != c.context_map(candidate) or data["cluster"] != c.CANDIDATE_CLUSTERS[candidate]:
            raise RuntimeError("cross-cluster contexts are forbidden")
    for item in result["vector_artifacts"]:
        c.verify_artifact(item["artifact"])
        for upstream in item["upstream"]:
            c.verify_artifact(upstream)
        with np.load(item["artifact"]["path"], allow_pickle=False) as arrays:
            if any(arrays[k].shape != (10, 512) or not np.isfinite(arrays[k]).all() for k in arrays.files):
                raise RuntimeError("graph vector shape/finite failure")
    return {"status": "STAGE4_FROZEN_EVALUATION_VERIFIED", "candidates": 3, "vector_artifacts": 18}


def main():
    a = c.parser(__doc__).parse_args()
    result = c.preflight(a.source_root) if a.preflight else verify(a.source_root) if a.verify else run(a.source_root, "smoke" if a.smoke else "full")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
