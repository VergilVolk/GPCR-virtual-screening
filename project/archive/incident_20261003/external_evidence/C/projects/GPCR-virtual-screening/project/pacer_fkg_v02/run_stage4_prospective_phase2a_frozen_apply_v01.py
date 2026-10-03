"""Stage4 Phase2a: fixed stride-5 blocks, historical frozen-state application.

No calibration, bandwidth selection or random-state generation is exposed.
--smoke uses one completed FULL Phase1 cache, never a fabricated MD input.
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
from project.pacer_fkg_v02 import run_stage4_prospective_phase1_bs256_v01 as p1

REPORT = c.RESULT_ROOT / "phase2a_frozen_apply"
RECEIPT = REPORT / "PHASE2A_FREEZE_RECEIPT_v01.json"


def construct(frames):
    """Exact Stage-B reductions with Stage4's explicit physical-time contract."""
    if frames.shape != (c.RAW_FRAMES, c.N_RESIDUES, c.BS_WIDTH):
        raise ValueError(f"expected raw BS256 [1000,270,256], got {frames.shape}")
    if not np.isfinite(frames).all():
        raise ValueError("non-finite BS256 input")
    sampled = np.asarray(frames[::c.TEMPORAL_STRIDE], dtype=np.float32)
    if sampled.shape != (c.ANALYSIS_FRAMES, c.N_RESIDUES, c.BS_WIDTH):
        raise ValueError("stride-5 shape contract failed")
    blocks = sampled.reshape(c.BLOCKS_PER_TRAJECTORY, c.BLOCK_FRAMES, c.N_RESIDUES, c.BS_WIDTH)
    mu = blocks.mean(axis=1, dtype=np.float64).astype(np.float32)
    delta = np.diff(blocks, axis=1)
    drift = delta.mean(axis=1, dtype=np.float64).astype(np.float32)
    rms = np.sqrt(np.mean(np.square(delta, dtype=np.float64), axis=1)).astype(np.float32)
    endpoint = (blocks[:, -1] - blocks[:, 0]) / np.float32(c.BLOCK_FRAMES - 1)
    if not np.allclose(drift, endpoint, rtol=2e-5, atol=2e-6):
        raise RuntimeError("SIGNED_DRIFT endpoint identity failed")
    return {"STATE_MOTION": np.concatenate((mu, rms), axis=-1).astype(np.float32, copy=False),
            "SIGNED_DRIFT": drift}


def paths(job, branch, mode):
    root = c.CACHE_ROOT / "phase2a_frozen_apply" / mode / job.system / f"replica_{job.replica:02d}"
    manifest = REPORT / "manifests" / mode / f"{job.key}__{branch}.json"
    return root / f"{branch}.npy", root / f"{branch}_RFF.npy", manifest


def expected_identity(job, branch, mode):
    return {**job.identity(), "branch": branch, "mode": mode,
            "input_manifest": c.artifact(p1.cache_paths(job, "full")[1]),
            "phase1_receipt": c.artifact(p1.FULL_RECEIPT),
            "historical_freeze_sha256": c.REQUIRED_FREEZE_SHA256,
            "temporal_contract": c.TEMPORAL_CONTRACT, "implementation": c.code_identity(),
            "prospective_replica_role": "FROZEN_STATE_EVALUATION_ONLY"}


def validate_pair(job, branch, mode):
    dp, rp, mp = paths(job, branch, mode)
    if not any(p.exists() for p in (dp, rp, mp)):
        return None
    if not all(p.is_file() for p in (dp, rp, mp)):
        raise RuntimeError(f"incomplete Phase2a triple; refusing overwrite: {job.key}/{branch}")
    m = c.load_json(mp)
    for k, v in expected_identity(job, branch, mode).items():
        if m.get(k) != v:
            raise RuntimeError(f"Phase2a identity mismatch {k}: {job.key}")
    for key, path, width in (("descriptor", dp, 512 if branch == "STATE_MOTION" else 256), ("rff", rp, 512)):
        if Path(m[key]["path"]).resolve() != path.resolve():
            raise RuntimeError("Phase2a artifact path differs")
        c.verify_artifact(m[key])
        a = np.load(path, mmap_mode="r", allow_pickle=False)
        try:
            if a.shape != (10, 270, width) or a.dtype != np.float32 or not np.isfinite(a).all():
                raise RuntimeError("Phase2a array shape/dtype/finite failure")
        finally:
            a._mmap.close()
    return m


def run(source, mode="full"):
    engine, anchor = c.verify_frozen_anchor()
    p1.verify(source)
    jobs = p1.receipt_jobs_unhashed(source)
    if mode == "smoke":
        jobs = jobs[:1]
    states = {b: engine.load_frozen_state(anchor, b) for b in c.BRANCHES}
    manifests = []
    for job in jobs:
        cached = {b: validate_pair(job, b, mode) for b in c.BRANCHES}
        if all(cached.values()):
            manifests += [c.artifact(paths(job, b, mode)[2]) for b in c.BRANCHES]
            continue
        cache_path, _ = p1.cache_paths(job, "full")
        frames = np.load(cache_path, mmap_mode="r", allow_pickle=False)
        try:
            descriptors = construct(frames)
        finally:
            frames._mmap.close()
        for branch in c.BRANCHES:
            dp, rp, mp = paths(job, branch, mode)
            if cached[branch] is None:
                mapped = engine.rff_frozen(engine.normalize_frozen(descriptors[branch], states[branch]), states[branch])
                if mapped.shape != (10, 270, 512) or not np.isfinite(mapped).all():
                    raise RuntimeError("historical RFF output contract failed")
                c.write_array(dp, descriptors[branch])
                c.write_array(rp, mapped)
                c.write_json(mp, {**expected_identity(job, branch, mode), "descriptor": c.artifact(dp),
                                  "rff": c.artifact(rp), "historical_anchor": anchor["receipt"],
                                  "forbidden_operations_performed": c.FORBIDDEN})
            manifests.append(c.artifact(mp))
        print(f"frozen-state applied {job.key}: 10 x 1-ns nominal blocks", flush=True)
    receipt = {"status": f"PHASE2A_{mode.upper()}_FROZEN", "job_count": len(jobs),
               "manifest_count": len(manifests), "manifests": manifests,
               "phase1_receipt": c.artifact(p1.FULL_RECEIPT), "historical_anchor": anchor["receipt"],
               "temporal_contract": c.TEMPORAL_CONTRACT, "forbidden_operations_performed": c.FORBIDDEN}
    c.write_json(RECEIPT if mode == "full" else REPORT / "SMOKE_RECEIPT_v01.json", receipt)
    return {"status": receipt["status"], "job_count": len(jobs), "manifest_count": len(manifests)}


def verify(source=c.SOURCE_ROOT):
    c.verify_frozen_anchor()
    p1.verify(source)
    r = c.load_json(RECEIPT)
    if r.get("status") != "PHASE2A_FULL_FROZEN" or r.get("job_count") != 36 or r.get("manifest_count") != 72:
        raise RuntimeError("Phase2a receipt count/status mismatch")
    if r.get("temporal_contract") != c.TEMPORAL_CONTRACT:
        raise RuntimeError("Phase2a temporal contract changed")
    c.verify_artifact(r["phase1_receipt"])
    records = {x["path"]: x for x in r["manifests"]}
    if len(records) != 72:
        raise RuntimeError("Phase2a receipt manifest set mismatch")
    for job in p1.receipt_jobs_unhashed(source):
        for branch in c.BRANCHES:
            mp = paths(job, branch, "full")[2]
            c.verify_artifact(records[str(mp.resolve())])
            if validate_pair(job, branch, "full") is None:
                raise RuntimeError("missing Phase2a artifact")
    return {"status": "PHASE2A_VERIFIED", "jobs": 36, "branches": 2, "receipt": c.artifact(RECEIPT)}


def main():
    a = c.parser(__doc__).parse_args()
    result = c.preflight(a.source_root) if a.preflight else verify(a.source_root) if a.verify else run(a.source_root, "smoke" if a.smoke else "full")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
