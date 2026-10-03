"""Stage4 Phase1: authenticate 36 inputs, extract frozen C1-BS256, no fitting.

Run --preflight during download. --smoke authenticates all inputs but embeds only
four frames of PACER0073/CP/R1. --run requires that smoke receipt to verify.
"""
from __future__ import annotations
import importlib.util
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from project.pacer_fkg_v02 import stage4_prospective_common_v01 as c

REPORT = c.RESULT_ROOT / "phase1_bs256"
INPUT_RECEIPT = c.RESULT_ROOT / "input_audit/INPUT_AUTHENTICATION_v01.json"
SMOKE_RECEIPT = REPORT / "SMOKE_RECEIPT_v01.json"
FULL_RECEIPT = REPORT / "PHASE1_FREEZE_RECEIPT_v01.json"
SMOKE_IDS = (0, 1, 2, 3)
DEFAULT_CHECKPOINT = REPO / "project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth"


def load_phase1_module():
    name = "_stage4_private_frozen_phase1"
    spec = importlib.util.spec_from_file_location(name, REPO / "project/pacer_fkg_v02/run_phase1_bs256.py")
    engine = importlib.util.module_from_spec(spec)
    sys.modules[name] = engine  # dataclasses resolves its defining module here
    spec.loader.exec_module(engine)
    engine.RESNAME_MAP = {**engine.RESNAME_MAP, **c.RESNAME_ALIASES}
    return engine


def mapped_topology(engine, job):
    # Only this new private instance receives the authenticated construct
    # selection. The historical imported module and topology stay unchanged.
    engine.SELECTION = job.selection
    return engine.topology_mapping(job)


def cache_paths(job, mode):
    base = c.CACHE_ROOT / "bs256" / mode / job.system / f"replica_{job.replica:02d}"
    return base / "C1_BS256.npy", REPORT / "manifests" / mode / f"{job.key}.json"


def identity(job, mode):
    frames = SMOKE_IDS if mode == "smoke" else tuple(range(c.RAW_FRAMES))
    return {**job.identity(), "mode": mode, "frame_ids": list(frames),
            "shape": [len(frames), c.N_RESIDUES, c.BS_WIDTH], "dtype": "float32",
            "input": job.input_record, "sequence": job.sequence, "selection": job.selection,
            "historical_freeze_sha256": c.REQUIRED_FREEZE_SHA256,
            "temporal_contract": c.TEMPORAL_CONTRACT, "implementation": c.code_identity(),
            "preprocessing": {"make_whole": False, "center": False, "align": False, "pbc_transform": False}}


def validate_pair(job, mode):
    import numpy as np
    path, manifest_path = cache_paths(job, mode)
    if not path.exists() and not manifest_path.exists():
        return None
    if not path.is_file() or not manifest_path.is_file():
        raise RuntimeError(f"incomplete cache/manifest pair; refusing overwrite: {job.key}")
    m = c.load_json(manifest_path)
    for k, v in identity(job, mode).items():
        if m.get(k) != v:
            raise RuntimeError(f"cache identity changed ({k}): {job.key}")
    if Path(m["cache"]["path"]).resolve() != path.resolve():
        raise RuntimeError("cache manifest points outside expected job path")
    c.verify_artifact(m["cache"])
    arr = np.load(path, mmap_mode="r", allow_pickle=False)
    try:
        if arr.shape != tuple(m["shape"]) or arr.dtype != np.float32 or not np.isfinite(arr).all():
            raise RuntimeError(f"cache shape/dtype/finite failure: {job.key}")
        gly = [aa == "G" for aa in job.sequence]
        if np.count_nonzero(arr[:, gly, 128:]) != 0:
            raise RuntimeError("nonzero Gly sidechain block")
    finally:
        arr._mmap.close()
    if m.get("replay", {}).get("passed") is not True or m["replay"].get("frame_ids") != list(SMOKE_IDS):
        raise RuntimeError("cache deterministic replay receipt failed")
    return m


def extract_one(job, mode, engine, model, device, model_record):
    import numpy as np
    cached = validate_pair(job, mode)
    if cached is not None:
        print(f"verified existing cache, skip {job.key}", flush=True)
        return cached
    path, manifest_path = cache_paths(job, mode)
    path = c.guard_output(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Any stale partial is evidence of an incomplete extraction, never erased
    # or treated as a valid resumable cache.
    if list(path.parent.glob(".C1_BS256.npy.partial.*")):
        raise RuntimeError(f"incomplete extraction temporary present: {path.parent}")
    tmp = c.guard_output(path.with_name(f".{path.name}.partial.{os.getpid()}"))
    u = out = None
    frame_ids = SMOKE_IDS if mode == "smoke" else tuple(range(c.RAW_FRAMES))
    try:
        u, atom_idx, atom_z, res_idx, bb, sc, records = mapped_topology(engine, job)
        out = np.lib.format.open_memmap(tmp, mode="w+", dtype=np.float32,
                                        shape=(len(frame_ids), c.N_RESIDUES, c.BS_WIDTH))
        for start in range(0, len(frame_ids), engine.BATCH_SIZE):
            ids = frame_ids[start:start + engine.BATCH_SIZE]
            out[start:start + len(ids)] = engine.infer_bs256(model, u, atom_idx, atom_z, res_idx, bb, sc, ids, device)
            out.flush()
            print(f"{job.key} {start + len(ids)}/{len(frame_ids)}", flush=True)
        if np.count_nonzero(out[:, [aa == 'G' for aa in job.sequence], 128:]) != 0:
            raise RuntimeError("Gly sidechain output nonzero")
        observed = np.asarray(out[:4]).copy()
        replay = engine.replay_metrics(observed, engine.infer_bs256(model, u, atom_idx, atom_z, res_idx, bb, sc, SMOKE_IDS, device))
        replay["frame_ids"] = list(SMOKE_IDS)
        if not replay["passed"]:
            raise RuntimeError("deterministic replay failed")
        out._mmap.close()
        out = None
        # Prove the source did not change during extraction. All hashes were
        # first authenticated before inference, including full coordinate scan.
        for k in ("trajectory", "topology", "state", "progress"):
            c.verify_artifact(job.input_record[k])
        os.link(tmp, path)
        manifest = {**identity(job, mode), "cache": c.artifact(path), "replay": replay,
                    "model_provenance": model_record, "atom14_heavy_atoms": int(len(atom_idx)),
                    "forbidden_operations_performed": c.FORBIDDEN}
        c.write_json(manifest_path, manifest)
        return manifest
    finally:
        if out is not None:
            out._mmap.close()
        if tmp.exists():
            tmp.unlink()  # only this invocation's unpublished temporary
        if u is not None:
            u.trajectory.close()


def jobs_from_receipt(source):
    r = c.load_json(INPUT_RECEIPT)
    if r.get("status") != "36_OF_36_ACCEPTED" or r.get("implementation") != c.code_identity() or r.get("temporal_contract") != c.TEMPORAL_CONTRACT:
        raise RuntimeError("input authentication receipt is not current/complete")
    if Path(r["source_root"]).resolve() != Path(source).resolve():
        raise RuntimeError("source-root differs from authenticated source")
    c.verify_artifact(r["job_matrix"])
    records = {(x["system"], x["replica"]): x for x in r["jobs"]}
    if len(r["jobs"]) != 36 or len(records) != 36:
        raise RuntimeError("input receipt must contain exactly 36 jobs")
    from dataclasses import replace
    jobs = []
    for job in c.expected_jobs():
        x = records[(job.system, job.replica)]
        if any(x[k] != v for k, v in job.identity().items()):
            raise RuntimeError("candidate/cluster/context/seed identity mismatch")
        paths = c.discover_inputs(job, source)
        for key in ("trajectory", "topology", "state", "progress"):
            if Path(x[key]["path"]).resolve() != paths[key].resolve():
                raise RuntimeError("authenticated source path changed")
            c.verify_artifact(x[key])
        c.verify_artifact(x["sanitized_topology"])
        jobs.append(replace(job, topology=Path(x["sanitized_topology"]["path"]), trajectory=paths["trajectory"],
                            sequence=x["receptor"]["sequence"], selection=x["selection"], input_record=x))
    return jobs


def verify(source=c.SOURCE_ROOT, mode="full", verify_sources=True):
    c.verify_frozen_anchor()
    path = SMOKE_RECEIPT if mode == "smoke" else FULL_RECEIPT
    r = c.load_json(path)
    expected_count = 1 if mode == "smoke" else 36
    if r.get("status") != f"PHASE1_{mode.upper()}_FROZEN" or r.get("cache_count") != expected_count:
        raise RuntimeError("Phase1 receipt contract failed")
    if r.get("historical_freeze_sha256") != c.REQUIRED_FREEZE_SHA256 or r.get("input_receipt_sha256") != c.sha256(INPUT_RECEIPT):
        raise RuntimeError("Phase1 upstream receipt changed")
    # Revalidation of source hashes is mandatory for standalone verification.
    jobs = jobs_from_receipt(source) if verify_sources else receipt_jobs_unhashed(source)
    if mode == "smoke":
        jobs = [j for j in jobs if (j.candidate, j.context, j.replica) == ("PACER0073", "CP", 1)]
    records = {x["path"]: x for x in r["manifests"]}
    if len(records) != expected_count:
        raise RuntimeError("Phase1 manifest inventory mismatch")
    for job in jobs:
        _, mp = cache_paths(job, mode)
        c.verify_artifact(records[str(mp.resolve())])
        if validate_pair(job, mode) is None:
            raise RuntimeError("missing Phase1 cache")
    return {"status": f"PHASE1_{mode.upper()}_VERIFIED", "cache_count": len(jobs), "receipt": c.artifact(path)}


def receipt_jobs_unhashed(source):
    # Used only immediately after fresh authenticate_inputs in the same run.
    from dataclasses import replace
    r = c.load_json(INPUT_RECEIPT)
    rec = {(x["system"], x["replica"]): x for x in r["jobs"]}
    return [replace(j, topology=Path(rec[(j.system, j.replica)]["sanitized_topology"]["path"]),
                    trajectory=Path(rec[(j.system, j.replica)]["trajectory"]["path"]),
                    sequence=rec[(j.system, j.replica)]["receptor"]["sequence"],
                    selection=rec[(j.system, j.replica)]["selection"], input_record=rec[(j.system, j.replica)])
            for j in c.expected_jobs()]


def run(source, mode, checkpoint, device, geom2vec_source):
    jobs, inputs = c.authenticate_inputs(source)
    if mode == "full":
        verify(source, "smoke", verify_sources=False)
    if mode == "smoke":
        jobs = [j for j in jobs if (j.candidate, j.context, j.replica) == ("PACER0073", "CP", 1)]
    engine = load_phase1_module()
    model, model_record = engine.load_model(checkpoint, device, geom2vec_source)
    manifests = []
    for job in jobs:
        extract_one(job, mode, engine, model, device, model_record)
        manifests.append(c.artifact(cache_paths(job, mode)[1]))
    receipt = {"status": f"PHASE1_{mode.upper()}_FROZEN", "cache_count": len(jobs),
               "historical_freeze_sha256": c.REQUIRED_FREEZE_SHA256,
               "input_receipt_sha256": c.sha256(INPUT_RECEIPT), "temporal_contract": c.TEMPORAL_CONTRACT,
               "manifests": manifests, "forbidden_operations_performed": c.FORBIDDEN,
               "prospective_replica_policy": "R1/R2/R3 evaluation only"}
    c.write_json(SMOKE_RECEIPT if mode == "smoke" else FULL_RECEIPT, receipt)
    return {"status": receipt["status"], "cache_count": len(jobs), "device": device}


def main():
    p = c.parser(__doc__)
    p.add_argument("--device", default="cuda")
    p.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    p.add_argument("--geom2vec-source", type=Path, default=Path("C:/projects/geom2vec-source"))
    p.add_argument("--mode", choices=("smoke", "full"), default="full", help="cache mode for --verify")
    a = p.parse_args()
    if a.preflight:
        result = c.preflight(a.source_root)
    elif a.verify:
        result = verify(a.source_root, a.mode)
    else:
        result = run(a.source_root, "smoke" if a.smoke else "full", a.checkpoint, a.device, a.geom2vec_source)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
