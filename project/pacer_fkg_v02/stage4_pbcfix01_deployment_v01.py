"""Deployment acceptance tools; no training, no modification of source assets."""
import argparse
import hashlib
import importlib.metadata
import json
import sys
from dataclasses import replace
from pathlib import Path

from project.pacer_fkg_v02 import stage4_pbcfix01_common_v01 as c

REPORT = c.RESULT_ROOT / "deployment_acceptance"


def assets(checkpoint, geom2vec_source, device):
    import torch
    from project.pacer_fkg_v02 import run_stage4_pbcfix01_phase1_bs256_v01 as p1
    _, anchor = c.verify_frozen_anchor()
    engine = p1.load_phase1_module()
    model, provenance = engine.load_model(checkpoint, device, geom2vec_source)
    versions = {name: importlib.metadata.version(name) for name in
                ("numpy", "MDAnalysis", "torch", "torch-geometric", "geom2vec")}
    result = {"status": "FROZEN_ASSETS_AND_MODEL_LOAD_PASS", "historical_anchor": anchor["receipt"],
              "model_provenance": provenance, "versions": versions, "python": sys.version,
              "device": device, "cuda_available": torch.cuda.is_available(), "implementation": c.code_identity()}
    implementation = hashlib.sha256(json.dumps(c.code_identity(), sort_keys=True).encode()).hexdigest()[:12]
    c.write_json(REPORT / "implementations" / implementation / f"ASSETS_{device}_v01.json", result)
    return model, engine, result


@c.exclusive_analysis
def topologies(source, checkpoint, geom2vec_source, device):
    import numpy as np
    import MDAnalysis as mda
    c.validate_matrix(source)
    c.verify_sums(source)
    model, engine, _ = assets(checkpoint, geom2vec_source, device)
    sequence = c.frozen_sequence()
    records = []
    for job in c.expected_jobs()[::3]:
        top = c.member(source, f"systems/{job.system}/minimized.pdb")
        build = c.member(source, f"systems/{job.system}/build_audit.json")
        audit = c.load_json(build)
        c.require(audit.get("status") == "COMPLETE" and audit.get("protocol_id") == c.PROTOCOL
                  and audit["outputs"]["minimized.pdb"] == c.sha256(top), "final topology build/hash mismatch")
        sanitized = c.sanitize_topology(top, job.system)
        u = mda.Universe(str(sanitized))
        try:
            selection, mapping = c.topology_selection(u, sequence)
            count = c.validate_selected_atoms(u, selection, sequence)
        finally:
            u.trajectory.close()
        engine.SELECTION = selection
        u, idx, z, resid, bb, sc, _ = engine.topology_mapping(replace(job, topology=sanitized, trajectory=sanitized, sequence=sequence))
        try:
            bs = engine.infer_bs256(model, u, idx, z, resid, bb, sc, [0], device)
            c.require(bs.shape == (1, 270, 256) and np.isfinite(bs).all(), "final topology encoder smoke failed")
        finally:
            u.trajectory.close()
        records.append({"system": job.system, "topology": c.artifact(top), "build_audit": c.artifact(build),
                        "atoms": len(u.atoms), "atom14_heavy_atoms": count, "receptor": mapping,
                        "encoder_shape": list(bs.shape), "encoder_finite": True})
        print(f"final topology + encoder {len(records)}/12 {job.system}", flush=True)
    result = {"status": "12_FINAL_MEMBRANE_TOPOLOGIES_ENCODER_PASS", "systems": records,
              "implementation": c.code_identity(), "forbidden_operations_performed": c.FORBIDDEN}
    c.write_json(REPORT / c.input_scope(source) / "FINAL_12_TOPOLOGIES_v01.json", result)
    return {"status": result["status"], "systems": len(records)}


@c.exclusive_analysis
def sample_encoder(source, job_id, checkpoint, geom2vec_source, device):
    import numpy as np
    from project.pacer_fkg_v02 import run_stage4_pbcfix01_inputs_v01 as inputs
    from project.pacer_fkg_v02 import run_stage4_pbcfix01_phase1_bs256_v01 as p1
    inputs.sample(source, job_id, 2)
    job = next(j for j in c.expected_jobs() if c.job_row(j)["job_id"] == job_id)
    sample_path = c.RESULT_ROOT / "input_audit/samples" / c.input_scope(source) / f"{job.key}__002.json"
    record = c.load_json(sample_path)["input"]
    job = replace(job, topology=Path(record["sanitized_topology"]["path"]),
                  trajectory=tuple(Path(x["path"]) for x in record["trajectory"]),
                  sequence=record["receptor"]["sequence"], selection=record["selection"], input_record=record)
    model, engine, _ = assets(checkpoint, geom2vec_source, device)
    u, idx, z, resid, bb, sc, _ = p1.mapped_topology(engine, job)
    try:
        bs = np.concatenate([engine.infer_bs256(model, u, idx, z, resid, bb, sc, list(range(i, i+4)), device)
                             for i in range(0, 20, 4)])
        replay_ids = [0, 9, 10, 19]
        replay = engine.replay_metrics(bs[replay_ids], engine.infer_bs256(model, u, idx, z, resid, bb, sc, replay_ids, device))
        c.require(replay["passed"], "real-segment boundary encoder replay failed")
    finally:
        u.trajectory.close()
    frozen, anchor = c.verify_frozen_anchor()
    # Four retained 50-ps samples cannot form a formal 20-frame Stage4 block.
    # Exercise the same frozen transforms/graph on this explicitly diagnostic
    # short window, without publishing Phase2a/2b scientific completion receipts.
    sampled = bs[::5]
    delta = np.diff(sampled, axis=0)
    descriptors = {"STATE_MOTION": np.concatenate((sampled.mean(axis=0, dtype=np.float64).astype(np.float32),
                    np.sqrt(np.square(delta, dtype=np.float64).mean(axis=0)).astype(np.float32)), axis=-1)[None],
                   "SIGNED_DRIFT": delta.mean(axis=0, dtype=np.float64).astype(np.float32)[None]}
    graph = c.load_json(c.FROZEN_ROOT / anchor["freeze"]["definitions"]["graph_definitions"]["path"])
    transition, regions = frozen.graph_transition(graph), frozen.region_indices(graph)
    diagnostics = {}
    for branch, values in descriptors.items():
        state = frozen.load_frozen_state(anchor, branch)
        rff = frozen.rff_frozen(frozen.normalize_frozen(values, state), state)
        diffused = frozen.diffuse(rff, transition)
        c.require(np.isfinite(diffused).all(), "real sample frozen numeric/graph path failed")
        diagnostics[branch] = {"rff_shape": list(rff.shape), "region_shapes":
                               {name: list(diffused[:, indices].mean(axis=1).shape) for name, indices in regions.items()}}
    c.verify_input_record(record)
    result = {"status": "REAL_TWO_SEGMENTS_ENCODER_AND_FROZEN_PATH_PASS", "job": job.identity(),
              "sample_receipt": c.artifact(sample_path), "encoder_shape": list(bs.shape),
              "replay": {**replay, "frame_ids": replay_ids}, "diagnostic_frozen_path": diagnostics,
              "short_window_is_formal_phase2_block": False, "four_context_contrasts_tested": False,
              "full_360ns_acceptance": False, "implementation": c.code_identity(),
              "forbidden_operations_performed": c.FORBIDDEN}
    c.write_json(REPORT / c.input_scope(source) / "REAL_TWO_SEGMENTS_ENCODER_v01.json", result)
    return {k: result[k] for k in ("status", "encoder_shape", "full_360ns_acceptance")}


def transfer_manifest(source, manifest, create=False):
    """Portable relative SHA256 ledger; publish outside the live campaign."""
    source, manifest = Path(source).resolve(), Path(manifest).resolve()
    c.require(not manifest.is_relative_to(source), "transfer manifest must be outside production/source tree")
    if create:
        c.validate_matrix(source)
        for job in c.expected_jobs():
            paths = c.discover_inputs(job, source)
            progress = c.load_json(paths["progress"])
            c.require(progress.get("status") == "COMPLETE" and progress.get("completed_steps") == c.FINAL_STEP,
                      "only generate final transfer manifest after all 36 tasks complete")
        files = sorted(p for p in source.rglob("*") if p.is_file()
                       and not any(".interrupted." in part for part in p.parts)
                       and p.suffix not in (".lock", ".log", ".tmp", ".pyc"))
        manifest.parent.mkdir(parents=True, exist_ok=True)
        temporary = manifest.with_name(manifest.name + ".partial")
        c.require(not manifest.exists(), "transfer manifest already exists")
        snapshot = [(p, p.stat().st_size, p.stat().st_mtime_ns) for p in files]
        with temporary.open("x", encoding="utf-8", newline="\n") as output:
            for path in files:
                c.require(not path.is_symlink() and path.resolve().is_relative_to(source), "unsafe transfer member")
                before = path.stat()
                digest = c.sha256(path)
                after = path.stat()
                c.require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), "source changed during hashing")
                output.write(f"{digest}  {path.relative_to(source).as_posix()}\n")
        for path, size, mtime in snapshot:
            c.require((path.stat().st_size, path.stat().st_mtime_ns) == (size, mtime), "source changed during transfer manifest generation")
        import os
        os.link(temporary, manifest)
        temporary.unlink()
    else:
        seen = set()
        for line in manifest.read_text(encoding="utf-8-sig").splitlines():
            digest, relative = line.split("  ", 1)
            c.require(relative not in seen, "duplicate transfer member")
            seen.add(relative)
            path = c.member(source, relative)
            c.require(c.sha256(path) == digest, f"transfer SHA256 mismatch: {relative}")
        for job in c.expected_jobs():
            paths = c.discover_inputs(job, source)
            for path in [paths["topology"], paths["progress"], *(d / name for d in paths["segments"] for name in
                         ("trajectory.dcd", "state.csv", "checkpoint.chk", "receipt.json"))]:
                c.require(path.relative_to(source).as_posix() in seen, "transfer manifest missing required production member")
    return {"status": "TRANSFER_MANIFEST_CREATED" if create else "TRANSFER_SHA256_PASS", "manifest": str(manifest)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("assets", "topologies", "sample-encoder", "create-sha256", "verify-sha256"))
    p.add_argument("--source-root", type=Path)
    p.add_argument("--manifest", type=Path)
    p.add_argument("--job-id", default=c.job_row(c.expected_jobs()[0])["job_id"])
    p.add_argument("--checkpoint", type=Path, default=c.REPO / "project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth")
    p.add_argument("--geom2vec-source", type=Path, default=Path("C:/projects/geom2vec-source"))
    p.add_argument("--device", default="cuda")
    a = p.parse_args()
    if a.action != "assets" and not a.source_root:
        p.error("--source-root is required")
    if a.action.endswith("sha256"):
        if not a.manifest:
            p.error("--manifest is required")
        result = transfer_manifest(a.source_root, a.manifest, a.action == "create-sha256")
    elif a.action == "assets":
        _, _, result = assets(a.checkpoint, a.geom2vec_source, a.device)
        result = {k: result[k] for k in ("status", "device", "versions", "cuda_available")}
    elif a.action == "topologies":
        result = topologies(a.source_root, a.checkpoint, a.geom2vec_source, a.device)
    else:
        result = sample_encoder(a.source_root, a.job_id, a.checkpoint, a.geom2vec_source, a.device)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
