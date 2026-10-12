"""Read-only PBCFIX01 contract, download preflight, and segment authentication."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from project.pacer_fkg_v02 import stage4_pbcfix01_common_v01 as c


@c.exclusive_analysis
def sample(source, job_id, count):
    c.require(1 <= count <= 100, "sample-segments must be between 1 and 100")
    rows = c.validate_matrix(source)
    assets = c.verify_sums(source)
    c.verify_frozen_anchor()
    matches = [j for j in c.expected_jobs() if c.job_row(j)["job_id"] == job_id]
    c.require(len(matches) == 1, "sample job-id must be one of the 36 matrix jobs")
    job = matches[0]
    identity, identity_assets = c.production_identity(job, source, rows[job_id])
    paths = c.discover_inputs(job, source, count)
    before = c.artifact(paths["progress"])
    progress = c.load_json(paths["progress"])
    c.require(progress.get("identity") == identity, "sample progress identity mismatch")
    c.require(progress.get("status") in ("RUNNING", "COMPLETE")
              and count * c.SEGMENT_STEPS <= progress.get("completed_steps", -1) <= c.FINAL_STEP
              and progress["completed_steps"] % c.SEGMENT_STEPS == 0
              and progress.get("completed_ns") == progress["completed_steps"] * .002 / 1000,
              "sample progress steps/time/status mismatch")
    record = c.inspect_segments(job, paths, identity, c.frozen_sequence())
    record["identity_assets"] = identity_assets
    c.verify_artifact(before)
    c.verify_input_record(record)
    for asset in assets:
        c.verify_artifact(asset)
    result = {"status": "SAMPLE_ONLY_NOT_FULL_ACCEPTANCE", "matrix_contract_jobs": 36,
              "segments_checked": count, "input": record, "implementation": c.code_identity(),
              "forbidden_operations_performed": c.FORBIDDEN}
    c.write_json(c.RESULT_ROOT / "input_audit/samples" / c.input_scope(source) / f"{job.key}__{count:03d}.json", result)
    return {"status": result["status"], "segments_checked": count, "frames": count * 10,
            "atom14_heavy_atoms": record["receptor"]["atom14_heavy_atoms"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for action in ("contract", "preflight", "sample", "run", "all", "recover-phase1"):
        group.add_argument(f"--{action}", action="store_true")
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--job-id")
    parser.add_argument("--sample-segments", type=int, default=2)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--checkpoint", type=Path, default=c.REPO / "project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth")
    parser.add_argument("--geom2vec-source", type=Path, default=Path("C:/projects/geom2vec-source"))
    args = parser.parse_args()
    if args.contract:
        rows = c.validate_matrix(args.source_root)
        assets = c.verify_sums(args.source_root)
        result = {"status": "CONTRACT_ONLY_36_JOBS_VALID", "jobs": len(rows),
                  "systems": 12, "total_ns": 360, "segments": 3600,
                  "package_assets_verified": len(assets), "coordinates_verified": False}
    elif args.preflight:
        result = c.preflight(args.source_root)
    elif args.sample:
        if not args.job_id:
            parser.error("--sample requires --job-id")
        result = sample(args.source_root, args.job_id, args.sample_segments)
    elif args.recover_phase1:
        from project.pacer_fkg_v02 import run_stage4_pbcfix01_phase1_bs256_v01 as p1
        result = p1.recover(args.source_root)
    elif args.all:
        from project.pacer_fkg_v02 import run_stage4_pbcfix01_phase1_bs256_v01 as p1
        from project.pacer_fkg_v02 import run_stage4_pbcfix01_phase2a_frozen_apply_v01 as p2
        from project.pacer_fkg_v02 import run_stage4_pbcfix01_phase2b_graph_region_v01 as p3
        p1.run(args.source_root, "smoke", args.checkpoint, args.device, args.geom2vec_source)
        p1.run(args.source_root, "full", args.checkpoint, args.device, args.geom2vec_source)
        p2.run(args.source_root)
        p3.run(args.source_root)
        result = p3.verify(args.source_root)
    else:
        _, receipt = c.authenticate_inputs(args.source_root)
        result = {"status": receipt["status"], "jobs": len(receipt["jobs"])}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
