#!/usr/bin/env python
"""PACER-DC 20 ns closure: frozen C1-BS256 extraction only.

Scope:
- LY2119620 no-probe / +probe
- R1/R2/R3
- 400 stored frames per trajectory
- frozen C1-BS256 encoder
- no calibration, normalization, RFF, graph diffusion, contrasts, or biology

This runner reuses the frozen PACER-FKG-v02 encoder implementation but writes
only into the close-loop result/cache namespace.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from project.pacer_fkg_v02.run_phase1_bs256 import (  # noqa: E402
    BACKBONE_ATOMS,
    BATCH_SIZE,
    BS_WIDTH,
    N_RESIDUES,
    RESTYPE_3TO1,
    infer_bs256,
    load_model,
    replay_metrics,
    topology_mapping,
    verify_phase0_and_frozen_candidate,
)

VERSION = "PACER_DC_CLOSE_LOOP_PHASE1_BS256_v01"

CLOSURE_ROOT = REPO / "project/results/pacer_dc_close_loop_20ns_v01"
PRODUCTION_ROOT = CLOSURE_ROOT / "production"
TOPOLOGY_ROOT = CLOSURE_ROOT / "topology_sanitized"

CACHE_ROOT = REPO / "project/cache/pacer_dc_close_loop_20ns_v01/track_b_fkg_v02/bs256"
REPORT_ROOT = CLOSURE_ROOT / "track_b_fkg_v02/phase1_bs256"

EXEC_MANIFEST = CLOSURE_ROOT / "CLOSE_LOOP_EXECUTION_MANIFEST_v01.json"
PREFLIGHT = CLOSURE_ROOT / "LY2119620_CLOSE_LOOP_PREFLIGHT_v01.json"
TRANSFER_MANIFEST = CLOSURE_ROOT / "TRANSFER_SHA256_v01.txt"
TOPOLOGY_AUDIT = TOPOLOGY_ROOT / "LY_TOPOLOGY_SANITIZATION_v01.json"

MAPPING = (
    REPO
    / "project/results/pacer_dc_four_context_v01/compound110/"
    / "G2_RESIDUE_MAPPING_v01.csv"
)

DEFAULT_CHECKPOINT = (
    REPO
    / "project/tools/geom2vec_checkpoints/"
    / "visnet_l6_h64_rbf64_r75.pth"
)

EXPECTED_EXEC_SHA = "d0333145061bc5ff6a4b52df33269cf1893c77472e2c839252c03c84e6a39b77"
EXPECTED_PREFLIGHT_SHA = "6b5d6c33ee7405c152ac0890c259da425091927fd87e02d311d9457c20313b3b"
EXPECTED_TRANSFER_SHA = "d0383b9fda1bf6de8d227b4bf26a623c4621218292b602812def9aed41799696"

SYSTEMS = (
    "LY2119620__candidate_no_probe",
    "LY2119620__candidate_probe",
)
REPLICAS = (1, 2, 3)
SEEDS = {1: 27101, 2: 38201, 3: 49301}

N_FRAMES = 400
REPLAY_FRAME_IDS = (0, 1, 2, 3)


@dataclass(frozen=True)
class ClosureJob:
    mode: str
    system: str
    replica: int
    topology: Path
    trajectory: Path
    trajectory_sha256: str
    topology_sha256: str
    sequence: str
    frame_ids: tuple[int, ...]

    @property
    def key(self) -> str:
        return f"{self.system}__replica_{self.replica:02d}"

    @property
    def output_shape(self) -> tuple[int, int, int]:
        return (len(self.frame_ids), N_RESIDUES, BS_WIDTH)

    @property
    def cache_path(self) -> Path:
        leaf = (
            "C1_BS256.frames_0000_0003.npy"
            if self.mode == "smoke"
            else "C1_BS256.npy"
        )
        return (
            CACHE_ROOT
            / self.mode
            / self.system
            / f"replica_{self.replica:02d}"
            / leaf
        )

    @property
    def manifest_path(self) -> Path:
        return (
            REPORT_ROOT
            / "manifests"
            / self.mode
            / f"{self.key}.json"
        )


def utc_now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.partial.{os.getpid()}")
    tmp.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    os.replace(tmp, path)


def close_memmap(array) -> None:
    mapping = getattr(array, "_mmap", None) if array is not None else None
    if mapping is not None:
        mapping.close()


def transfer_hashes() -> dict[str, str]:
    records = {}
    for line in TRANSFER_MANIFEST.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, path = line.split(None, 1)
        records[path.strip().replace("\\", "/")] = digest
    return records


def frozen_sequence() -> str:
    with MAPPING.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != N_RESIDUES:
        raise RuntimeError(
            f"frozen mapping has {len(rows)} residues, expected {N_RESIDUES}"
        )

    sequence = []
    for row in rows:
        resname = row["resname"].strip()
        aa = RESTYPE_3TO1.get(resname)
        if aa is None:
            raise RuntimeError(f"unsupported residue in frozen mapping: {resname}")
        sequence.append(aa)

    return "".join(sequence)


def verify_closure_anchor() -> dict:
    if sha256(EXEC_MANIFEST) != EXPECTED_EXEC_SHA:
        raise RuntimeError("closure execution manifest SHA mismatch")

    if sha256(PREFLIGHT) != EXPECTED_PREFLIGHT_SHA:
        raise RuntimeError("LY preflight SHA mismatch")

    if sha256(TRANSFER_MANIFEST) != EXPECTED_TRANSFER_SHA:
        raise RuntimeError("transfer manifest SHA mismatch")

    execution = load_json(EXEC_MANIFEST)
    preflight = load_json(PREFLIGHT)
    topology = load_json(TOPOLOGY_AUDIT)

    if execution.get("status") != "STAGE0_PASS":
        raise RuntimeError("Stage 0 is not PASS")

    if preflight.get("status") != "LY2119620_CLOSE_LOOP_PREFLIGHT_PASS":
        raise RuntimeError("LY preflight is not PASS")

    new_production = preflight["new_production"]

    if new_production["trajectory_ns"] != 20:
        raise RuntimeError("closure trajectory length changed")

    if new_production["trajectory_ps"] != 50:
        raise RuntimeError("closure frame spacing changed")

    if new_production["frames_per_trajectory"] != N_FRAMES:
        raise RuntimeError("closure frame count changed")

    topo_records = {
        item["system"]: item
        for item in topology["records"]
    }

    if set(topo_records) != set(SYSTEMS):
        raise RuntimeError("sanitized topology audit system set mismatch")

    return {
        "execution_manifest": {
            "path": str(EXEC_MANIFEST),
            "sha256": EXPECTED_EXEC_SHA,
        },
        "preflight": {
            "path": str(PREFLIGHT),
            "sha256": EXPECTED_PREFLIGHT_SHA,
        },
        "transfer_manifest": {
            "path": str(TRANSFER_MANIFEST),
            "sha256": EXPECTED_TRANSFER_SHA,
        },
        "topology_audit": {
            "path": str(TOPOLOGY_AUDIT),
            "sha256": sha256(TOPOLOGY_AUDIT),
        },
        "topology_records": topo_records,
    }


def build_jobs(mode: str, anchor: dict) -> list[ClosureJob]:
    hashes = transfer_hashes()
    sequence = frozen_sequence()
    jobs = []

    for system in SYSTEMS:
        topo_record = anchor["topology_records"][system]

        topology = TOPOLOGY_ROOT / f"{system}_minimized_no_conect.pdb"

        if not topology.is_file():
            raise FileNotFoundError(topology)

        topology_hash = sha256(topology)

        if topology_hash != topo_record["sanitized_sha256"]:
            raise RuntimeError(
                f"sanitized topology SHA mismatch: {system}"
            )

        for replica in REPLICAS:
            replica_dir = (
                PRODUCTION_ROOT
                / system
                / f"replica_{replica:02d}"
            )

            trajectory = replica_dir / "trajectory.dcd"
            progress_path = replica_dir / "progress.json"

            progress = load_json(progress_path)

            if progress.get("status") != "complete":
                raise RuntimeError(f"incomplete trajectory: {system} R{replica}")

            if float(progress.get("completed_ns", -1)) != 20.0:
                raise RuntimeError(f"wrong completed_ns: {system} R{replica}")

            if int(progress.get("seed", -1)) != SEEDS[replica]:
                raise RuntimeError(f"seed mismatch: {system} R{replica}")

            relative = trajectory.relative_to(REPO).as_posix()

            expected_hash = hashes.get(relative)

            if expected_hash is None:
                raise RuntimeError(
                    f"trajectory absent from transfer manifest: {relative}"
                )

            actual_hash = sha256(trajectory)

            if actual_hash != expected_hash:
                raise RuntimeError(
                    f"trajectory SHA mismatch: {system} R{replica}"
                )

            frame_ids = (
                REPLAY_FRAME_IDS
                if mode == "smoke"
                else tuple(range(N_FRAMES))
            )

            jobs.append(
                ClosureJob(
                    mode=mode,
                    system=system,
                    replica=replica,
                    topology=topology,
                    trajectory=trajectory,
                    trajectory_sha256=actual_hash,
                    topology_sha256=topology_hash,
                    sequence=sequence,
                    frame_ids=tuple(frame_ids),
                )
            )

    if mode == "smoke":
        # One representative trajectory only.
        jobs = [jobs[0]]

    expected = 1 if mode == "smoke" else 6
    if len(jobs) != expected:
        raise RuntimeError(
            f"job count {len(jobs)} != expected {expected}"
        )

    return jobs


def source_record(job: ClosureJob) -> dict:
    return {
        "trajectory": {
            "path": str(job.trajectory),
            "sha256": job.trajectory_sha256,
            "bytes": job.trajectory.stat().st_size,
        },
        "topology": {
            "path": str(job.topology),
            "sha256": job.topology_sha256,
            "bytes": job.topology.stat().st_size,
        },
        "frame_ids": list(job.frame_ids),
    }


def extract_job(
    job: ClosureJob,
    model,
    device: str,
    model_record: dict,
) -> dict:

    target = job.cache_path
    manifest_path = job.manifest_path

    if target.exists() or manifest_path.exists():
        raise FileExistsError(
            f"refusing overwrite for {job.key}: "
            f"{target} / {manifest_path}"
        )

    target.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    temporary = target.with_name(
        f".{target.name}.partial.{os.getpid()}"
    )

    universe = None
    output = None
    cache = None

    started = time.perf_counter()

    try:
        (
            universe,
            atom_indices,
            atom_z,
            residue_index,
            backbone_mask,
            sidechain_mask,
            atom_records,
        ) = topology_mapping(job)

        if len(universe.trajectory) != N_FRAMES:
            raise RuntimeError(
                f"{job.key}: expected {N_FRAMES} frames, "
                f"got {len(universe.trajectory)}"
            )

        output = np.lib.format.open_memmap(
            temporary,
            mode="w+",
            dtype=np.float32,
            shape=job.output_shape,
        )

        total = len(job.frame_ids)

        for start in range(0, total, BATCH_SIZE):
            frame_ids = job.frame_ids[start:start + BATCH_SIZE]

            values = infer_bs256(
                model,
                universe,
                atom_indices,
                atom_z,
                residue_index,
                backbone_mask,
                sidechain_mask,
                frame_ids,
                device,
            )

            output[start:start + len(frame_ids)] = values
            output.flush()

            completed = start + len(frame_ids)
            elapsed = time.perf_counter() - started

            print(
                f"[{job.key}] "
                f"{completed:4d}/{total} "
                f"({100*completed/total:6.2f}%) "
                f"elapsed={elapsed:.1f}s",
                flush=True,
            )

        gly = np.asarray([aa == "G" for aa in job.sequence])

        if np.count_nonzero(output[:, gly, 128:]) != 0:
            raise RuntimeError(
                f"Gly side-chain block nonzero: {job.key}"
            )

        close_memmap(output)
        output = None

        replay_ids = tuple(
            x for x in REPLAY_FRAME_IDS
            if x in job.frame_ids
        )

        cache = np.load(
            temporary,
            mmap_mode="r",
            allow_pickle=False,
        )

        positions = [
            job.frame_ids.index(frame)
            for frame in replay_ids
        ]

        observed = np.asarray(cache[positions])

        replayed = infer_bs256(
            model,
            universe,
            atom_indices,
            atom_z,
            residue_index,
            backbone_mask,
            sidechain_mask,
            replay_ids,
            device,
        )

        replay = replay_metrics(
            observed,
            replayed,
        )
        replay["frame_ids"] = list(replay_ids)

        if not replay["passed"]:
            raise RuntimeError(
                f"deterministic replay failed: {job.key}"
            )

        close_memmap(cache)
        cache = None

        cache_hash = sha256(temporary)

        os.replace(temporary, target)

        manifest = {
            "schema":
                "pacer_dc.close_loop_20ns.bs256_cache_manifest.v1",

            "version":
                VERSION,

            "created_at":
                utc_now(),

            "mode":
                job.mode,

            "system":
                job.system,

            "replica":
                job.replica,

            "shape":
                list(job.output_shape),

            "dtype":
                "float32",

            "sequence":
                job.sequence,

            "selection":
                "chainID E and protein",

            "preprocessing": {
                "topology":
                    "CONECT-stripped copy of authenticated minimized.pdb",
                "make_whole": False,
                "center": False,
                "align": False,
                "pbc_transform": False,
            },

            "source":
                source_record(job),

            "cache": {
                "path": str(target),
                "sha256": cache_hash,
                "bytes": target.stat().st_size,
            },

            "model_provenance":
                model_record,

            "atom_mapping": {
                "expected_heavy_atoms":
                    int(len(atom_indices)),
                "backbone_atoms":
                    list(BACKBONE_ATOMS),
                "sidechain_present_residues":
                    int(
                        sum(
                            item["sidechain_present"]
                            for item in atom_records
                        )
                    ),
                "glycine_count":
                    job.sequence.count("G"),
            },

            "replay":
                replay,

            "elapsed_seconds":
                time.perf_counter() - started,

            "output_scope":
                "C1-BS256 only; no calibration, normalization, "
                "RFF, graph diffusion, contrasts, PACER-FKG "
                "scores, or biological interpretation.",
        }

        atomic_json(
            manifest_path,
            manifest,
        )

        return {
            "status": "PASS",
            "job": job.key,
            "cache": str(target),
            "cache_sha256": cache_hash,
            "shape": list(job.output_shape),
            "replay": replay,
            "elapsed_seconds": manifest["elapsed_seconds"],
        }

    except Exception:
        close_memmap(cache)
        close_memmap(output)

        if temporary.exists():
            temporary.unlink()

        raise

    finally:
        if universe is not None:
            universe.trajectory.close()


def run(
    mode: str,
    checkpoint: Path,
    device: str,
    geom2vec_source: Path,
) -> dict:

    closure_anchor = verify_closure_anchor()

    # Reuse the historical frozen encoder protection gate.
    frozen_encoder = verify_phase0_and_frozen_candidate()

    jobs = build_jobs(
        mode,
        closure_anchor,
    )

    print(
        f"Closure anchor PASS; jobs={len(jobs)}",
        flush=True,
    )

    model, model_record = load_model(
        checkpoint,
        device,
        geom2vec_source,
    )

    results = []

    for index, job in enumerate(jobs, 1):
        print(
            f"=== {mode.upper()} {index}/{len(jobs)}: "
            f"{job.key} ===",
            flush=True,
        )

        results.append(
            extract_job(
                job,
                model,
                device,
                model_record,
            )
        )

    return {
        "schema":
            "pacer_dc.close_loop_20ns.phase1_bs256_run.v1",

        "version":
            VERSION,

        "created_at":
            utc_now(),

        "status":
            "PASS",

        "mode":
            mode,

        "device":
            device,

        "batch_size":
            BATCH_SIZE,

        "closure_anchor":
            closure_anchor,

        "frozen_encoder":
            frozen_encoder,

        "jobs":
            results,

        "forbidden_downstream_operations_performed":
            [],
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    modes = parser.add_mutually_exclusive_group(
        required=True
    )

    modes.add_argument(
        "--smoke",
        action="store_true",
        help="extract frames 0-3 of LY no-probe R1 only",
    )

    modes.add_argument(
        "--full",
        action="store_true",
        help="extract all six authenticated 400-frame LY trajectories",
    )

    parser.add_argument(
        "--device",
        default="cuda",
    )

    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=DEFAULT_CHECKPOINT,
    )

    parser.add_argument(
        "--geom2vec-source",
        type=Path,
        default=Path(r"C:\projects\geom2vec-source"),
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    mode = (
        "smoke"
        if args.smoke
        else "full"
    )

    started = time.perf_counter()

    try:
        report = run(
            mode,
            args.checkpoint,
            args.device,
            args.geom2vec_source,
        )

        report["elapsed_seconds"] = (
            time.perf_counter() - started
        )

        REPORT_ROOT.mkdir(
            parents=True,
            exist_ok=True,
        )

        report_path = (
            REPORT_ROOT
            / f"RUN_{mode.upper()}_v01.json"
        )

        if report_path.exists():
            raise FileExistsError(
                f"refusing overwrite: {report_path}"
            )

        atomic_json(
            report_path,
            report,
        )

        print(
            json.dumps(
                {
                    "status": "PASS",
                    "mode": mode,
                    "report": str(report_path),
                    "jobs": len(report["jobs"]),
                },
                indent=2,
            ),
            flush=True,
        )

        return 0

    except Exception as exc:
        failure = {
            "schema":
                "pacer_dc.close_loop_20ns.phase1_bs256_failure.v1",
            "version":
                VERSION,
            "created_at":
                utc_now(),
            "status":
                "FAIL",
            "mode":
                mode,
            "error_type":
                type(exc).__name__,
            "error":
                str(exc),
            "traceback":
                traceback.format_exc(),
        }

        REPORT_ROOT.mkdir(
            parents=True,
            exist_ok=True,
        )

        failure_path = (
            REPORT_ROOT
            / f"RUN_{mode.upper()}_FAILURE_v01.json"
        )

        atomic_json(
            failure_path,
            failure,
        )

        print(
            json.dumps(
                failure,
                indent=2,
            ),
            flush=True,
        )

        return 1


if __name__ == "__main__":
    raise SystemExit(main())