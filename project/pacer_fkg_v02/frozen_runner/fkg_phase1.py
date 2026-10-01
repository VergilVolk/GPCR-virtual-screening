"""Generic Phase-1 driver: frozen C1-BS256 extraction for any candidate.

All encoder numerics come from the frozen `run_phase1_bs256` module.  This file
only supplies job construction from a run spec and the receipt/manifest plumbing.
"""

from __future__ import annotations

import csv
import os
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from project.pacer_fkg_v02.frozen_runner import engines, paths as pathlib_policy
from project.pacer_fkg_v02.frozen_runner.io_utils import check_npy
from project.pacer_fkg_v02.frozen_runner.spec import RunSpec, sha256_file

VERSION = "PACER_FKG_V02_FROZEN_RUNNER_PHASE1_v01"
MAPPING_RELPATH = "project/results/pacer_dc_four_context_v01/compound110/G2_RESIDUE_MAPPING_v01.csv"
DEFAULT_CHECKPOINT_RELPATH = "project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth"


def _now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat()


def _atomic_json(path: Path, payload: Any) -> None:
    tmp = path.with_name(f".{path.name}.partial.{os.getpid()}")
    tmp.write_text(__import__("json").dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def frozen_sequence(repo_root: Path, phase1) -> str:
    """Reconstruct the 270-residue receptor sequence from the frozen mapping."""
    mapping = repo_root / MAPPING_RELPATH
    with mapping.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != phase1.N_RESIDUES:
        raise RuntimeError(f"frozen residue mapping count mismatch: {len(rows)} != {phase1.N_RESIDUES}")
    seq = []
    for row in rows:
        aa = phase1.RESTYPE_3TO1.get(str(row["resname"]).strip())
        if aa is None:
            raise RuntimeError(f"unsupported residue {row['resname']}")
        seq.append(aa)
    return "".join(seq)


@dataclass(frozen=True)
class CandidateJob:
    """Duck-types the frozen `run_phase1_bs256.Job` attribute surface."""

    system: str
    replica: int
    topology: Path
    trajectory: Path
    trajectory_sha256: str
    topology_sha256: str
    sequence: str
    n_frames: int

    @property
    def key(self) -> str:
        return f"{self.system}__replica_{self.replica:02d}"

    @property
    def frame_ids(self) -> tuple[int, ...]:
        return tuple(range(self.n_frames))


def build_jobs(spec: RunSpec, sequence: str, repo_root: Path) -> list[CandidateJob]:
    jobs: list[CandidateJob] = []
    for system in spec.candidate_systems:
        topology = spec.paths.topology_root / f"{system}_minimized_no_conect.pdb"
        if not topology.is_file():
            raise FileNotFoundError(f"sanitized topology missing for {system}: {topology}")
        for replica in spec.replicas:
            directory = spec.paths.production_root / system / f"replica_{replica:02d}"
            trajectory = directory / "trajectory.dcd"
            progress = directory / "progress.json"
            for required in (trajectory, progress):
                if not required.is_file():
                    raise FileNotFoundError(required)
            payload = __import__("json").loads(progress.read_text(encoding="utf-8-sig"))
            if payload.get("status") != "complete" or float(payload.get("completed_ns", -1)) != spec.n_frames / 20.0:
                raise RuntimeError(f"incomplete production trajectory {system} R{replica}")
            if int(payload.get("seed", -1)) != spec.seeds[replica]:
                raise RuntimeError(f"seed mismatch {system} R{replica}: {payload.get('seed')} != {spec.seeds[replica]}")
            jobs.append(CandidateJob(system, replica, topology, trajectory,
                                     sha256_file(trajectory), sha256_file(topology), sequence, spec.n_frames))
    return jobs


def run(spec: RunSpec, *, repo_root: Path, device: str = "cuda",
        checkpoint: Path | None = None, geom2vec_source: Path | None = None,
        mapping_override: Path | None = None) -> dict[str, Any]:
    phase1 = engines.load_phase1_module()
    engines.assert_frozen_numerics()

    root = repo_root.resolve()
    receipt_path = pathlib_policy.guard_new_output_noclobber(
        spec.paths.report_root / "phase1" / "PHASE1_FREEZE_RECEIPT_v01.json", root, spec.paths.allowed_output_roots)

    sequence = frozen_sequence(root, phase1)
    jobs = build_jobs(spec, sequence, root)
    expected = len(spec.candidate_systems) * len(spec.replicas)
    if len(jobs) != expected:
        raise RuntimeError(f"built {len(jobs)} jobs, expected {expected}")

    ckpt = Path(checkpoint) if checkpoint else (root / DEFAULT_CHECKPOINT_RELPATH)
    source = Path(geom2vec_source) if geom2vec_source else Path(r"C:\projects\geom2vec-source")
    model, model_record = phase1.load_model(ckpt, device, source)

    records = []
    for index, job in enumerate(jobs, 1):
        cache = pathlib_policy.guard_new_output_noclobber(
            spec.paths.phase1_cache_root / job.system / f"replica_{job.replica:02d}" / "C1_BS256.npy",
            root, spec.paths.allowed_output_roots)
        manifest_path = pathlib_policy.guard_new_output_noclobber(
            spec.paths.report_root / "phase1" / "manifests" / f"{job.key}.json", root, spec.paths.allowed_output_roots)
        cache.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        tmp = cache.with_name(f".{cache.name}.partial.{os.getpid()}")
        universe = out = loaded = None
        try:
            universe, atom_indices, atom_z, residue_index, bb, sc, mapping_records = phase1.topology_mapping(job)
            out = np.lib.format.open_memmap(tmp, mode="w+", dtype=np.float32,
                                            shape=(job.n_frames, phase1.N_RESIDUES, phase1.BS_WIDTH))
            for start in range(0, job.n_frames, phase1.BATCH_SIZE):
                ids = job.frame_ids[start:start + phase1.BATCH_SIZE]
                out[start:start + len(ids)] = phase1.infer_bs256(
                    model, universe, atom_indices, atom_z, residue_index, bb, sc, ids, device)
                out.flush()
            glycine = np.asarray([aa == "G" for aa in job.sequence])
            if np.count_nonzero(out[:, glycine, 128:]) != 0:
                raise RuntimeError("Gly side-chain output is not exactly zero")
            out._mmap.close(); out = None
            loaded = np.load(tmp, mmap_mode="r", allow_pickle=False)
            observed = np.asarray(loaded[list(phase1.REPLAY_FRAME_IDS)])
            replayed = phase1.infer_bs256(model, universe, atom_indices, atom_z, residue_index, bb, sc,
                                          phase1.REPLAY_FRAME_IDS, device)
            replay = phase1.replay_metrics(observed, replayed)
            replay["frame_ids"] = list(phase1.REPLAY_FRAME_IDS)
            if not replay["passed"]:
                raise RuntimeError(f"deterministic replay failed for {job.key}")
            loaded._mmap.close(); loaded = None
            digest = sha256_file(tmp)
            os.replace(tmp, cache)
            manifest = {
                "schema": "pacer.final.frozen_runner.phase1_bs256.v1", "version": VERSION,
                "created_at": _now(), "run_id": spec.run_id, "candidate_id": spec.molecule.candidate_id,
                "system": job.system, "replica": job.replica, "seed": spec.seeds[job.replica],
                "shape": [job.n_frames, phase1.N_RESIDUES, phase1.BS_WIDTH], "dtype": "float32",
                "sequence": job.sequence, "selection": phase1.SELECTION,
                "preprocessing": {"make_whole": False, "center": False, "align": False, "pbc_transform": False},
                "historical_fkg_v02_freeze_sha256": engines.REQUIRED_FREEZE_SHA256,
                "trajectory": {"path": str(job.trajectory), "sha256": job.trajectory_sha256,
                               "bytes": job.trajectory.stat().st_size},
                "topology": {"path": str(job.topology), "sha256": job.topology_sha256,
                             "bytes": job.topology.stat().st_size},
                "cache": {"path": str(cache), "sha256": digest, "bytes": cache.stat().st_size},
                "model_provenance": model_record,
                "atom_mapping": {"heavy_atoms": int(len(atom_indices)),
                                 "backbone_atoms": list(phase1.BACKBONE_ATOMS),
                                 "sidechain_present_residues": int(sum(r["sidechain_present"] for r in mapping_records)),
                                 "glycine_count": job.sequence.count("G")},
                "replay": replay, "elapsed_seconds": time.perf_counter() - started,
                "output_scope": "C1-BS256 only; no calibration and no biological metrics",
            }
            _atomic_json(manifest_path, manifest)
            records.append(manifest)
        finally:
            for handle in (loaded, out):
                mmap = getattr(handle, "_mmap", None) if handle is not None else None
                if mmap is not None:
                    mmap.close()
            if tmp.exists():
                tmp.unlink()
            if universe is not None:
                universe.trajectory.close()

    receipt = {
        "schema": "pacer.final.frozen_runner.phase1_freeze_receipt.v1", "version": VERSION,
        "created_at": _now(), "status": "PHASE1_FROZEN", "run_id": spec.run_id,
        "candidate_id": spec.molecule.candidate_id,
        "historical_fkg_v02_freeze_sha256": engines.REQUIRED_FREEZE_SHA256,
        "run_spec": {"path": str(spec.source_path), "sha256": sha256_file(spec.source_path)},
        "geometry": {"n_frames": spec.n_frames, "block_frames": spec.block_frames,
                     "blocks_per_trajectory": spec.blocks_per_trajectory},
        "seeds": {str(k): v for k, v in spec.seeds.items()},
        "cache_count": len(records),
        "caches": [{"system": m["system"], "replica": m["replica"], "trajectory": m["trajectory"],
                    "topology": m["topology"], "cache": m["cache"], "replay": m["replay"]} for m in records],
        "forbidden_operations_performed": {"calibration": False, "normalization_refit": False,
                                           "bandwidth_refit": False, "rff_refit": False, "graph_refit": False,
                                           "region_refit": False, "threshold_selection": False,
                                           "outcome_driven_tuning": False},
    }
    _atomic_json(receipt_path, receipt)
    return {"status": "PHASE1_FROZEN", "receipt": str(receipt_path),
            "receipt_sha256": sha256_file(receipt_path), "cache_count": len(records)}


def verify(spec: RunSpec, *, repo_root: Path) -> dict[str, Any]:
    import json

    phase1 = engines.load_phase1_module()
    engines.assert_frozen_numerics()
    root = repo_root.resolve()
    receipt_path = spec.paths.report_root / "phase1" / "PHASE1_FREEZE_RECEIPT_v01.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8-sig"))
    if receipt.get("status") != "PHASE1_FROZEN":
        raise RuntimeError("invalid Phase-1 receipt status")
    expected = len(spec.candidate_systems) * len(spec.replicas)
    if int(receipt.get("cache_count", -1)) != expected:
        raise RuntimeError(f"Phase-1 receipt cache_count {receipt.get('cache_count')} != {expected}")
    for item in receipt["caches"]:
        cache = Path(item["cache"]["path"])
        if not cache.is_file() or sha256_file(cache) != item["cache"]["sha256"]:
            raise RuntimeError(f"Phase-1 cache verification failed: {cache}")
        try:
            check_npy(cache, (spec.n_frames, phase1.N_RESIDUES, phase1.BS_WIDTH))
        except ValueError as exc:
            raise RuntimeError(f"Phase-1 cache contract failed: {cache} ({exc})") from exc
    return {"status": "PHASE1_VERIFIED", "receipt_sha256": sha256_file(receipt_path), "cache_count": expected}
