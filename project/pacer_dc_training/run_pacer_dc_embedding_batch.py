#!/usr/bin/env python3
"""Run audited PACER-DC atom14 and frozen OneProt-MD stages in batches.

Run from the Windows canonical repository root in the ``oneprot-md`` environment.
This wrapper never runs production MD, edits raw trajectory windows, or trains a model.
It delegates preprocessing and inference to the existing, versioned extractor.

Examples:
    python project/pacer_dc_training/run_pacer_dc_embedding_batch.py --candidate compound110 --replicas 3
    python project/pacer_dc_training/run_pacer_dc_embedding_batch.py --candidate compound110 --replicas 2 3 --dry-run
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

CONTEXTS = ("candidate_probe", "candidate_no_probe", "probe_only", "apo")


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def check_raw(folder: Path, replica: int, window: int) -> dict:
    manifest = folder / "trajectory_windows.json"
    if not manifest.is_file():
        raise RuntimeError(f"Missing raw-window manifest: {manifest}")
    rows = json.loads(manifest.read_text(encoding="utf-8"))["rows"]
    if len(rows) != 4 or {r["context"] for r in rows} != set(CONTEXTS):
        raise RuntimeError(f"Four contexts missing or duplicated: {manifest}")
    for row in rows:
        c = row["context"]
        if not (
            row["status"] == "ok"
            and row["replica_id"] == replica
            and row["window_id"] == window
            and row["start_frame"] == window * 100
            and row["end_frame"] == window * 100 + 99
            and row["frame_count"] == 100
            and row["receptor_residues"] == 270
            and row["receptor_atoms"] == 4349
        ):
            raise RuntimeError(f"Bad raw-window manifest row: {replica=} {window=} {c=}")
        raw_path = folder / row["raw_path"]
        if not raw_path.is_file():
            raise RuntimeError(f"Missing raw NPZ: {raw_path}")
        with np.load(raw_path, allow_pickle=False) as raw:
            coords = raw["coords"]
            if coords.shape != (100, 4349, 3) or not np.isfinite(coords).all():
                raise RuntimeError(f"Invalid raw coordinates: {raw_path}")
        if replica == 3 and c == "apo" and not row["trajectory_path"].replace("\\", "/").endswith(
            "apo/replica_03/trajectory_corrected.dcd"
        ):
            raise RuntimeError("R3 apo must come from corrected DCD, not the original")
    return {r["context"]: r for r in rows}


def stage_paths(folder: Path, window: int, stage: str) -> list[Path]:
    if stage == "atom14":
        return [folder / "atom14" / f"{c}_w{window:03d}{suffix}"
                for c in CONTEXTS for suffix in (".npy", ".csv")]
    return ([folder / "embeddings" / f"{c}_w{window:03d}.npy" for c in CONTEXTS]
            + [folder / "embeddings" / f"dPAM_w{window:03d}.npy",
               folder / "unit_summary.json", folder / "four_context_long.csv",
               folder / "four_context_differential.csv"])


def validate_atom14(folder: Path, window: int) -> None:
    for c in CONTEXTS:
        arr = np.load(folder / "atom14" / f"{c}_w{window:03d}.npy", mmap_mode="r")
        if arr.shape != (100, 270, 14, 3):
            raise RuntimeError(f"Unexpected atom14 shape for {c}: {arr.shape}")
        if not (folder / "atom14" / f"{c}_w{window:03d}.csv").is_file():
            raise RuntimeError(f"Missing seqres CSV for {c}")


def validate_embed(folder: Path, window: int, replica: int) -> dict:
    zs = {}
    for c in CONTEXTS:
        path = folder / "embeddings" / f"{c}_w{window:03d}.npy"
        z = np.load(path, allow_pickle=False)
        if z.shape != (1024,) or not np.isfinite(z).all() or not np.isclose(np.linalg.norm(z), 1.0, atol=1e-4):
            raise RuntimeError(f"Invalid embedding: {path}")
        zs[c] = z
    d = np.load(folder / "embeddings" / f"dPAM_w{window:03d}.npy", allow_pickle=False)
    if d.shape != (1024,) or not np.isfinite(d).all() or not np.allclose(
        d, zs["candidate_probe"] - zs["probe_only"], rtol=1e-6, atol=1e-7
    ):
        raise RuntimeError(f"Invalid dPAM for R{replica} W{window:03d}")
    summary = json.loads((folder / "unit_summary.json").read_text(encoding="utf-8"))
    if (summary.get("completeness") != "complete" or
            summary.get("replica_id") != replica or
            summary.get("window_id") != window or
            set(summary.get("contexts_embedded", [])) != set(CONTEXTS)):
        raise RuntimeError(f"Incomplete or mismatched unit summary: R{replica} W{window:03d}")
    return {"replica": replica, "window": window,
            "dpam_l2": float(np.linalg.norm(d)),
            "checkpoint_sha256": summary.get("checkpoint_sha256", "")}


def run_stage(extractor: Path, stage: str, candidate: str, replica: int,
              window: int, folder: Path, logdir: Path) -> None:
    command = [sys.executable, str(extractor), "--stage", stage,
               "--candidate", candidate, "--replica", str(replica),
               "--window", str(window), "--window-frames", "100",
               "--contexts", ",".join(CONTEXTS), "--out-dir", str(folder)]
    logdir.mkdir(parents=True, exist_ok=True)
    log = logdir / f"R{replica}_W{window:03d}_{stage}.log"
    print("RUN:", " ".join(command), flush=True)
    with log.open("w", encoding="utf-8") as out:
        out.write("COMMAND: " + " ".join(command) + "\n")
        out.write("START_UTC: " + datetime.now(timezone.utc).isoformat() + "\n")
        out.flush()
        proc = subprocess.run(command, stdout=out, stderr=subprocess.STDOUT, check=False)
        out.write(f"\nEXIT_CODE: {proc.returncode}\n")
    if proc.returncode != 0:
        raise RuntimeError(f"{stage} failed. See {log}")
    print("STAGE_DONE:", log)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidate", default="compound110")
    parser.add_argument("--replicas", type=int, nargs="+", default=[3])
    parser.add_argument("--windows", type=int, nargs="+", default=list(range(5)))
    parser.add_argument("--results-root", type=Path,
                        default=Path("project/results/pacer_dc_four_context_v01"))
    parser.add_argument("--extractor", type=Path,
                        default=Path("project/pacer_dc_training/extract_pacer_dc_four_context_embeddings.py"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if len(set(args.replicas)) != len(args.replicas) or len(set(args.windows)) != len(args.windows):
        parser.error("duplicate replicas/windows are not supported")
    if any(r < 1 for r in args.replicas) or any(w not in range(5) for w in args.windows):
        parser.error("replicas must be positive; windows must be within 0..4")
    if not args.extractor.is_file():
        parser.error(f"Extractor not found: {args.extractor}")

    root = args.results_root / args.candidate
    logdir = root / "batch_logs"
    complete = []
    for replica in args.replicas:
        for window in args.windows:
            folder = root / f"replica_{replica:02d}" / f"window_{window:03d}"
            check_raw(folder, replica, window)
            for stage in ("atom14", "embed"):
                paths = stage_paths(folder, window, stage)
                existing = [p.exists() for p in paths]
                if all(existing):
                    if stage == "atom14":
                        validate_atom14(folder, window)
                    else:
                        result = validate_embed(folder, window, replica)
                        complete.append(result)
                    print(f"SKIP_VERIFIED: R{replica} W{window:03d} {stage}", flush=True)
                    continue
                if any(existing):
                    raise RuntimeError(
                        f"Partial {stage} outputs in {folder}. Stop to audit; "
                        "no automatic overwriting is allowed."
                    )
                if args.dry_run:
                    print(f"WOULD_RUN: R{replica} W{window:03d} {stage}", flush=True)
                    continue
                run_stage(args.extractor, stage, args.candidate, replica, window, folder, logdir)
                if not all(p.is_file() for p in paths):
                    raise RuntimeError(f"Missing outputs after {stage}: {folder}")
                if stage == "atom14":
                    validate_atom14(folder, window)
                else:
                    result = validate_embed(folder, window, replica)
                    complete.append(result)
                print(f"PASS: R{replica} W{window:03d} {stage}", flush=True)

    if args.dry_run:
        print("DRY_RUN_COMPLETE: no extraction performed")
        return
    # This is a batch report, not evidence of statistical independence or biological activity.
    checksums = {r["checkpoint_sha256"] for r in complete if r["checkpoint_sha256"]}
    if len(checksums) > 1:
        raise RuntimeError(f"Checkpoint mismatch among completed units: {sorted(checksums)}")
    logdir.mkdir(parents=True, exist_ok=True)
    report = logdir / (f"{args.candidate}_R{'-'.join(map(str,args.replicas))}_batch_report.json")
    payload = {"candidate": args.candidate, "replicas": args.replicas,
               "windows": args.windows, "units_verified": complete,
               "checkpoint_sha256": next(iter(checksums), ""),
               "generated_utc": datetime.now(timezone.utc).isoformat(),
               "note": "Descriptive frozen-encoder output. Does not establish biological mechanism or training eligibility."}
    report.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"ALL_REQUESTED_WINDOWS_PASS: {len(complete)}")
    print(f"REPORT: {report}")


if __name__ == "__main__":
    main()
