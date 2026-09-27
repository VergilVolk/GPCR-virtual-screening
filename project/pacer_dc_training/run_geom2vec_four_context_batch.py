#!/usr/bin/env python
"""Batch frozen Geom2Vec extraction for existing PACER-DC atom14 windows."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from extract_geom2vec_atom14 import (
    create_model_from_checkpoint,
    embed,
    flatten_atom14,
    install_cpu_radius_graph_fallback,
    read_sequence,
    sha256,
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True,
                   help="candidate root containing replica_*/window_*/atom14")
    p.add_argument("--output-root", type=Path, required=True)
    p.add_argument("--checkpoint", type=Path, required=True)
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--batch-size", type=int, default=1)
    p.add_argument("--device", default="cuda")
    p.add_argument("--replicas", nargs="*", type=int)
    p.add_argument("--limit-files", type=int)
    args = p.parse_args()

    files = sorted(args.input_root.glob("replica_*/window_*/atom14/*_w*.npy"))
    if args.replicas:
        allowed = {f"replica_{r:02d}" for r in args.replicas}
        files = [path for path in files if any(part in allowed for part in path.parts)]
    if args.limit_files:
        files = files[: args.limit_files]
    if not files:
        raise SystemExit(f"no atom14 windows under {args.input_root}")

    fallback = install_cpu_radius_graph_fallback()
    model = create_model_from_checkpoint(str(args.checkpoint), device=args.device)
    rows = []
    started_all = time.time()
    for number, atom_path in enumerate(files, 1):
        rel = atom_path.relative_to(args.input_root)
        output = args.output_root / rel.with_suffix(".geom2vec.npz")
        output.parent.mkdir(parents=True, exist_ok=True)
        sequence = read_sequence(atom_path.with_suffix(".csv"))
        coords = np.load(atom_path, mmap_mode="r")
        frame_ids = np.arange(0, coords.shape[0], args.stride, dtype=np.int64)
        xyz, atom_z, ridx = flatten_atom14(np.asarray(coords[frame_ids]), sequence)
        started = time.time()
        features = embed(model, xyz, atom_z, ridx, len(sequence), args.batch_size, args.device)
        elapsed = time.time() - started
        np.savez_compressed(
            output,
            residue_features=features,
            global_features=features.mean(axis=1),
            frame_ids=frame_ids,
            sequence=np.asarray(sequence),
        )
        row = {
            "input": str(rel).replace("\\", "/"),
            "output": str(output.relative_to(args.output_root)).replace("\\", "/"),
            "input_sha256": sha256(atom_path),
            "frames": int(len(frame_ids)),
            "shape": list(features.shape),
            "finite": bool(np.isfinite(features).all()),
            "seconds": round(elapsed, 3),
            "status": "ok",
        }
        rows.append(row)
        print(f"[{number}/{len(files)}] {rel}: {features.shape}, {elapsed:.2f}s", flush=True)

    audit = {
        "method": "frozen Geom2Vec ViSNet batch extraction",
        "evidence_level": "encoder_feasibility_pilot",
        "checkpoint_sha256": sha256(args.checkpoint),
        "device": args.device,
        "stride": args.stride,
        "batch_size": args.batch_size,
        "cpu_radius_graph_fallback": fallback,
        "files_total": len(rows),
        "files_ok": sum(row["status"] == "ok" for row in rows),
        "elapsed_seconds": round(time.time() - started_all, 3),
        "rows": rows,
        "claim_boundary": "Embedding extraction only; qualification requires independent-replica audit.",
    }
    args.output_root.mkdir(parents=True, exist_ok=True)
    (args.output_root / "batch_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
