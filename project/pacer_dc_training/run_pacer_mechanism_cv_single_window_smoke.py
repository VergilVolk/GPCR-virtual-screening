#!/usr/bin/env python3
"""Single-window real-data smoke test for PACER-MCV (never a replica gate)."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from pacer_mechanism_cv import CONTEXTS, extract_frame_features, factorial_contrast, load_config, load_mapping, load_sequence


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--atom14-dir", type=Path, required=True)
    p.add_argument("--mapping", type=Path, required=True)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--window", type=int, default=0)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    cfg, mapping = load_config(args.config), load_mapping(args.mapping)
    means, provenance, names = {}, [], None
    for context in CONTEXTS:
        stem = f"{context}_w{args.window:03d}"
        npy, seq_path = args.atom14_dir / f"{stem}.npy", args.atom14_dir / f"{stem}.csv"
        arr = np.load(npy, mmap_mode="r"); sequence = load_sequence(seq_path)
        local_names, trace = extract_frame_features(np.asarray(arr), sequence, mapping, cfg)
        if names is None: names = local_names
        if names != local_names: raise ValueError("feature order mismatch")
        means[context] = trace.mean(axis=0)
        provenance.append({"context": context, "path": str(npy), "sha256": sha256(npy), "frames": int(arr.shape[0])})
    rows, axes = [], {}
    for axis, signs in cfg["axes"].items():
        vector = factorial_contrast(means, signs); axes[axis] = float(np.linalg.norm(vector))
        for name, value in zip(names, vector): rows.append({"axis": axis, "feature": name, "raw_contrast": float(value)})
    args.output.mkdir(parents=True)
    with (args.output / "PACER_MCV_R1_FEATURES.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    report = {"status": "PASS_SOFTWARE_SMOKE", "evidence_level": "single_replica_single_window_smoke",
              "feature_count": len(names), "axis_raw_norms": axes, "provenance": provenance,
              "claim_boundary": "Execution and finite physical features only; no replication, significance, synergy or PAM claim."}
    (args.output / "PACER_MCV_R1_SMOKE_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__": main()
