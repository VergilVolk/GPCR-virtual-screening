#!/usr/bin/env python
"""Matched-window R2/R3 qualification for frozen Geom2Vec PACER-DC features."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


CONTEXTS = ("apo", "probe_only", "candidate_no_probe", "candidate_probe")


def cosine(a, b):
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / denom) if denom else float("nan")


def load_regions(path, n_residues):
    raw = json.loads(path.read_text(encoding="utf-8"))["regions"]
    out = {name: np.asarray([x["embedding_index"] for x in rows], int)
           for name, rows in raw.items() if rows}
    out["global"] = np.arange(n_residues)
    return out


def locate(root, replica, window, context):
    pattern = f"replica_{replica:02d}/window_{window:03d}/atom14/{context}_w{window:03d}.geom2vec.npz"
    path = root / pattern
    if not path.exists():
        raise FileNotFoundError(path)
    with np.load(path) as z:
        return z["residue_features"].mean(axis=0).astype(np.float64)


def bootstrap_median(values, seed=20260927, draws=10000):
    values = np.asarray(values, float)
    rng = np.random.default_rng(seed)
    sampled = rng.choice(values, size=(draws, len(values)), replace=True)
    medians = np.median(sampled, axis=1)
    return [float(x) for x in np.quantile(medians, [0.025, 0.975])]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-root", type=Path, required=True)
    p.add_argument("--region-map", type=Path, required=True)
    p.add_argument("--replicas", nargs=2, type=int, default=(2, 3))
    p.add_argument("--windows", nargs="+", type=int, default=range(5))
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    probe = locate(args.embedding_root, args.replicas[0], args.windows[0], "apo")
    regions = load_regions(args.region_map, probe.shape[0])
    vectors = {}
    for replica in args.replicas:
        for window in args.windows:
            z = {c: locate(args.embedding_root, replica, window, c) for c in CONTEXTS}
            for region, idx in regions.items():
                pooled = {c: z[c][idx].mean(axis=0) for c in CONTEXTS}
                vectors[(replica, window, region, "dAGO")] = pooled["candidate_no_probe"] - pooled["apo"]
                vectors[(replica, window, region, "dINT")] = (
                    pooled["candidate_probe"] - pooled["probe_only"]
                    - pooled["candidate_no_probe"] + pooled["apo"]
                )

    rows = []
    for region in regions:
        for axis in ("dAGO", "dINT"):
            cosines = [cosine(vectors[(args.replicas[0], w, region, axis)],
                              vectors[(args.replicas[1], w, region, axis)])
                       for w in args.windows]
            ci = bootstrap_median(cosines)
            rows.append({
                "region": region,
                "axis": axis,
                "matched_window_cosines": cosines,
                "median_cosine": float(np.median(cosines)),
                "bootstrap_window_ci95": ci,
                "positive_windows": int(np.sum(np.asarray(cosines) > 0)),
                "pilot_gate": bool(np.median(cosines) > 0.5 and ci[0] > 0),
                "independence_warning": "Windows are contiguous augmentations, not independent replicas.",
            })

    functional = {"ACh_pocket", "compound110_pocket", "ECV_anchor", "W435_gate", "species_probe"}
    functional_pass = [r for r in rows if r["region"] in functional and r["axis"] == "dINT" and r["pilot_gate"]]
    distal = next(r for r in rows if r["region"] == "distal_control" and r["axis"] == "dINT")
    report = {
        "method": "matched-window Geom2Vec R2/R3 four-context audit",
        "evidence_level": "two_replica_encoder_qualification",
        "replicas": args.replicas,
        "windows": args.windows,
        "rows": rows,
        "decision": {
            "pilot_pass": bool(functional_pass and not distal["pilot_gate"]),
            "functional_regions_passing_dINT": [r["region"] for r in functional_pass],
            "distal_dINT_gate": distal["pilot_gate"],
            "full_replacement_authorized": False,
            "reason": "A second pharmacological control (LY2119620) is still required even if the two-replica pilot passes.",
        },
        "claim_boundary": "Two-replica matched-window qualification; not a potency or PAM efficacy benchmark.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report["decision"], indent=2))


if __name__ == "__main__":
    main()
