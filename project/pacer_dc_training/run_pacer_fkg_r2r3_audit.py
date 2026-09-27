#!/usr/bin/env python
"""Run the frozen PACER-FKG analysis over an R2/R3 Geom2Vec archive.

This wrapper exists so the GPU teammate can evaluate the already-extracted
40 NPZ files with one command.  It does not rerun MD or the encoder.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np


AXES = ("synergy_interaction", "intrinsic_agonism", "conditional_pam_effect")
EXCLUDED_REGIONS = {"distal_control", "stable_core_control"}


def rankdata(values):
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks


def spearman(a, b):
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(rankdata(np.asarray(a)), rankdata(np.asarray(b)))[0, 1])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--graph", type=Path, required=True)
    p.add_argument("--output-root", type=Path, required=True)
    p.add_argument("--replicas", nargs="+", type=int, default=[2, 3])
    p.add_argument("--windows", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    args = p.parse_args()

    analyzer = Path(__file__).with_name("pacer_factorial_kernel_graph.py")
    args.output_root.mkdir(parents=True, exist_ok=True)
    records = []
    for replica in args.replicas:
        for window in args.windows:
            embedding_root = args.input_root / f"replica_{replica:02d}" / f"window_{window:03d}" / "atom14"
            expected = [embedding_root / f"{c}_w{window:03d}.geom2vec.npz" for c in
                        ("apo", "probe_only", "candidate_no_probe", "candidate_probe")]
            missing = [str(path) for path in expected if not path.exists()]
            if missing:
                raise FileNotFoundError("missing Geom2Vec inputs:\n" + "\n".join(missing))
            output = args.output_root / f"R{replica}_W{window:03d}.fkg.json"
            command = [
                sys.executable, str(analyzer),
                "--embedding-root", str(embedding_root),
                "--pattern", f"{{context}}_w{window:03d}.geom2vec.npz",
                "--graph", str(args.graph),
                "--output", str(output),
                "--preprocess", "global_residual_channel_robust",
            ]
            subprocess.run(command, check=True)
            result = json.loads(output.read_text(encoding="utf-8"))
            for row in result["region_scores"]:
                records.append({"replica": replica, "window": window, **row})

    stable = {(r["replica"], r["window"], axis): r[f"{axis}_kernel_u2"]
              for r in records if r["region"] == "stable_core_control" for axis in AXES}
    summaries = []
    regions = sorted({r["region"] for r in records})
    for region in regions:
        for axis in AXES:
            by_replica = {}
            for replica in args.replicas:
                rows = sorted((r for r in records if r["region"] == region and r["replica"] == replica),
                              key=lambda r: r["window"])
                deltas = [r[f"{axis}_kernel_u2"] - stable[(replica, r["window"], axis)] for r in rows]
                by_replica[str(replica)] = {
                    "window_u2_minus_stable": deltas,
                    "median_u2_minus_stable": float(np.median(deltas)),
                    "positive_windows": int(np.sum(np.asarray(deltas) > 0)),
                }
            first = by_replica[str(args.replicas[0])]["window_u2_minus_stable"]
            second = by_replica[str(args.replicas[1])]["window_u2_minus_stable"] if len(args.replicas) == 2 else []
            rho = spearman(first, second) if second else None
            gate = (
                region not in EXCLUDED_REGIONS
                and all(v["median_u2_minus_stable"] > 0 for v in by_replica.values())
                and all(v["positive_windows"] >= 4 for v in by_replica.values())
                and rho is not None and rho >= 0
            )
            summaries.append({
                "region": region, "axis": axis, "replicas": by_replica,
                "matched_window_spearman": rho,
                "preregistered_specificity_gate": gate,
            })

    report = {
        "method": "PACER-FKG R2/R3 frozen batch audit v01",
        "evidence_level": "two_replica_method_qualification",
        "input_root": str(args.input_root),
        "graph": str(args.graph),
        "replicas": args.replicas,
        "windows": args.windows,
        "window_results": records,
        "cross_replica_summary": summaries,
        "passing_synergy_regions": [s["region"] for s in summaries
                                    if s["axis"] == "synergy_interaction" and s["preregistered_specificity_gate"]],
        "claim_boundary": (
            "Two replicas and contiguous windows qualify an estimator only. "
            "They do not establish PAM efficacy, potency, convergence or generalization."
        ),
    }
    audit = args.output_root / "PACER_FKG_R2R3_AUDIT.json"
    audit.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({
        "audit": str(audit),
        "passing_synergy_regions": report["passing_synergy_regions"],
        "windows_completed": len(args.replicas) * len(args.windows),
    }, indent=2))


if __name__ == "__main__":
    main()
