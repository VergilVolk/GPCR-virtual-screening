#!/usr/bin/env python
"""Descriptive OneProt-MD and differential stability across replicas/windows.

Read-only over existing unit directories. No training, no AUC, no PAM score.
PCA/tICA/VAMP are intentionally excluded until a train-only shared basis exists.

    python audit_pacer_dc_four_context_stability.py \
      --units <unit_dir> [<unit_dir> ...] \
      --output project/results/pacer_dc_four_context_v01/stability_report.json
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

def cos(a, b) -> float:
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--units", nargs="+", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    records = {}
    differentials = {}
    for unit in args.units:
        long_csv = unit / "four_context_long.csv"
        if not long_csv.is_file():
            continue
        with long_csv.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                if row["status"] != "ok" or not row["embedding_path"]:
                    continue
                key = f"{row['candidate_id']}/r{row['replica_id']}/w{int(row['start_frame']):04d}/{row['context']}"
                emb = np.load(unit / row["embedding_path"])
                records[key] = {
                    "unit": str(unit), "context": row["context"],
                    "replica": row["replica_id"], "start_frame": row["start_frame"],
                    "candidate": row["candidate_id"], "one_prot": emb,
                    "embedding_norm": float(np.linalg.norm(emb)),
                    "finite": bool(np.isfinite(emb).all()),
                    "nonzero": int((emb != 0).sum()),
                    "embedding_dim": int(emb.size),
                }
        diff_csv = unit / "four_context_differential.csv"
        if diff_csv.is_file():
            with diff_csv.open(encoding="utf-8") as fh:
                row = next(csv.DictReader(fh), None)
            if row and row["completeness"] == "complete":
                key = (f"{row['candidate_id']}/r{row['replica_id']}/"
                       f"w{int(row['start_frame']):04d}")
                differentials[key] = {
                    "candidate": row["candidate_id"],
                    "replica": row["replica_id"],
                    "start_frame": row["start_frame"],
                    "dPAM": np.load(unit / row["d_PAM_path"]),
                    "dAGO": np.load(unit / row["d_AGO_path"]),
                }

    print("=== available context embeddings (descriptive) ===")
    for k in sorted(records):
        r = records[k]
        print(f"  {k:<46} z={r['embedding_dim']} norm={r['embedding_norm']:.6f} "
              f"finite={r['finite']} nonzero={r['nonzero']}")

    # pairs: same context, same window, different replica  -> replica-to-replica
    #        same context, same replica, different window  -> window-to-window reference
    pairs = []
    keys = sorted(records)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            ra, rb = records[a], records[b]
            if ra["context"] != rb["context"]:
                continue
            same_window = ra["start_frame"] == rb["start_frame"]
            same_replica = ra["replica"] == rb["replica"]
            kind = ("replica-to-replica" if same_window and not same_replica else
                    "window-to-window" if same_replica and not same_window else None)
            if kind is None:
                continue
            pairs.append({
                "kind": kind, "context": ra["context"], "a": a, "b": b,
                "one_prot_cosine": round(cos(ra["one_prot"], rb["one_prot"]), 6),
                "one_prot_l2": round(float(np.linalg.norm(ra["one_prot"] - rb["one_prot"])), 6),
            })

    print("\n=== descriptive stability (cosine of the representation vectors) ===")
    hdr = f"  {'kind':<20}{'context':<10}{'pair':<44}{'OneProt':<10}"
    print(hdr)
    for p in sorted(pairs, key=lambda x: (x["kind"], x["a"])):
        print(f"  {p['kind']:<20}{p['context']:<10}{p['a'] + ' vs ' + p['b']:<44}"
              f"{p['one_prot_cosine']:<10}")

    diff_pairs = []
    diff_keys = sorted(differentials)
    for i, a in enumerate(diff_keys):
        for b in diff_keys[i + 1:]:
            ra, rb = differentials[a], differentials[b]
            if ra["candidate"] != rb["candidate"]:
                continue
            same_window = ra["start_frame"] == rb["start_frame"]
            same_replica = ra["replica"] == rb["replica"]
            kind = ("replica-to-replica" if same_window and not same_replica else
                    "window-to-window" if same_replica and not same_window else None)
            if kind:
                diff_pairs.append({
                    "kind": kind, "a": a, "b": b,
                    "dPAM_cosine": round(cos(ra["dPAM"], rb["dPAM"]), 6),
                    "dAGO_cosine": round(cos(ra["dAGO"], rb["dAGO"]), 6),
                    "dPAM_l2": round(float(np.linalg.norm(ra["dPAM"] - rb["dPAM"])), 6),
                    "dAGO_l2": round(float(np.linalg.norm(ra["dAGO"] - rb["dAGO"])), 6),
                })

    report = {
        "units": [str(u) for u in args.units],
        "context_embeddings": {k: {kk: vv for kk, vv in v.items() if not isinstance(vv, np.ndarray)}
                               for k, v in records.items()},
        "pairs": pairs,
        "four_context_differential_stability": {
            "computable": bool(diff_pairs),
            "n_complete_units": len(differentials),
            "pairs": diff_pairs,
            "reason": ("" if diff_pairs else
                       "need at least two matched windows or replicas for one candidate"),
        },
        "sanity_checks": {
            "d_PAM_finite_nonzero": all(np.isfinite(r["dPAM"]).all() and
                                         np.linalg.norm(r["dPAM"]) > 0
                                         for r in differentials.values()),
            "d_AGO_finite_nonzero": all(np.isfinite(r["dAGO"]).all() and
                                         np.linalg.norm(r["dAGO"]) > 0
                                         for r in differentials.values()),
            "available_context_norm_is_one": all(abs(r["embedding_norm"] - 1.0) < 1e-5
                                                 for r in records.values()),
            "available_context_finite_nonzero": all(r["finite"] and r["nonzero"] > 0
                                                    for r in records.values()),
        },
        "claim_boundary": (
            "Descriptive representation stability only. No AUC, no PAM classification, "
            "no biological interpretation; windows/frames are not independent samples."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
