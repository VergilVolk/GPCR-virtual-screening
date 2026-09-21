#!/usr/bin/env python
"""Phase 4b/5/6 analysis: baseline representations on the identical window, plus
descriptive stability of the context embedding across replicas and windows.

Read-only over existing unit directories. No training, no AUC, no PAM score.

    python audit_pacer_dc_four_context_stability.py \
      --units <unit_dir> [<unit_dir> ...] \
      --output project/results/pacer_dc_four_context_v01/stability_report.json
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path

import numpy as np

EXTRACTOR = Path(__file__).resolve().parent / "extract_pacer_dc_four_context_embeddings.py"


def load_module():
    spec = importlib.util.spec_from_file_location("fce", EXTRACTOR)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cos(a, b) -> float:
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--units", nargs="+", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    fce = load_module()

    records = {}
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
                arr = np.load(unit / row["atom14_path"])
                ca = arr[:, :, 1, :].reshape(arr.shape[0], -1).astype(np.float64)
                reps = fce.baseline_representations(ca)
                records[key] = {
                    "unit": str(unit), "context": row["context"],
                    "replica": row["replica_id"], "start_frame": row["start_frame"],
                    "one_prot": emb, "pca": reps["pca"], "tica": reps["tica"],
                    "vamp": reps["vamp"], "vamp2_score": reps["vamp2_score"],
                    "embedding_norm": float(np.linalg.norm(emb)),
                    "finite": bool(np.isfinite(emb).all()),
                    "nonzero": int((emb != 0).sum()),
                    "embedding_dim": int(emb.size),
                }

    print("=== available context embeddings (descriptive) ===")
    for k in sorted(records):
        r = records[k]
        print(f"  {k:<46} z={r['embedding_dim']} norm={r['embedding_norm']:.6f} "
              f"finite={r['finite']} nonzero={r['nonzero']} vamp2={r['vamp2_score']:.4g}")

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
                "pca_cosine": round(cos(ra["pca"], rb["pca"]), 6),
                "tica_cosine": round(cos(ra["tica"], rb["tica"]), 6),
                "vamp_cosine": round(cos(ra["vamp"], rb["vamp"]), 6),
            })

    print("\n=== descriptive stability (cosine of the representation vectors) ===")
    hdr = f"  {'kind':<20}{'context':<10}{'pair':<44}{'OneProt':<10}{'PCA':<9}{'tICA':<9}{'VAMP':<9}"
    print(hdr)
    for p in sorted(pairs, key=lambda x: (x["kind"], x["a"])):
        print(f"  {p['kind']:<20}{p['context']:<10}{p['a'] + ' vs ' + p['b']:<44}"
              f"{p['one_prot_cosine']:<10}{p['pca_cosine']:<9}{p['tica_cosine']:<9}{p['vamp_cosine']:<9}")

    report = {
        "units": [str(u) for u in args.units],
        "context_embeddings": {k: {kk: vv for kk, vv in v.items() if not isinstance(vv, np.ndarray)}
                               for k, v in records.items()},
        "pairs": pairs,
        "four_context_differential_stability": {
            "computable": False,
            "reason": ("no candidate has all four contexts in the production tree; "
                       "d_PAM / d_AGO stability cannot be computed"),
        },
        "sanity_checks": {
            "z_CA_ne_z_A": "not computable (candidate_probe / probe_only absent)",
            "z_C_ne_z_0": "not computable (candidate_no_probe absent)",
            "d_PAM_finite_nonzero": "not computed (unit incomplete)",
            "d_AGO_finite_nonzero": "not computed (unit incomplete)",
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
