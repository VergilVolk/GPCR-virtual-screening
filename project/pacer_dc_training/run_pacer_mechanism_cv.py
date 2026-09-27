#!/usr/bin/env python3
"""Run PACER-MCV on four-context atom14 windows and audit R2/R3 replication."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from pacer_mechanism_cv import (
    CONTEXTS, cosine, extract_frame_features, factorial_contrast, load_config,
    load_mapping, load_sequence, robust_location_scale, summarize_trace,
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def find_input(root: Path, replica: int, window: int, context: str) -> tuple[Path, Path]:
    folder = root / f"replica_{replica:02d}" / f"window_{window:03d}" / "atom14"
    npy = folder / f"{context}_w{window:03d}.npy"
    csv_path = folder / f"{context}_w{window:03d}.csv"
    if not npy.is_file() or not csv_path.is_file():
        raise FileNotFoundError(f"missing atom14/sequence pair: {npy}, {csv_path}")
    return npy, csv_path


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("no rows")
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def bootstrap_cosine(window_vectors: dict[int, np.ndarray], draws: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    n2, n3 = len(window_vectors[2]), len(window_vectors[3])
    values = []
    for _ in range(draws):
        # Independent window resampling: window IDs across stochastic replicas are not paired observations.
        a = window_vectors[2][rng.integers(0, n2, n2)].mean(axis=0)
        b = window_vectors[3][rng.integers(0, n3, n3)].mean(axis=0)
        c = cosine(a, b)
        if c is not None:
            values.append(c)
    q = np.quantile(values, [0.025, 0.5, 0.975]) if values else [np.nan] * 3
    return {"q025": float(q[0]), "median": float(q[1]), "q975": float(q[2]), "draws": len(values)}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--mapping", type=Path, required=True)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--replicas", nargs="+", type=int, default=[2, 3])
    p.add_argument("--windows", nargs="+", type=int, default=list(range(5)))
    p.add_argument("--bootstrap", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20260927)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    if args.replicas != [2, 3]:
        raise ValueError("v0.1 frozen audit requires --replicas 2 3")

    config = load_config(args.config)
    mapping = load_mapping(args.mapping)
    raw: dict[int, dict[int, dict[str, dict[str, np.ndarray]]]] = {}
    provenance, feature_names = [], None
    for replica in args.replicas:
        raw[replica] = {}
        for window in args.windows:
            raw[replica][window] = {}
            expected_sequence = None
            for context in CONTEXTS:
                npy, seq = find_input(args.input_root, replica, window, context)
                sequence = load_sequence(seq)
                if expected_sequence is None: expected_sequence = sequence
                if sequence != expected_sequence: raise ValueError(f"sequence mismatch R{replica} W{window}")
                arr = np.load(npy, mmap_mode="r")
                names, trace = extract_frame_features(np.asarray(arr), sequence, mapping, config)
                if feature_names is None: feature_names = names
                if names != feature_names: raise ValueError("feature order changed")
                raw[replica][window][context] = summarize_trace(trace)
                provenance.append({"replica": replica, "window": window, "context": context,
                                   "path": str(npy), "sha256": sha256(npy), "frames": int(arr.shape[0])})

    # Each statistic has different units/distribution.  Calibrate each one only
    # on R2 and freeze it before R3; never mix Angstrom distances with lag-1.
    statistics = ("mean", "std", "q10", "q90", "lag1")
    scale_by_stat = {}
    for stat in statistics:
        calibration = np.stack([raw[2][w][c][stat] for w in args.windows for c in CONTEXTS])
        _, scale_by_stat[stat] = robust_location_scale(calibration)
    mechanism_mask = np.asarray(["distal_control" not in name for name in feature_names])
    if mechanism_mask.sum() < 2:
        raise ValueError("mechanism feature panel is empty")
    window_rows, summaries, feature_rows = [], [], []
    for stat in statistics:
        for axis, signs in config["axes"].items():
            by_rep = {}
            for replica in args.replicas:
                vectors = []
                for window in args.windows:
                    contexts = {c: raw[replica][window][c][stat] for c in CONTEXTS}
                    vector = factorial_contrast(contexts, signs) / scale_by_stat[stat]
                    vectors.append(vector)
                    for name, value in zip(feature_names, vector):
                        window_rows.append({"replica": replica, "window": window, "statistic": stat,
                                            "axis": axis, "feature": name, "scaled_contrast": float(value)})
                by_rep[replica] = np.stack(vectors)
            pooled2_all, pooled3_all = by_rep[2].mean(axis=0), by_rep[3].mean(axis=0)
            pooled2, pooled3 = pooled2_all[mechanism_mask], pooled3_all[mechanism_mask]
            sign = np.sign(pooled2) == np.sign(pooled3)
            active = (np.abs(pooled2) + np.abs(pooled3)) > 1e-8
            mechanism_windows = {rep: by_rep[rep][:, mechanism_mask] for rep in args.replicas}
            summaries.append({
                "statistic": stat, "axis": axis,
                "cross_replica_cosine": cosine(pooled2, pooled3),
                "feature_sign_agreement": float(sign[active].mean()) if active.any() else None,
                "r2_vector_norm": float(np.linalg.norm(pooled2)),
                "r3_vector_norm": float(np.linalg.norm(pooled3)),
                "mechanism_feature_count": int(mechanism_mask.sum()),
                "distal_control_excluded_from_primary_vector": True,
                "bootstrap_cosine": bootstrap_cosine(mechanism_windows, args.bootstrap, args.seed + len(summaries)),
            })
            for index, name in enumerate(feature_names):
                feature_rows.append({
                    "statistic": stat, "axis": axis, "feature": name,
                    "r2_pooled_scaled_contrast": float(pooled2_all[index]),
                    "r3_pooled_scaled_contrast": float(pooled3_all[index]),
                    "same_direction": bool(np.sign(pooled2_all[index]) == np.sign(pooled3_all[index])),
                    "minimum_absolute_effect": float(min(abs(pooled2_all[index]), abs(pooled3_all[index]))),
                    "included_in_primary_mechanism_vector": bool(mechanism_mask[index]),
                })

    primary = next(x for x in summaries if x["statistic"] == "mean" and x["axis"] == "synergy_interaction")
    # Frozen conservative signal rule; this is a replication screen, not a pharmacology test.
    primary["mechanism_signal_gate"] = bool(
        primary["cross_replica_cosine"] is not None
        and primary["cross_replica_cosine"] >= 0.50
        and primary["feature_sign_agreement"] >= 0.65
        and primary["bootstrap_cosine"]["q025"] > 0.0
    )
    report = {
        "method": config["model_id"], "status": "COMPLETED_DESCRIPTIVE",
        "evidence_level": "two_replica_physical_mechanism_screen",
        "feature_count": len(feature_names), "features": feature_names,
        "calibration": "separate R2-only robust scale per statistic; R3 held out; distal control excluded from primary vector",
        "primary_endpoint": "mean/synergy_interaction",
        "primary_result": primary, "all_endpoints": summaries,
        "provenance": {
            "config": str(args.config), "config_sha256": sha256(args.config),
            "mapping": str(args.mapping), "mapping_sha256": sha256(args.mapping),
            "input_files": provenance, "input_count": len(provenance),
        },
        "claim_boundary": config["claim_boundary"],
        "limitations": [
            "Only one candidate molecule; cannot estimate PAM classification performance.",
            "Two replicas and short windows; bootstrap is descriptive and not a molecular-level p-value.",
            "Receptor-only atom14 input cannot measure ACh or candidate pose stability.",
        ],
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "PACER_MCV_AUDIT_v01.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_csv(args.output / "PACER_MCV_WINDOW_FEATURES_v01.csv", window_rows)
    write_csv(args.output / "PACER_MCV_FEATURE_REPLICATION_v01.csv", feature_rows)
    print(json.dumps({"status": report["status"], "feature_count": len(feature_names),
                      "primary_result": primary, "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
