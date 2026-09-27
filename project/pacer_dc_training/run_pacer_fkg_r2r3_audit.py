#!/usr/bin/env python
"""Run the frozen PACER-FKG analysis over an R2/R3 Geom2Vec archive.

This wrapper exists so the GPU teammate can evaluate the already-extracted
40 NPZ files with one command.  It does not rerun MD or the encoder.
"""

from __future__ import annotations

import argparse
import hashlib
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


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def preflight(args):
    """Read-only input, sequence, graph and output checks; writes no results."""
    if not args.input_root.is_dir():
        raise NotADirectoryError(args.input_root)
    if not args.graph.is_file():
        raise FileNotFoundError(args.graph)
    if args.output_root.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing FKG output directory: {args.output_root}"
        )
    if (args.input_root.resolve() == args.output_root.resolve()
            or args.input_root.resolve() in args.output_root.resolve().parents):
        raise ValueError("FKG output must not be inside the input Embedding archive")
    if len(args.replicas) != 2 or len(set(args.replicas)) != 2:
        raise ValueError("R2/R3 audit requires two distinct replica IDs")
    if not args.windows or len(set(args.windows)) != len(args.windows):
        raise ValueError("Windows must be nonempty and unique")

    graph = json.loads(args.graph.read_text(encoding="utf-8"))
    nodes = graph["nodes"]
    edges = graph["edges"]
    if graph.get("version") != "M4_MULTISTRUCTURE_GRAPH_v01":
        raise ValueError("Unexpected graph version")
    if len(nodes) != 270 or len(edges) != 1304:
        raise ValueError("Frozen graph must contain 270 nodes and 1304 edges")
    if any(node.get("embedding_index") != index for index, node in enumerate(nodes)):
        raise ValueError("Graph node order is inconsistent with embedding indices")
    if any(not (0 <= edge["source"] < len(nodes) and 0 <= edge["target"] < len(nodes))
           for edge in edges):
        raise ValueError("Out-of-range graph edge")
    if any(nodes[member["embedding_index"]]["site"] != member["site"]
           for members in graph["regions"].values() for member in members):
        raise ValueError("Region residue mapping does not match graph node order")
    if "stable_core_control" not in graph["regions"]:
        raise ValueError("Frozen stable-core control is absent")

    manifest = []
    baseline_seq = None
    contexts = ("apo", "probe_only", "candidate_no_probe", "candidate_probe")
    for replica in args.replicas:
        for window in args.windows:
            root = args.input_root / f"replica_{replica:02d}" / f"window_{window:03d}" / "atom14"
            ids_in_window = None
            for context in contexts:
                path = root / f"{context}_w{window:03d}.geom2vec.npz"
                if not path.is_file():
                    raise FileNotFoundError(f"Missing FKG input: {path}")
                with np.load(path, allow_pickle=False) as data:
                    for key in ("residue_features", "sequence", "frame_ids"):
                        if key not in data:
                            raise ValueError(f"{path}: missing NPZ field {key}")
                    features = data["residue_features"]
                    sequence = str(data["sequence"].item())
                    ids = data["frame_ids"]
                    if features.shape != (args.expected_frames, len(nodes), 128):
                        raise ValueError(f"{path}: unexpected feature shape {features.shape}")
                    if not np.isfinite(features).all():
                        raise ValueError(f"{path}: nonfinite residue features")
                    if len(sequence) != len(nodes):
                        raise ValueError(f"{path}: sequence length mismatch")
                    if any(sequence[i] != node["site"][0] for i, node in enumerate(nodes)):
                        raise ValueError(f"{path}: sequence disagrees with frozen graph")
                    if baseline_seq is None:
                        baseline_seq = sequence
                    elif baseline_seq != sequence:
                        raise ValueError(f"{path}: inconsistent receptor sequence")
                    if ids.shape != (args.expected_frames,) or np.any(np.diff(ids) <= 0):
                        raise ValueError(f"{path}: missing/nonmonotonic frame IDs")
                    if ids_in_window is None:
                        ids_in_window = ids
                    elif not np.array_equal(ids_in_window, ids):
                        raise ValueError(f"{path}: within-window context frame IDs differ")
                manifest.append({
                    "path": str(path.relative_to(args.input_root)).replace("\\\\", "/"),
                    "sha256": file_sha256(path),
                })
    return {
        "graph_sha256": file_sha256(args.graph),
        "graph_version": graph["version"],
        "input_files": manifest,
        "n_inputs": len(manifest),
        "expected_frames": args.expected_frames,
        "sequence": baseline_seq,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--graph", type=Path, required=True)
    p.add_argument("--output-root", type=Path, required=True)
    p.add_argument("--replicas", nargs="+", type=int, default=[2, 3])
    p.add_argument("--windows", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    p.add_argument("--expected-frames", type=int, default=100)
    p.add_argument("--preflight-only", action="store_true")
    args = p.parse_args()

    provenance = preflight(args)
    if args.preflight_only:
        print(json.dumps({
            "preflight": "PASS",
            "inputs_checked": provenance["n_inputs"],
            "graph_sha256": provenance["graph_sha256"],
            "output_root_untouched": str(args.output_root),
        }, indent=2))
        return

    analyzer = Path(__file__).with_name("pacer_factorial_kernel_graph.py")
    args.output_root.mkdir(parents=True, exist_ok=False)
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
    stable_graph = {(r["replica"], r["window"], axis): r[f"{axis}_graph_diffused_mean"]
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
                graph_deltas = [r[f"{axis}_graph_diffused_mean"] -
                                stable_graph[(replica, r["window"], axis)] for r in rows]
                by_replica[str(replica)] = {
                    "window_u2_minus_stable": deltas,
                    "median_u2_minus_stable": float(np.median(deltas)),
                    "positive_windows": int(np.sum(np.asarray(deltas) > 0)),
                    "window_graph_diffused_minus_stable": graph_deltas,
                    "median_graph_diffused_minus_stable": float(np.median(graph_deltas)),
                }
            first = by_replica[str(args.replicas[0])]["window_u2_minus_stable"]
            second = by_replica[str(args.replicas[1])]["window_u2_minus_stable"] if len(args.replicas) == 2 else []
            rho = spearman(first, second) if second else None
            first_graph = by_replica[str(args.replicas[0])]["window_graph_diffused_minus_stable"]
            second_graph = (by_replica[str(args.replicas[1])]["window_graph_diffused_minus_stable"]
                            if len(args.replicas) == 2 else [])
            graph_rho = spearman(first_graph, second_graph) if second_graph else None
            gate = (
                region not in EXCLUDED_REGIONS
                and all(v["median_u2_minus_stable"] > 0 for v in by_replica.values())
                and all(v["positive_windows"] >= 4 for v in by_replica.values())
                and rho is not None and rho >= 0
            )
            summaries.append({
                "region": region, "axis": axis, "replicas": by_replica,
                "matched_window_spearman": rho,
                "matched_window_graph_diffused_spearman_descriptive": graph_rho,
                "preregistered_specificity_gate": gate,
            })

    report = {
        "method": "PACER-FKG R2/R3 frozen batch audit v01",
        "evidence_level": "two_replica_method_qualification",
        "input_root": str(args.input_root),
        "graph": str(args.graph),
        "replicas": args.replicas,
        "windows": args.windows,
        "provenance": {
            **provenance,
            "batch_script_sha256": file_sha256(Path(__file__)),
            "estimator_script_sha256": file_sha256(analyzer),
            "preprocess": "global_residual_channel_robust",
            "bandwidth": "per-call median heuristic",
            "graph_diffusion": {"alpha": 0.65, "steps": 20},
        },
        "gate_interpretation": (
            "The existing gate compares non-diffused region kernel squared estimates "
            "against stable core and matched-window scalar Spearman, not RKHS directional "
            "agreement and not graph-specific improvement. Independent prior registration "
            "of the exact numerical thresholds has not been verified."
        ),
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
