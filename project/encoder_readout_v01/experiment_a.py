"""Standalone Experiment A: authenticated cached features -> unscaled FP64 R/G.

No encoder imports, fitted transforms, kernels, contrasts, or statistical gates.
Run with -B to avoid bytecode artifacts. See README.md for the data contract.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np

VERSION = "encoder_readout_A_v01"
WORKSPACE = Path(__file__).resolve().parents[2]
MANIFEST = WORKSPACE / "project/results/pacer_dc_fkg_R2R3_v01/PACER_FKG_R2R3_AUDIT.json"
MANIFEST_SHA = "b7dbcb11f07b543694d85a680125c25b1fac2edf16be92ad0f9ff84f174bd90d"
GRAPH = WORKSPACE / "project/results/pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json"
BATCH = WORKSPACE / "project/results/pacer_dc_geom2vec_R2R3_full_v01/batch_audit.json"
BATCH_SHA = "5849b2c4e9c0e098a8480aea75bc121f35d10d4ba05b541e3c8406fa924621b5"
OUTPUT_BASE = WORKSPACE / "project/results/encoder_readout_A_v01"
CONTEXTS = ("apo", "probe_only", "candidate_no_probe", "candidate_probe")
# Numerical reconstruction bound only; not a biological/statistical gate.
EPS_MULTIPLIER = 128


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json_verified(path, expected):
    raw = Path(path).read_bytes()
    if sha256_bytes(raw) != expected:
        raise ValueError(f"Frozen JSON SHA256 mismatch: {path}")
    return json.loads(raw)


def expected_paths():
    return {
        f"replica_{r:02d}/window_{w:03d}/atom14/{c}_w{w:03d}.geom2vec.npz"
        for r in (2, 3) for w in range(5) for c in CONTEXTS
    }


def region_indices(graph):
    nodes = graph["nodes"]
    if any(n["embedding_index"] != i for i, n in enumerate(nodes)):
        raise ValueError("Noncanonical graph node order")
    regions = {}
    for name, members in graph["regions"].items():
        indices = [m["embedding_index"] for m in members]
        if not indices or len(set(indices)) != len(indices):
            raise ValueError(f"Empty or duplicate region membership: {name}")
        for m in members:
            i = m["embedding_index"]
            if not isinstance(i, int) or not 0 <= i < len(nodes):
                raise ValueError(f"Invalid residue index: {name}")
            if m["site"] != nodes[i]["site"]:
                raise ValueError(f"Residue identity mismatch: {name}")
        regions[name] = np.asarray(indices, dtype=np.int64)
    return regions


def load_contract():
    manifest = read_json_verified(MANIFEST, MANIFEST_SHA)
    provenance = manifest["provenance"]
    rows = provenance["input_files"]
    if len(rows) != 40 or {r["path"] for r in rows} != expected_paths():
        raise ValueError("Manifest must contain the exact 40 unique archived paths")
    graph = read_json_verified(GRAPH, provenance["graph_sha256"])
    if graph["version"] != "M4_MULTISTRUCTURE_GRAPH_v01" or len(graph["nodes"]) != 270:
        raise ValueError("Unexpected frozen graph")
    regions = region_indices(graph)
    if len(regions) != 9:
        raise ValueError("Expected nine frozen regions")
    sequence = provenance["sequence"]
    if "".join(n["site"][0] for n in graph["nodes"]) != sequence:
        raise ValueError("Graph and archived sequence disagree")
    batch = read_json_verified(BATCH, BATCH_SHA)
    if batch["files_total"] != 40 or batch["stride"] != 1:
        raise ValueError("Unexpected archived extraction metadata")
    return manifest, graph, regions, batch


def read_input(root, row, sequence, shape=(100, 270, 128)):
    """Hash and parse the same bytes; never write to the input archive."""
    root = Path(root).resolve()
    path = (root / row["path"]).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Input path escapes input root")
    if not path.is_file():
        raise FileNotFoundError(path)
    raw = path.read_bytes()
    if sha256_bytes(raw) != row["sha256"]:
        raise ValueError(f"NPZ SHA256 mismatch: {row['path']}")
    with np.load(io.BytesIO(raw), allow_pickle=False) as data:
        required = {"residue_features", "global_features", "sequence", "frame_ids"}
        if not required.issubset(data.files):
            raise ValueError(f"Missing NPZ fields: {row['path']}")
        arrays = {name: data[name] for name in data.files}
    # Reject nonfinite numbers in any numeric field, including optional metadata.
    for name, value in arrays.items():
        if np.issubdtype(value.dtype, np.number) and not np.isfinite(value).all():
            raise ValueError(f"Nonfinite {name}: {row['path']}")
    h, g = arrays["residue_features"], arrays["global_features"]
    if h.shape != shape or h.dtype != np.dtype("float32"):
        raise ValueError(f"Unexpected feature shape/dtype: {h.shape}/{h.dtype}")
    if g.shape != (shape[0], shape[2]) or g.dtype != np.dtype("float32"):
        raise ValueError("Unexpected archived global feature shape/dtype")
    seq = arrays["sequence"]
    if seq.shape != () or seq.dtype.kind != "U" or str(seq.item()) != sequence:
        raise ValueError("Sequence metadata mismatch")
    ids = arrays["frame_ids"]
    if ids.dtype != np.dtype("int64") or not np.array_equal(ids, np.arange(shape[0], dtype=np.int64)):
        raise ValueError("Expected stride-one local int64 frame IDs 0..T-1")
    if "residue_index" in arrays:
        ri = arrays["residue_index"]
        if ri.dtype != np.dtype("int64") or not np.array_equal(ri, np.arange(shape[1])):
            raise ValueError("Residue index metadata mismatch")
    # Preserve original storage layout for this historical FP32 replay.
    if not np.array_equal(g, h.mean(axis=1)):
        raise ValueError("Archived global features disagree with FP32 residue mean")
    return arrays


def preflight(root):
    manifest, graph, regions, batch = load_contract()
    records = []
    for row in manifest["provenance"]["input_files"]:
        arrays = read_input(root, row, manifest["provenance"]["sequence"])
        records.append({**row, "fields": {
            k: {"shape": list(v.shape), "dtype": str(v.dtype)} for k, v in arrays.items()
        }})
    return manifest, graph, regions, batch, records


def decompose(h, regions):
    """All arithmetic is FP64 on C-contiguous input; no scaling/fitting/casting back."""
    if h.ndim != 3 or min(h.shape) < 1 or h.dtype.kind != "f" or not np.isfinite(h).all():
        raise ValueError("Expected a finite, nonempty floating (frames,residues,channels) tensor")
    x = np.array(h, dtype=np.float64, order="C", copy=True)
    global_mean = x.mean(axis=1, dtype=np.float64)
    residual = x - global_mean[:, None, :]
    local, original = [], []
    for name, ix in regions.items():
        ix = np.asarray(ix)
        if ix.ndim != 1 or ix.dtype.kind not in "iu" or not len(ix):
            raise ValueError(f"Invalid region indices: {name}")
        if np.any(ix < 0) or np.any(ix >= x.shape[1]) or len(np.unique(ix)) != len(ix):
            raise ValueError(f"Out-of-range or duplicate region indices: {name}")
        local.append(residual[:, ix, :].mean(axis=1, dtype=np.float64))
        original.append(x[:, ix, :].mean(axis=1, dtype=np.float64))
    if not local:
        raise ValueError("At least one region is required")
    return {
        "global_mean": global_mean,
        "local_residual": residual,
        "region_local_mean": np.stack(local, axis=1),
        "region_original_mean": np.stack(original, axis=1),
    }


def error_metrics(actual, expected):
    delta = np.asarray(actual, dtype=np.float64) - np.asarray(expected, dtype=np.float64)
    denominator = float(np.linalg.norm(expected))
    return {"max_abs": float(np.max(np.abs(delta))),
            "rmse": float(np.sqrt(np.mean(delta * delta, dtype=np.float64))),
            "relative_l2": float(np.linalg.norm(delta) / denominator) if denominator else None}


def verify(h, branches, regions):
    x = np.array(h, dtype=np.float64, order="C", copy=True)
    g, r = branches["global_mean"], branches["local_residual"]
    expected_shapes = {"global_mean": (x.shape[0], x.shape[2]),
                       "local_residual": x.shape,
                       "region_local_mean": (x.shape[0], len(regions), x.shape[2]),
                       "region_original_mean": (x.shape[0], len(regions), x.shape[2])}
    for key, shape in expected_shapes.items():
        a = branches[key]
        if a.shape != shape or a.dtype != np.float64 or not np.isfinite(a).all():
            raise ValueError(f"Invalid output tensor: {key}")
    bound = EPS_MULTIPLIER * np.finfo(np.float64).eps * max(1.0, float(np.max(np.abs(x))))
    checks = {"residue_reconstruction": error_metrics(r + g[:, None, :], x),
              "global_mean": error_metrics(g, x.mean(axis=1, dtype=np.float64)),
              "residual_zero_mean": error_metrics(r.mean(axis=1, dtype=np.float64), np.zeros_like(g))}
    per_region = {}
    for j, (name, ix) in enumerate(regions.items()):
        original = x[:, ix, :].mean(axis=1, dtype=np.float64)
        reconstructed = branches["region_local_mean"][:, j, :] + g
        per_region[name] = error_metrics(reconstructed, original)
        checks[f"stored_region_original/{name}"] = error_metrics(branches["region_original_mean"][:, j, :], original)
    if any(m["max_abs"] > bound for m in [*checks.values(), *per_region.values()]):
        raise ValueError("FP64 reconstruction bound exceeded")
    repeated = decompose(h, regions)
    repeat = {k: branches[k].tobytes(order="C") == v.tobytes(order="C") for k, v in repeated.items()}
    if not all(repeat.values()):
        raise ValueError("Numerical repeatability failed")
    return {"absolute_error_bound": bound, "checks": checks, "regions": per_region,
            "repeat_bitwise_equal": repeat}


def check_output_root(output, input_root):
    output = Path(output).resolve()
    base = OUTPUT_BASE.resolve()
    # Also reject a symlink/junction that redirects the versioned base outside the worktree.
    if not base.is_relative_to(WORKSPACE) or output == base or not output.is_relative_to(base):
        raise ValueError(f"Output must be a new run directory beneath {OUTPUT_BASE}")
    source = Path(input_root).resolve()
    if output == source or output.is_relative_to(source) or source.is_relative_to(output):
        raise ValueError("Input/output paths overlap")
    if output.exists():
        raise FileExistsError(f"Refusing existing output directory: {output}")
    return output


def git_text(*args):
    return subprocess.check_output(["git", *args], cwd=WORKSPACE, text=True).strip()


def run(input_root, output_root=None):
    output = check_output_root(output_root, input_root) if output_root is not None else None
    manifest, graph, regions, batch, records = preflight(input_root)
    membership_count = sum(len(ix) for ix in regions.values())
    summary = {"version": VERSION, "inputs_verified": len(records),
               "region_counts": {k: len(v) for k, v in regions.items()},
               "total_region_memberships": membership_count,
               "unique_region_residues": len(set(np.concatenate(list(regions.values())).tolist())),
               "membership_resolution": "109 memberships; prior 119 was an arithmetic error, not a map change."}
    if output is None:
        return {"status": "PREFLIGHT_PASS", **summary}

    # No output directory is created until every input and all frozen contracts pass.
    output.mkdir(parents=True, exist_ok=False)
    names = np.asarray(list(regions))
    indices = np.concatenate(list(regions.values()))
    offsets = np.asarray([0, *np.cumsum([len(ix) for ix in regions.values()])], dtype=np.int64)
    sequence = manifest["provenance"]["sequence"]
    report = {"status": "RUNNING", **summary, "started_utc": datetime.now(timezone.utc).isoformat(),
              "input_root": str(Path(input_root).resolve()), "output_root": str(output),
              "arithmetic": {"input": "float32", "working_and_storage": "float64",
                             "layout": "C-contiguous working H", "reduction_axis": 1,
                             "normalization": "none; unscaled arithmetic mean and residual",
                             "error_bound": "128 * eps(float64) * max(1, max(abs(H)))",
                             "repeatability": "two independent decompositions; array bytes, not ZIP bytes"},
              "frame_contract": "Original local IDs preserved; no absolute times or source-DCD offsets inferred.",
              "source_npz_manifest": {"path": str(MANIFEST), "sha256": MANIFEST_SHA},
              "graph": {"path": str(GRAPH), "sha256": file_sha256(GRAPH), "version": graph["version"],
                        "nodes": graph["nodes"], "regions": graph["regions"]},
              "historical_extraction": {"path": str(BATCH), "sha256": BATCH_SHA, "audit": batch,
                                        "limitation": "Archived provenance, not verification of encoder weight coverage or installed-source equivalence."},
              "source_files": [{"path": str(p.relative_to(WORKSPACE)), "sha256": file_sha256(p)}
                               for p in sorted(Path(__file__).parent.rglob("*"))
                               if p.is_file() and p.suffix in {".py", ".md"}],
              "runtime": {"python": sys.version, "executable": sys.executable, "numpy": np.__version__,
                          "platform": platform.platform(), "machine": platform.machine(),
                          "thread_environment": {k: os.environ.get(k) for k in
                                                 ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")}},
              "git": {"head": git_text("rev-parse", "HEAD"), "branch": git_text("branch", "--show-current"),
                      "status_before_outputs": git_text("status", "--short")},
              "scope": "Cached-feature reconstruction only; no encoder, fit, contrasts, kernels, or statistical gates.",
              "files": []}
    for number, record in enumerate(records, 1):
        # Reauthenticate bytes before processing to detect changes since preflight.
        arrays = read_input(input_root, record, sequence)
        h = arrays["residue_features"]
        branches = decompose(h, regions)
        metrics = verify(h, branches, regions)
        relative = Path(record["path"]).with_suffix(".experiment_a.npz")
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = {**branches, "frame_ids": arrays["frame_ids"], "sequence": arrays["sequence"],
                   "residue_index": np.arange(270, dtype=np.int64),
                   "residue_sites": np.asarray([n["site"] for n in graph["nodes"]]),
                   "region_names": names, "region_indices": indices, "region_offsets": offsets,
                   "source_relative_path": np.asarray(record["path"]), "source_sha256": np.asarray(record["sha256"]),
                   "graph_sha256": np.asarray(report["graph"]["sha256"]), "version": np.asarray(VERSION)}
        with target.open("xb") as stream:
            np.savez_compressed(stream, **payload)
        with np.load(target, allow_pickle=False) as stored:
            if set(stored.files) != set(payload) or any(
                stored[k].dtype != v.dtype or stored[k].shape != v.shape or
                stored[k].tobytes(order="C") != v.tobytes(order="C") for k, v in payload.items()
            ):
                raise ValueError(f"Storage round-trip mismatch: {target}")
        report["files"].append({**record, "output": relative.as_posix(), "output_sha256": file_sha256(target),
                                "output_bytes": target.stat().st_size, "storage_roundtrip_bitwise_equal": True,
                                "numerical_validation": metrics})
        print(f"[{number}/40] {record['path']}: reconstruction and repeatability PASS", flush=True)
    report["summary"] = {
        "residue_reconstruction_max_abs": max(r["numerical_validation"]["checks"]["residue_reconstruction"]["max_abs"] for r in report["files"]),
        "region_reconstruction_max_abs": max(m["max_abs"] for r in report["files"] for m in r["numerical_validation"]["regions"].values()),
        "residual_mean_max_abs": max(r["numerical_validation"]["checks"]["residual_zero_mean"]["max_abs"] for r in report["files"]),
        "all_repeat_bitwise_equal": True, "all_storage_roundtrips_bitwise_equal": True,
        "total_output_npz_bytes": sum(r["output_bytes"] for r in report["files"])}
    report["status"] = "PASS"
    report["completed_utc"] = datetime.now(timezone.utc).isoformat()
    with (output / "validation_report.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return {"status": "PASS", **summary, **report["summary"], "report": str(output / "validation_report.json")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    if args.preflight_only and args.output_root is not None:
        parser.error("--preflight-only must not specify --output-root")
    if not args.preflight_only and args.output_root is None:
        parser.error("validation requires --output-root; otherwise use --preflight-only")
    print(json.dumps(run(args.input_root, args.output_root), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
