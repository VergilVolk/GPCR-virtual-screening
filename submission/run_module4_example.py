"""Recompute one real four-context MD block using unchanged frozen methods."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
import sys
import time
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parent
VENDOR = ROOT / "vendor/geom2vec"
# Resolve the exact bundled upstream source before historical imports.
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(VENDOR / "src"))

import MDAnalysis as mda
import numpy as np
import torch
from project.pacer_fkg_v02 import stage4_prospective_common_v01 as c
from project.pacer_fkg_v02 import run_stage4_prospective_phase1_bs256_v01 as p1
from project.pacer_fkg_v02 import run_stage4_prospective_phase2b_graph_region_v01 as p2b
from project.encoder_intermediate_v01.core import verify_and_strict_load
from src.examples.extract_four_context_sample import verify as verify_sample

CONTEXTS = dict(zip(c.ALIASES, (
    "cluster0__apo", "cluster0__probe_only",
    "PACER0073__candidate_no_probe", "PACER0073__candidate_probe")))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_vendor():
    import geom2vec
    manifest = json.loads((VENDOR / "PINNED_SOURCE_MANIFEST.json").read_text())
    if manifest["commit"] != "371d642ec1061664f16e49fcac702d07fc8d0b51":
        raise RuntimeError("Upstream pin changed")
    if Path(geom2vec.__file__).resolve() != (VENDOR / "src/geom2vec/__init__.py").resolve():
        raise RuntimeError("Unexpected geom2vec import source")
    for item in manifest["files"]:
        path = VENDOR / item["path"]
        if path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            raise RuntimeError(f"Pinned upstream source changed: {path}")
    return {"commit": manifest["commit"], "file_count": len(manifest["files"]),
            "manifest_sha256": sha256(VENDOR / "PINNED_SOURCE_MANIFEST.json")}


def sample_descriptors(frames):
    """Stage4 reductions on its already-strided first block; no second stride."""
    if frames.shape != (20, 270, 256) or not np.isfinite(frames).all():
        raise ValueError("Expected one finite frozen-stride BS256 block")
    blocks = frames.reshape(1, 20, 270, 256)
    mu = blocks.mean(axis=1, dtype=np.float64).astype(np.float32)
    delta = np.diff(blocks, axis=1)
    drift = delta.mean(axis=1, dtype=np.float64).astype(np.float32)
    rms = np.sqrt(np.mean(np.square(delta, dtype=np.float64), axis=1)).astype(np.float32)
    endpoint = (blocks[:, -1] - blocks[:, 0]) / np.float32(19)
    if not np.allclose(drift, endpoint, rtol=2e-5, atol=2e-6):
        raise RuntimeError("SIGNED_DRIFT endpoint identity failed")
    return {"STATE_MOTION": np.concatenate((mu, rms), axis=-1), "SIGNED_DRIFT": drift}


def run(sample_root, output, asset_root, device):
    start = time.perf_counter()
    if output.exists():
        raise FileExistsError(f"Use a new output directory: {output}")
    c.FROZEN_ROOT = asset_root.resolve()
    vendor = verify_vendor()
    verify_sample(sample_root)
    manifest_path = sample_root / "manifest.json"
    manifest_hash = sha256(manifest_path)
    sample = json.loads(manifest_path.read_text())
    records = {r["system"]: r for r in sample["systems"]}
    if (sample["candidate"] != "PACER0073" or set(records) != set(CONTEXTS.values())
            or any(r["replica"] != 1 or r["seed"] != 27101
                   or r["raw_frame_ids"] != list(range(0, 100, 5)) for r in records.values())):
        raise ValueError("Sample four-context identity changed")
    engine, anchor = c.verify_frozen_anchor()
    graph = c.load_json(asset_root / anchor["freeze"]["definitions"]["graph_definitions"]["path"])
    sequence = "".join(c.AA[n["resname"]] for n in graph["nodes"]) if "resname" in graph["nodes"][0] else "".join(n["site"][0] for n in graph["nodes"])
    checkpoint = ROOT / "project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth"
    model, model_report = verify_and_strict_load(checkpoint)
    if not model_report["complete_strict_coverage"]:
        raise RuntimeError("Frozen checkpoint strict-load failed")
    if device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable; use --device cpu")
    # Require the official compiled backend. No alternate neighbor algorithm.
    import torch_cluster
    from torch_geometric.nn import radius_graph
    radius_graph(torch.zeros((2, 3), device=device), r=1.0,
                 batch=torch.zeros(2, dtype=torch.long, device=device), max_num_neighbors=1)
    model = model.to(device).eval()
    extraction = p1.load_phase1_module()
    states = {b: engine.load_frozen_state(anchor, b) for b in c.BRANCHES}
    arrays, mapped, receptor_records, replay_records = {}, {b: {} for b in c.BRANCHES}, {}, {}
    for alias, system in CONTEXTS.items():
        record = records[system]
        top = sample_root / record["sample_files"]["topology"]["path"]
        traj = sample_root / record["sample_files"]["trajectory"]["path"]
        universe = mda.Universe(str(top), str(traj))
        try:
            selection, mapping = c.topology_selection(universe, sequence)
        finally:
            universe.trajectory.close()
        extraction.SELECTION = selection
        extraction.RESNAME_MAP = {**extraction.RESNAME_MAP, **c.RESNAME_ALIASES}
        job = SimpleNamespace(topology=top, trajectory=traj, sequence=sequence, key=system)
        u, atoms, atom_z, res_idx, bb, sc, _ = extraction.topology_mapping(job)
        try:
            frames = extraction.infer_bs256(model, u, atoms, atom_z, res_idx, bb, sc, list(range(20)), device)
            replay = extraction.infer_bs256(model, u, atoms, atom_z, res_idx, bb, sc, list(range(4)), device)
        finally:
            u.trajectory.close()
        replay_record = extraction.replay_metrics(frames[:4], replay)
        if not replay_record["passed"]:
            raise RuntimeError("Frozen numerical replay criterion failed")
        if np.count_nonzero(frames[:, [aa == "G" for aa in sequence], 128:]):
            raise RuntimeError("Gly side-chain feature contract failed")
        arrays[f"BS256__{alias}"] = frames
        receptor_records[alias], replay_records[alias] = mapping, replay_record
        for branch, descriptor in sample_descriptors(frames).items():
            arrays[f"descriptor__{branch}__{alias}"] = descriptor
            value = engine.rff_frozen(engine.normalize_frozen(descriptor, states[branch]), states[branch])
            if value.shape != (1, 270, 512) or not np.isfinite(value).all():
                raise RuntimeError("Frozen RFF contract failed")
            arrays[f"RFF__{branch}__{alias}"] = value
            mapped[branch][alias] = value
        print(f"Recomputed frozen BS256 and RFF: {alias} / {system}", flush=True)
    transition, regions = engine.graph_transition(graph), engine.region_indices(graph)
    rows = []
    for branch, contexts in mapped.items():
        for name, formula in p2b.all_contrasts(engine).items():
            diffused = engine.diffuse(p2b.contrast(contexts, name, engine), transition)
            arrays[f"graph__{branch}__{name}"] = diffused
            for region, indices in regions.items():
                vector = diffused[:, indices, :].mean(axis=1, dtype=np.float64).astype(np.float32)
                arrays[f"region__{branch}__{name}__{region}"] = vector
                rows.append({"candidate": "PACER0073", "replica": 1, "branch": branch,
                             "contrast": name, "formula": json.dumps(formula, sort_keys=True),
                             "region": region, "residues": len(indices), "blocks": 1,
                             "vector_norm": float(np.linalg.norm(vector[0])),
                             "interpretation": "single-block calculation example; no PAM classification"})
    verify_sample(sample_root)
    if sha256(manifest_path) != manifest_hash:
        raise RuntimeError("Sample manifest changed during run")
    c.verify_frozen_anchor()
    verify_vendor()
    output.mkdir(parents=True)
    np.savez_compressed(output / "computed_arrays.npz", **arrays)
    with (output / "module4_example_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    runtime = {name: importlib.metadata.version(name) for name in
               ("numpy", "scipy", "MDAnalysis", "torch", "torch-geometric", "torch-cluster", "torch-scatter")}
    report = {"status": "MODULE4_SAMPLE_CORE_COMPUTATION_PASS", "contexts": 4,
              "frames_per_context": 20, "blocks_per_context": 1, "replicas": [1],
              "frozen_receptor_residues": 270, "BS256_width": 256, "RFF_width": 512,
              "graph_alpha": engine.ALPHA, "graph_steps": engine.DIFFUSION_STEPS,
              "sample_manifest_sha256": manifest_hash, "vendor": vendor,
              "frozen_manifest_sha256": c.REQUIRED_FREEZE_SHA256,
              "checkpoint_sha256": model_report["actual_sha256"],
              "receptor_mapping": receptor_records, "numerical_replay": replay_records,
              "runtime": runtime, "python": platform.python_version(), "platform": platform.platform(),
              "neighbor_backend": "official torch_cluster.radius_graph through PyG; max_num_neighbors=32",
              "cross_device_numerical_equivalence_established": False,
              "device": device, "device_name": torch.cuda.get_device_name() if device.startswith("cuda") else "CPU",
              "elapsed_seconds": time.perf_counter() - start,
              "outputs": {p.name: sha256(p) for p in output.iterdir() if p.is_file()},
              "forbidden_operations_performed": c.FORBIDDEN,
              "claim_boundary": "Actual trajectory-to-FKG computation on one R1 block only; no replication agreement, full Stage4 result or PAM efficacy claim."}
    (output / "run_receipt.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "output": str(output), "rows": len(rows)}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", type=Path, default=ROOT / "data/examples/module4_pacer0073_r1")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "run/module4_example")
    parser.add_argument("--asset-root", type=Path, default=ROOT)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args()
    run(args.sample_dir, args.output_dir, args.asset_root, args.device)


if __name__ == "__main__":
    main()
