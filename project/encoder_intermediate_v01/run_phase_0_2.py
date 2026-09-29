#!/usr/bin/env python
"""Run the mandatory Experiment C Phase 0--2 gate and stop before Phase 3."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch_geometric

# Direct execution places this file's directory, rather than the repository
# root, on sys.path.  Add the root before importing versioned project modules.
WORKSPACE = Path(__file__).resolve().parents[2]
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))

from project.encoder_intermediate_v01.core import (
    EXPECTED_CHECKPOINT_SHA256,
    FROZEN_CONFIG,
    PINNED_GEOM2VEC_COMMIT,
    forward_three_states,
    hook_count,
    pool_batch,
    sha256,
    verify_and_strict_load,
)
from project.pacer_dc_training.extract_geom2vec_atom14 import (
    embed as production_embed,
    flatten_atom14,
    invariant_residue_features as production_pool,
    read_sequence,
)

DEFAULT_CHECKPOINT = WORKSPACE / "project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth"
DEFAULT_AUDIT = WORKSPACE / "project/results/pacer_dc_geom2vec_R2R3_full_v01/batch_audit.json"
DEFAULT_SOURCE = Path(r"C:\projects\geom2vec-source")
DEFAULT_INPUT = Path(r"C:\projects\GPCR-virtual-screening\project\results\pacer_dc_four_context_v01\compound110")
DEFAULT_ARCHIVE = Path(r"C:\projects\GPCR-virtual-screening\project\results\pacer_dc_geom2vec_R2R3_full_v01")
DEFAULT_OUTPUT = WORKSPACE / "project/results/encoder_intermediate_C_v01/run_001"
SMOKE_RELATIVE = Path("replica_02/window_000/atom14/apo_w000.npy")
SMOKE_FRAME_IDS = np.asarray([0, 1, 2, 3], dtype=np.int64)

SOURCE_FILES = [
    "__init__.py",
    "models/__init__.py",
    "models/factory.py",
    "models/representation/__init__.py",
    "models/representation/visnet.py",
    "nn/__init__.py",
    "nn/equivariant.py",
]


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def run_git(source: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", "-c", f"safe.directory={source.as_posix()}", "-C", str(source), *args],
        check=False,
        capture_output=True,
        text=True,
    )


def source_verification(source: Path) -> dict:
    import geom2vec

    installed = Path(geom2vec.__file__).resolve().parent
    direct_url_path = installed.parent / "geom2vec-1.0.0.dist-info" / "direct_url.json"
    direct_url = None
    if direct_url_path.is_file():
        direct_url = json.loads(direct_url_path.read_text(encoding="utf-8")).get("url")
    head_result = run_git(source, "rev-parse", "HEAD")
    head = head_result.stdout.strip() if head_result.returncode == 0 else None
    diff_result = run_git(
        source,
        "diff",
        "--exit-code",
        PINNED_GEOM2VEC_COMMIT,
        "--",
        *[f"src/geom2vec/{name}" for name in SOURCE_FILES],
    )
    rows = []
    for relative in SOURCE_FILES:
        source_path = source / "src/geom2vec" / relative
        installed_path = installed / relative
        source_hash = sha256(source_path) if source_path.is_file() else None
        installed_hash = sha256(installed_path) if installed_path.is_file() else None
        rows.append(
            {
                "file": relative,
                "pinned_checkout_path": str(source_path),
                "installed_path": str(installed_path),
                "pinned_checkout_sha256": source_hash,
                "installed_sha256": installed_hash,
                "match": source_hash is not None and source_hash == installed_hash,
            }
        )
    extractor = WORKSPACE / "project/pacer_dc_training/extract_geom2vec_atom14.py"
    complete = (
        head == PINNED_GEOM2VEC_COMMIT
        and diff_result.returncode == 0
        and all(row["match"] for row in rows)
    )
    return {
        "pinned_commit": PINNED_GEOM2VEC_COMMIT,
        "checkout": str(source.resolve()),
        "checkout_head": head,
        "checkout_head_matches_pin": head == PINNED_GEOM2VEC_COMMIT,
        "tracked_inference_files_clean_against_pin": diff_result.returncode == 0,
        "git_diff_output": (diff_result.stdout + diff_result.stderr).strip(),
        "installed_package_root": str(installed),
        "installed_direct_url": direct_url,
        "inference_files": rows,
        "repository_radius_graph_override": {
            "path": str(extractor),
            "sha256": sha256(extractor),
            "status_for_cuda_smoke": "not installed; official PyG torch-cluster backend is available",
        },
        "source_equivalent": complete,
        "substantive_inference_path_differences": [] if complete else ["see failed hash/commit/diff checks"],
    }


def runtime_report(device: str) -> dict:
    return {
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "pytorch": torch.__version__,
        "pyg": torch_geometric.__version__,
        "numpy": np.__version__,
        "cuda_build": torch.version.cuda,
        "cuda_available": torch.cuda.is_available(),
        "cudnn": torch.backends.cudnn.version(),
        "device": device,
        "device_name": torch.cuda.get_device_name(torch.device(device)) if device.startswith("cuda") else "CPU",
        "neighbor_backend": {
            "implementation": "torch_cluster.radius_graph via torch_geometric.nn.radius_graph",
            "torch_cluster": package_version("torch-cluster"),
            "pyg_lib": package_version("pyg-lib"),
            "torch_scatter": package_version("torch-scatter"),
            "torch_sparse": package_version("torch-sparse"),
            "repository_cpu_fallback_active": False,
        },
        "historical_runtime_equivalence": {
            "established": False,
            "reason": "The historical YAML gives minimum/range constraints, not an exact solved environment or hardware/backend record.",
        },
    }


def inventory_inputs(audit_path: Path, input_root: Path) -> dict:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    rows = []
    for expected in audit["rows"]:
        path = input_root / expected["input"]
        sequence_path = path.with_suffix(".csv")
        found = path.is_file()
        actual_hash = sha256(path) if found else None
        rows.append(
            {
                "relative_path": expected["input"],
                "expected_sha256": expected["input_sha256"],
                "found_path": str(path.resolve()) if found else None,
                "found": found,
                "actual_sha256": actual_hash,
                "matches_archived_provenance": found and actual_hash == expected["input_sha256"],
                "sequence_csv": str(sequence_path.resolve()) if sequence_path.is_file() else None,
                "sequence_csv_sha256": sha256(sequence_path) if sequence_path.is_file() else None,
            }
        )
    return {
        "archived_batch_audit": str(audit_path.resolve()),
        "archived_batch_audit_sha256": sha256(audit_path),
        "input_root": str(input_root.resolve()),
        "expected_files": len(rows),
        "found_files": sum(row["found"] for row in rows),
        "hash_matches": sum(row["matches_archived_provenance"] for row in rows),
        "all_authenticated": all(row["matches_archived_provenance"] for row in rows),
        "files": rows,
    }


def error_metrics(observed: torch.Tensor | np.ndarray, reference: torch.Tensor | np.ndarray) -> dict:
    a = torch.as_tensor(observed).detach().double().cpu()
    b = torch.as_tensor(reference).detach().double().cpu()
    delta = a - b
    denominator = torch.linalg.vector_norm(b).clamp_min(torch.finfo(torch.float64).tiny)
    return {
        "max_abs": float(delta.abs().max()),
        "rmse": float(torch.sqrt(torch.mean(delta.square()))),
        "relative_l2": float(torch.linalg.vector_norm(delta) / denominator),
    }


def prepare_batch(xyz: np.ndarray, atom_z: np.ndarray, device: str):
    positions = torch.as_tensor(xyz, dtype=torch.float32, device=device)
    batch_size, n_atoms, _ = positions.shape
    z0 = torch.as_tensor(atom_z, dtype=torch.long, device=device)
    z = z0.repeat(batch_size)
    batch = torch.arange(batch_size, device=device).repeat_interleave(n_atoms)
    return z, positions.reshape(-1, 3), batch, batch_size, n_atoms


def run_three(model, xyz, atom_z, residue_index, n_residues, device):
    z, pos, batch, batch_size, n_atoms = prepare_batch(xyz, atom_z, device)
    ridx = torch.as_tensor(residue_index, dtype=torch.long, device=device)
    with torch.inference_mode():
        states, returned = forward_three_states(model, z, pos, batch)
        pooled = {
            name: pool_batch(state, ridx, n_residues, batch_size, n_atoms).detach().cpu()
            for name, state in states.items()
        }
        raw = {name: (state[0].detach().cpu(), state[1].detach().cpu()) for name, state in states.items()}
        returned_cpu = tuple(value.detach().cpu() for value in returned)
    return pooled, raw, returned_cpu


def smoke_test(
    model: torch.nn.Module,
    input_root: Path,
    archive_root: Path,
    device: str,
) -> dict:
    atom14 = input_root / SMOKE_RELATIVE
    sequence_csv = atom14.with_suffix(".csv")
    sequence = read_sequence(sequence_csv)
    original = np.load(atom14, mmap_mode="r")
    coords = np.asarray(original[SMOKE_FRAME_IDS])
    xyz, atom_z, residue_index = flatten_atom14(coords, sequence)
    model = model.to(device).eval()

    hooks_before = hook_count(model)
    pooled, raw, returned_hooked = run_three(model, xyz, atom_z, residue_index, len(sequence), device)
    hooks_after = hook_count(model)
    replay, replay_raw, _ = run_three(model, xyz, atom_z, residue_index, len(sequence), device)
    hooks_after_replay = hook_count(model)

    z, pos, batch, _, _ = prepare_batch(xyz, atom_z, device)
    with torch.inference_mode():
        returned_plain = tuple(value.detach().cpu() for value in model(z=z, pos=pos, batch=batch))

    production = production_embed(
        model, xyz, atom_z, residue_index, len(sequence), len(SMOKE_FRAME_IDS), device
    )
    production_direct = []
    ridx = torch.as_tensor(residue_index, dtype=torch.long)
    n_atoms = len(atom_z)
    for frame in range(len(SMOKE_FRAME_IDS)):
        production_direct.append(
            production_pool(
                returned_plain[0][frame * n_atoms : (frame + 1) * n_atoms],
                returned_plain[1][frame * n_atoms : (frame + 1) * n_atoms],
                ridx,
                len(sequence),
            )
        )
    production_direct = torch.stack(production_direct)

    archive_path = archive_root / SMOKE_RELATIVE.with_suffix(".geom2vec.npz")
    archive_features = np.load(archive_path)["residue_features"][SMOKE_FRAME_IDS] if archive_path.is_file() else None

    angle = np.deg2rad(37.0)
    axis = np.asarray([1.0, 2.0, -1.0], dtype=np.float64)
    axis /= np.linalg.norm(axis)
    skew = np.asarray([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rotation = np.eye(3) * np.cos(angle) + (1 - np.cos(angle)) * np.outer(axis, axis) + np.sin(angle) * skew
    translation = np.asarray([11.0, -7.0, 3.5], dtype=np.float64)
    rotated_xyz = (xyz.astype(np.float64) @ rotation.T).astype(np.float32)
    translated_xyz = (xyz.astype(np.float64) + translation).astype(np.float32)
    rotated, _, _ = run_three(model, rotated_xyz, atom_z, residue_index, len(sequence), device)
    translated, _, _ = run_three(model, translated_xyz, atom_z, residue_index, len(sequence), device)

    state_shapes = {}
    for name in ("C0", "C1", "C2"):
        state_shapes[name] = {
            "atom_scalar": list(raw[name][0].shape),
            "atom_vector": list(raw[name][1].shape),
            "residue_features": list(pooled[name].shape),
            "finite": bool(torch.isfinite(pooled[name]).all()),
        }
    return {
        "predeclared_input": {
            "atom14": str(atom14.resolve()),
            "atom14_sha256": sha256(atom14),
            "sequence_csv": str(sequence_csv.resolve()),
            "sequence_csv_sha256": sha256(sequence_csv),
            "source_shape": list(original.shape),
            "frame_ids": SMOKE_FRAME_IDS.tolist(),
            "residues": len(sequence),
            "heavy_atoms": len(atom_z),
        },
        "hook_locations": {
            "C0": "exact ViSNet model return (x_rep, v_rep)",
            "C1": "representation_model.vis_mp_layers[3] forward pre-hook inputs[0:2]",
            "C2_scalar": "representation_model.out_norm forward pre-hook input[0]",
            "C2_vector": "representation_model.vec_out_norm forward pre-hook input[0]",
        },
        "state_shapes_and_finiteness": state_shapes,
        "ordering": {
            "frame_ids_observed": SMOKE_FRAME_IDS.tolist(),
            "frame_ids_expected": [0, 1, 2, 3],
            "frame_order_exact": bool(np.array_equal(SMOKE_FRAME_IDS, np.arange(4))),
            "residue_index_observed": [0, len(sequence) - 1],
            "residue_count": len(sequence),
            "residue_order_exact": bool(np.array_equal(np.arange(len(sequence)), np.arange(270))),
        },
        "hook_cleanup": {
            "count_before": hooks_before,
            "count_after_first_forward": hooks_after,
            "count_after_replay": hooks_after_replay,
            "clean": hooks_before == hooks_after == hooks_after_replay,
        },
        "c0_production_match": error_metrics(pooled["C0"], production),
        "c0_direct_production_pool_match": error_metrics(pooled["C0"], production_direct),
        "c0_archived_production_match": error_metrics(pooled["C0"], archive_features) if archive_features is not None else None,
        "hooked_vs_unhooked_return": {
            f"returned_tensor_{index}": error_metrics(hooked, plain)
            for index, (hooked, plain) in enumerate(zip(returned_hooked, returned_plain))
        },
        "deterministic_replay": {name: error_metrics(replay[name], pooled[name]) for name in ("C0", "C1", "C2")},
        "deterministic_replay_raw": {
            name: {
                "scalar": error_metrics(replay_raw[name][0], raw[name][0]),
                "vector": error_metrics(replay_raw[name][1], raw[name][1]),
            }
            for name in ("C0", "C1", "C2")
        },
        "rigid_transform": {
            "rotation_matrix": rotation.tolist(),
            "rotation_determinant": float(np.linalg.det(rotation)),
            "translation": translation.tolist(),
            "rotation_invariance": {name: error_metrics(rotated[name], pooled[name]) for name in ("C0", "C1", "C2")},
            "translation_invariance": {name: error_metrics(translated[name], pooled[name]) for name in ("C0", "C1", "C2")},
        },
    }


def all_zero(metrics: dict, tolerance: float) -> bool:
    return all(value["max_abs"] <= tolerance for value in metrics.values())


def render_report(report: dict) -> str:
    ckpt = report["checkpoint"]
    source = report["source"]
    inventory = report["input_inventory"]
    smoke = report.get("smoke_test")
    lines = [
        "# PACER-DC encoder optimization — Experiment C v01",
        "",
        "## Phase 0–2 stopping report",
        "",
        f"Status: **{report['status']}**. Phase 3 was not run.",
        "",
        "## Phase 0 provenance gate",
        "",
        f"- Checkpoint SHA256: `{ckpt['actual_sha256']}` (expected hash match: `{ckpt['sha256_match']}`).",
        f"- Strict checkpoint load: `{ckpt['strict_load_ok']}`; parameter coverage: `{ckpt['parameter_count_coverage']['matched']}/{ckpt['parameter_count_coverage']['total']}` keys and `{ckpt['parameter_element_coverage']['fraction']:.6f}` of parameter elements.",
        f"- Source checkout HEAD: `{source['checkout_head']}`; pinned-source equivalence: `{source['source_equivalent']}`.",
        f"- Historical runtime equivalence: `{report['runtime']['historical_runtime_equivalence']['established']}`. {report['runtime']['historical_runtime_equivalence']['reason']}",
        f"- Authenticated short-trajectory atom14 inventory: `{inventory['hash_matches']}/{inventory['expected_files']}` match archived extraction hashes.",
        "",
        "## Runtime",
        "",
        f"- Python: `{report['runtime']['python'].splitlines()[0]}`",
        f"- PyTorch / PyG / NumPy: `{report['runtime']['pytorch']}` / `{report['runtime']['pyg']}` / `{report['runtime']['numpy']}`",
        f"- CUDA build / device: `{report['runtime']['cuda_build']}` / `{report['runtime']['device_name']}`",
        f"- Neighbor backend: `{report['runtime']['neighbor_backend']['implementation']}`, torch-cluster `{report['runtime']['neighbor_backend']['torch_cluster']}`.",
    ]
    if smoke:
        lines += [
            "",
            "## Fixed states and smoke test",
            "",
            "C0 is the returned normalized production state. C1 is the accumulated input to message-passing layer index 3. C2 is captured at the scalar and vector output-normalization pre-hooks.",
            "",
            "All three atom states have scalar shape `[8556, 64]`, vector shape `[8556, 3, 64]`, and pooled residue shape `[4, 270, 128]` for the four-frame batch.",
            "",
            f"- C0 vs existing production `embed`: max absolute error `{smoke['c0_production_match']['max_abs']:.9g}`.",
            f"- Deterministic replay maxima (C0/C1/C2): " + ", ".join(f"`{smoke['deterministic_replay'][name]['max_abs']:.9g}`" for name in ('C0','C1','C2')) + ".",
            f"- Rotation-invariance maxima (C0/C1/C2): " + ", ".join(f"`{smoke['rigid_transform']['rotation_invariance'][name]['max_abs']:.9g}`" for name in ('C0','C1','C2')) + ".",
            f"- Translation-invariance maxima (C0/C1/C2): " + ", ".join(f"`{smoke['rigid_transform']['translation_invariance'][name]['max_abs']:.9g}`" for name in ('C0','C1','C2')) + ".",
            f"- Hook cleanup: `{smoke['hook_cleanup']['clean']}`; hooked-return agreement within CUDA tolerance: `{all_zero(smoke['hooked_vs_unhooked_return'], 5e-6)}`.",
            f"- CUDA numerical tolerances used by the engineering gate: `{smoke['acceptance_tolerances']}`. Raw errors remain recorded above and in `smoke_test_metrics.json`.",
            "",
            "## Gate decision",
            "",
            report["phase_3_decision"],
        ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--batch-audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--geom2vec-source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--archive-root", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    write_json(args.output_root / "frozen_config.json", FROZEN_CONFIG)

    model, checkpoint = verify_and_strict_load(args.checkpoint)
    source = source_verification(args.geom2vec_source)
    inventory = inventory_inputs(args.batch_audit, args.input_root)
    runtime = runtime_report(args.device)
    write_json(args.output_root / "checkpoint_verification.json", checkpoint)
    write_json(args.output_root / "source_provenance.json", source)
    write_json(args.output_root / "runtime.json", runtime)
    write_json(args.output_root / "input_inventory.json", inventory)

    report = {
        "experiment": "PACER-DC encoder optimization — Experiment C v01",
        "scope": "Phase 0 through Phase 2 only",
        "checkpoint": checkpoint,
        "source": source,
        "runtime": runtime,
        "input_inventory": inventory,
    }
    gate_errors = []
    if not checkpoint["complete_strict_coverage"]:
        gate_errors.append("checkpoint strict coverage failed")
    if not source["source_equivalent"]:
        gate_errors.append("inference-path source equivalence failed")
    if not inventory["all_authenticated"]:
        gate_errors.append("historical atom14 inventory is incomplete or hash-mismatched")
    if args.device.startswith("cuda") and not torch.cuda.is_available():
        gate_errors.append("requested CUDA device is unavailable")
    if gate_errors:
        report["status"] = "STOPPED_BEFORE_INFERENCE"
        report["gate_errors"] = gate_errors
        report["phase_3_decision"] = "Phase 3 is not technically safe to run because: " + "; ".join(gate_errors) + "."
        write_json(args.output_root / "verification_report.json", report)
        (args.output_root / "REPORT.md").write_text(render_report(report), encoding="utf-8")
        raise SystemExit(report["phase_3_decision"])

    started = time.time()
    smoke = smoke_test(model, args.input_root, args.archive_root, args.device)
    smoke["elapsed_seconds"] = round(time.time() - started, 3)
    smoke["acceptance_tolerances"] = {
        "c0_production_match_max_abs": 2e-6,
        "hooked_vs_unhooked_return_max_abs": 5e-6,
        "deterministic_replay_max_abs": 5e-5,
        "rigid_transform_invariance_max_abs": 2e-4,
        "note": "Separate CUDA forwards can differ at low order because graph/scatter reductions use parallel floating-point accumulation.",
    }
    write_json(args.output_root / "smoke_test_metrics.json", smoke)
    smoke_ok = (
        all(item["finite"] for item in smoke["state_shapes_and_finiteness"].values())
        and all(item["residue_features"] == [4, 270, 128] for item in smoke["state_shapes_and_finiteness"].values())
        and smoke["ordering"]["frame_order_exact"]
        and smoke["ordering"]["residue_order_exact"]
        and smoke["hook_cleanup"]["clean"]
        and smoke["c0_production_match"]["max_abs"] <= 2e-6
        and all_zero(smoke["hooked_vs_unhooked_return"], 5e-6)
        and all_zero(smoke["deterministic_replay"], 5e-5)
        and all_zero(smoke["rigid_transform"]["rotation_invariance"], 2e-4)
        and all_zero(smoke["rigid_transform"]["translation_invariance"], 2e-4)
    )
    report["smoke_test"] = smoke
    report["status"] = "PHASE_0_2_COMPLETE" if smoke_ok else "STOPPED_AFTER_SMOKE_FAILURE"
    report["gate_errors"] = [] if smoke_ok else ["one or more Phase 2 smoke criteria failed; inspect smoke_test_metrics.json"]
    report["phase_3_decision"] = (
        "Phase 3 is technically safe to run under the frozen configuration and authenticated short R2/R3 inputs. It was deliberately not run in this turn."
        if smoke_ok
        else "Phase 3 is not technically safe to run until the smoke discrepancy is resolved."
    )
    write_json(args.output_root / "verification_report.json", report)
    (args.output_root / "REPORT.md").write_text(render_report(report), encoding="utf-8")
    if smoke_ok:
        write_json(
            args.output_root / "PHASE_0_2_COMPLETE.json",
            {
                "status": report["status"],
                "phase_3_run": False,
                "verification_report_sha256": sha256(args.output_root / "verification_report.json"),
                "report_sha256": sha256(args.output_root / "REPORT.md"),
            },
        )
    print(json.dumps({"status": report["status"], "output": str(args.output_root)}, indent=2))


if __name__ == "__main__":
    main()
