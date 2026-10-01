from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np


SYSTEMS = {
    "apo": "apo",
    "probe_only": "probe",
    "compound110__candidate_no_probe": "candidate",
    "compound110__candidate_probe": "candidate_probe",
}

REPLICAS = (1, 2, 3)

EXPECTED_CONTRACT_SHA256 = (
    "8d2e7d49e8f5970485db62cc5669002459c2e0f0d490275316e9f2e4c2e36d77"
)
EXPECTED_G2C_AUDIT_SHA256 = (
    "1559d58293196ee2506cfb63a286ea04192832632675b580366b4a9b1253bfff"
)
EXPECTED_REFERENCE_SHA256 = (
    "27e4d3b466fbabbba0b8c7042844c3565ddea1958f3008fcb9c369000ac6fc6f"
)

EXPECTED_WIDTH = 16.476339519583565
EXPECTED_BASE_SEED = 271828
EXPECTED_RFF_SEED = 272599
EXPECTED_RFF_DIM = 256
EXPECTED_BLOCK_FRAMES = 20
EXPECTED_BLOCKS = 50
EXPECTED_FRAMES = 1000
EXPECTED_RESIDUES = 270
EXPECTED_CHANNELS = 128


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def array_sha256(a: np.ndarray) -> str:
    a = np.ascontiguousarray(a)
    h = hashlib.sha256()
    h.update(str(a.dtype).encode("ascii"))
    h.update(str(tuple(a.shape)).encode("ascii"))
    h.update(a.tobytes(order="C"))
    return h.hexdigest()


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(msg)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_rff_basis(width: float, seed: int, dim: int):
    # Exact G2-C generation semantics:
    # rng=np.random.default_rng(args.seed+771)
    # weights=rng.normal(0,1/width,size=(128,args.rff_dim)).astype(np.float32)
    # bias=rng.uniform(0,2*np.pi,size=args.rff_dim).astype(np.float32)
    rng = np.random.default_rng(seed)
    weights = rng.normal(
        0.0,
        1.0 / width,
        size=(EXPECTED_CHANNELS, dim),
    ).astype(np.float32)
    bias = rng.uniform(
        0.0,
        2.0 * np.pi,
        size=dim,
    ).astype(np.float32)
    return weights, bias


def audit(args):
    contract_path = args.contract.resolve()
    feature_manifest_path = args.feature_manifest.resolve()
    g2c_path = args.g2c_audit.resolve()
    reference_path = args.reference.resolve()

    require(contract_path.exists(), f"MISSING contract: {contract_path}")
    require(feature_manifest_path.exists(), f"MISSING feature manifest: {feature_manifest_path}")
    require(g2c_path.exists(), f"MISSING G2-C audit: {g2c_path}")
    require(reference_path.exists(), f"MISSING reference: {reference_path}")

    require(
        sha256(contract_path) == EXPECTED_CONTRACT_SHA256,
        "CONTRACT_SHA256_MISMATCH",
    )
    require(
        sha256(g2c_path) == EXPECTED_G2C_AUDIT_SHA256,
        "G2C_AUDIT_SHA256_MISMATCH",
    )
    require(
        sha256(reference_path) == EXPECTED_REFERENCE_SHA256,
        "ERRATUM_REFERENCE_SHA256_MISMATCH",
    )

    contract = load_json(contract_path)
    manifest = load_json(feature_manifest_path)
    g2c = load_json(g2c_path)

    require(
        contract.get("contract_version") == "PACER_FKG_LONG_MD_v02",
        "WRONG_CONTRACT_VERSION",
    )
    require(
        contract.get("biological_contrasts_computed_at_contract_creation") is False,
        "CONTRACT_ALREADY_CONTAINS_CONTRASTS",
    )

    temporal = contract["temporal_contract"]
    kernel = contract["common_kernel"]
    legacy = contract["legacy_policy"]
    reporting = contract["cross_replica_reporting"]

    require(
        temporal["block_frames"] == EXPECTED_BLOCK_FRAMES,
        "BLOCK_FRAMES_MISMATCH",
    )
    require(
        temporal["blocks_per_trajectory"] == EXPECTED_BLOCKS,
        "BLOCK_COUNT_MISMATCH",
    )
    require(
        temporal["overlap"] is False,
        "OVERLAPPING_BLOCKS_FORBIDDEN",
    )
    require(
        temporal["blocks_are_independent_biological_replicates"] is False,
        "BLOCKS_MUST_NOT_BE_BIOLOGICAL_REPLICATES",
    )
    require(
        temporal["independent_replication_unit"] == "replica_seed",
        "WRONG_INDEPENDENT_REPLICATION_UNIT",
    )

    require(
        int(kernel["rff_dimension"]) == EXPECTED_RFF_DIM,
        "RFF_DIM_MISMATCH",
    )
    require(
        int(kernel["base_seed"]) == EXPECTED_BASE_SEED,
        "BASE_SEED_MISMATCH",
    )
    require(
        int(kernel["rff_rng_seed_offset"]) == 771,
        "RFF_SEED_OFFSET_MISMATCH",
    )
    require(
        int(kernel["actual_rff_rng_seed"]) == EXPECTED_RFF_SEED,
        "ACTUAL_RFF_SEED_MISMATCH",
    )
    require(
        abs(float(kernel["shared_node_bandwidth"]) - EXPECTED_WIDTH) <= 1e-12,
        "SHARED_BANDWIDTH_MISMATCH",
    )
    require(
        kernel["bandwidth_refit_on_long_md"] is False,
        "LONG_MD_BANDWIDTH_REFIT_FORBIDDEN",
    )
    require(
        kernel["rff_basis_refit_on_long_md"] is False,
        "LONG_MD_RFF_REFIT_FORBIDDEN",
    )
    require(
        kernel["common_rff_basis_required_across_all_contexts_replicas_blocks"] is True,
        "COMMON_RFF_BASIS_NOT_REQUIRED",
    )

    require(
        legacy["kernel_u2_scalar_gate"] == "REVOKED",
        "LEGACY_KERNEL_U2_GATE_NOT_REVOKED",
    )
    require(
        reporting["new_qualification_gate"] == "NOT_DEFINED",
        "UNDECLARED_NEW_GATE_PRESENT",
    )
    require(
        reporting["automatic_functional_gate"] is False,
        "AUTOMATIC_FUNCTIONAL_GATE_FORBIDDEN",
    )

    cal = g2c["calibration"]
    require(
        abs(float(cal["shared_node_bandwidth"]) - EXPECTED_WIDTH) <= 1e-12,
        "G2C_WIDTH_DOES_NOT_MATCH_CONTRACT",
    )
    require(
        int(cal["rff_seed"]) == EXPECTED_RFF_SEED,
        "G2C_RFF_SEED_DOES_NOT_MATCH_CONTRACT",
    )
    require(
        int(cal["rff_dim"]) == EXPECTED_RFF_DIM,
        "G2C_RFF_DIM_DOES_NOT_MATCH_CONTRACT",
    )

    require(manifest.get("status") == "PASS", "FEATURE_MANIFEST_NOT_PASS")
    require(manifest.get("n_trajectories") == 12, "EXPECTED_12_FEATURE_ARTIFACTS")
    require(
        manifest.get("total_encoded_frames") == 12000,
        "EXPECTED_12000_ENCODED_FRAMES",
    )

    entries = manifest["entries"]
    require(len(entries) == 12, "FEATURE_ENTRY_COUNT_MISMATCH")

    index = {}
    for e in entries:
        key = (e["system"], int(e["replica"]))
        require(key not in index, f"DUPLICATE_FEATURE_ENTRY: {key}")
        index[key] = e

    expected_keys = {
        (system, replica)
        for system in SYSTEMS
        for replica in REPLICAS
    }
    require(set(index) == expected_keys, "FOUR_CONTEXT_THREE_REPLICA_MATRIX_INCOMPLETE")

    print(
        "system                              R  feature_sha  "
        "shape             blocks block_frames finite order status"
    )

    feature_hashes = {}

    for system in SYSTEMS:
        for replica in REPLICAS:
            e = index[(system, replica)]
            p = Path(e["path"])

            # Manifest paths are normally repo-relative.
            if not p.is_absolute():
                p = (args.repo_root / p).resolve()

            require(p.exists(), f"MISSING_FEATURE: {p}")

            actual_sha = sha256(p)
            require(
                actual_sha == e["sha256"],
                f"FEATURE_SHA256_MISMATCH: {system} R{replica}",
            )

            with np.load(p, allow_pickle=False) as z:
                f = z["residue_features"]
                frame_ids = z["frame_ids"]
                residue_index = z["residue_index"]

                shape_ok = f.shape == (
                    EXPECTED_FRAMES,
                    EXPECTED_RESIDUES,
                    EXPECTED_CHANNELS,
                )
                finite = bool(np.isfinite(f).all())
                order = bool(
                    np.array_equal(frame_ids, np.arange(EXPECTED_FRAMES))
                    and np.array_equal(
                        residue_index,
                        np.arange(EXPECTED_RESIDUES),
                    )
                )

                require(shape_ok, f"FEATURE_SHAPE_FAIL: {system} R{replica}")
                require(
                    f.dtype == np.float32,
                    f"FEATURE_DTYPE_FAIL: {system} R{replica}",
                )
                require(finite, f"NONFINITE_FEATURES: {system} R{replica}")
                require(order, f"FEATURE_ORDER_FAIL: {system} R{replica}")

                # Only reshape to audit the temporal contract.
                # No context subtraction and no contrast is computed here.
                blocked = f.reshape(
                    EXPECTED_BLOCKS,
                    EXPECTED_BLOCK_FRAMES,
                    EXPECTED_RESIDUES,
                    EXPECTED_CHANNELS,
                )
                require(
                    blocked.shape == (50, 20, 270, 128),
                    f"BLOCK_RESHAPE_FAIL: {system} R{replica}",
                )

            feature_hashes[f"{system}:R{replica}"] = actual_sha

            print(
                f"{system:35s} R{replica} "
                f"{actual_sha[:12]} "
                f"{str((1000,270,128)):17s} "
                f"{EXPECTED_BLOCKS:6d} {EXPECTED_BLOCK_FRAMES:12d} "
                f"{str(finite):6s} {str(order):5s} PASS"
            )

    weights, bias = build_rff_basis(
        EXPECTED_WIDTH,
        EXPECTED_RFF_SEED,
        EXPECTED_RFF_DIM,
    )

    require(weights.shape == (128, 256), "RFF_WEIGHT_SHAPE_FAIL")
    require(bias.shape == (256,), "RFF_BIAS_SHAPE_FAIL")
    require(weights.dtype == np.float32, "RFF_WEIGHT_DTYPE_FAIL")
    require(bias.dtype == np.float32, "RFF_BIAS_DTYPE_FAIL")
    require(np.isfinite(weights).all(), "RFF_WEIGHT_NONFINITE")
    require(np.isfinite(bias).all(), "RFF_BIAS_NONFINITE")

    print()
    print("FKG_LONG_MD_AUDIT_ONLY_PASS")
    print("contexts", len(SYSTEMS))
    print("paired_replicas", len(REPLICAS))
    print("feature_artifacts", len(entries))
    print("frames_per_trajectory", EXPECTED_FRAMES)
    print("blocks_per_trajectory", EXPECTED_BLOCKS)
    print("block_frames", EXPECTED_BLOCK_FRAMES)
    print("independent_replication_unit replica_seed")
    print("shared_node_bandwidth", repr(EXPECTED_WIDTH))
    print("base_seed", EXPECTED_BASE_SEED)
    print("actual_rff_rng_seed", EXPECTED_RFF_SEED)
    print("rff_dimension", EXPECTED_RFF_DIM)
    print("rff_weights_sha256", array_sha256(weights))
    print("rff_bias_sha256", array_sha256(bias))
    print("bandwidth_refit_on_long_md False")
    print("rff_basis_refit_on_long_md False")
    print("kernel_u2_gate REVOKED")
    print("new_gate NOT_DEFINED")
    print("biological_contrasts_computed False")


def cli():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
    )
    p.add_argument(
        "--contract",
        type=Path,
        required=True,
    )
    p.add_argument(
        "--feature-manifest",
        type=Path,
        required=True,
    )
    p.add_argument(
        "--g2c-audit",
        type=Path,
        required=True,
    )
    p.add_argument(
        "--reference",
        type=Path,
        required=True,
    )
    p.add_argument(
        "--audit-only",
        action="store_true",
        required=True,
    )
    args = p.parse_args()

    if not args.audit_only:
        raise SystemExit("THIS_VERSION_SUPPORTS_AUDIT_ONLY")

    audit(args)


if __name__ == "__main__":
    try:
        cli()
    except Exception as exc:
        print(f"FKG_LONG_MD_AUDIT_ONLY_FAIL: {type(exc).__name__}: {exc}")
        sys.exit(1)
