"""Fail-closed run-spec loader for prospective frozen-method candidate runs.

The run spec is the ONLY place where molecule identity, replica seeds, frame
geometry, asset locations and output roots may vary.  Every frozen scientific
quantity is validated against the values imported from the frozen v02 modules.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from project.pacer_fkg_v02.frozen_runner import engines

SCHEMA = "pacer.final.pipeline_run_spec.v1"

# Context-matching contract: a prospective candidate MUST reuse the historical
# paired seed group so its Delta_* contrasts are comparable with the frozen
# apo / probe_only shared controls.
FROZEN_SEEDS = {1: 27101, 2: 38201, 3: 49301}
FROZEN_REPLICAS = (1, 2, 3)
FROZEN_BLOCK_FRAMES = 20
ALLOWED_N_FRAMES = (400, 1000)  # 20 ns and 50 ns frozen observation windows
EXPECTED_CONTEXT_KEYS = ("A", "P", "C", "CP")
EXPECTED_PANEL_KEYS = ("C", "CP")

FORBIDDEN_SCOPE_FLAGS = (
    "pam_classifier_claim",
    "potency_predictor_claim",
    "pam_probability_output",
    "magnitude_threshold_defined",
    "outcome_driven_rule",
    "calibration_performed",
    "model_training_performed",
    "adapter_fitting_performed",
    "score_fusion_invented",
)


class RunSpecError(ValueError):
    """Raised when a run spec is missing, malformed or violates the freeze."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunSpecError(message)


@dataclass(frozen=True)
class Molecule:
    candidate_id: str
    canonical_smiles: str
    role: str
    chemotype: str | None = None
    provenance: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Paths:
    topology_root: Path
    production_root: Path
    phase1_cache_root: Path
    phase2a_cache_root: Path
    report_root: Path
    allowed_output_roots: tuple[Path, ...]


@dataclass(frozen=True)
class RunSpec:
    schema: str
    run_id: str
    created_at: str
    molecule: Molecule
    scope: dict[str, Any]
    contexts: dict[str, str]
    replicas: tuple[int, ...]
    seeds: dict[int, int]
    n_frames: int
    block_frames: int
    blocks_per_trajectory: int
    reference_panel: dict[str, str]
    paths: Paths
    inputs: tuple[dict[str, Any], ...]
    calibration: dict[str, Any]
    raw: dict[str, Any]
    source_path: Path

    # -- derived -----------------------------------------------------------
    @property
    def candidate_systems(self) -> tuple[str, ...]:
        return (self.contexts["C"], self.contexts["CP"])

    @property
    def all_systems(self) -> tuple[str, ...]:
        return (self.contexts["A"], self.contexts["P"], self.contexts["C"], self.contexts["CP"])

    @property
    def reference_systems(self) -> dict[str, str]:
        return dict(self.reference_panel)


def _load_json(path: Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def load_run_spec(path: Path | str, repo_root: Path | None = None) -> RunSpec:
    """Load and fully validate a run spec.  Raises `RunSpecError` on any breach."""
    source = Path(path).resolve()
    if not source.is_file():
        raise RunSpecError(f"run spec not found: {source}")
    repo = (repo_root or Path(__file__).resolve().parents[3]).resolve()
    raw = _load_json(source)

    _require(raw.get("schema") == SCHEMA, f"unsupported run-spec schema: {raw.get('schema')!r}")
    run_id = str(raw.get("run_id") or "")
    _require(bool(run_id), "run_id is required")

    # ---- scope -----------------------------------------------------------
    scope = raw.get("scope")
    _require(isinstance(scope, dict), "scope object is required")
    _require(scope.get("prospective_frozen_method_candidate_evaluation") is True,
             "scope.prospective_frozen_method_candidate_evaluation must be true")
    for flag in FORBIDDEN_SCOPE_FLAGS:
        _require(scope.get(flag, False) is False, f"scope.{flag} must be false for a frozen apply-only run")

    # ---- molecule --------------------------------------------------------
    mol_raw = raw.get("molecule")
    _require(isinstance(mol_raw, dict), "molecule object is required")
    for key in ("candidate_id", "canonical_smiles", "role"):
        _require(bool(mol_raw.get(key)), f"molecule.{key} is required")
    molecule = Molecule(
        candidate_id=str(mol_raw["candidate_id"]),
        canonical_smiles=str(mol_raw["canonical_smiles"]),
        role=str(mol_raw["role"]),
        chemotype=mol_raw.get("chemotype"),
        provenance=dict(mol_raw.get("provenance") or {}),
    )
    _require(molecule.role in {
        "prospective_computational_hypothesis",
        "known_PAM_positive_control",
        "allosteric_agonist_specificity_control",
        "strict_experimental_inactive_binding_to_be_tested",
        "same_source_PAM_positive_control",
        "regression_replay_control",
    }, f"unrecognised molecule.role: {molecule.role}")

    # ---- contexts --------------------------------------------------------
    contexts = raw.get("contexts")
    _require(isinstance(contexts, dict), "contexts object is required")
    _require(tuple(sorted(contexts)) == tuple(sorted(EXPECTED_CONTEXT_KEYS)),
             f"contexts must define exactly {EXPECTED_CONTEXT_KEYS}")
    for key in EXPECTED_CONTEXT_KEYS:
        _require(isinstance(contexts[key], str) and contexts[key], f"contexts.{key} must be a non-empty string")
    values = [contexts[k] for k in EXPECTED_CONTEXT_KEYS]
    _require(len(set(values)) == len(values), "context system names must be distinct")

    # ---- replicas / seeds (frozen) --------------------------------------
    replicas = tuple(int(r) for r in raw.get("replicas", []))
    _require(replicas == FROZEN_REPLICAS, f"replicas must be exactly {FROZEN_REPLICAS}, got {replicas}")
    seeds_raw = raw.get("seeds") or {}
    seeds = {int(k): int(v) for k, v in seeds_raw.items()}
    _require(seeds == FROZEN_SEEDS,
             f"seeds must equal the frozen paired seed group {FROZEN_SEEDS} for context matching, got {seeds}")

    # ---- frame geometry --------------------------------------------------
    frames = raw.get("frames")
    _require(isinstance(frames, dict), "frames object is required")
    n_frames = int(frames.get("n_frames", -1))
    block_frames = int(frames.get("block_frames", -1))
    blocks = int(frames.get("blocks_per_trajectory", -1))
    _require(n_frames in ALLOWED_N_FRAMES, f"n_frames must be one of {ALLOWED_N_FRAMES}, got {n_frames}")
    _require(block_frames == FROZEN_BLOCK_FRAMES, f"block_frames must be {FROZEN_BLOCK_FRAMES}, got {block_frames}")
    _require(blocks * block_frames == n_frames,
             f"blocks_per_trajectory * block_frames must equal n_frames ({blocks}*{block_frames} != {n_frames})")

    # ---- reference panel -------------------------------------------------
    panel = raw.get("reference_panel")
    _require(isinstance(panel, dict), "reference_panel object is required")
    _require(tuple(sorted(panel)) == tuple(sorted(EXPECTED_PANEL_KEYS)),
             f"reference_panel must define exactly {EXPECTED_PANEL_KEYS}")
    for key in EXPECTED_PANEL_KEYS:
        _require(isinstance(panel[key], str) and panel[key], f"reference_panel.{key} must be a non-empty string")

    # ---- calibration policy ---------------------------------------------
    calibration = dict(raw.get("calibration") or {})
    usage = str(calibration.get("usage", ""))
    _require(usage in {"forbidden", "verify_only"},
             f"calibration.usage must be 'forbidden' or 'verify_only', got {usage!r}")
    _require(calibration.get("r2_fitting_allowed", False) is False, "calibration.r2_fitting_allowed must be false")

    # ---- paths -----------------------------------------------------------
    paths_raw = raw.get("paths")
    _require(isinstance(paths_raw, dict), "paths object is required")

    def _resolve(value: Any) -> Path:
        p = Path(str(value))
        return p if p.is_absolute() else (repo / p)

    allowed = tuple(_resolve(x) for x in (paths_raw.get("allowed_output_roots") or []))
    _require(len(allowed) >= 1, "paths.allowed_output_roots must list at least one root")
    paths = Paths(
        topology_root=_resolve(paths_raw["topology_root"]),
        production_root=_resolve(paths_raw["production_root"]),
        phase1_cache_root=_resolve(paths_raw["phase1_cache_root"]),
        phase2a_cache_root=_resolve(paths_raw["phase2a_cache_root"]),
        report_root=_resolve(paths_raw["report_root"]),
        allowed_output_roots=allowed,
    )

    # ---- inputs ----------------------------------------------------------
    inputs = tuple(dict(x) for x in (raw.get("inputs") or []))
    for item in inputs:
        _require(bool(item.get("role")), "every input needs a role")
        _require(bool(item.get("path")), "every input needs a path")
        value = str(item.get("sha256", ""))
        _require(len(value) == 64 and all(c in "0123456789abcdef" for c in value.lower()),
                 f"input {item.get('path')} needs a 64-hex sha256")

    # ---- frozen anchor present ------------------------------------------
    anchor = raw.get("frozen_anchor") or {}
    _require(anchor.get("freeze_sha256") == engines.REQUIRED_FREEZE_SHA256,
             f"frozen_anchor.freeze_sha256 must equal {engines.REQUIRED_FREEZE_SHA256}")

    return RunSpec(
        schema=SCHEMA, run_id=run_id, created_at=str(raw.get("created_at", "")),
        molecule=molecule, scope=scope, contexts={k: str(contexts[k]) for k in EXPECTED_CONTEXT_KEYS},
        replicas=replicas, seeds=seeds, n_frames=n_frames, block_frames=block_frames,
        blocks_per_trajectory=blocks, reference_panel={k: str(panel[k]) for k in EXPECTED_PANEL_KEYS},
        paths=paths, inputs=inputs, calibration=calibration, raw=raw, source_path=source,
    )


def authenticate_inputs(spec: RunSpec, repo_root: Path | None = None) -> dict[str, Any]:
    """Recompute every declared input hash.  Returns a per-input record."""
    repo = (repo_root or Path(__file__).resolve().parents[3]).resolve()
    records: list[dict[str, Any]] = []
    for item in spec.inputs:
        p = Path(str(item["path"]))
        full = p if p.is_absolute() else (repo / p)
        if not full.is_file():
            raise RunSpecError(f"declared input is missing: {full}")
        actual = sha256_file(full)
        expected = str(item["sha256"]).lower()
        records.append({"role": item["role"], "path": str(full), "expected": expected,
                        "actual": actual, "match": actual == expected})
    mismatched = [r for r in records if not r["match"]]
    if mismatched:
        raise RunSpecError("input hash mismatch: " + ", ".join(r["path"] for r in mismatched))
    return {"inputs": records, "count": len(records), "all_match": True}


def assert_frozen_spec(spec: RunSpec) -> None:
    """Re-assert the imported frozen numerics and the spec's frozen fields."""
    engines.assert_frozen_numerics()
    _require(spec.block_frames == engines.P2_BLOCK_FRAMES, "block_frames drifted from the frozen Phase-2 value")
    _require(engines.N_RESIDUES == engines.P2_N_RESIDUES, "frozen residue count is inconsistent between modules")
    _require(engines.RFF_FEATURES == engines.P2_RFF_FEATURES, "frozen RFF feature count is inconsistent")
