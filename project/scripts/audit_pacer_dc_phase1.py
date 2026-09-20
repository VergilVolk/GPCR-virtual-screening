#!/usr/bin/env python
"""Audit PACER-DC readiness without manufacturing missing evidence.

The audit separates simulation, evidence scoring, and supervised-learning
gates. Trajectory frames never count as independent functional systems.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ALLOWED_CONTEXTS = {"candidate_probe", "candidate_no_probe", "probe_only", "apo"}
LIGAND_CONTEXTS = {"candidate_probe", "candidate_no_probe"}
CONTROL_CONTEXTS = {"probe_only", "apo"}


def audit_manifest(manifest: pd.DataFrame) -> dict[str, object]:
    required = {
        "run_id", "candidate_id", "context", "replicate_id", "role",
        "paired_seed_group", "status",
    }
    missing_columns = sorted(required - set(manifest.columns))
    errors: list[str] = []
    warnings: list[str] = []
    if missing_columns:
        return {"valid": False, "errors": [f"missing columns: {missing_columns}"]}
    if manifest.run_id.duplicated().any():
        errors.append("run_id is not unique")
    unknown = sorted(set(manifest.context.astype(str)) - ALLOWED_CONTEXTS)
    if unknown:
        errors.append(f"unknown contexts: {unknown}")

    shared = manifest[manifest.candidate_id == "SHARED_CONTROL"]
    shared_contexts = set(shared.context.astype(str))
    if shared_contexts != CONTROL_CONTEXTS:
        errors.append(
            f"shared controls must contain {sorted(CONTROL_CONTEXTS)}, got {sorted(shared_contexts)}"
        )
    shared_replicas = {
        context: set(shared.loc[shared.context == context, "replicate_id"].astype(str))
        for context in CONTROL_CONTEXTS
    }
    if shared_replicas.get("probe_only") != shared_replicas.get("apo"):
        errors.append("probe_only and apo shared controls do not have matched replicas")

    candidates = sorted(set(manifest.candidate_id.astype(str)) - {"SHARED_CONTROL"})
    candidate_audit: dict[str, object] = {}
    roles: dict[str, str] = {}
    for candidate in candidates:
        own = manifest[manifest.candidate_id.astype(str) == candidate]
        contexts = set(own.context.astype(str))
        reps = {
            context: set(own.loc[own.context == context, "replicate_id"].astype(str))
            for context in LIGAND_CONTEXTS
        }
        own_errors: list[str] = []
        if contexts != LIGAND_CONTEXTS:
            own_errors.append(f"requires {sorted(LIGAND_CONTEXTS)}, got {sorted(contexts)}")
        if reps.get("candidate_probe") != reps.get("candidate_no_probe"):
            own_errors.append("ligand contexts do not have matched replicas")
        if reps.get("candidate_probe") != shared_replicas.get("probe_only"):
            own_errors.append("ligand and shared-control replica ids are not matched")
        if len(reps.get("candidate_probe", set())) < 3:
            own_errors.append("fewer than three paired replicas per ligand context")
        unique_roles = sorted(set(own.role.astype(str)))
        if len(unique_roles) != 1:
            own_errors.append(f"candidate has inconsistent roles: {unique_roles}")
        else:
            roles[candidate] = unique_roles[0]
        candidate_audit[candidate] = {
            "contexts": sorted(contexts),
            "replicate_count_per_ligand_context": {key: len(value) for key, value in reps.items()},
            "errors": own_errors,
        }
        errors.extend(f"{candidate}: {message}" for message in own_errors)

    statuses = manifest.status.astype(str).value_counts().to_dict()
    if statuses.get("not_started", 0) == len(manifest):
        warnings.append("all manifest runs are not_started")
    positive_controls = [c for c, role in roles.items() if "positive_control" in role.lower()]
    agonist_controls = [c for c, role in roles.items() if "agonist" in role.lower()]
    neutral_controls = [
        c for c, role in roles.items()
        if "neutral" in role.lower() or "no_function" in role.lower()
    ]
    return {
        "valid": not errors,
        "row_count": int(len(manifest)),
        "candidate_count_excluding_shared_control": len(candidates),
        "statuses": statuses,
        "positive_control_candidates": positive_controls,
        "allosteric_agonist_controls": agonist_controls,
        "neutral_or_no_function_controls": neutral_controls,
        "candidate_details": candidate_audit,
        "errors": errors,
        "warnings": warnings,
    }


def audit_evidence(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"exists": False, "complete_candidate_count": 0, "warning": "evidence file is absent"}
    evidence = pd.read_csv(path)
    required = {
        "candidate_id", "context", "replicate_id", "metric", "value",
        "evidence_level", "source_id",
    }
    missing = sorted(required - set(evidence.columns))
    if missing:
        return {"exists": True, "valid": False, "errors": [f"missing columns: {missing}"]}
    coverage = (
        evidence.groupby(["candidate_id", "context", "metric"], dropna=False)
        .agg(
            replica_count=("replicate_id", "nunique"),
            evidence_levels=("evidence_level", lambda x: ";".join(sorted(set(map(str, x))))),
            source_count=("source_id", "nunique"),
        )
        .reset_index()
    )
    candidates = sorted(set(evidence.candidate_id.astype(str)) - {"SHARED_CONTROL"})
    required_pairs = {
        ("candidate_probe", "coupling_coordinate"),
        ("probe_only", "coupling_coordinate"),
        ("candidate_probe", "orthosteric_pose_rmsd_A"),
        ("probe_only", "orthosteric_pose_rmsd_A"),
        ("candidate_no_probe", "coupling_coordinate"),
        ("apo", "coupling_coordinate"),
        ("candidate_probe", "binding_compatibility"),
        ("candidate_no_probe", "binding_compatibility"),
    }
    complete: list[str] = []
    candidate_missing: dict[str, list[str]] = {}
    shared = evidence[evidence.candidate_id.astype(str) == "SHARED_CONTROL"]
    for candidate in candidates:
        own = evidence[evidence.candidate_id.astype(str) == candidate]
        own_trajectory = own[own.evidence_level.astype(str) == "trajectory"]
        shared_trajectory = shared[shared.evidence_level.astype(str) == "trajectory"]
        available = set(zip(own_trajectory.context.astype(str), own_trajectory.metric.astype(str)))
        available |= set(zip(shared_trajectory.context.astype(str), shared_trajectory.metric.astype(str)))
        missing_pairs = sorted(required_pairs - available)
        candidate_missing[candidate] = [f"{a}/{b}" for a, b in missing_pairs]
        if not missing_pairs:
            complete.append(candidate)
    return {
        "exists": True,
        "valid": True,
        "row_count": int(len(evidence)),
        "candidate_count": len(candidates),
        "complete_candidate_count": len(complete),
        "complete_candidates": complete,
        "missing_by_candidate": candidate_missing,
        "coverage": coverage.to_dict(orient="records"),
    }


def render_report(audit: dict[str, object]) -> str:
    manifest = audit["manifest"]
    evidence = audit["evidence"]
    environment = audit["environment"]
    lines = [
        "# PACER-DC 第一阶段闸门报告", "",
        f"- 清单结构有效：**{manifest.get('valid', False)}**",
        f"- MD 生产环境就绪：**{environment.get('production_ready', False)}**",
        f"- 完整四上下文候选数：**{evidence.get('complete_candidate_count', 0)}**",
        f"- 监督式轨迹模型可开发：**{audit['gates']['supervised_representation_learning_ready']}**",
        f"- PACER-DC 候选排序可执行：**{audit['gates']['candidate_scoring_ready']}**",
        "", "## 结论", "",
    ]
    if audit["gates"]["simulation_ready"]:
        lines.append("MD 清单与环境均通过，可以启动冻结的试验。")
    else:
        lines.append("当前不能启动可重复的生产 MD；先解决环境或清单中的阻断项。")
    if not audit["gates"]["supervised_representation_learning_ready"]:
        lines.append("当前独立功能类别和完整上下文不足，不能训练或评价 Triplet/PAM 分类器；轨迹帧数不计作独立样本数。")
    if not audit["gates"]["candidate_scoring_ready"]:
        lines.append("当前没有候选具备完整四上下文证据，不输出 PAM 概率或候选优劣结论。")
    lines.extend(["", "## 阻断项", ""])
    lines.extend(f"- {item}" for item in (audit["blockers"] or ["无"]))
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=ROOT / "config" / "pacer_dc_pilot_md_manifest.csv")
    parser.add_argument("--evidence", type=Path, default=ROOT / "data" / "benchmarks" / "pacer_dc_reference_evidence.csv")
    parser.add_argument("--environment-audit", type=Path, default=ROOT / "results" / "pacer_dc_md_environment_audit.json")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results" / "pacer_dc_phase1_audit")
    args = parser.parse_args()

    manifest = audit_manifest(pd.read_csv(args.manifest))
    evidence = audit_evidence(args.evidence)
    if args.environment_audit.exists():
        environment = json.loads(args.environment_audit.read_text(encoding="utf-8"))
    else:
        environment = {"production_ready": False, "blocking_reason": "environment audit is absent"}
    roles_ok = (
        len(manifest.get("positive_control_candidates", [])) >= 3
        and len(manifest.get("allosteric_agonist_controls", [])) >= 1
        and len(manifest.get("neutral_or_no_function_controls", [])) >= 1
    )
    complete_count = int(evidence.get("complete_candidate_count", 0))
    gates = {
        "simulation_ready": bool(manifest.get("valid")) and bool(environment.get("production_ready")),
        "candidate_scoring_ready": complete_count > 0,
        "supervised_representation_learning_ready": roles_ok and complete_count >= 5,
    }
    blockers: list[str] = list(manifest.get("errors", []))
    if not environment.get("production_ready", False):
        blockers.append(str(environment.get("blocking_reason", "MD environment is not production-ready")))
    if complete_count == 0:
        blockers.append("no candidate has all required matched-context evidence")
    if len(manifest.get("positive_control_candidates", [])) < 3:
        blockers.append("fewer than three independent PAM positive-control chemotypes are in the MD ledger")
    if not manifest.get("neutral_or_no_function_controls", []):
        blockers.append("no confirmed binding-but-non-PAM control is in the MD ledger")
    audit = {
        "method": "PACER-DC-v1", "manifest": manifest, "evidence": evidence,
        "environment": environment, "gates": gates, "blockers": blockers,
        "claim_boundary": "This audit reports readiness only; it is not model performance or PAM confirmation.",
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    (args.outdir / "REPORT.md").write_text(render_report(audit), encoding="utf-8")
    print(json.dumps({"gates": gates, "blockers": blockers}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
