from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from pacer_dc_score import _validate, pareto_front, score_candidate  # noqa: E402


def test_complete_vector_and_directionality() -> None:
    rows = []
    for replica in ("r1", "r2", "r3"):
        rows.extend(
            [
                ["X", "candidate_probe", replica, "coupling_coordinate", 1.2, "trajectory", "test"],
                ["X", "probe_only", replica, "coupling_coordinate", 0.4, "trajectory", "test"],
                ["X", "candidate_probe", replica, "orthosteric_pose_rmsd_A", 1.0, "trajectory", "test"],
                ["X", "probe_only", replica, "orthosteric_pose_rmsd_A", 2.0, "trajectory", "test"],
                ["X", "candidate_no_probe", replica, "coupling_coordinate", 0.2, "trajectory", "test"],
                ["X", "apo", replica, "coupling_coordinate", 0.1, "trajectory", "test"],
                ["X", "candidate_probe", replica, "binding_compatibility", 0.8, "trajectory", "test"],
                ["X", "candidate_no_probe", replica, "binding_compatibility", 0.6, "trajectory", "test"],
            ]
        )
    df = pd.DataFrame(
        rows,
        columns=[
            "candidate_id", "context", "replicate_id", "metric", "value",
            "evidence_level", "source_id",
        ],
    )
    _validate(df)
    result = score_candidate("X", df, seed=1, iterations=100)
    assert result["evidence_complete"] is True
    assert abs(result["CoupledShift"] - 0.8) < 1e-12
    assert abs(result["OrthostericStabilization_A"] - 1.0) < 1e-12
    assert abs(result["IntrinsicActivationRisk"] - 0.1) < 1e-12
    assert abs(result["BindingCompatibility"] - 0.7) < 1e-12
    assert result["pam_probability_allowed"] is False


def test_incomplete_candidate_cannot_enter_pareto() -> None:
    scored = pd.DataFrame(
        [
            {
                "candidate_id": "partial",
                "CoupledShift": 99.0,
                "OrthostericStabilization_A": 99.0,
                "IntrinsicActivationRisk": -99.0,
                "BindingCompatibility": 99.0,
                "evidence_complete": False,
            },
            {
                "candidate_id": "complete",
                "CoupledShift": 1.0,
                "OrthostericStabilization_A": 1.0,
                "IntrinsicActivationRisk": 0.1,
                "BindingCompatibility": 1.0,
                "evidence_complete": True,
            },
        ]
    )
    front = pareto_front(scored)
    assert not bool(front.iloc[0])
    assert bool(front.iloc[1])


def test_shared_controls_are_reused_without_becoming_candidate() -> None:
    rows = [
        ["X", "candidate_probe", "r1", "coupling_coordinate", 1.0, "trajectory", "test"],
        ["X", "candidate_probe", "r1", "orthosteric_pose_rmsd_A", 1.0, "trajectory", "test"],
        ["X", "candidate_probe", "r1", "binding_compatibility", 0.8, "trajectory", "test"],
        ["X", "candidate_no_probe", "r1", "coupling_coordinate", 0.2, "trajectory", "test"],
        ["X", "candidate_no_probe", "r1", "binding_compatibility", 0.7, "trajectory", "test"],
        ["SHARED_CONTROL", "probe_only", "r1", "coupling_coordinate", 0.4, "trajectory", "test"],
        ["SHARED_CONTROL", "probe_only", "r1", "orthosteric_pose_rmsd_A", 2.0, "trajectory", "test"],
        ["SHARED_CONTROL", "apo", "r1", "coupling_coordinate", 0.1, "trajectory", "test"],
    ]
    df = pd.DataFrame(
        rows,
        columns=[
            "candidate_id", "context", "replicate_id", "metric", "value",
            "evidence_level", "source_id",
        ],
    )
    result = score_candidate("X", df, seed=1, iterations=20)
    assert result["evidence_complete"] is True
    assert abs(result["CoupledShift"] - 0.6) < 1e-12
    assert abs(result["IntrinsicActivationRisk"] - 0.1) < 1e-12


def test_static_evidence_cannot_complete_dynamic_vector() -> None:
    rows = []
    for context, metric, value in [
        ("candidate_probe", "coupling_coordinate", 1.0),
        ("probe_only", "coupling_coordinate", 0.4),
        ("candidate_probe", "orthosteric_pose_rmsd_A", 1.0),
        ("probe_only", "orthosteric_pose_rmsd_A", 2.0),
        ("candidate_no_probe", "coupling_coordinate", 0.2),
        ("apo", "coupling_coordinate", 0.1),
        ("candidate_probe", "binding_compatibility", 0.8),
        ("candidate_no_probe", "binding_compatibility", 0.7),
    ]:
        rows.append(["X", context, "r1", metric, value, "static_structure", "test"])
    df = pd.DataFrame(
        rows,
        columns=[
            "candidate_id", "context", "replicate_id", "metric", "value",
            "evidence_level", "source_id",
        ],
    )
    result = score_candidate("X", df, seed=1, iterations=20)
    assert result["evidence_complete"] is False
    assert "not_trajectory" in result["missing_requirements"]
