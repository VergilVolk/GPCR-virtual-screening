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


def test_pareto_front_domination_requires_pairwise_comparison() -> None:
    """Domination is evaluated pairwise, so a multi-eligible case is required.

    ``test_incomplete_candidate_cannot_enter_pareto`` holds a single eligible
    row, so its inner comparison loop never executes and cannot detect a
    regression in the row-vector handling.  pandas 3.x returns read-only row
    views, so negating the intrinsic-risk axis in place fails unless the row
    vector is copied first.
    """
    scored = pd.DataFrame(
        [
            # A dominates B on every axis.
            {"candidate_id": "A", "CoupledShift": 2.0, "OrthostericStabilization_A": 2.0,
             "IntrinsicActivationRisk": 0.10, "BindingCompatibility": 2.0,
             "evidence_complete": True},
            {"candidate_id": "B", "CoupledShift": 1.0, "OrthostericStabilization_A": 1.0,
             "IntrinsicActivationRisk": 0.20, "BindingCompatibility": 1.0,
             "evidence_complete": True},
            # C trades coupling against stability and compatibility, but carries
            # the lowest intrinsic risk, so A and C do not dominate each other.
            {"candidate_id": "C", "CoupledShift": 3.0, "OrthostericStabilization_A": 0.5,
             "IntrinsicActivationRisk": 0.05, "BindingCompatibility": 1.5,
             "evidence_complete": True},
            # D matches A except for a worse intrinsic risk, so A dominates it.
            {"candidate_id": "D", "CoupledShift": 2.0, "OrthostericStabilization_A": 2.0,
             "IntrinsicActivationRisk": 0.30, "BindingCompatibility": 2.0,
             "evidence_complete": True},
        ]
    )
    front = pareto_front(scored)
    assert sorted(scored.loc[front, "candidate_id"]) == ["A", "C"]


def test_intrinsic_activation_risk_is_minimized_not_maximized() -> None:
    """Lower IntrinsicActivationRisk must win; the axis is negated internally."""
    scored = pd.DataFrame(
        [
            {"candidate_id": "low_risk", "CoupledShift": 1.0,
             "OrthostericStabilization_A": 1.0, "IntrinsicActivationRisk": 0.10,
             "BindingCompatibility": 1.0, "evidence_complete": True},
            {"candidate_id": "high_risk", "CoupledShift": 1.0,
             "OrthostericStabilization_A": 1.0, "IntrinsicActivationRisk": 0.90,
             "BindingCompatibility": 1.0, "evidence_complete": True},
        ]
    )
    front = pareto_front(scored)
    assert sorted(scored.loc[front, "candidate_id"]) == ["low_risk"]


def test_ineligible_candidate_is_excluded_from_every_comparison() -> None:
    """An ineligible row must neither enter the front nor dominate others."""
    scored = pd.DataFrame(
        [
            # Dominates everything numerically, but carries no complete evidence.
            {"candidate_id": "ineligible", "CoupledShift": 99.0,
             "OrthostericStabilization_A": 99.0, "IntrinsicActivationRisk": -99.0,
             "BindingCompatibility": 99.0, "evidence_complete": False},
            {"candidate_id": "A", "CoupledShift": 2.0, "OrthostericStabilization_A": 2.0,
             "IntrinsicActivationRisk": 0.10, "BindingCompatibility": 2.0,
             "evidence_complete": True},
            {"candidate_id": "C", "CoupledShift": 3.0, "OrthostericStabilization_A": 0.5,
             "IntrinsicActivationRisk": 0.05, "BindingCompatibility": 1.5,
             "evidence_complete": True},
        ]
    )
    front = pareto_front(scored)
    assert sorted(scored.loc[front, "candidate_id"]) == ["A", "C"]
