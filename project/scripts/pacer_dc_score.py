#!/usr/bin/env python
"""Score matched dual-context M4 trajectories without inventing a PAM probability.

Input is a long-form evidence CSV.  PACER-DC reports a vector of mechanistically
separate endpoints and performs Pareto sorting only when every required context
is present.  Replicas, not trajectory frames, are the independent units.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


CONTEXTS = {"candidate_probe", "candidate_no_probe", "probe_only", "apo"}
METRICS = {
    "coupling_coordinate",
    "orthosteric_pose_rmsd_A",
    "binding_compatibility",
}
REQUIRED = {
    "CoupledShift": (
        ("candidate_probe", "coupling_coordinate"),
        ("probe_only", "coupling_coordinate"),
    ),
    "OrthostericStabilization": (
        ("candidate_probe", "orthosteric_pose_rmsd_A"),
        ("probe_only", "orthosteric_pose_rmsd_A"),
    ),
    "IntrinsicActivationRisk": (
        ("candidate_no_probe", "coupling_coordinate"),
        ("apo", "coupling_coordinate"),
    ),
    "BindingCompatibility": (
        ("candidate_probe", "binding_compatibility"),
        ("candidate_no_probe", "binding_compatibility"),
    ),
}


def _validate(df: pd.DataFrame) -> None:
    required_columns = {
        "candidate_id",
        "context",
        "replicate_id",
        "metric",
        "value",
        "evidence_level",
        "source_id",
    }
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    bad_context = sorted(set(df.context.astype(str)) - CONTEXTS)
    bad_metric = sorted(set(df.metric.astype(str)) - METRICS)
    if bad_context:
        raise ValueError(f"Unknown contexts: {bad_context}")
    if bad_metric:
        raise ValueError(f"Unknown metrics: {bad_metric}")
    if not np.isfinite(pd.to_numeric(df.value, errors="coerce")).all():
        raise ValueError("Every evidence row must have a finite numeric value")
    duplicates = df.duplicated(
        ["candidate_id", "context", "replicate_id", "metric"], keep=False
    )
    if duplicates.any():
        cols = ["candidate_id", "context", "replicate_id", "metric"]
        raise ValueError("Duplicate evidence keys:\n" + df.loc[duplicates, cols].to_string(index=False))


def _values(df: pd.DataFrame, context: str, metric: str) -> np.ndarray:
    rows = df[(df.context == context) & (df.metric == metric)]
    return rows.value.to_numpy(dtype=float)


def _delta(df: pd.DataFrame, left: tuple[str, str], right: tuple[str, str]) -> float:
    a = _values(df, *left)
    b = _values(df, *right)
    if not len(a) or not len(b):
        return np.nan
    return float(a.mean() - b.mean())


def _paired_values(
    df: pd.DataFrame, left: tuple[str, str], right: tuple[str, str]
) -> np.ndarray:
    a = df[(df.context == left[0]) & (df.metric == left[1])][
        ["replicate_id", "value"]
    ].rename(columns={"value": "left"})
    b = df[(df.context == right[0]) & (df.metric == right[1])][
        ["replicate_id", "value"]
    ].rename(columns={"value": "right"})
    paired = a.merge(b, on="replicate_id", how="inner")
    return (paired.left - paired.right).to_numpy(dtype=float)


def _bootstrap_delta(
    df: pd.DataFrame,
    left: tuple[str, str],
    right: tuple[str, str],
    rng: np.random.Generator,
    iterations: int,
) -> tuple[float, float]:
    paired = _paired_values(df, left, right)
    if len(paired) >= 2:
        draws = np.empty(iterations, dtype=float)
        for i in range(iterations):
            draws[i] = rng.choice(paired, len(paired), replace=True).mean()
        lo, hi = np.quantile(draws, [0.025, 0.975])
        return float(lo), float(hi)
    a = _values(df, *left)
    b = _values(df, *right)
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    draws = np.empty(iterations, dtype=float)
    for i in range(iterations):
        draws[i] = rng.choice(a, len(a), replace=True).mean() - rng.choice(
            b, len(b), replace=True
        ).mean()
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return float(lo), float(hi)


def _has(df: pd.DataFrame, context: str, metric: str) -> bool:
    return bool(((df.context == context) & (df.metric == metric)).any())


def _has_trajectory(df: pd.DataFrame, context: str, metric: str) -> bool:
    return bool(
        (
            (df.context == context)
            & (df.metric == metric)
            & (df.evidence_level == "trajectory")
        ).any()
    )


def score_candidate(
    candidate: str, df: pd.DataFrame, seed: int, iterations: int
) -> dict[str, object]:
    own = df[df.candidate_id == candidate].copy()
    shared = df[
        (df.candidate_id == "SHARED_CONTROL")
        & df.context.isin(["probe_only", "apo"])
    ].copy()
    # Candidate-specific controls take precedence; otherwise reuse the matched
    # shared control simulation rather than rerunning an identical system.
    for context in ("probe_only", "apo"):
        own_metrics = set(own.loc[own.context == context, "metric"])
        add = shared[(shared.context == context) & ~shared.metric.isin(own_metrics)]
        own = pd.concat([own, add], ignore_index=True)
    cdf = own
    rng = np.random.default_rng(seed)

    coupled = _delta(
        cdf,
        ("candidate_probe", "coupling_coordinate"),
        ("probe_only", "coupling_coordinate"),
    )
    coupled_ci = _bootstrap_delta(
        cdf,
        ("candidate_probe", "coupling_coordinate"),
        ("probe_only", "coupling_coordinate"),
        rng,
        iterations,
    )
    # Lower pose RMSD is better, so reverse the subtraction.
    stabilization = _delta(
        cdf,
        ("probe_only", "orthosteric_pose_rmsd_A"),
        ("candidate_probe", "orthosteric_pose_rmsd_A"),
    )
    stabilization_ci = _bootstrap_delta(
        cdf,
        ("probe_only", "orthosteric_pose_rmsd_A"),
        ("candidate_probe", "orthosteric_pose_rmsd_A"),
        rng,
        iterations,
    )
    intrinsic = _delta(
        cdf,
        ("candidate_no_probe", "coupling_coordinate"),
        ("apo", "coupling_coordinate"),
    )
    intrinsic_ci = _bootstrap_delta(
        cdf,
        ("candidate_no_probe", "coupling_coordinate"),
        ("apo", "coupling_coordinate"),
        rng,
        iterations,
    )

    binding_rows = cdf[
        (cdf.metric == "binding_compatibility")
        & cdf.context.isin(["candidate_probe", "candidate_no_probe"])
    ]
    binding = float(binding_rows.value.mean()) if len(binding_rows) else np.nan

    missing: list[str] = []
    for endpoint, requirements in REQUIRED.items():
        for context, metric in requirements:
            if not _has(cdf, context, metric):
                missing.append(f"{endpoint}:{context}/{metric}:absent")
            elif not _has_trajectory(cdf, context, metric):
                missing.append(f"{endpoint}:{context}/{metric}:not_trajectory")

    trajectory_contexts = sorted(
        set(cdf.loc[cdf.evidence_level == "trajectory", "context"].astype(str))
    )
    complete = not missing
    return {
        "candidate_id": candidate,
        "CoupledShift": coupled,
        "CoupledShift_CI95_low": coupled_ci[0],
        "CoupledShift_CI95_high": coupled_ci[1],
        "OrthostericStabilization_A": stabilization,
        "OrthostericStabilization_CI95_low": stabilization_ci[0],
        "OrthostericStabilization_CI95_high": stabilization_ci[1],
        "IntrinsicActivationRisk": intrinsic,
        "IntrinsicActivationRisk_CI95_low": intrinsic_ci[0],
        "IntrinsicActivationRisk_CI95_high": intrinsic_ci[1],
        "BindingCompatibility": binding,
        "n_evidence_rows": int(len(cdf)),
        "n_paired_coupling_replicas": int(
            len(
                _paired_values(
                    cdf,
                    ("candidate_probe", "coupling_coordinate"),
                    ("probe_only", "coupling_coordinate"),
                )
            )
        ),
        "n_paired_stability_replicas": int(
            len(
                _paired_values(
                    cdf,
                    ("probe_only", "orthosteric_pose_rmsd_A"),
                    ("candidate_probe", "orthosteric_pose_rmsd_A"),
                )
            )
        ),
        "n_paired_intrinsic_replicas": int(
            len(
                _paired_values(
                    cdf,
                    ("candidate_no_probe", "coupling_coordinate"),
                    ("apo", "coupling_coordinate"),
                )
            )
        ),
        "trajectory_contexts": ";".join(trajectory_contexts),
        "missing_requirements": ";".join(missing),
        "evidence_complete": complete,
        "pam_probability_allowed": False,
        "interpretation": (
            "complete_mechanistic_vector_not_experimental_confirmation"
            if complete
            else "partial_evidence_do_not_classify_or_rank_as_PAM"
        ),
    }


def pareto_front(scored: pd.DataFrame) -> pd.Series:
    """Return nondominated rows; higher is better except intrinsic risk."""
    eligible = scored.evidence_complete.astype(bool)
    front = pd.Series(False, index=scored.index)
    cols = [
        "CoupledShift",
        "OrthostericStabilization_A",
        "IntrinsicActivationRisk",
        "BindingCompatibility",
    ]
    for i in scored.index[eligible]:
        x = scored.loc[i, cols].to_numpy(float, copy=True)
        x[2] *= -1.0
        dominated = False
        for j in scored.index[eligible & (scored.index != i)]:
            y = scored.loc[j, cols].to_numpy(float, copy=True)
            y[2] *= -1.0
            if np.all(y >= x) and np.any(y > x):
                dominated = True
                break
        front.loc[i] = not dominated
    return front


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=Path, required=True, nargs="+",
        help="One or more compatible replica-evidence CSV files",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--bootstrap", type=int, default=5000)
    args = parser.parse_args()

    df = pd.concat([pd.read_csv(path) for path in args.input], ignore_index=True)
    _validate(df)
    df["value"] = df.value.astype(float)
    candidates = sorted(
        c for c in df.candidate_id.astype(str).unique() if c != "SHARED_CONTROL"
    )
    rows = [
        score_candidate(c, df, args.seed + i, args.bootstrap)
        for i, c in enumerate(candidates)
    ]
    scored = pd.DataFrame(rows)
    scored["pareto_eligible"] = scored.evidence_complete.astype(bool)
    scored["pareto_front"] = pareto_front(scored)

    args.outdir.mkdir(parents=True, exist_ok=True)
    scored.to_csv(args.outdir / "pacer_dc_evidence_vectors.csv", index=False)
    audit = {
        "method": "PACER-DC",
        "independent_unit": "trajectory replica or independent structure; never frame",
        "uncertainty_design": "paired replica bootstrap when replicate_id matches; otherwise unpaired replica bootstrap",
        "candidate_count": int(len(scored)),
        "complete_candidate_count": int(scored.evidence_complete.sum()),
        "pareto_candidate_count": int(scored.pareto_front.sum()),
        "pam_probability_emitted": False,
        "claim_boundary": (
            "PACER-DC is a mechanistic evidence vector. Only complete matched-context "
            "candidates enter Pareto comparison; wet functional assays are required to confirm PAM activity."
        ),
    }
    (args.outdir / "audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
