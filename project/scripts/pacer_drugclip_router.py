#!/usr/bin/env python3
"""Deploy the frozen PACER DrugCLIP generation router.

The router consumes scores produced for the *same target--molecule rows* by:

1. the 2023 DrugCLIP backbone with the GPCR-LOTO adapter; and
2. the 2026 Science DrugCLIP backbone with the 13-target family-aug adapter.

It does not load either neural network.  That separation is intentional: model
inference and score routing are independently auditable stages.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "pair_id", "target", "score_2023_gpcr_loto", "score_2026_13t_famaug"
}
M4_ALIASES = frozenset({"M4", "M4R", "CHRM4"})


@dataclass(frozen=True)
class RouterConfig:
    """Frozen deployment configuration.

    Weights apply only to non-M4 targets. M4 always routes to the 2023 GPCR
    adapter because the 2026 family-aug branch showed M4 negative transfer.
    """

    old_weight: float = 0.5
    new_weight: float = 0.5
    old_model_id: str = "drugclip2023_gpcr_loto"
    new_model_id: str = "drugclip2026_13t_family_aug"
    protocol_id: str = "pacer_drugclip_m4_safe_router_v01"

    def validate(self) -> None:
        if self.old_weight < 0 or self.new_weight < 0:
            raise ValueError("Router weights must be non-negative")
        if not np.isclose(self.old_weight + self.new_weight, 1.0):
            raise ValueError("Router weights must sum to 1")

    @property
    def is_frozen_v01(self) -> bool:
        return np.isclose(self.old_weight, 0.5) and np.isclose(self.new_weight, 0.5)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_input(frame: pd.DataFrame) -> None:
    missing = REQUIRED_COLUMNS - set(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if frame.empty:
        raise ValueError("Input table is empty")
    if frame.pair_id.isna().any() or frame.pair_id.duplicated().any():
        raise ValueError("pair_id must be non-null and unique")
    if frame.target.isna().any():
        raise ValueError("target must be non-null")
    score_columns = ["score_2023_gpcr_loto", "score_2026_13t_famaug"]
    numeric = frame[score_columns].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
        raise ValueError("Both model score columns must contain finite numeric values")
    target_norm = frame.target.astype(str).str.upper()
    routing_target = np.where(target_norm.isin(M4_ALIASES), "M4R", target_norm)
    counts = pd.Series(routing_target).value_counts()
    if (counts < 2).any():
        bad = counts[counts < 2].index.tolist()
        raise ValueError(f"Percentile routing requires at least two candidates per target: {bad}")


def _percentile_rank(frame: pd.DataFrame, column: str) -> pd.Series:
    group_column = "_routing_target" if "_routing_target" in frame else "target"
    return frame.groupby(group_column, sort=False)[column].rank(method="average", pct=True)


def route_scores(frame: pd.DataFrame, config: RouterConfig | None = None) -> pd.DataFrame:
    """Route two model-score columns into one deployable ranking.

    Parameters
    ----------
    frame:
        One row per target--molecule pair. Required columns are ``pair_id``,
        ``target``, ``score_2023_gpcr_loto`` and
        ``score_2026_13t_famaug``. Additional columns are preserved.
    config:
        Frozen v01 defaults to 0.5/0.5 rank fusion outside M4.

    Returns
    -------
    pandas.DataFrame
        Original columns plus per-model percentile ranks, routed score,
        target-local final rank, route name and claim boundary.
    """
    config = config or RouterConfig()
    config.validate()
    _validate_input(frame)
    out = frame.copy()
    target_norm = out.target.astype(str).str.upper()
    out["_routing_target"] = np.where(target_norm.isin(M4_ALIASES), "M4R", target_norm)
    out["rankpct_2023_gpcr_loto"] = _percentile_rank(out, "score_2023_gpcr_loto")
    out["rankpct_2026_13t_famaug"] = _percentile_rank(out, "score_2026_13t_famaug")
    is_m4 = target_norm.isin(M4_ALIASES)
    fused = (config.old_weight * out.rankpct_2023_gpcr_loto
             + config.new_weight * out.rankpct_2026_13t_famaug)
    out["pacer_binding_score"] = np.where(is_m4, out.rankpct_2023_gpcr_loto, fused)
    out["route"] = np.where(is_m4, "m4_2023_gpcr_loto", "gpcr_fixed_rank_fusion")
    out["final_rank"] = out.groupby("_routing_target", sort=False).pacer_binding_score.rank(
        method="first", ascending=False
    ).astype(int)
    out["router_protocol"] = config.protocol_id
    out["frozen_protocol"] = bool(config.is_frozen_v01)
    out["eligible_for_pam_claim"] = False
    out["claim_boundary"] = "binding-candidate retrieval only; not PAM function or potency"
    out = out.drop(columns="_routing_target")
    return out.sort_values(["target", "final_rank", "pair_id"], kind="mergesort").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="CSV containing pair_id, target and both model scores")
    parser.add_argument("--output", type=Path, required=True, help="Routed CSV")
    parser.add_argument("--audit", type=Path, help="Audit JSON; defaults beside output")
    parser.add_argument("--old-weight", type=float, default=0.5)
    parser.add_argument("--new-weight", type=float, default=0.5)
    args = parser.parse_args()
    config = RouterConfig(old_weight=args.old_weight, new_weight=args.new_weight)
    frame = pd.read_csv(args.input)
    routed = route_scores(frame, config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    routed.to_csv(args.output, index=False)
    audit_path = args.audit or args.output.with_suffix(".audit.json")
    audit = {
        "protocol": asdict(config),
        "frozen_protocol": config.is_frozen_v01,
        "input": str(args.input),
        "input_sha256": sha256(args.input),
        "output": str(args.output),
        "output_sha256": sha256(args.output),
        "rows": int(len(routed)),
        "targets": {str(k): int(v) for k, v in routed.target.value_counts().sort_index().items()},
        "routes": {str(k): int(v) for k, v in routed.route.value_counts().sort_index().items()},
        "claim_boundary": "Binding-candidate retrieval only; PACER-FKG or experiment is required for a PAM-function claim.",
    }
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
