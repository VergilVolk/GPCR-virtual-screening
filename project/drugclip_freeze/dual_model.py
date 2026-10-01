"""Frozen dual-model decision rule.

    decision = OLD_TOP50  INTERSECT  NEW_TOP25

Model A (official 2023 weights) is the M4 decision line; Model B (family-
augmented fine-tune) is a second opinion.  Score fusion was explicitly evaluated
and rejected upstream, so no fused score may be produced or consumed here.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable

OLD_TOP = 50
NEW_TOP = 25
FUSION_COLUMNS_FORBIDDEN = ("fixed_rank_fusion", "fusion_rank", "fused_score", "blend", "stacker")


class DualModelError(ValueError):
    """Raised when the frozen dual-model rule would be violated."""


def decide(records: Iterable[dict[str, Any]], *, old_top: int = OLD_TOP, new_top: int = NEW_TOP) -> dict[str, Any]:
    rows = []
    for record in records:
        old_rank = record.get("old_rank")
        new_rank = record.get("new_rank")
        if old_rank is None or new_rank is None:
            raise DualModelError(f"record lacks old_rank/new_rank: {record}")
        rows.append({"candidate_id": str(record["candidate_id"]),
                     "old_rank": int(old_rank), "new_rank": int(new_rank)})
    if not rows:
        raise DualModelError("no candidate records supplied")
    dual = sorted((r for r in rows if r["old_rank"] <= old_top and r["new_rank"] <= new_top),
                  key=lambda r: (r["old_rank"], r["new_rank"]))
    return {
        "schema": "pacer.drugclip_freeze.dual_model_decision.v1",
        "rule": f"old_top{old_top} INTERSECT new_top{new_top}",
        "old_top": old_top,
        "new_top": new_top,
        "n_input": len(rows),
        "n_dual_top": len(dual),
        "dual_top_candidates": dual,
        "fusion_used": False,
        "model_a_role": "M4 decision line (official frozen weights)",
        "model_b_role": "second opinion (frozen projection-only apply)",
        "claim_boundary": ("Dual-model intersection is a prioritization device only; it is not a "
                           "PAM probability, not a fused score and not experimental confirmation."),
    }


def assert_no_fusion(columns: Iterable[str]) -> None:
    offenders = sorted(set(columns) & set(FUSION_COLUMNS_FORBIDDEN))
    if offenders:
        raise DualModelError(f"fused-score columns are forbidden in a frozen dual-model artifact: {offenders}")


def verify_committed_shortlist(path: Path, *, old_top: int = OLD_TOP, new_top: int = NEW_TOP) -> dict[str, Any]:
    """Re-derive the frozen rule against the committed v02 dual-model shortlist."""
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        columns = reader.fieldnames or []
    # The decision consumes ONLY the two frozen ranks.  Descriptive fusion columns may
    # legitimately appear on a benchmark artifact; they must never enter the decision.
    decision_inputs = [c for c in columns if c in ("old_rank", "new_rank", "candidate_id")]
    assert_no_fusion([c for c in columns if c in ("old_rank", "new_rank")])

    checked, violations = [], []
    for row in rows:
        old_rank = int(float(row["drugclip_rank"]))
        new_rank = int(float(row["new_rank"]))
        expected_dual = old_rank <= old_top and new_rank <= new_top
        merge = str(row.get("_merge", "")).strip().lower()
        observed_dual = merge == "both"
        record = {"candidate_id": row["candidate_id"], "old_rank": old_rank, "new_rank": new_rank,
                  "_merge": merge, "expected_dual": expected_dual, "observed_dual": observed_dual}
        checked.append(record)
        if expected_dual != observed_dual:
            violations.append(record)
    return {
        "shortlist_path": str(path),
        "n_rows": len(rows),
        "columns_checked": ["drugclip_rank", "new_rank", "_merge"],
        "old_top": old_top,
        "new_top": new_top,
        "records": checked,
        "violations": violations,
        "rule_reproduced": not violations,
        "note": ("The committed shortlist is the frozen output of the dual-model rule; the "
                 "cross-model fused columns it also carries are descriptive only and are NOT used "
                 "by the decision."),
    }
