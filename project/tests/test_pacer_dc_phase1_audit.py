from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_pacer_dc_phase1 import audit_manifest  # noqa: E402


def _manifest() -> pd.DataFrame:
    rows = []
    for context in ("probe_only", "apo"):
        for replicate in ("r1", "r2", "r3"):
            rows.append([f"shared-{context}-{replicate}", "SHARED_CONTROL", context, replicate, "shared_context_control", replicate[-1], "not_started"])
    for context in ("candidate_probe", "candidate_no_probe"):
        for replicate in ("r1", "r2", "r3"):
            rows.append([f"pam-{context}-{replicate}", "PAM1", context, replicate, "known_PAM_positive_control", replicate[-1], "not_started"])
    return pd.DataFrame(rows, columns=["run_id", "candidate_id", "context", "replicate_id", "role", "paired_seed_group", "status"])


def test_valid_balanced_manifest() -> None:
    result = audit_manifest(_manifest())
    assert result["valid"] is True
    assert result["candidate_count_excluding_shared_control"] == 1


def test_missing_context_is_rejected() -> None:
    manifest = _manifest()
    manifest = manifest[~((manifest.candidate_id == "PAM1") & (manifest.context == "candidate_no_probe"))]
    result = audit_manifest(manifest)
    assert result["valid"] is False
    assert any("PAM1" in message for message in result["errors"])
