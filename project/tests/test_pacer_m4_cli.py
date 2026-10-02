from __future__ import annotations

from pathlib import Path

import pytest

from pacer_m4.runner import run_stage, stage_capabilities
from pacer_m4.stages import STAGES


ROOT = Path(__file__).resolve().parents[2]


def test_registry_covers_all_scientific_stage_families() -> None:
    categories = {stage.category for stage in STAGES.values()}
    assert {"binding", "docking", "potency", "selection", "four_context_md", "dynamic_analysis"} <= categories


def test_registered_scripts_exist() -> None:
    report = stage_capabilities(ROOT)
    assert all(stage["script_exists"] for stage in report["stages"])


def test_dry_run_records_claim_boundary_without_execution(tmp_path: Path) -> None:
    receipt = run_stage(
        "drugclip-route",
        ["--input", "missing.csv", "--output", "unused.csv"],
        root=ROOT,
        dry_run=True,
        receipt_path=tmp_path / "receipt.json",
    )
    assert receipt["status"] == "dry_run"
    assert "not PAM" in str(receipt["claim_boundary"])
    assert (tmp_path / "receipt.json").is_file()


def test_unknown_stage_is_rejected() -> None:
    with pytest.raises(KeyError, match="Unknown stage"):
        run_stage("does-not-exist", root=ROOT, dry_run=True)
