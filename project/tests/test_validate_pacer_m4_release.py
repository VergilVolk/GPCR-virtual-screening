from __future__ import annotations

import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_pacer_m4_release.py"
SPEC = importlib.util.spec_from_file_location("release_validator", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


def _manifest(path: Path) -> None:
    path.write_text(
        json.dumps(
            {
                "release_id": "test",
                "required": [{"path": "required.txt", "role": "required"}],
                "optional": [{"path": "optional.bin", "role": "optional"}],
                "claim_boundary": "test only",
            }
        ),
        encoding="utf-8",
    )


def test_missing_optional_does_not_invalidate_source_release(tmp_path: Path) -> None:
    (tmp_path / "required.txt").write_text("ok", encoding="utf-8")
    _manifest(tmp_path / "manifest.json")
    result = MOD.validate_release(tmp_path, Path("manifest.json"))
    assert result["valid"] is True
    assert result["missing_optional"] == ["optional.bin"]


def test_missing_required_is_rejected(tmp_path: Path) -> None:
    _manifest(tmp_path / "manifest.json")
    result = MOD.validate_release(tmp_path, Path("manifest.json"))
    assert result["valid"] is False
    assert result["missing_required"] == ["required.txt"]


def test_strict_optional_requires_optional_assets(tmp_path: Path) -> None:
    (tmp_path / "required.txt").write_text("ok", encoding="utf-8")
    _manifest(tmp_path / "manifest.json")
    result = MOD.validate_release(tmp_path, Path("manifest.json"), strict_optional=True)
    assert result["valid"] is False
