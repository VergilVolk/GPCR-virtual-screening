#!/usr/bin/env python3
"""Validate the reviewer-facing PACER-M4 repository contract.

The checker is intentionally dependency-free.  It verifies that every required
entry point in the release manifest exists, records file hashes, reports
optional assets separately, and warns about large tracked files.  Missing
optional trajectories or checkpoints do not make the source release invalid.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


DEFAULT_MANIFEST = Path("project/config/pacer_m4_release_manifest_v01.json")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _entry_status(root: Path, entry: dict[str, str]) -> dict[str, Any]:
    path = root / entry["path"]
    result: dict[str, Any] = {**entry, "exists": path.exists()}
    if path.is_file():
        result.update({"kind": "file", "bytes": path.stat().st_size, "sha256": sha256(path)})
    elif path.is_dir():
        files = [item for item in path.rglob("*") if item.is_file()]
        result.update({"kind": "directory", "files": len(files)})
    else:
        result["kind"] = "missing"
    return result


def _tracked_large_files(root: Path, threshold_mb: float) -> list[dict[str, Any]]:
    try:
        completed = subprocess.run(
            ["git", "ls-files", "-z"], cwd=root, check=True, capture_output=True
        )
    except (OSError, subprocess.CalledProcessError):
        return []
    threshold = int(threshold_mb * 1024 * 1024)
    large: list[dict[str, Any]] = []
    for raw in completed.stdout.split(b"\0"):
        if not raw:
            continue
        relative = raw.decode("utf-8", errors="surrogateescape")
        path = root / relative
        if path.is_file() and path.stat().st_size >= threshold:
            large.append({"path": relative, "bytes": path.stat().st_size})
    return sorted(large, key=lambda item: item["bytes"], reverse=True)


def validate_release(
    root: Path,
    manifest_path: Path,
    *,
    strict_optional: bool = False,
    large_file_mb: float = 25.0,
) -> dict[str, Any]:
    manifest = json.loads((root / manifest_path).read_text(encoding="utf-8"))
    required = [_entry_status(root, entry) for entry in manifest["required"]]
    optional = [_entry_status(root, entry) for entry in manifest.get("optional", [])]
    missing_required = [entry["path"] for entry in required if not entry["exists"]]
    missing_optional = [entry["path"] for entry in optional if not entry["exists"]]
    valid = not missing_required and (not strict_optional or not missing_optional)
    return {
        "release_id": manifest["release_id"],
        "manifest": str(manifest_path),
        "valid": valid,
        "strict_optional": strict_optional,
        "required": required,
        "optional": optional,
        "missing_required": missing_required,
        "missing_optional": missing_optional,
        "tracked_files_at_least_mb": large_file_mb,
        "tracked_large_files": _tracked_large_files(root, large_file_mb),
        "claim_boundary": manifest["claim_boundary"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--strict-optional", action="store_true")
    parser.add_argument("--large-file-mb", type=float, default=25.0)
    args = parser.parse_args()

    report = validate_release(
        args.root.resolve(),
        args.manifest,
        strict_optional=args.strict_optional,
        large_file_mb=args.large_file_mb,
    )
    text = json.dumps(report, indent=2, ensure_ascii=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    raise SystemExit(0 if report["valid"] else 1)


if __name__ == "__main__":
    main()
