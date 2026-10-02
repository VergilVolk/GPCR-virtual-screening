"""Safe subprocess runner for registered PACER-M4 stages."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from .stages import STAGES


def repository_root(explicit: Path | None = None) -> Path:
    if explicit is not None:
        return explicit.resolve()
    configured = os.environ.get("PACER_M4_ROOT")
    if configured:
        return Path(configured).resolve()
    return Path(__file__).resolve().parents[1]


def git_commit(root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def stage_capabilities(root: Path) -> dict[str, object]:
    return {
        "repository_root": str(root),
        "git_commit": git_commit(root),
        "stages": [STAGES[key].to_dict(root) for key in sorted(STAGES)],
    }


def run_stage(
    stage_id: str,
    stage_args: Sequence[str] = (),
    *,
    root: Path | None = None,
    dry_run: bool = False,
    receipt_path: Path | None = None,
    log_path: Path | None = None,
) -> dict[str, object]:
    if stage_id not in STAGES:
        raise KeyError(f"Unknown stage '{stage_id}'. Available: {', '.join(sorted(STAGES))}")
    root = repository_root(root)
    spec = STAGES[stage_id]
    script = root / spec.script
    if not script.is_file():
        raise FileNotFoundError(f"Registered stage script is missing: {script}")
    arguments = list(stage_args)
    if arguments and arguments[0] == "--":
        arguments = arguments[1:]
    command = [sys.executable, str(script), *arguments]
    started = datetime.now(timezone.utc)
    receipt: dict[str, object] = {
        "stage_id": stage_id,
        "category": spec.category,
        "command": command,
        "repository_root": str(root),
        "git_commit": git_commit(root),
        "started_utc": started.isoformat(),
        "dry_run": dry_run,
        "claim_boundary": spec.claim_boundary,
    }
    if dry_run:
        receipt.update({"status": "dry_run", "returncode": None})
    else:
        destination = None
        handle = None
        if log_path:
            destination = log_path if log_path.is_absolute() else root / log_path
            destination.parent.mkdir(parents=True, exist_ok=True)
            handle = destination.open("w", encoding="utf-8")
        try:
            completed = subprocess.run(
                command,
                cwd=root,
                check=False,
                stdout=handle,
                stderr=subprocess.STDOUT if handle else None,
                text=True if handle else None,
            )
        finally:
            if handle:
                handle.close()
        receipt.update(
            {
                "status": "completed" if completed.returncode == 0 else "failed",
                "returncode": completed.returncode,
                "finished_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
        if destination:
            receipt["log"] = str(destination)
    if receipt_path:
        destination = receipt_path if receipt_path.is_absolute() else root / receipt_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    return receipt
