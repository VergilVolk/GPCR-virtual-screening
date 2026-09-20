# -*- coding: utf-8 -*-
"""Reproducible end-to-end M4R GaMD ensemble validation and candidate evidence run."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

P = Path(__file__).resolve().parents[1]
PY = sys.executable


def run(script: str, *args: str) -> None:
    cmd = [PY, str(P / "scripts" / script), *args]
    print("\nRUN", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=P.parent, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--skip-preparation", action="store_true")
    ap.add_argument("--skip-validation-docking", action="store_true")
    ap.add_argument("--skip-candidate-docking", action="store_true")
    args = ap.parse_args()
    workers = str(args.workers)

    if not args.skip_preparation:
        run("prepare_m4_gamd_ensemble.py")
        run("build_m4_gamd_validation_set.py")
    if not args.skip_validation_docking:
        run("dock_m4_gamd_ensemble.py", "--set", "validation", "--clusters", "0-9",
            "--workers", workers)
    run("analyze_m4_gamd_ensemble.py")
    if not args.skip_candidate_docking:
        run("dock_m4_gamd_ensemble.py", "--set", "candidates", "--clusters", "0-9",
            "--workers", workers)
    run("rank_m4_gamd_candidates.py")


if __name__ == "__main__":
    main()
