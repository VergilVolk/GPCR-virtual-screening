# -*- coding: utf-8 -*-
"""One-command reproduction of the PACER-FS evidence package."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


P = Path(__file__).resolve().parents[1]


def run(name):
    print("\n>>", name, flush=True)
    subprocess.run([sys.executable, str(P / "scripts" / name)], cwd=P.parent, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-audit", action="store_true", help="include slow nested negative-control audit")
    args = parser.parse_args()
    for script in [
        "run_pacer_fs_baselines.py", "analyze_pacer_fs_v01.py",
        "run_pacer_series_normalized.py", "run_pacer_centered_model_suite.py",
        "run_pacer_centering_robustness.py", "audit_pacer_centered_domain.py",
        "run_pacer_centered_temporal.py", "run_pacer_temporal_fewshot.py",
        "run_pacer_dual_channel_loso.py", "finalize_pacer_fs.py",
    ]:
        run(script)
    if args.full_audit:
        run("run_pacer_fs_nested.py")
    print("\nPACER-FS evidence package complete: project/results/pacer_fs_final/")


if __name__ == "__main__":
    main()
