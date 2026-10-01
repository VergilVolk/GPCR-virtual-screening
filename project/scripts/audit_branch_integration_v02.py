#!/usr/bin/env python3
"""Verify that every frozen project line is an ancestor of the integration HEAD."""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


EXPECTED = {
    "oneprot_g0": "e04949f6",
    "encoder_optimization": "0a4b4f23",
    "fkg_long_md": "07ce8ae1",
    "geom2vec_pilot": "14d86967",
    "close_loop_20ns": "6d12bcef",
    "cm00734_analysis": "853aaa12",
    "cm00734_stage_b_synthesis": "5e86149b",
    "drugclip_handoff": "c19b7443",
    "project_benchmark_v01": "a2e01fe1",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = {}
    for name, commit in EXPECTED.items():
        present = subprocess.run(["git", "cat-file", "-e", f"{commit}^{{commit}}"], cwd=args.repo_root).returncode == 0
        included = present and subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"], cwd=args.repo_root).returncode == 0
        rows[name] = {"commit": commit, "object_present": present, "ancestor_of_head": included}
    status = "PASS" if all(row["ancestor_of_head"] for row in rows.values()) else "FAIL"
    result = {"status": status, "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=args.repo_root, text=True).strip(),
              "frozen_lines": rows}
    output = args.output or args.repo_root / "project/results/project_wide_integration_benchmark_v02/branch_integration_audit_v02.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
