#!/usr/bin/env python3
"""One-command, fail-closed integration benchmark for the PACER project."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run(command: list[str], cwd: Path) -> dict:
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    record = {"command": command, "returncode": completed.returncode,
              "stdout_tail": completed.stdout[-4000:], "stderr_tail": completed.stderr[-4000:]}
    if completed.returncode:
        raise RuntimeError(json.dumps(record, ensure_ascii=False, indent=2))
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--artifact-root", type=Path, required=True,
                        help="Checkout containing git-ignored raw predictions")
    parser.add_argument("--old-ranking", type=Path, required=True)
    parser.add_argument("--famaug-scores", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=20000)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    out = root / "project/results/project_wide_integration_benchmark_v02"
    py = sys.executable
    steps = []
    steps.append(run([py, "project/scripts/recompute_binding_benchmark_v02.py",
                      "--repo-root", str(root), "--artifact-root", str(args.artifact_root.resolve()),
                      "--bootstrap", str(args.bootstrap)], root))
    steps.append(run([py, "project/scripts/audit_candidate_handoff_v02.py",
                      "--old-ranking", str(args.old_ranking.resolve()),
                      "--famaug-scores", str(args.famaug_scores.resolve()),
                      "--output-dir", str(out)], root))
    steps.append(run([py, "project/scripts/build_benchmark_registry_v02.py", "--repo-root", str(root)], root))
    steps.append(run([py, "project/scripts/make_project_benchmark_figures_v02.py", "--repo-root", str(root)], root))
    steps.append(run([py, "project/scripts/audit_branch_integration_v02.py", "--repo-root", str(root)], root))

    required = [out / "binding_recomputed_v02.json", out / "candidate_handoff_audit_v02.json",
                out / "benchmark_method_registry_v02.csv", out / "figures/fig_v02_01_binding_fair_bars.png",
                out / "figures/fig_v02_06_candidate_handoff_scatter.png", out / "branch_integration_audit_v02.json"]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    audit = json.loads((out / "candidate_handoff_audit_v02.json").read_text(encoding="utf-8"))
    if missing or audit["status"] != "PASS" or audit["matched_ids_and_smiles"] != 200:
        raise RuntimeError(f"integration gate failed: missing={missing}, candidate_audit={audit['status']}")
    receipt = {
        "status": "PASS", "version": "v02", "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "commands": steps,
        "required_artifacts": [{"path": str(path.relative_to(root)), "sha256": digest(path)} for path in required],
        "claim_boundary": "Integration success is reproducibility, not PAM confirmation or universal SOTA.",
    }
    receipt_path = out / "integration_run_receipt_v02.json"
    receipt_path.write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "receipt": str(receipt_path)}, indent=2))


if __name__ == "__main__":
    main()
