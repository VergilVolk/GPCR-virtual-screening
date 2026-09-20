"""Reproduce the frozen external and post-hoc multidomain PACER evidence."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "project" / "scripts"
STEPS = [
    "build_suven_patent_external.py",
    "audit_suven_patent_external.py",
    "run_pacer_suven_external.py",
    "analyze_suven_factor_sar.py",
    "evaluate_factor_sar_nested.py",
    "evaluate_multidomain_sar.py",
    "run_multidomain_vu_sequential.py",
    "analyze_generation4_stereo.py",
    "build_monash_ly2033298_allostery.py",
    "run_pacer_monash_allostery.py",
    "build_m4_mechanistic_multidomain.py",
    "run_pacer_mechanistic_multidomain.py",
    "run_pacer_mechanistic_fewshot.py",
    "run_pacer_assay_anchor_robustness.py",
    "build_acadia_2025_external.py",
    "run_pacer_acadia_external.py",
    "run_pacer_acadia_anchor_robustness.py",
    "analyze_acadia_anchor_failure.py",
    "analyze_acadia_pam_agonism_tradeoff.py",
    "evaluate_pacer_adaptive_anchors.py",
    "run_pacer_us20260055116_external.py",
    "analyze_us20260055116_functional_divergence.py",
    "evaluate_pacer_context_kernel.py",
    "build_candidate_function_space_panel.py",
    "plot_pacer_function_space_v3.py",
]


def main() -> None:
    for step in STEPS:
        print(f"\n=== {step} ===", flush=True)
        subprocess.run([sys.executable, str(SCRIPTS / step)], cwd=ROOT, check=True)
    print("\nPACER external frontier reproduction complete.")


if __name__ == "__main__":
    main()
