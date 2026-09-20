# -*- coding: utf-8 -*-
"""One-command reproduction of PACER functional-frontier and biology analyses."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

P = Path(__file__).resolve().parents[1]


def run(name):
    print("\n>>", name, flush=True)
    subprocess.run([sys.executable, str(P / "scripts" / name)], cwd=P.parent, check=True)


def main():
    run("run_pacer_delta_residual.py")
    run("run_pacer_meta_metric.py")
    run("analyze_residue_potency_associations.py")
    run("dock_candidates_7trq_mechanism.py")
    run("make_pacer_biology_figure.py")
    run("validate_pacer_release.py")
    print("\nFunctional-frontier reproduction complete.")


if __name__ == "__main__":
    main()
