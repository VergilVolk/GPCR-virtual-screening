# -*- coding: utf-8 -*-
"""One-command reproduction of Thompson–Miao M4 metrics and PACER-XR."""
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--recluster", action="store_true")
    args = ap.parse_args()
    run("reproduce_thompson_miao_m4_metrics.py")
    run("run_cross_target_rank_fusion.py")
    run("audit_pacer_xr_functional_semantics.py")
    audit = P / "results" / "m4_gamd_public_recluster_v01" / "audit.json"
    if args.recluster or not audit.exists():
        run("download_m4_gamd_figshare.py")
        run("recluster_m4_gamd_public_stride.py")
        run("compare_m4_cluster_representatives.py")
    else:
        print("\n>> verified public-stride clustering cache reused:", audit)
    print("\nReproduction complete. Read project/docs/THOMPSON_MIAO_2026_EXACT_REPRODUCTION.md")


if __name__ == "__main__":
    main()
