# -*- coding: utf-8 -*-
"""Reproduce the PACER-M4 v1 evidence and candidate pipeline.

Expensive docking is resumable.  Existing frozen outputs are reused unless
--force-docking is supplied.
"""
from __future__ import annotations
import argparse,subprocess,sys
from pathlib import Path
P=Path(__file__).resolve().parents[1]
def run(script,*args):
 cmd=[sys.executable,str(P/"scripts"/script),*map(str,args)];print("\n>>",script,flush=True);subprocess.run(cmd,cwd=P.parent,check=True)
def main():
 a=argparse.ArgumentParser();a.add_argument("--force-docking",action="store_true");a.add_argument("--workers",type=int,default=8);x=a.parse_args()
 # Audits and strict retrospective tests.
 run("run_lightgbm_loso.py");run("run_lightgbm_v2_loso.py");run("evaluate_local_domain_algorithm.py");run("run_temporal_benchmark.py");run("reproduce_m4r_enrichment.py")
 # Generation audit and bounded candidate pool.
 run("compare_generators.py");run("build_candidate_portfolio.py")
 features=P/"results"/"pacer_candidates_v01"/"7trs_docking"/"features.csv"
 if x.force_docking or not features.exists():run("dock_candidate_portfolio.py","--workers",x.workers,"--exhaustiveness",4)
 else:print("\n>> docking cache reused:",features)
 run("select_pareto_candidates.py");run("calibrate_candidate_intervals.py");run("export_candidate_sdf.py");run("make_pacer_final_figures.py");run("validate_pacer_release.py")
 # Freeze the prospective assay batch and relative-only pre-assay predictions.
 run("build_pacer_assay_closed_loop.py");run("apply_pacer_fs_candidates.py")
 print("\nPACER-M4 v1 complete. Read project/PACER_M4_FINAL_REPORT.md and results/pacer_candidates_v01/final/.")
if __name__=="__main__":main()
