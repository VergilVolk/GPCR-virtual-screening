#!/usr/bin/env python3
"""Aggregate 3-seed CGDA-on-Science-2026 LOSO runs: mean macro metrics and
paired target-level bootstrap CIs for cgda - reference and cgda - pocket."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--runs", type=Path, nargs="+", required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=20000)
    p.add_argument("--seed", type=int, default=20260926)
    a = p.parse_args()
    runs = [json.load(open(r, encoding="utf-8")) for r in a.runs]
    targets = runs[0]["targets"] if "targets" in runs[0] else sorted(runs[0]["folds"])
    assert all((r["targets"] if "targets" in r else sorted(r["folds"])) == targets for r in runs)
    names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]

    def per_run_macro(r):
        return {m: {n: float(np.mean([r["folds"][t][m][n] for t in targets])) for n in names}
                for m in ("pocket", "reference", "cgda")}

    macro = {m: {n: float(np.mean([per_run_macro(r)[m][n] for r in runs])) for n in names}
             for m in ("pocket", "reference", "cgda")}
    rng = np.random.default_rng(a.seed)
    deltas = {}
    for cmp in ("reference", "pocket"):
        # per (run, target, metric) matrix
        mat = np.stack([[[r["folds"][t]["cgda"][n] - r["folds"][t][cmp][n] for n in names]
                         for t in targets] for r in runs])  # runs x targets x metrics
        mean_per_run = mat.mean(axis=1)                      # runs x metrics
        boot = np.stack([mean_per_run[rng.integers(len(runs), size=len(runs))].mean(axis=0)
                         for _ in range(a.bootstrap)])       # bootstrap x metrics
        deltas[cmp] = {n: [float(np.quantile(boot[:, i], 0.025)),
                           float(mean_per_run[:, i].mean()),
                           float(np.quantile(boot[:, i], 0.975))]
                       for i, n in enumerate(names)}
    report = {"protocol": "mean of 3 independent seeds; paired target bootstrap over seeds",
              "n_targets": len(targets), "n_runs": len(runs),
              "macro": macro, "delta_95ci": deltas,
              "claim_boundary": "Retrospective LOSO; reference-assisted binding screen only."}
    a.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"macro": macro, "delta_95ci": deltas}, indent=2))

if __name__ == "__main__":
    main()
