#!/usr/bin/env python3
"""Pre-specified 50:50 target-rank fusion of ECFP and PACER structural scores."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from finetune_drugclip_gpcr_screening import binary_metrics


def evaluate(frame: pd.DataFrame, column: str, targets: list[str]) -> dict:
    per = {}
    for target in targets:
        keep = frame.target.eq(target).to_numpy()
        per[target] = binary_metrics(frame.loc[keep, "label"].to_numpy(int), frame.loc[keep, column].to_numpy())
    names = ["roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"]
    return {"macro": {name: float(np.mean([per[t][name] for t in targets])) for name in names},
            "per_target": per}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ecfp", type=Path, required=True)
    parser.add_argument("--pacer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10000)
    args = parser.parse_args()
    ecfp = pd.read_csv(args.ecfp)
    pacer = pd.read_csv(args.pacer)
    columns = ["pair_id", "target", "label"]
    frame = ecfp[columns + ["ecfp_loto"]].merge(
        pacer[columns + ["science2026_raw", "pacer_cgm"]], on=columns, validate="one_to_one"
    )
    for source in ["ecfp_loto", "science2026_raw", "pacer_cgm"]:
        frame[source + "_rank"] = frame.groupby("target")[source].rank(method="average", pct=True)
    frame["fixed_fusion"] = 0.5 * frame.ecfp_loto_rank + 0.5 * frame.pacer_cgm_rank
    targets = sorted(frame.target.unique())
    reports = {name: evaluate(frame, name, targets) for name in
               ["science2026_raw", "pacer_cgm", "ecfp_loto", "fixed_fusion"]}
    metrics = list(reports["fixed_fusion"]["macro"])
    rng = np.random.default_rng(20260927)
    draws = rng.integers(0, len(targets), size=(args.bootstrap, len(targets)))
    ci = {}
    for baseline in ["science2026_raw", "pacer_cgm", "ecfp_loto"]:
        ci[baseline] = {}
        for metric in metrics:
            delta = np.asarray([
                reports["fixed_fusion"]["per_target"][target][metric]
                - reports[baseline]["per_target"][target][metric] for target in targets
            ])
            ci[baseline][metric] = list(map(float, np.quantile(delta[draws].mean(1), [0.025, 0.5, 0.975])))
    report = {
        "method": "Pre-specified equal-weight within-target rank fusion",
        "weight_search": False,
        "metrics": reports,
        "delta_fixed_fusion_vs_baseline_target_bootstrap_95ci": ci,
        "claim_boundary": "Development-set complementarity test; fixed weight, not nested-selected and not SOTA evidence.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    frame.to_csv(args.output.with_suffix(".predictions.csv"), index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
