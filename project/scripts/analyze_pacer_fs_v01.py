# -*- coding: utf-8 -*-
"""Paired uncertainty analysis for PACER-FS v0.1.

All comparisons use exactly the same query compounds within an episode.  The
bootstrap resamples medicinal-chemistry series, not individual molecules, so
large series cannot dominate the confidence interval.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


P = Path(__file__).resolve().parents[1]
IN = P / "results" / "pacer_fs_v01" / "episodes.csv"
OUT = P / "results" / "pacer_fs_v01"
N_BOOT = 50_000
SEED = 20260829


def paired_table(df: pd.DataFrame, shot: int, policy: str, challenger: str) -> pd.DataFrame:
    sub = df[(df.shot == shot) & (df.policy == policy) & df.method.isin(["Base", challenger])]
    keys = ["group", "repeat"]
    wide = sub.pivot(index=keys, columns="method", values=["Spearman", "MAE"]).dropna().reset_index()
    wide.columns = ["_".join(x).rstrip("_") if isinstance(x, tuple) else x for x in wide.columns]
    return wide


def aggregate_repeats(wide: pd.DataFrame, challenger: str) -> pd.DataFrame:
    cols = [f"Spearman_{challenger}", "Spearman_Base", f"MAE_{challenger}", "MAE_Base"]
    return wide.groupby("group", as_index=False)[cols].mean()


def bootstrap_delta(per_series: pd.DataFrame, challenger: str, rng: np.random.Generator) -> dict:
    dr = (per_series[f"Spearman_{challenger}"] - per_series["Spearman_Base"]).to_numpy(float)
    dm = (per_series[f"MAE_{challenger}"] - per_series["MAE_Base"]).to_numpy(float)
    draw = rng.integers(0, len(per_series), size=(N_BOOT, len(per_series)))
    br = dr[draw].mean(axis=1)
    bm = dm[draw].mean(axis=1)
    return {
        "n_series": int(len(per_series)),
        "delta_macro_spearman": float(dr.mean()),
        "delta_macro_spearman_ci95": [float(x) for x in np.quantile(br, [0.025, 0.975])],
        "probability_delta_spearman_gt_0": float((br > 0).mean()),
        "delta_macro_mae": float(dm.mean()),
        "delta_macro_mae_ci95": [float(x) for x in np.quantile(bm, [0.025, 0.975])],
        "probability_delta_mae_lt_0": float((bm < 0).mean()),
    }


def main() -> None:
    df = pd.read_csv(IN)
    rng = np.random.default_rng(SEED)
    results = {}
    all_rows = []

    for policy in ["random", "predicted_span"]:
        for shot in [1, 3, 5]:
            for challenger in ["KernelResidual", "TanimotoGP", "DeltaSARHybrid"]:
                wide = paired_table(df, shot, policy, challenger)
                per = aggregate_repeats(wide, challenger)
                key = f"{policy}_{shot}shot_{challenger}"
                results[key] = bootstrap_delta(per, challenger, rng)
                per = per.rename(columns={
                    f"Spearman_{challenger}": "Spearman_method",
                    f"MAE_{challenger}": "MAE_method",
                })
                per.insert(1, "policy", policy)
                per.insert(2, "shot", shot)
                per.insert(3, "method", challenger)
                per["delta_Spearman"] = per.Spearman_method - per.Spearman_Base
                per["delta_MAE"] = per.MAE_method - per.MAE_Base
                all_rows.append(per)

    per_series = pd.concat(all_rows, ignore_index=True)
    per_series.to_csv(OUT / "paired_per_series.csv", index=False)

    # Development-set comparison of deterministic predicted-span anchors with
    # the mean over 20 random anchor selections.  This is diagnostic, not an
    # independent test because the policy was discovered on these series.
    method = "DeltaSARHybrid"
    for shot in [1, 3, 5]:
        a = aggregate_repeats(paired_table(df, shot, "predicted_span", method), method)
        b = aggregate_repeats(paired_table(df, shot, "random", method), method)
        z = a.merge(b, on="group", suffixes=("_span", "_random"))
        dr = (z[f"Spearman_{method}_span"] - z[f"Spearman_{method}_random"]).to_numpy()
        dm = (z[f"MAE_{method}_span"] - z[f"MAE_{method}_random"]).to_numpy()
        draw = rng.integers(0, len(z), size=(N_BOOT, len(z)))
        results[f"predicted_span_vs_random_{shot}shot_{method}"] = {
            "status": "development_diagnostic_not_blind_test",
            "n_series": int(len(z)),
            "delta_macro_spearman": float(dr.mean()),
            "delta_macro_spearman_ci95": [float(x) for x in np.quantile(dr[draw].mean(1), [0.025, 0.975])],
            "delta_macro_mae": float(dm.mean()),
            "delta_macro_mae_ci95": [float(x) for x in np.quantile(dm[draw].mean(1), [0.025, 0.975])],
        }

    (OUT / "paired_bootstrap.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    focus = results["predicted_span_3shot_DeltaSARHybrid"]
    print(json.dumps(focus, indent=2))
    print("\nPer-series 3-shot predicted-span DeltaSARHybrid:")
    print(per_series[(per_series.policy == "predicted_span") & (per_series.shot == 3) &
                     (per_series.method == "DeltaSARHybrid")].to_string(index=False))


if __name__ == "__main__":
    main()
