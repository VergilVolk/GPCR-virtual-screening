"""Diagnose why the frozen chemistry-diverse Acadia anchors failed."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "results" / "pacer_external_acadia_anchor_robustness_v01" / "repeats.csv"
OUT = PROJECT / "results" / "pacer_external_acadia_anchor_failure_v01"
METHOD = "PACER-FS-DeltaSARHybrid-Centered"
FEATURES = [
    "anchor_pEC50_span", "anchor_base_prediction_span",
    "anchor_knn_prediction_span", "anchor_residual_sd",
    "anchor_pair_tanimoto_mean", "anchor_pair_tanimoto_max",
    "mean_query_max_anchor_tanimoto",
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(SOURCE)
    random = d[(d.policy == "random") & (d.method == METHOD)].copy()
    deterministic = d[(d.policy == "deterministic_diverse") & (d.method == METHOD)].iloc[0]
    rows = []
    for feature in FEATURES:
        rho_rank = float(spearmanr(random[feature], random.Spearman).statistic)
        rho_mae = float(spearmanr(random[feature], random.MAE).statistic)
        value = float(deterministic[feature])
        percentile = float((random[feature] <= value).mean())
        rows.append({
            "feature": feature, "spearman_with_query_rank": rho_rank,
            "spearman_with_query_MAE": rho_mae,
            "deterministic_anchor_value": value,
            "deterministic_percentile_among_random": percentile,
        })
    diagnostics = pd.DataFrame(rows).sort_values(
        "spearman_with_query_rank", key=lambda x: x.abs(), ascending=False)
    diagnostics.to_csv(OUT / "feature_diagnostics.csv", index=False)
    payload = {
        "status": "post_hoc_failure_diagnosis",
        "method": METHOD, "n_random_anchor_triples": len(random),
        "deterministic_anchor_ids": deterministic.anchor_ids,
        "deterministic_anchor_pEC50_span": float(deterministic.anchor_pEC50_span),
        "random_median_anchor_pEC50_span": float(random.anchor_pEC50_span.median()),
        "span_vs_rank_spearman": float(spearmanr(
            random.anchor_pEC50_span, random.Spearman).statistic),
        "interpretation": "Chemistry diversity did not guarantee functional-response bracketing; this motivates sequential functional anchoring.",
        "claim_boundary": "Post-hoc diagnosis on Acadia; cannot validate a new anchor policy.",
    }
    (OUT / "audit.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lines = ["# Why the frozen Acadia anchors failed", "",
             f"The deterministic chemistry-diverse anchors spanned only {payload['deterministic_anchor_pEC50_span']:.3f} log units, versus a random-triplet median of {payload['random_median_anchor_pEC50_span']:.3f}.",
             f"Across 1000 random triples, response span correlated with query ranking at rho={payload['span_vs_rank_spearman']:.3f}.", "",
             "Chemical coverage was already saturated, so adding more structural diversity did not solve functional bracketing. The next method must choose later anchors sequentially after observing earlier functional measurements.", "",
             diagnostics.to_markdown(index=False, floatfmt=".3f"), "",
             "This is a post-hoc mechanism diagnosis, not independent support for a replacement policy."]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print(diagnostics.to_string(index=False))


if __name__ == "__main__":
    main()
