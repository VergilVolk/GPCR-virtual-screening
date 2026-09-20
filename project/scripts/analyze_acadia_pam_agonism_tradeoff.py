"""Functional SAR trade-off analysis for the frozen Acadia M4 campaign."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "external_acadia_2025" / "external_molecule_benchmark.csv"
OUT = PROJECT / "results" / "acadia_pam_agonism_tradeoff_v01"
SEED = 20260830
N_BOOT = 20_000

MATCHED = [(12, 30), (5, 32), (6, 31), (9, 46), (10, 47),
           (11, 39), (16, 43), (17, 45), (37, 38), (27, 49)]


def core_class(smiles: str) -> str:
    # Patent series differs chiefly at the fused core 5-position.  This explicit
    # assignment was checked against Table 1 names and avoids label-based grouping.
    return "2,5-dimethyl" if "2,5-dimethyl" in smiles else "unknown"


def correlation_ci(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    estimate = float(spearmanr(x, y).statistic)
    rng = np.random.default_rng(SEED)
    values = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(x), len(x))
        value = spearmanr(x[idx], y[idx]).statistic
        if np.isfinite(value):
            values.append(value)
    return {"n": len(x), "rho": estimate,
            "ci95": [float(v) for v in np.quantile(values, [.025, .975])],
            "bootstrap_replicates_valid": len(values)}


def pareto_front(frame: pd.DataFrame) -> pd.Series:
    # Larger is better for potency and PAM RE; smaller is better for agonist RE.
    values = frame[["pEC50", "hM4_PAM_RE_pct", "hM4_agonist_RE_pct"]].to_numpy(float)
    front = np.ones(len(frame), dtype=bool)
    for i, row in enumerate(values):
        dominated = ((values[:, 0] >= row[0]) & (values[:, 1] >= row[1]) &
                     (values[:, 2] <= row[2]) &
                     ((values[:, 0] > row[0]) | (values[:, 1] > row[1]) |
                      (values[:, 2] < row[2]))).any()
        front[i] = not dominated
    return pd.Series(front, index=frame.index)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(DATA)
    exact = d[d.eligible_exact_potency.astype(bool)].copy()

    # Core labels are taken from patent systematic names, not inferred from outcomes.
    exact["core_group"] = np.select(
        [exact.patent_name.str.contains("2,5-dimethyl", case=False, na=False),
         exact.patent_name.str.contains("2-ethyl", case=False, na=False)],
        ["2,5-dimethyl", "2-ethyl"], default="2-methyl")
    core = exact.groupby("core_group").agg(
        n=("compound_id", "size"), median_pEC50=("pEC50", "median"),
        median_PAM_RE=("hM4_PAM_RE_pct", "median"),
        median_agonist_RE=("hM4_agonist_RE_pct", "median"),
        fraction_agonist_RE_ge_50=("intrinsic_agonist_re_ge_50", "mean"),
    ).reset_index().sort_values("median_pEC50", ascending=False)
    core.to_csv(OUT / "core_groups.csv", index=False)

    example_to_row = {}
    for _, row in exact.iterrows():
        for example in str(row.source_example_ids).split("|"):
            example_to_row[int(example)] = row
    pair_rows = []
    for low, high in MATCHED:
        if low not in example_to_row or high not in example_to_row:
            continue
        a, b = example_to_row[low], example_to_row[high]
        pair_rows.append({
            "lower_example": low, "higher_example": high,
            "lower_id": a.compound_id, "higher_id": b.compound_id,
            "delta_pEC50": b.pEC50 - a.pEC50,
            "fold_potency": 10 ** (b.pEC50 - a.pEC50),
            "delta_PAM_RE": b.hM4_PAM_RE_pct - a.hM4_PAM_RE_pct,
            "delta_intrinsic_agonist_RE": b.hM4_agonist_RE_pct - a.hM4_agonist_RE_pct,
        })
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(OUT / "matched_pairs.csv", index=False)

    exact["pareto_potency_pamre_lowagonism"] = pareto_front(exact)
    exact["M4_over_M2_EC50_ratio"] = exact.hM2_PAM_EC50_nM / exact.hM4_PAM_EC50_nM
    pareto = exact[exact.pareto_potency_pamre_lowagonism].copy()
    keep = ["compound_id", "example_id", "pEC50", "hM4_PAM_EC50_nM",
            "hM4_PAM_RE_pct", "hM4_agonist_RE_pct", "hM2_PAM_EC50_nM",
            "M4_over_M2_EC50_ratio", "core_group", "canonical_smiles"]
    pareto[keep].sort_values(["hM4_agonist_RE_pct", "pEC50"], ascending=[True, False]).to_csv(
        OUT / "pareto_candidates.csv", index=False)
    balanced = exact[(exact.hM4_PAM_EC50_nM <= 60) &
                     (exact.hM4_PAM_RE_pct >= 80) &
                     (exact.hM4_agonist_RE_pct < 40)].copy()
    balanced[keep].sort_values(["hM4_agonist_RE_pct", "pEC50"], ascending=[True, False]).to_csv(
        OUT / "balanced_reference_controls.csv", index=False)

    correlations = {
        "pEC50_vs_intrinsic_agonist_RE": correlation_ci(exact.pEC50, exact.hM4_agonist_RE_pct),
        "PAM_RE_vs_intrinsic_agonist_RE": correlation_ci(exact.hM4_PAM_RE_pct, exact.hM4_agonist_RE_pct),
    }
    summary = {
        "dataset": "M4_ACADIA_WO2025122811_FUNCTIONAL_HIGHCONF_V01",
        "n_exact_potency": len(exact), "n_matched_pairs": len(pairs),
        "n_pareto": len(pareto), "n_balanced_reference_controls": len(balanced),
        "balanced_control_rule_post_hoc": "M4 PAM EC50 <=60 nM, PAM RE >=80%, intrinsic agonist RE <40%",
        "correlations": correlations,
        "matched_pair_fraction_potency_and_agonism_increase": float(
            ((pairs.delta_pEC50 > 0) & (pairs.delta_intrinsic_agonist_RE > 0)).mean()),
        "matched_pair_potency_increase_count": int((pairs.delta_pEC50 > 0).sum()),
        "matched_pair_potency_sign_test_two_sided_p": float(
            binomtest(int((pairs.delta_pEC50 > 0).sum()), len(pairs), .5).pvalue),
        "matched_pair_median_fold_potency": float(pairs.fold_potency.median()),
        "matched_pair_agonism_increase_count": int((pairs.delta_intrinsic_agonist_RE > 0).sum()),
        "matched_pair_agonism_sign_test_two_sided_p": float(
            binomtest(int((pairs.delta_intrinsic_agonist_RE > 0).sum()), len(pairs), .5).pvalue),
        "matched_pair_median_intrinsic_agonist_RE_delta": float(
            pairs.delta_intrinsic_agonist_RE.median()),
        "claim_boundary": "Within-patent SAR association; not universal causality and not prospective validation.",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    c1 = correlations["pEC50_vs_intrinsic_agonist_RE"]
    c2 = correlations["PAM_RE_vs_intrinsic_agonist_RE"]
    lines = ["# Acadia M4 PAM potency–agonism trade-off", "",
             f"Among {len(exact)} exact PAM curves, pEC50 correlated with intrinsic agonist RE (rho={c1['rho']:.3f}, bootstrap 95% CI [{c1['ci95'][0]:.3f}, {c1['ci95'][1]:.3f}]).",
             f"PAM RE also correlated with intrinsic agonist RE (rho={c2['rho']:.3f}, 95% CI [{c2['ci95'][0]:.3f}, {c2['ci95'][1]:.3f}]).", "",
             "This establishes a series-specific optimization trade-off: potency alone is an incomplete objective. A deployable ranker should jointly optimize PAM potency/efficacy, low intrinsic agonism, and subtype selectivity.", "",
             "## Core groups", "", core.to_markdown(index=False, floatfmt=".3f"), "",
             f"Across {len(pairs)} structure-matched 2-methyl to 2,5-dimethyl comparisons, potency increased in {summary['matched_pair_potency_increase_count']}/{len(pairs)} pairs (median {summary['matched_pair_median_fold_potency']:.1f}-fold; two-sided sign-test p={summary['matched_pair_potency_sign_test_two_sided_p']:.4f}). Intrinsic agonist RE increased in {summary['matched_pair_agonism_increase_count']}/{len(pairs)} pairs (median +{summary['matched_pair_median_intrinsic_agonist_RE_delta']:.1f} percentage points; p={summary['matched_pair_agonism_sign_test_two_sided_p']:.3f}).", "",
             f"A descriptive balanced-control rule retained {len(balanced)} known examples (EC50 <=60 nM, PAM RE >=80%, agonist RE <40%). This threshold was defined after viewing the campaign and is not a validated classifier.", "",
             "The Pareto and balanced-control files contain known patent controls, not newly discovered compounds."]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(pareto[keep].sort_values(["hM4_agonist_RE_pct", "pEC50"], ascending=[True, False]).to_string(index=False))


if __name__ == "__main__":
    main()
