# -*- coding: utf-8 -*-
"""Audit PACER-XR against assay-specific M4 PAM labels without label inflation.

This is deliberately a semantic stress test, not a promotable functional
benchmark.  The exact-SMILES overlap contains no A-tier confirmed inactive
controls; all overlapping negatives are B-tier database-text records and are
also labelled active in the ASD benchmark used to construct the docking set.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


P = Path(__file__).resolve().parents[1]
CROSSWALK = P / "results" / "m4_gamd_ensemble" / "miao_vs_functional_label_crosswalk.csv"
PREDICTIONS = P / "results" / "pacer_xr_v01" / "m4_blind_test_predictions.csv"
FUNCTIONAL = P / "data" / "benchmarks" / "m4_pam_v1" / "pam_vs_inactive.csv"
OUT = P / "results" / "pacer_xr_functional_semantics_v01"
OUT.mkdir(parents=True, exist_ok=True)

METHODS = (
    "PACER_XR_Cascade1pct",
    "PACER_XR",
    "Glide_BEmin",
    "Vina_BEmin",
)
LEAD_NEIGHBOURS = {
    "PACER0076": "CM00240",
    "PACER0057": "CM00202",
    "PACER0026": "CM00594",
}


def bootstrap_auc(y: np.ndarray, score: np.ndarray, seed: int = 829, n_boot: int = 2000):
    """Stratified bootstrap CI; descriptive only for this non-independent set."""
    rng = np.random.default_rng(seed)
    pos = np.flatnonzero(y == 1)
    neg = np.flatnonzero(y == 0)
    values = []
    for _ in range(n_boot):
        idx = np.r_[rng.choice(pos, len(pos), replace=True), rng.choice(neg, len(neg), replace=True)]
        values.append(roc_auc_score(y[idx], score[idx]))
    return [float(x) for x in np.quantile(values, [0.025, 0.975])]


def main() -> None:
    cross = pd.read_csv(CROSSWALK)
    pred = pd.read_csv(PREDICTIONS)
    functional = pd.read_csv(FUNCTIONAL)

    # A molecule can occur more than once in the official library.  Aggregate
    # rank scores by exact SMILES instead of selecting an arbitrary row or the
    # most favourable occurrence.
    pred_mol = pred.groupby("SMILES", as_index=False).agg(
        official_rows=("ligand_id", "size"),
        official_ids_in_xr=("ligand_id", lambda x: ";".join(map(str, x))),
        **{method: (method, "mean") for method in METHODS},
    )
    overlap = cross.merge(pred_mol, left_on="canonical_smiles", right_on="SMILES", how="inner")
    overlap["negative_is_A_tier"] = overlap.negative_evidence_tier.fillna("").str.startswith("A_")
    overlap["negative_is_B_tier"] = overlap.negative_evidence_tier.fillna("").str.startswith("B_")
    overlap.to_csv(OUT / "overlap_predictions.csv", index=False)

    y = overlap.target.to_numpy(int)
    metric_rows = []
    for method in METHODS:
        s = overlap[method].to_numpy(float)
        metric_rows.append({
            "method": method,
            "n_molecules": int(len(overlap)),
            "n_positive": int(y.sum()),
            "n_negative": int((y == 0).sum()),
            "ROC_AUC_apparent": float(roc_auc_score(y, s)),
            "ROC_AUC_bootstrap_CI_descriptive": json.dumps(bootstrap_auc(y, s)),
            "PR_AUC_apparent": float(average_precision_score(y, s)),
            "promotable_functional_metric": False,
        })
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "apparent_metrics_nonpromotable.csv", index=False)

    # Audit whether our lead chemotype is actually represented near the top of
    # the broad-AM benchmark.  Scores are rank percentiles where one is best.
    func_by_id = functional.set_index("canonical_molecule_id")
    lead_rows = []
    for lead, neighbour in LEAD_NEIGHBOURS.items():
        smiles = func_by_id.loc[neighbour, "canonical_smiles"]
        hit = pred_mol[pred_mol.SMILES == smiles]
        row = {
            "candidate_id": lead,
            "nearest_known_functional_PAM": neighbour,
            "mapped_to_official_XR": bool(len(hit)),
        }
        if len(hit):
            h = hit.iloc[0]
            row["official_ids"] = h.official_ids_in_xr
            for method in METHODS:
                row[f"{method}_rank_percentile"] = float(h[method])
                row[f"{method}_top_fraction_pct"] = float((1.0 - h[method]) * 100.0)
        lead_rows.append(row)
    leads = pd.DataFrame(lead_rows)
    leads.to_csv(OUT / "lead_chemotype_transfer_audit.csv", index=False)

    neg = overlap[overlap.target == 0]
    audit = {
        "analysis_role": "semantic_stress_test_only",
        "exact_smiles_crosswalk_rows": int(len(cross)),
        "mapped_unique_molecules": int(len(overlap)),
        "unmapped_unique_molecules": int(len(cross) - len(overlap)),
        "functional_positive": int((overlap.target == 1).sum()),
        "functional_negative": int((overlap.target == 0).sum()),
        "A_tier_confirmed_negative": int(neg.negative_is_A_tier.sum()),
        "B_tier_database_text_negative": int(neg.negative_is_B_tier.sum()),
        "negative_also_ASD_active": int((neg.official_max_label == 1).sum()),
        "negative_also_ASD_active_fraction": float((neg.official_max_label == 1).mean()),
        "class_prevalence": float(overlap.target.mean()),
        "aggregation": "mean rank score across duplicate official rows sharing exact SMILES",
        "primary_functional_benchmark": (
            "A-tier 27 PAM + 27 primary-confirmed inactive set in m4_gamd_ensemble; "
            "cluster00/BEmin ROC AUC 0.499 and BEavg 0.457"
        ),
        "claim_boundary": (
            "The apparent overlap AUC cannot establish functional PAM prediction because no A-tier "
            "negative is present, all negatives are lower-tier database-text records, and those same "
            "structures are labelled active in the ASD benchmark."
        ),
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")

    m = metrics.set_index("method")
    c = m.loc["PACER_XR_Cascade1pct"]
    v = m.loc["Vina_BEmin"]
    g = m.loc["Glide_BEmin"]
    lead_lines = "\n".join(
        f"- {r.candidate_id} nearest known PAM {r.nearest_known_functional_PAM}: "
        f"cascade top {r.PACER_XR_Cascade1pct_top_fraction_pct:.2f}%, "
        f"XR top {r.PACER_XR_top_fraction_pct:.2f}%, "
        f"Glide-BEmin top {r.Glide_BEmin_top_fraction_pct:.2f}%."
        for r in leads.itertuples()
    )
    report = f"""# PACER-XR functional-semantic stress test

## Result

Exact-SMILES mapping retained **{len(overlap)}/{len(cross)}** functional molecules: {int(y.sum())} positives and {int((y == 0).sum())} negatives. Molecule-level mean-rank aggregation gives:

- PACER-XR Cascade apparent ROC AUC **{c.ROC_AUC_apparent:.3f}**, PR AUC **{c.PR_AUC_apparent:.3f}**;
- PACER-XR apparent ROC AUC **{m.loc['PACER_XR'].ROC_AUC_apparent:.3f}**;
- Vina-BEmin apparent ROC AUC **{v.ROC_AUC_apparent:.3f}**;
- Glide-BEmin apparent ROC AUC **{g.ROC_AUC_apparent:.3f}**.

These numbers are **not a valid functional SOTA claim**. All {len(neg)} negatives are B-tier database-text “inactive” records, none is an A-tier primary-table-confirmed inactive, and all {audit['negative_also_ASD_active']} are labelled active at least once in the ASD-derived docking benchmark. The positive prevalence is {overlap.target.mean():.1%}. The apparent AUC therefore measures a mixture of provenance and label semantics, not clean PAM efficacy discrimination.

## Lead-domain transfer

{lead_lines}

The broad-AM benchmark does not place every nearest neighbour of our lead series near the top. Therefore PACER-XR's public-benchmark improvement cannot be transferred to PACER0076/0057/0026 without computing the missing official score channels and, ultimately, running a probe-dependent functional assay.

## Promotable conclusion

PACER-XR Cascade is a deterministic improvement on the authors' public **broad allosteric-modulator enrichment** benchmark. It is not yet a validated functional-PAM classifier. The primary clean functional test remains the balanced A-tier 27+27 panel, where static/ensemble docking is approximately random. This separation is the scientific reason PACER combines broad-AM retrieval, few-shot within-series potency modelling, receptor-state evidence, and an explicit assay gate instead of calling docking score “PAM efficacy”.
"""
    (OUT / "VALIDATION_REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(audit, indent=2))
    print(metrics.to_string(index=False))
    print(leads.to_string(index=False))


if __name__ == "__main__":
    main()
