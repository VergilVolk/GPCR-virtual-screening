# -*- coding: utf-8 -*-
"""Apply the frozen PACER-FS relative-SAR head to the 24 candidate hypotheses.

Before functional results exist, this emits only within-reference-series ranks.
After three QC-passing lead-series anchor EC50 values are entered in the round-1
template, it estimates one series offset and emits calibrated pEC50 hypotheses.
No offset is transferred to another chemotype/reference series.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor

from run_pacer_fs_baselines import DATA, features


P = Path(__file__).resolve().parents[1]
CANDIDATES = P / "results" / "pacer_candidates_v01" / "final" / "final_candidate_hypotheses.csv"
ROUND1 = P / "results" / "pacer_assay_closed_loop_v01" / "round1_results_template.csv"
AUDIT = P / "results" / "pacer_assay_closed_loop_v01" / "selection_audit.json"
OUT = P / "results" / "pacer_assay_closed_loop_v01"


def qc_true(value):
    return str(value).strip().lower() in {"1", "true", "yes", "pass", "passed", "y"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=ROUND1)
    args = ap.parse_args()

    train = pd.read_csv(DATA).reset_index(drop=True)
    cand = pd.read_csv(CANDIDATES).reset_index(drop=True)
    series_map = dict(zip(train.canonical_molecule_id, train.source_component))
    cand["reference_series"] = cand.nearest_known_id.map(series_map)

    all_smiles = pd.concat([train.canonical_smiles, cand.canonical_smiles], ignore_index=True)
    X, D, _ = features(all_smiles)
    n = len(train)
    y = train.pEC50.to_numpy(float)
    g = train.source_component.astype(str).to_numpy()
    centered = np.zeros(n, float)
    for group in set(g):
        idx = np.where(g == group)[0]
        centered[idx] = y[idx] - y[idx].mean()
    sel = np.argsort(X[:n].var(0))[-1024:]
    xx = np.c_[X[:, sel], D]
    model = LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1, random_state=42,
        verbosity=-1, n_jobs=6,
    ).fit(xx[:n], centered)
    cand["pacer_relative_sar"] = model.predict(xx[n:])
    cand["within_series_rank"] = cand.groupby("reference_series")["pacer_relative_sar"].rank(
        method="first", ascending=False
    ).astype(int)
    cand["calibrated_pEC50"] = np.nan
    cand["calibration_status"] = "relative_only_no_functional_anchors"

    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    lead = audit["lead_reference_series"]
    offset = None
    valid = []
    if args.results.exists():
        r = pd.read_csv(args.results, dtype=str).fillna("")
        anchors = r[r.assay_role == "calibration_anchor"]
        for _, row in anchors.iterrows():
            try:
                ec50 = float(row.pam_ec50_nM)
            except (TypeError, ValueError):
                continue
            if ec50 > 0 and qc_true(row.qc_pass):
                p_ec50 = 9.0 - np.log10(ec50)
                rel = float(cand.loc[cand.candidate_id == row.candidate_id, "pacer_relative_sar"].iloc[0])
                valid.append((row.candidate_id, p_ec50 - rel))
    if len(valid) == 3:
        offset = float(np.mean([x[1] for x in valid]))
        mask = cand.reference_series == lead
        cand.loc[mask, "calibrated_pEC50"] = cand.loc[mask, "pacer_relative_sar"] + offset
        cand.loc[mask, "calibration_status"] = "three_anchor_series_calibrated_hypothesis"
        cand.loc[~mask, "calibration_status"] = "different_chemotype_offset_not_transferable"

    cols = [
        "candidate_id", "canonical_smiles", "reference_series", "portfolio_role",
        "pacer_relative_sar", "within_series_rank", "calibrated_pEC50",
        "calibration_status", "potency_reference_lower90", "potency_reference_upper90",
        "strict_inactive_risk_ref", "max_tanimoto",
    ]
    cand[cols].sort_values(["reference_series", "within_series_rank"]).to_csv(
        OUT / "candidate_pacer_fs_predictions.csv", index=False
    )
    payload = {
        "lead_reference_series": lead,
        "valid_anchor_count": len(valid),
        "anchor_ids_used": [x[0] for x in valid],
        "series_offset": offset,
        "status": "calibrated" if offset is not None else "awaiting_three_qc_passing_anchor_results",
        "boundary": "Relative SAR is not PAM identity; calibration never crosses reference series.",
    }
    (OUT / "candidate_calibration_audit.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
