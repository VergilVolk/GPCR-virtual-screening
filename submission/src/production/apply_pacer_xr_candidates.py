# -*- coding: utf-8 -*-
"""Create and apply the frozen six-channel PACER-XR candidate score contract.

No score is imputed.  New molecules must be processed with the same receptors,
software and score definitions as the Thompson-Miao public benchmark.  Raw
scores are converted to empirical rank percentiles against the frozen official
M4 reference distributions, then fused without fitting M4 labels.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


P = Path(__file__).resolve().parents[1]
ROOT = P / "tools" / "gpcr-am-ensemble-docking" / "docking_scores" / "M4R"
FINAL = P / "results" / "pacer_candidates_v01" / "final" / "final_candidate_hypotheses.csv"
OUT = P / "results" / "pacer_xr_candidate_transfer_v01"
TEMPLATE = OUT / "candidate_six_channel_raw_scores.csv"

FILES = {
    "Glide_PDB": (ROOT / "PDB" / "M4R_PDB_Glide_scores.csv", "glide_gscore"),
    "Glide_BEmin": (ROOT / "Ensemble" / "M4R_Ensemble_Glide_BEmin_ranked.csv", "BE_min"),
    "Glide_BEavg": (ROOT / "Ensemble" / "M4R_Ensemble_Glide_BEavg_ranked.csv", "BE_avg"),
    "Vina_PDB": (ROOT / "PDB" / "M4R_PDB_Vina_scores.csv", "vina_score"),
    "Vina_BEmin": (ROOT / "Ensemble" / "M4R_Ensemble_Vina_BEmin_ranked.csv", "BE_min"),
    "Vina_BEavg": (ROOT / "Ensemble" / "M4R_Ensemble_Vina_BEavg_ranked.csv", "BE_avg"),
}


def make_template() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(FINAL, usecols=["candidate_id", "canonical_smiles"])
    for name in FILES:
        candidates[f"{name}_raw"] = np.nan
    candidates["score_protocol_status"] = "awaiting_protocol_matched_six_channel_scores"
    candidates.to_csv(TEMPLATE, index=False)
    contract = {
        "required_channels": list(FILES),
        "missing_score_policy": "fail; no imputation and no partial XR claim",
        "direction": "lower raw energy is better for every released channel",
        "normalization": "empirical mid-rank percentile against each frozen official M4 distribution",
        "cascade": "official Glide-BEmin top-1% threshold, then six-channel equal-rank mean",
        "candidate_claim": "ranking hypothesis only; not functional PAM probability",
        "commercial_dependency": "three Glide channels require a protocol-matched Schrodinger Glide run",
    }
    (OUT / "SCORE_CONTRACT.json").write_text(json.dumps(contract, indent=2), encoding="utf-8")
    print(f"Wrote blank score contract: {TEMPLATE}")


def empirical_rank(raw: pd.Series, reference: np.ndarray) -> np.ndarray:
    """One is best; mid-rank handles exact ties deterministically."""
    ref = np.sort(reference[np.isfinite(reference)])
    x = raw.to_numpy(float)
    less = np.searchsorted(ref, x, side="left")
    right = np.searchsorted(ref, x, side="right")
    return 1.0 - (less + 0.5 * (right - less)) / len(ref)


def apply_scores(path: Path) -> None:
    scores = pd.read_csv(path)
    required = ["candidate_id", "canonical_smiles"] + [f"{x}_raw" for x in FILES]
    missing_cols = [x for x in required if x not in scores]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    if scores.candidate_id.duplicated().any():
        raise ValueError("candidate_id must be unique")
    raw_cols = [f"{x}_raw" for x in FILES]
    scores[raw_cols] = scores[raw_cols].apply(pd.to_numeric, errors="coerce")
    missing = scores[raw_cols].isna().sum()
    if missing.any():
        raise ValueError(
            "PACER-XR requires all six protocol-matched scores; no imputation is allowed. "
            f"Missing counts: {missing[missing > 0].to_dict()}"
        )

    reference_n = {}
    for name, (file, col) in FILES.items():
        ref = pd.read_csv(file, usecols=[col])[col].to_numpy(float)
        reference_n[name] = int(np.isfinite(ref).sum())
        scores[f"{name}_rank"] = empirical_rank(scores[f"{name}_raw"], ref)
    scores["PACER_XR"] = scores[[f"{x}_rank" for x in FILES]].mean(axis=1)

    glide_ref = pd.read_csv(FILES["Glide_BEmin"][0], usecols=[FILES["Glide_BEmin"][1]])[
        FILES["Glide_BEmin"][1]
    ].dropna().to_numpy(float)
    top1_threshold = float(np.quantile(glide_ref, 0.01, method="higher"))
    scores["cascade_head_top1pct"] = scores.Glide_BEmin_raw <= top1_threshold
    scores["cascade_stage"] = np.where(scores.cascade_head_top1pct, 1, 2)
    # Within stage 1 preserve Glide-BEmin; within stage 2 sort by XR.
    scores["cascade_stage_score"] = np.where(
        scores.cascade_head_top1pct, scores.Glide_BEmin_rank, scores.PACER_XR
    )
    scores = scores.sort_values(
        ["cascade_stage", "cascade_stage_score"],
        ascending=[True, False],
        kind="stable",
    ).reset_index(drop=True)
    scores["candidate_cascade_rank"] = np.arange(1, len(scores) + 1)
    scores["claim_boundary"] = "broad-AM retrieval hypothesis; functional PAM assay required"
    scores.to_csv(OUT / "candidate_pacer_xr_ranking.csv", index=False)
    audit = {
        "n_candidates": int(len(scores)),
        "all_six_channels_complete": True,
        "reference_n": reference_n,
        "Glide_BEmin_official_top1pct_raw_threshold": top1_threshold,
        "cascade_head_candidates": int(scores.cascade_head_top1pct.sum()),
        "functional_PAM_probability": False,
        "claim_boundary": "This ranking cannot establish PAM efficacy or cooperativity.",
    }
    (OUT / "application_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(scores[["candidate_cascade_rank", "candidate_id", "cascade_stage", "PACER_XR"]].to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, help="Completed six-channel CSV; defaults to the template")
    ap.add_argument("--create-template", action="store_true")
    args = ap.parse_args()
    if args.create_template or not TEMPLATE.exists():
        make_template()
    if args.input:
        apply_scores(args.input)


if __name__ == "__main__":
    main()
