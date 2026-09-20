# -*- coding: utf-8 -*-
"""Join GaMD ensemble evidence to PACER candidates and compute an honest Pareto set."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analyze_m4_gamd_ensemble import aggregate

P = Path(__file__).resolve().parents[1]
ROOT = P / "results" / "m4_gamd_ensemble"
SOURCE = P / "results" / "pacer_candidates_v01" / "final" / "final_candidate_hypotheses.csv"
FRAME = ROOT / "candidates_docking_protocol_v127_ph7_rank1" / "framewise_scores.csv"


def pareto_front(df: pd.DataFrame) -> np.ndarray:
    """Minimize all columns after converting favorable directions."""
    x = np.column_stack([
        df.BE_min.to_numpy(float),
        -df.pocket_coverage_mean.to_numpy(float),
        -df.potency_reference_lower90.to_numpy(float),
        df.strict_inactive_risk_ref.to_numpy(float),
    ])
    keep = np.ones(len(x), dtype=bool)
    for i in range(len(x)):
        dominated = np.any(np.all(x <= x[i], axis=1) & np.any(x < x[i], axis=1))
        keep[i] = not dominated
    return keep


def main() -> None:
    candidates = pd.read_csv(SOURCE)
    frame = pd.read_csv(FRAME)
    frame = frame[(frame.error.isna() | frame.error.eq("")) & frame.vina_affinity.notna()]
    dynamic = aggregate(frame)
    out = candidates.merge(dynamic, left_on="candidate_id", right_on="molecule_id", validate="one_to_one")
    out = out[out.complete_ensemble].copy()
    out["dynamic_pareto"] = pareto_front(out).astype(int)
    out["claim"] = "computational PAM hypothesis; dynamic binding evidence is not functional efficacy"
    out = out.sort_values(["dynamic_pareto", "potency_reference_lower90", "BE_min"],
                          ascending=[False, False, True])
    keep = [
        "candidate_id", "canonical_smiles", "portfolio_role", "generator", "dynamic_pareto",
        "potency_reference_lower90", "potency_reference_upper90", "strict_inactive_risk_ref",
        "max_tanimoto", "cluster00_vina", "BE_min", "BE_avg", "best_cluster", "vina_sd",
        "vina_range", "pocket_coverage_mean", "pocket_coverage_min", "claim",
    ]
    out[keep].to_csv(ROOT / "candidate_ensemble_evidence.csv", index=False)
    rho, p = spearmanr(out.vina_affinity, out.BE_min, nan_policy="omit")
    audit = {
        "n_complete": int(len(out)), "n_dynamic_pareto": int(out.dynamic_pareto.sum()),
        "single_7TRS_vs_GaMD_BEmin_spearman": {"rho": float(rho), "p": float(p)},
        "ranking_policy": "Pareto only; no fitted weighted sum on prospective candidates.",
        "objectives": ["minimize BE_min", "maximize pocket coverage", "maximize lower90 potency reference",
                       "minimize strict inactive risk"],
        "claim_boundary": "No candidate is called a confirmed PAM without ternary functional assay.",
    }
    (ROOT / "candidate_ensemble_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
