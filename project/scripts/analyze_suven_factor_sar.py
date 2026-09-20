"""Diagnose the independent Suven matrix as headgroup x receptor-side core SAR."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.preprocessing import OneHotEncoder


P = Path(__file__).resolve().parents[1]
DATA = P / "data" / "benchmarks" / "m4_pam_v1" / "external_suven_2025" / "external_patent_potency.csv"
PRED = P / "results" / "pacer_external_suven_v01" / "predictions.csv"
OUT = P / "results" / "pacer_suven_factor_sar_v01"


def headgroup(name: str) -> str:
    rules = [
        ("Dihydro-benzo[1,4]dioxin-6", "benzodioxin6"),
        ("Dihydro-isobenzofuran-5", "isobenzofuran5"),
        ("Dihydro-benzofuran-5", "benzofuran5"),
        ("Dihydro-benzofuran-6", "benzofuran6"),
        ("6-Methoxy-pyridin-3", "methoxypyridin3"),
        ("4-Methoxy-phenoxy", "methoxyphenyl4"),
        ("6-Fluoro-pyridin-3", "fluoropyridin3"),
        ("6-methyl-pyridin-3", "methylpyridin3"),
        ("benzonitrile", "benzonitrile3"),
    ]
    for token, label in rules:
        if token.lower() in name.lower():
            return label
    raise ValueError(f"Unclassified headgroup: {name}")


def core_state(name: str) -> tuple[str, str]:
    lower = name.lower()
    family = (
        "triazolo43" if "triazolo[4,3-b]" in lower else
        "triazolo15" if "triazolo[1,5-b]" in lower else
        "imidazo12" if "imidazo[1,2-b]" in lower else "unknown"
    )
    substitutions = [
        "3-methoxymethyl-7,8-dimethyl", "3-difluoromethyl", "3,7,8-trimethyl",
        "3,7-dimethyl", "2,7,8-trimethyl", "2,7-dimethyl", "7,8-dimethyl", "7-methyl",
    ]
    state = next((token for token in substitutions if token in lower), "unresolved")
    return family, f"{family}:{state}"


def rho(y, prediction):
    value = spearmanr(y, prediction).statistic
    return float(value) if np.isfinite(value) else 0.0


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(PRED)
    frame["headgroup"] = frame.name.map(headgroup)
    states = frame.name.map(core_state)
    frame["core_family"] = [item[0] for item in states]
    frame["core_state"] = [item[1] for item in states]

    # Label-held-out two-way additive model. This is diagnostic only and is not
    # an external result because the patent labels motivated this analysis.
    design = frame[["headgroup", "core_state"]]
    target = frame.cre_luc_pEC50.to_numpy(float)
    transformer = ColumnTransformer(
        [("factor", OneHotEncoder(handle_unknown="ignore"), ["headgroup", "core_state"])],
        remainder="drop",
    )
    encoded = transformer.fit_transform(design)
    additive = cross_val_predict(Ridge(alpha=1.0), encoded, target, cv=LeaveOneOut())
    frame["FactorSAR_additive_LOO"] = additive

    head_summary = frame.groupby("headgroup").agg(
        n=("compound_id", "size"), mean_pEC50=("cre_luc_pEC50", "mean"),
        sd_pEC50=("cre_luc_pEC50", "std"), median_EC50_nM=("cre_luc_ec50_nM", "median"),
    ).reset_index().sort_values("mean_pEC50", ascending=False)
    core_summary = frame.groupby(["core_family", "core_state"]).agg(
        n=("compound_id", "size"), mean_pEC50=("cre_luc_pEC50", "mean"),
        sd_pEC50=("cre_luc_pEC50", "std"), median_EC50_nM=("cre_luc_ec50_nM", "median"),
    ).reset_index().sort_values("mean_pEC50", ascending=False)

    # Matched headgroup effects across core states expose sign reversals.
    matched = []
    for head, subset in frame.groupby("headgroup"):
        rows = list(subset.itertuples())
        for i, left in enumerate(rows):
            for right in rows[i + 1:]:
                matched.append({
                    "headgroup": head,
                    "left_id": left.compound_id,
                    "right_id": right.compound_id,
                    "left_core": left.core_state,
                    "right_core": right.core_state,
                    "delta_pEC50_right_minus_left": right.cre_luc_pEC50 - left.cre_luc_pEC50,
                })
    matched_frame = pd.DataFrame(matched)

    frame.to_csv(OUT / "factorized_compounds.csv", index=False)
    head_summary.to_csv(OUT / "headgroup_summary.csv", index=False)
    core_summary.to_csv(OUT / "core_state_summary.csv", index=False)
    matched_frame.to_csv(OUT / "matched_core_transitions.csv", index=False)
    audit = {
        "status": "post_hoc_factor_sar_diagnostic",
        "n": len(frame),
        "n_headgroups": int(frame.headgroup.nunique()),
        "n_core_states": int(frame.core_state.nunique()),
        "n_core_families": int(frame.core_family.nunique()),
        "factor_sar_additive_loo_spearman": rho(target, additive),
        "factor_sar_additive_loo_mae": float(mean_absolute_error(target, additive)),
        "absolute_qsar_spearman": rho(target, frame["Absolute-QSAR"].to_numpy(float)),
        "warning": "Post-hoc developmental analysis; not independent external validation.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    print("\nHeadgroups:\n", head_summary.to_string(index=False))
    print("\nCore states:\n", core_summary.to_string(index=False))


if __name__ == "__main__":
    main()
