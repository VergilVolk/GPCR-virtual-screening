# -*- coding: utf-8 -*-
"""Dock 24 PACER candidates to 7TRQ and annotate preregistered residue hypotheses."""
from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

import dock_potency_coupling_features as docking


P = Path(__file__).resolve().parents[1]
CANDIDATES = P / "results" / "pacer_candidates_v01" / "final" / "final_candidate_hypotheses.csv"
REFERENCE = P / "results" / "pacer_structure_loso_v01" / "docking_features.csv"
OUT = P / "results" / "pacer_candidate_7trq_mechanism_v01"
POSES = OUT / "poses"


def favorable_percentile(reference, value, lower_is_better):
    ref = np.asarray(reference, float)
    return float(np.mean(ref >= value) if lower_is_better else np.mean(ref <= value))


def main():
    OUT.mkdir(parents=True, exist_ok=True); POSES.mkdir(parents=True, exist_ok=True)
    docking.TMP = POSES
    candidates = pd.read_csv(CANDIDATES)
    tasks = [(r.candidate_id, r.canonical_smiles, 42, 4) for r in candidates.itertuples()]
    rows = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(docking.dock_one, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result(); rows.append(row)
            print(row["canonical_molecule_id"], row.get("error", ""), flush=True)
    result = pd.DataFrame(rows).rename(columns={"canonical_molecule_id": "candidate_id"})
    result = candidates.merge(result, on="candidate_id", how="left")

    reference = pd.read_csv(REFERENCE)
    reference = reference[
        reference.canonical_molecule_id.astype(str).str.startswith("CM")
        & reference.error.fillna("").eq("")
    ].drop_duplicates("canonical_molecule_id")
    reference["res184_pairs_per_heavy"] = reference.res184_pairs / reference.n_heavy_atoms
    result["res184_pairs_per_heavy"] = result.res184_pairs / result.n_heavy_atoms
    result["Y89_proximity_favorable_percentile"] = [
        favorable_percentile(reference.res89_min_A, x, True) for x in result.res89_min_A
    ]
    result["Q184_contact_favorable_percentile"] = [
        favorable_percentile(reference.res184_pairs_per_heavy, x, False)
        for x in result.res184_pairs_per_heavy
    ]
    result["D432_proximity_favorable_percentile"] = [
        favorable_percentile(reference.res432_min_A, x, True) for x in result.res432_min_A
    ]
    cols = ["Y89_proximity_favorable_percentile", "Q184_contact_favorable_percentile",
            "D432_proximity_favorable_percentile"]
    result["mechanism_checks_above_median"] = (result[cols] >= .5).sum(axis=1)
    result["all_three_geometry_checks"] = result.mechanism_checks_above_median.eq(3)
    result["annotation_only"] = True
    result["claim_boundary"] = (
        "7TRQ pose consistency; associations did not survive physchem-adjusted global FDR"
    )
    result = result.sort_values(
        ["all_three_geometry_checks", "mechanism_checks_above_median", "pocket_residue_coverage"],
        ascending=[False, False, False],
    )
    result.to_csv(OUT / "candidate_mechanism_annotations.csv", index=False)
    audit = {
        "n_candidates": int(len(result)),
        "docking_failures": int(result.error.fillna("").ne("").sum()),
        "all_three_geometry_checks": int(result.all_three_geometry_checks.sum()),
        "use_as_ranking_score": False,
        "reason": "physchem-adjusted residue associations did not pass global FDR",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    show = ["candidate_id", *cols, "mechanism_checks_above_median", "all_three_geometry_checks"]
    report = "# PACER candidate 7TRQ mechanism annotation\n\n"
    report += result[show].to_markdown(index=False, floatfmt=".3f")
    report += (
        "\n\nThese are geometry-consistency annotations, not a PAM score. "
        "Y89/Q184/D432 associations passed the primary series-controlled analysis but not the "
        "physchem-adjusted global FDR, so they cannot override PACER-FS or functional assays.\n"
    )
    (OUT / "VALIDATION_REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(audit, indent=2)); print(result[show].to_string(index=False))


if __name__ == "__main__":
    main()
