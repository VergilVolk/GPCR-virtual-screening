# -*- coding: utf-8 -*-
"""Audit semantic overlap between ASD AM labels and assay-specific M4 PAM labels."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

P = Path(__file__).resolve().parents[1]
STRICT = P / "data" / "benchmarks" / "m4_pam_v1" / "pam_vs_inactive.csv"
RAW = P.parent / "CHRM4_PAM_Modeling_Handoff_v1.0" / "data" / "modeling_pam_vs_inactive.csv"
MIAO = P / "data" / "miao2026_m4r" / "M4R_Ensemble_Vina_framewise_scores.csv"
OUT = P / "results" / "m4_gamd_ensemble"


def main() -> None:
    strict = pd.read_csv(STRICT)
    raw = pd.read_csv(RAW, low_memory=False)
    miao = pd.read_csv(MIAO, usecols=["ligand_id", "SMILES", "is_active"])
    grouped = miao.groupby("SMILES").agg(
        official_ids=("ligand_id", lambda x: ";".join(map(str, x))),
        official_min_label=("is_active", "min"),
        official_max_label=("is_active", "max"),
        official_occurrences=("is_active", "size"),
    ).reset_index()
    overlap = strict.merge(grouped, left_on="canonical_smiles", right_on="SMILES", how="inner")
    overlap["official_status"] = overlap.apply(
        lambda r: "both_active_and_decoy" if r.official_min_label == 0 and r.official_max_label == 1
        else ("active_only" if r.official_max_label == 1 else "decoy_only"), axis=1)
    evidence = raw[["canonical_molecule_id", "assay_id", "readout", "receptor_species",
                    "negative_evidence_tier", "negative_evidence_type", "notes"]].drop_duplicates("canonical_molecule_id")
    overlap = overlap.merge(evidence, on="canonical_molecule_id", how="left")
    overlap["semantic_discordance"] = (
        ((overlap.target == 0) & (overlap.official_max_label == 1)) |
        ((overlap.target == 1) & (overlap.official_max_label == 0))
    )
    overlap.to_csv(OUT / "miao_vs_functional_label_crosswalk.csv", index=False)

    neg = overlap[overlap.target == 0]
    pos = overlap[overlap.target == 1]
    summary = {
        "strict_molecules": int(len(strict)),
        "exact_smiles_overlap_unique": int(len(overlap)),
        "overlap_positive": int(len(pos)),
        "overlap_negative": int(len(neg)),
        "negative_tagged_active_at_least_once": int(((neg.official_max_label == 1)).sum()),
        "negative_active_discordance_rate": float((neg.official_max_label == 1).mean()) if len(neg) else None,
        "positive_decoy_only": int((pos.official_max_label == 0).sum()),
        "appears_as_both_active_and_decoy": int((overlap.official_status == "both_active_and_decoy").sum()),
        "discordance_by_negative_tier": overlap[(overlap.target == 0) & overlap.semantic_discordance]
            .negative_evidence_tier.value_counts(dropna=False).to_dict(),
        "interpretation": (
            "ASD active means broad allosteric-modulator curation, while the strict endpoint is "
            "assay-specific functional potentiation. The labels are not interchangeable."
        ),
        "caution": (
            "Most discordant negatives are database-text tier and need primary-table verification; "
            "this is label-semantic discordance, not by itself proof that either source is wrong."
        ),
    }
    (OUT / "miao_functional_label_overlap_audit.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    report = f"""# ASD allosteric label vs M4 functional PAM label audit

- Strict assay-specific set: {summary['strict_molecules']} molecules.
- Exact-SMILES overlap with the Miao/ASD docking library: {summary['exact_smiles_overlap_unique']}.
- Among {summary['overlap_negative']} overlapping strict functional negatives, {summary['negative_tagged_active_at_least_once']} ({summary['negative_active_discordance_rate']:.1%}) are tagged as active at least once in the ASD benchmark.
- Among {summary['overlap_positive']} overlapping functional positives, {summary['positive_decoy_only']} occur only as decoys.
- {summary['appears_as_both_active_and_decoy']} structures appear in both official active and decoy rows.

## Meaning

The two labels answer different questions. ASD membership is a broad literature-level allosteric-modulator label; our strict target is ACh-probe, human-M4 functional potentiation under a specified assay. Therefore the public benchmark can validate docking enrichment against its own curation, but it is not an uncontaminated ground truth for PAM efficacy.

Most negative discordances come from database-text “Not Active” records whose primary per-compound tables were not manually verified. They remain an important warning, not a verdict that ASD or ChEMBL is definitively wrong. The primary conclusion must use only A-tier primary-source-confirmed inactive compounds.
"""
    (OUT / "MIAO_FUNCTIONAL_LABEL_AUDIT.md").write_text(report, encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
