"""Build an honest multi-endpoint functional-risk and assay panel for PACER candidates.

The calculated endpoint estimates are weak, similarity-weighted priors.  They are
used to allocate experiments, never to label a generated molecule as a PAM.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"
OUT = RESULTS / "pacer_candidate_function_space_v01"
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fp(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return FPGEN.GetFingerprint(mol)


def weighted_prior(query_fp, reference: pd.DataFrame, endpoint: str) -> dict:
    usable = reference[pd.to_numeric(reference[endpoint], errors="coerce").notna()].copy()
    values = pd.to_numeric(usable[endpoint]).to_numpy(float)
    similarities = np.asarray(DataStructs.BulkTanimotoSimilarity(
        query_fp, [fp(x) for x in usable.canonical_smiles]), float)
    order = np.argsort(similarities)[::-1][: min(5, len(similarities))]
    weights = np.maximum(similarities[order], 1e-3) ** 4
    mean = float(np.average(values[order], weights=weights))
    sd = float(np.sqrt(np.average((values[order] - mean) ** 2, weights=weights)))
    return {
        f"prior_{endpoint}": mean,
        f"prior_{endpoint}_neighbor_sd": sd,
        f"prior_{endpoint}_max_tanimoto": float(similarities[order[0]]),
        f"prior_{endpoint}_nearest_id": str(usable.iloc[order[0]].compound_id),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(
        RESULTS / "pacer_assay_closed_loop_v01" / "round1_assay_batch.csv")
    candidates = candidates[candidates.is_control.eq(0)].copy().reset_index(drop=True)
    acadia = pd.read_csv(PROJECT / "data" / "benchmarks" / "m4_pam_v1" /
                         "external_acadia_2025" / "external_molecule_benchmark.csv")

    priors = []
    for row in candidates.itertuples(index=False):
        query_fp = fp(row.canonical_smiles)
        item = {"candidate_id": row.candidate_id}
        for endpoint in ["pEC50", "hM4_PAM_RE_pct", "hM4_agonist_RE_pct"]:
            item.update(weighted_prior(query_fp, acadia, endpoint))
        priors.append(item)
    frame = candidates.merge(pd.DataFrame(priors), on="candidate_id", how="left")

    # A deliberately conservative domain flag.  No transferred functional prior
    # is considered decision-grade below 0.50 similarity to its nearest labelled
    # molecule; uncertainty then increases experimental priority rather than
    # becoming a confident negative prediction.
    frame["functional_prior_domain"] = np.where(
        frame.prior_hM4_agonist_RE_pct_max_tanimoto >= .50,
        "local_reference", "extrapolative_not_decision_grade")
    width = frame.potency_reference_upper90 - frame.potency_reference_lower90
    width_scaled = (width - width.min()) / max(width.max() - width.min(), 1e-9)
    novelty = 1 - frame.max_tanimoto.clip(0, 1)
    dynamic = (frame.vina_sd - frame.vina_sd.min()) / max(
        frame.vina_sd.max() - frame.vina_sd.min(), 1e-9)
    exploration = frame.portfolio_role.eq("exploratory_hypothesis").astype(float)
    frame["information_gain_priority"] = (
        .35 * width_scaled + .30 * novelty + .20 * dynamic + .15 * exploration)
    frame["functional_identity_status"] = "unknown_requires_m4_pam_assay"
    frame["intrinsic_agonism_status"] = "unknown_requires_candidate_alone_assay"
    frame["m2_selectivity_status"] = "unknown_requires_matched_m2_assay"
    frame = frame.sort_values("information_gain_priority", ascending=False)
    frame.to_csv(OUT / "candidate_functional_risk.csv", index=False)

    endpoints = [
        ("M4_PAM_ACh_EC20", "primary", "EC50_nM; Emax_pct", "Does ACh-potentiating activity exist?"),
        ("M4_candidate_alone", "primary", "EC50_nM; Emax_pct", "Intrinsic agonism liability"),
        ("M2_PAM_ACh_EC20", "primary", "EC50_nM; Emax_pct", "M4 over M2 functional selectivity"),
        ("M4_PAM_ACh_EC80", "secondary", "EC50_nM; Emax_pct", "Probe-window/cooperativity dependence"),
        ("M4_GTPgammaS_confirmation", "secondary_top_hits", "EC50_nM; Emax_pct", "Orthogonal proximal G-protein confirmation"),
    ]
    panel = []
    for row in frame.itertuples(index=False):
        for endpoint, tier, fields, question in endpoints:
            panel.append({"candidate_id": row.candidate_id,
                          "information_gain_priority": row.information_gain_priority,
                          "assay_endpoint": endpoint, "assay_tier": tier,
                          "required_fields": fields, "scientific_question": question,
                          "result_status": "not_measured"})
    pd.DataFrame(panel).to_csv(OUT / "minimal_function_space_assay_panel.csv", index=False)

    audit = {
        "n_candidates": int(len(frame)),
        "primary_endpoints_per_candidate": 3,
        "secondary_endpoints": 2,
        "transferred_priors_are_pam_labels": False,
        "decision_grade_similarity_threshold": 0.50,
        "n_candidates_with_local_functional_prior": int(
            frame.functional_prior_domain.eq("local_reference").sum()),
        "claim_boundary": (
            "Similarity-weighted Acadia values are experimental-design priors only; "
            "PAM identity, agonism, and selectivity remain unmeasured."),
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = [
        "# PACER candidate function-space panel", "",
        "The generated candidates are computational hypotheses, not predicted PAM confirmations.", "",
        frame[["candidate_id", "portfolio_role", "information_gain_priority",
               "prior_hM4_agonist_RE_pct", "prior_hM4_agonist_RE_pct_max_tanimoto",
               "functional_prior_domain"]].to_markdown(index=False, floatfmt=".3f"), "",
        "Every candidate receives matched M4 PAM, candidate-alone M4, and M2 PAM assays. "
        "ACh EC80 and GTPgammaS are confirmatory endpoints, not replacements for the primary panel.", "",
        audit["claim_boundary"],
    ]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(frame[["candidate_id", "information_gain_priority",
                 "prior_hM4_agonist_RE_pct_max_tanimoto",
                 "functional_prior_domain"]].to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
