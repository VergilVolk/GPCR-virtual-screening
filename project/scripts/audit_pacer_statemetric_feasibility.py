"""Quantify whether existing M4 data can support functional hard-negative metric learning."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "benchmarks" / "m4_pam_v1"
OUT = PROJECT / "results" / "pacer_statemetric_feasibility_v01"
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(smiles)
    return FPGEN.GetFingerprint(mol)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    potency = pd.read_csv(ROOT / "potency_molecules.csv").reset_index(drop=True)
    fps = [fp(x) for x in potency.canonical_smiles]
    positives, negatives = defaultdict(list), defaultdict(list)
    pair_rows = []
    for i in range(len(potency)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
        for offset, similarity in enumerate(sims, start=i + 1):
            j = offset
            same_series = potency.iloc[i].source_component == potency.iloc[j].source_component
            if not same_series or similarity < .50:
                continue
            delta = abs(float(potency.iloc[i].pEC50 - potency.iloc[j].pEC50))
            kind = None
            if delta <= .30:
                kind = "positive"
                positives[i].append(j); positives[j].append(i)
            elif similarity >= .55 and delta >= 1.0:
                kind = "hard_negative"
                negatives[i].append(j); negatives[j].append(i)
            if kind:
                pair_rows.append({"mol_i": potency.iloc[i].canonical_molecule_id,
                                  "mol_j": potency.iloc[j].canonical_molecule_id,
                                  "source_component": potency.iloc[i].source_component,
                                  "tanimoto": similarity, "abs_delta_pEC50": delta,
                                  "pair_role": kind})
    pair_frame = pd.DataFrame(pair_rows)
    anchors = []
    for i, row in potency.iterrows():
        n_pos, n_neg = len(positives[i]), len(negatives[i])
        anchors.append({"canonical_molecule_id": row.canonical_molecule_id,
                        "source_component": row.source_component,
                        "n_positive_partners": n_pos,
                        "n_hard_negative_partners": n_neg,
                        "raw_triplet_combinations": n_pos * n_neg})
    anchor_frame = pd.DataFrame(anchors)
    usable = anchor_frame[(anchor_frame.n_positive_partners > 0) &
                          (anchor_frame.n_hard_negative_partners > 0)]

    classification = pd.read_csv(ROOT / "pam_vs_inactive.csv").reset_index(drop=True)
    pos = classification[classification.target == 1]
    neg = classification[classification.target == 0]
    neg_fps = [fp(x) for x in neg.canonical_smiles]
    cross_pairs = []
    for row in pos.itertuples(index=False):
        query = fp(row.canonical_smiles)
        for j, similarity in enumerate(DataStructs.BulkTanimotoSimilarity(query, neg_fps)):
            if similarity >= .55:
                cross_pairs.append({"pam_id": row.canonical_molecule_id,
                                    "inactive_id": neg.iloc[j].canonical_molecule_id,
                                    "tanimoto": similarity,
                                    "pam_source": row.source_component,
                                    "inactive_source": neg.iloc[j].source_component})
    cross_frame = pd.DataFrame(cross_pairs)

    traj_files = sorted((PROJECT / "data" / "m4_gamd_figshare_33283491").glob(
        "run*-GaMD-stride1000ps.nc"))
    monash = pd.read_csv(ROOT / "external_monash_ly2033298" / "external_allostery.csv")
    mechanistic = pd.read_csv(ROOT / "mechanistic_chembl" / "mechanistic_multidomain.csv")
    us_summary = json.loads((PROJECT / "results" /
        "us20260055116_functional_divergence_v01" / "summary.json").read_text())

    by_series = usable.groupby("source_component").size().rename("usable_anchor_count").reset_index()
    by_series.to_csv(OUT / "usable_anchors_by_series.csv", index=False)
    pair_frame.to_csv(OUT / "potency_metric_pairs.csv", index=False)
    anchor_frame.to_csv(OUT / "anchor_triplet_capacity.csv", index=False)
    cross_frame.to_csv(OUT / "pam_inactive_hard_pairs.csv", index=False)

    audit = {
        "definitions": {
            "positive": "same source_component, ECFP4 Tanimoto >=0.50, abs delta pEC50 <=0.30",
            "potency_hard_negative": "same source_component, ECFP4 Tanimoto >=0.55, abs delta pEC50 >=1.0",
            "pam_inactive_hard_negative": "PAM vs experimental inactive, ECFP4 Tanimoto >=0.55",
        },
        "potency_molecules": len(potency),
        "positive_pairs": int((pair_frame.pair_role == "positive").sum()),
        "potency_hard_negative_pairs": int((pair_frame.pair_role == "hard_negative").sum()),
        "anchors_with_both_partner_types": int(len(usable)),
        "series_with_usable_anchors": int(usable.source_component.nunique()),
        "raw_triplet_combinations": int(usable.raw_triplet_combinations.sum()),
        "pam_inactive_hard_pairs": int(len(cross_frame)),
        "public_gamd_replicas": len(traj_files),
        "public_gamd_frames_known_from_reproduction": 3000,
        "public_gamd_unique_pam_chemotypes": 1,
        "monash_molecules": len(monash),
        "monash_complete_operational_triplets": int(monash[["functional_pKB", "log_tauB", "log_alpha_beta"]].notna().all(axis=1).sum()),
        "monash_explicit_inactive": int(monash.functional_status.eq("inactive").sum()),
        "mechanistic_molecules": len(mechanistic),
        "mechanistic_complete_camp_triplets": int(mechanistic[["pKB", "cAMP_log_alpha_beta", "cAMP_log_tauB"]].notna().all(axis=1).sum()),
        "us_local_pairs": us_summary["n_local_pairs"],
        "us_local_ge2bin_cliffs": us_summary["n_local_pairs_with_ge2_bin_endpoint_cliff"],
    }
    audit["molecular_triplet_go"] = bool(
        audit["anchors_with_both_partner_types"] >= 30 and
        audit["series_with_usable_anchors"] >= 5 and
        audit["pam_inactive_hard_pairs"] >= 10)
    audit["dynamic_statemil_go"] = False
    audit["dynamic_blocker"] = (
        "Six replicas cover one PAM chemotype and no matched ACh-only/inactive/weak-PAM controls; "
        "only unsupervised representation and replica-consistency pilots are currently valid.")
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# PACER-StateMetric data feasibility audit", "", "```json",
              json.dumps(audit, indent=2), "```", "",
              "Molecular triplet combinations are correlated augmentations, not independent molecule count. "
              "Triplets must be constructed after source/series splitting.", "",
              "Dynamic StateMIL remains No-Go for functional training until matched control trajectories and multiple PAM chemotypes exist."]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
