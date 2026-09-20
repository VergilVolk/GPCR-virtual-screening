# -*- coding: utf-8 -*-
"""Freeze an assay-ready, two-round PACER-M4 candidate validation loop.

Selection is label-free and auditable.  It intentionally avoids a scalar
"PAM score": potency, inactive-risk, novelty, and dynamic-pocket evidence are
kept as separate columns.  The first round contains calibration, exploitation,
and exploration molecules plus primary-source functional controls.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


P = Path(__file__).resolve().parents[1]
FINAL = P / "results" / "pacer_candidates_v01" / "final" / "final_candidate_hypotheses.csv"
DYNAMIC = P / "results" / "m4_gamd_ensemble" / "candidate_ensemble_evidence.csv"
CONTROLS = P / "results" / "m4_gamd_ensemble" / "validation_set_primary.csv"
POTENCY = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
OUT = P / "results" / "pacer_assay_closed_loop_v01"
OUT.mkdir(parents=True, exist_ok=True)


def fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def sim(a, b):
    return DataStructs.TanimotoSimilarity(a, b)


def pick_quantile_anchors(df, fps, k=3):
    """Predicted-span selection with a diversity penalty and no labels."""
    values = df.potency_reference_lower90.to_numpy(float)
    targets = np.quantile(values, [.25, .50, .75])
    scale = max(np.ptp(values), .1)
    chosen = []
    for target in targets[:k]:
        best = None
        for i in range(len(df)):
            if i in chosen:
                continue
            distance = abs(values[i] - target) / scale
            redundancy = max([sim(fps[i], fps[j]) for j in chosen], default=0.0)
            score = distance + .35 * redundancy
            key = (score, float(df.iloc[i].strict_inactive_risk_ref), df.iloc[i].candidate_id)
            if best is None or key < best[0]:
                best = (key, i)
        chosen.append(best[1])
    return chosen


def pick_diverse_anchors(indices, fps, k=3):
    """Medoid plus max-min ECFP diversity; label-free and deterministic."""
    chosen = [max(indices, key=lambda i: np.mean([sim(fps[i], fps[j]) for j in indices]))]
    while len(chosen) < k:
        remaining = [i for i in indices if i not in chosen]
        chosen.append(max(remaining, key=lambda i: min(1.0 - sim(fps[i], fps[j]) for j in chosen)))
    return chosen


def maxmin_pick(df, fps, pool, chosen, n, utility):
    picked = []
    for _ in range(n):
        best = None
        for i in pool:
            if i in chosen or i in picked:
                continue
            novelty = 1.0 - max([sim(fps[i], fps[j]) for j in chosen + picked], default=0.0)
            value = float(utility(df.iloc[i])) + .40 * novelty
            key = (value, -float(df.iloc[i].strict_inactive_risk_ref), df.iloc[i].candidate_id)
            if best is None or key > best[0]:
                best = (key, i)
        if best is not None:
            picked.append(best[1])
    return picked


def main():
    final = pd.read_csv(FINAL)
    dynamic = pd.read_csv(DYNAMIC)
    keep = [c for c in dynamic.columns if c not in final.columns or c == "candidate_id"]
    d = final.merge(dynamic[keep], on="candidate_id", how="left")
    potency = pd.read_csv(POTENCY)[["canonical_molecule_id", "source_component"]]
    series_map = dict(zip(potency.canonical_molecule_id, potency.source_component))
    d["reference_series"] = d.nearest_known_id.map(series_map)
    d = d.sort_values("candidate_id").reset_index(drop=True)
    fps = [fp(s) for s in d.canonical_smiles]

    # PACER-FS calibration is only valid inside one coherent medicinal-
    # chemistry series.  Never share an offset across unrelated chemotypes.
    lead_series = d.reference_series.value_counts().index[0]
    lead_pool = [i for i in range(len(d)) if d.iloc[i].reference_series == lead_series]
    anchors = pick_diverse_anchors(lead_pool, fps)
    exploit_pool = [
        i for i in lead_pool if d.iloc[i].portfolio_role == "local_exploitation"
    ]
    exploit = maxmin_pick(
        d, fps, exploit_pool, anchors, 3,
        lambda r: float(r.potency_reference_lower90) - 1.5 * float(r.strict_inactive_risk_ref),
    )
    explore_pool = [i for i in range(len(d)) if i not in lead_pool]
    explore = maxmin_pick(
        d, fps, explore_pool, anchors + exploit, 2,
        lambda r: (1.0 - float(r.max_tanimoto)) + .25 * float(r.pocket_coverage_mean),
    )

    roles = {i: "calibration_anchor" for i in anchors}
    roles.update({i: "exploitation_test" for i in exploit})
    roles.update({i: "exploration_test" for i in explore})
    selected = []
    for order, i in enumerate(anchors + exploit + explore, 1):
        r = d.iloc[i].to_dict()
        r.update({
            "round": 1, "round1_order": order, "assay_role": roles[i],
            "wet_status": "not_tested", "is_control": 0,
            "decision_claim": "computational hypothesis; PAM identity requires functional assay",
        })
        selected.append(r)

    # Two property-matched primary-source functional pairs are assay controls.
    controls = pd.read_csv(CONTROLS)
    for pair in (0, 1):
        for _, c in controls[controls.match_pair == pair].sort_values("target", ascending=False).iterrows():
            selected.append({
                "candidate_id": c.canonical_molecule_id,
                "canonical_smiles": c.canonical_smiles,
                "round": 1,
                "round1_order": len(selected) + 1,
                "assay_role": "positive_control" if int(c.target) == 1 else "inactive_control",
                "wet_status": "historical_primary_source_label_not_retested",
                "is_control": 1,
                "match_pair": int(pair),
                "decision_claim": "historical control; identity and purity must be verified before use",
            })

    round1 = pd.DataFrame(selected)
    core = [
        "round", "round1_order", "candidate_id", "canonical_smiles", "assay_role",
        "is_control", "wet_status", "portfolio_role", "dynamic_pareto",
        "potency_reference_lower90", "potency_reference_upper90",
        "strict_inactive_risk_ref", "max_tanimoto", "cluster00_vina", "BE_avg",
        "vina_sd", "pocket_coverage_mean", "decision_claim",
    ]
    for col in core:
        if col not in round1:
            round1[col] = np.nan
    round1[core].to_csv(OUT / "round1_assay_batch.csv", index=False)

    template = round1[round1.is_control == 0][["candidate_id", "canonical_smiles", "assay_role"]].copy()
    for col in [
        "assay_probe", "assay_readout",
        "compound_purity_pct", "ach_ec20_potentiation_pct", "pam_ec50_nM",
        "max_potentiation_fold", "agonism_alone_pct", "alpha_beta_estimate",
        "functional_pKB", "functional_pKB_ci_low", "functional_pKB_ci_high",
        "log_alpha_beta", "log_alpha_beta_ci_low", "log_alpha_beta_ci_high",
        "log_tauB", "log_tauB_ci_low", "log_tauB_ci_high", "operational_model_qc",
        "m1_counter_pct", "m2_counter_pct", "m3_counter_pct", "m5_counter_pct",
        "cytotoxicity_pct", "replicate_n", "qc_pass", "round2_decision", "notes",
    ]:
        template[col] = ""
    template.to_csv(OUT / "round1_results_template.csv", index=False)

    remaining = d[~d.candidate_id.isin(round1.candidate_id)].copy()
    remaining.insert(0, "round2_priority_status", "await_round1_anchor_calibration")
    remaining.to_csv(OUT / "round2_candidate_pool.csv", index=False)

    audit = {
        "method": "diverse three-anchor calibration plus diversity-constrained exploitation/exploration",
        "n_candidates_round1": 8,
        "n_controls_round1": 4,
        "n_remaining_round2": int(len(remaining)),
        "lead_reference_series": lead_series,
        "lead_series_size": len(lead_pool),
        "anchor_ids": [d.iloc[i].candidate_id for i in anchors],
        "exploitation_ids": [d.iloc[i].candidate_id for i in exploit],
        "exploration_ids": [d.iloc[i].candidate_id for i in explore],
        "claim_boundary": (
            "This is an assay-ready selection, not wet validation or confirmed PAM discovery. "
            "The three-anchor offset applies only to the lead reference series; other chemotypes remain uncalibrated."
        ),
    }
    (OUT / "selection_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
