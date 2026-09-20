# -*- coding: utf-8 -*-
"""Series-controlled residue/contact association analysis for M4 PAM potency."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.stats import rankdata, spearmanr


P = Path(__file__).resolve().parents[1]
DATA = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
STATE_FILES = {
    "7TRQ": P / "results" / "pacer_structure_loso_v01" / "docking_features.csv",
    "7TRP": P / "results" / "pacer_structure_7trp_v01" / "docking_features.csv",
    "7TRS": P / "results" / "pacer_structure_7trs_v01" / "docking_features.csv",
}
OUT = P / "results" / "m4_residue_potency_association_v01"
OUT.mkdir(parents=True, exist_ok=True)
RESIDUES = (89, 92, 93, 96, 184, 186, 190, 423, 432, 433, 435, 436, 439)
N_PERM = 5000
N_BOOT = 20000
SEED = 830


def load_state(name, path):
    frame = pd.read_csv(path)
    frame = frame[frame.error.fillna("").eq("")].drop_duplicates("canonical_molecule_id")
    keep = ["canonical_molecule_id", "n_heavy_atoms"]
    for residue in RESIDUES:
        keep += [f"res{residue}_min_A", f"res{residue}_pairs"]
    frame = frame[keep].copy()
    out = frame[["canonical_molecule_id"]].copy()
    for residue in RESIDUES:
        out[f"{name}|res{residue}|min_A"] = frame[f"res{residue}_min_A"]
        out[f"{name}|res{residue}|pairs_per_heavy"] = (
            frame[f"res{residue}_pairs"] / frame.n_heavy_atoms.clip(lower=1)
        )
    return out


def bh_fdr(pvalues):
    p = np.asarray(pvalues, float); n = len(p); order = np.argsort(p)
    q = np.empty(n, float); ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    q[order] = ranked
    return q


def zrank(x):
    r = rankdata(x, method="average").astype(float)
    return (r - r.mean()) / r.std(ddof=1)


def main():
    potency = pd.read_csv(DATA)
    states = {name: load_state(name, path) for name, path in STATE_FILES.items()}
    merged = potency[["canonical_molecule_id", "canonical_smiles", "pEC50", "source_component"]].copy()
    for frame in states.values():
        merged = merged.merge(frame, on="canonical_molecule_id", how="inner")
    merged["pEC50_centered"] = merged.pEC50 - merged.groupby("source_component").pEC50.transform("mean")

    # State-response contrasts use identical residue feature definitions.
    for residue in RESIDUES:
        for kind in ("min_A", "pairs_per_heavy"):
            merged[f"AChContrast|res{residue}|{kind}"] = (
                merged[f"7TRS|res{residue}|{kind}"]
                - .5 * (merged[f"7TRQ|res{residue}|{kind}"] + merged[f"7TRP|res{residue}|{kind}"])
            )

    features = [c for c in merged if "|res" in c]
    valid = [c for c in features if merged[c].nunique() > 2 and merged[c].notna().all()]
    xrank = np.column_stack([zrank(merged[c].to_numpy(float)) for c in valid])
    yrank = zrank(merged.pEC50_centered.to_numpy(float))
    observed = yrank @ xrank / (len(merged) - 1)

    descriptor_functions = [
        Descriptors.MolWt, Descriptors.MolLogP, Descriptors.NumHDonors,
        Descriptors.NumHAcceptors, Descriptors.TPSA, Descriptors.NumRotatableBonds,
        Descriptors.NumAromaticRings, Descriptors.FractionCSP3, Descriptors.HeavyAtomCount,
    ]
    descriptor_matrix = np.asarray([
        [fn(Chem.MolFromSmiles(s)) for fn in descriptor_functions]
        for s in merged.canonical_smiles
    ], float)
    descriptor_matrix = (descriptor_matrix - descriptor_matrix.mean(0)) / np.where(
        descriptor_matrix.std(0) > 1e-8, descriptor_matrix.std(0), 1
    )
    group_dummy = pd.get_dummies(merged.source_component, drop_first=False).to_numpy(float)
    nuisance = np.c_[np.ones(len(merged)), group_dummy, descriptor_matrix]
    raw_x = merged[valid].to_numpy(float)
    adjusted_x = raw_x - nuisance @ np.linalg.lstsq(nuisance, raw_x, rcond=None)[0]
    adjusted_xrank = np.column_stack([zrank(adjusted_x[:, j]) for j in range(adjusted_x.shape[1])])
    adjusted_observed = yrank @ adjusted_xrank / (len(merged) - 1)

    rng = np.random.default_rng(SEED)
    group_indices = [np.flatnonzero(merged.source_component.to_numpy() == g)
                     for g in sorted(merged.source_component.unique())]
    exceed = np.zeros(len(valid), int)
    adjusted_exceed = np.zeros(len(valid), int)
    for _ in range(N_PERM):
        perm = np.arange(len(merged))
        for idx in group_indices:
            perm[idx] = rng.permutation(idx)
        stat = yrank[perm] @ xrank / (len(merged) - 1)
        adjusted_stat = yrank[perm] @ adjusted_xrank / (len(merged) - 1)
        exceed += np.abs(stat) >= np.abs(observed)
        adjusted_exceed += np.abs(adjusted_stat) >= np.abs(adjusted_observed)
    pvalue = (exceed + 1) / (N_PERM + 1)
    qvalue = bh_fdr(pvalue)
    adjusted_pvalue = (adjusted_exceed + 1) / (N_PERM + 1)
    adjusted_qvalue = bh_fdr(adjusted_pvalue)

    rows = []
    eval_groups = [g for g, n in merged.source_component.value_counts().items() if n >= 12]
    for j, feature in enumerate(valid):
        per_group = []
        for group in eval_groups:
            d = merged[merged.source_component == group]
            rho = spearmanr(d.pEC50, d[feature]).statistic
            per_group.append(float(rho) if np.isfinite(rho) else 0.0)
        arr = np.asarray(per_group)
        boot = arr[rng.integers(0, len(arr), size=(N_BOOT, len(arr)))].mean(1)
        rows.append({
            "feature": feature, "pooled_centered_Spearman": float(observed[j]),
            "within_series_permutation_p": float(pvalue[j]), "BH_q": float(qvalue[j]),
            "physchem_adjusted_Spearman": float(adjusted_observed[j]),
            "physchem_adjusted_permutation_p": float(adjusted_pvalue[j]),
            "physchem_adjusted_BH_q": float(adjusted_qvalue[j]),
            "macro_series_Spearman": float(arr.mean()),
            "macro_series_CI_low": float(np.quantile(boot, .025)),
            "macro_series_CI_high": float(np.quantile(boot, .975)),
            "positive_series": int((arr > 0).sum()), "negative_series": int((arr < 0).sum()),
            "n_evaluable_series": len(arr),
        })
    result = pd.DataFrame(rows).sort_values(["BH_q", "within_series_permutation_p"])
    result.to_csv(OUT / "residue_associations.csv", index=False)
    merged.to_csv(OUT / "analysis_matrix.csv", index=False)

    significant = result[(result.BH_q < .05) &
                         (result[["positive_series", "negative_series"]].max(axis=1) >= 7)]
    adjusted_significant = result[result.physchem_adjusted_BH_q < .05]
    top = result.head(15)
    audit = {
        "n_molecules_complete_three_state": int(len(merged)),
        "n_series_total": int(merged.source_component.nunique()),
        "n_evaluable_series_ge12": int(len(eval_groups)),
        "n_hypotheses": int(len(result)), "within_series_permutations": N_PERM,
        "significant_consistent_features": int(len(significant)),
        "physchem_adjusted_significant_features": int(len(adjusted_significant)),
        "preregistered_conclusion": (
            "cross_series_static_residue_association_supported" if len(significant)
            else "no_universal_single_residue_static_potency_marker"
        ),
        "claim_boundary": "association analysis; no causal or PAM-identity claim",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = "# M4 residue-potency association result\n\n"
    report += f"Complete three-state molecules: {len(merged)}; tested hypotheses: {len(result)}; "
    report += f"FDR-significant and direction-consistent features: {len(significant)}.\n\n"
    report += "## Top preregistered tests\n\n"
    report += top.to_markdown(index=False, floatfmt=".4f") + "\n\n"
    if len(significant):
        report += "The listed significant features support a cross-series static residue association, not causality.\n"
    else:
        report += (
            "No single residue/contact feature survived the preregistered FDR and direction-consistency gate. "
            "The supported conclusion is that M4 PAM potency is not encoded by one universal static contact; "
            "distributed, chemotype-dependent and dynamic coupling remains the testable mechanism.\n"
        )
    (OUT / "VALIDATION_REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(audit, indent=2)); print(top.to_string(index=False))


if __name__ == "__main__":
    main()
