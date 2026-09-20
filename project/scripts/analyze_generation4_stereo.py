# -*- coding: utf-8 -*-
"""Compare 2D/chiral QSAR and GaMD ensemble on Generation-4 stereoisomers."""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from scipy.stats import spearmanr

from analyze_m4_gamd_ensemble import aggregate
from run_pacer_fs_baselines import DESC, features


P = Path(__file__).resolve().parents[1]
EXT = P / "data" / "benchmarks" / "m4_pam_v1" / "external_2026_vu6025733"
ROOT = P / "results" / "m4_gamd_ensemble"
DOCK = ROOT / "generation4_stereo_docking_protocol_v127_ph7_rank1" / "framewise_scores.csv"
OUT = P / "results" / "pacer_generation4_stereo_v01"


def chiral_features(smiles):
    fps, desc = [], []
    for value in smiles:
        mol = Chem.MolFromSmiles(value)
        bit = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048, useChirality=True)
        arr = np.zeros(2048, np.float32); DataStructs.ConvertToNumpyArray(bit, arr)
        fps.append(arr); desc.append([f(mol) for f in DESC])
    return np.asarray(fps, np.float32), np.asarray(desc, np.float32)


def fit_predict(X, D, y, n_train, seed):
    train = np.arange(n_train); test = np.arange(n_train, len(X))
    selected = np.argsort(X[train].var(0))[-1024:]
    xx = np.c_[X[:, selected], D]
    model = LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1,
        random_state=seed, verbosity=-1, n_jobs=6,
    )
    return model.fit(xx[train], y).predict(xx[test])


def rho(y, p):
    value = spearmanr(y, p).statistic
    return float(value) if np.isfinite(value) else 0.0


def exact_permutation_p(y, p):
    observed = abs(rho(y, p)); values = []
    for order in itertools.permutations(range(len(p))):
        values.append(abs(rho(y, p[np.asarray(order)])))
    return float(np.mean(np.asarray(values) >= observed - 1e-12))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    hist = pd.read_csv(P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv")
    table1 = pd.read_csv(EXT / "external_potency.csv")
    gen3 = pd.read_csv(EXT / "generation3_potency.csv")
    gen3_exact = gen3[gen3.censoring.eq("exact")]
    gen4 = pd.read_csv(EXT / "generation4_stereo_potency.csv")
    train = pd.concat([
        hist[["canonical_smiles", "pEC50"]], table1[["canonical_smiles", "pEC50"]],
        gen3_exact[["canonical_smiles", "pEC50"]],
    ], ignore_index=True)
    all_smiles = pd.concat([train.canonical_smiles, gen4.canonical_smiles], ignore_index=True)
    X2, D2, _ = features(all_smiles)
    Xc, Dc = chiral_features(all_smiles)
    pred_2d = fit_predict(X2, D2, train.pEC50.to_numpy(float), len(train), 201)
    pred_chiral = fit_predict(Xc, Dc, train.pEC50.to_numpy(float), len(train), 202)

    frame = pd.read_csv(DOCK)
    frame = frame[(frame.error.isna() | frame.error.eq("")) & frame.vina_affinity.notna()]
    dock = aggregate(frame).rename(columns={"molecule_id": "compound_id"})
    merged = gen4.merge(dock, on="compound_id", how="left", validate="one_to_one")
    merged["Cumulative-2D-QSAR"] = pred_2d
    merged["Cumulative-Chiral-QSAR"] = pred_chiral
    merged.to_csv(OUT / "predictions.csv", index=False)

    primary = merged[merged.primary_stereochemical_query.astype(bool)].copy()
    methods = {
        "Cumulative-2D-QSAR": primary["Cumulative-2D-QSAR"].to_numpy(float),
        "Cumulative-Chiral-QSAR": primary["Cumulative-Chiral-QSAR"].to_numpy(float),
        "Single-cluster00": -primary.cluster00_vina.to_numpy(float),
        "GaMD-BEmin": -primary.BE_min.to_numpy(float),
        "GaMD-BEavg": -primary.BE_avg.to_numpy(float),
        "GaMD-vina-mean": -primary.vina_mean.to_numpy(float),
    }
    y = primary.pEC50.to_numpy(float)
    rows = []
    for name, pred in methods.items():
        rows.append({
            "method": name, "n": len(y), "Spearman": rho(y, pred),
            "exact_permutation_p_two_sided": exact_permutation_p(y, pred),
            "prediction_range": float(np.ptp(pred)),
        })
    metrics = pd.DataFrame(rows).sort_values("Spearman", ascending=False)
    metrics.to_csv(OUT / "metrics.csv", index=False)

    isotope_pairs = [("33h", "33i"), ("33j", "33k"), ("33l", "33m"),
                     ("33n", "33o"), ("33p", "33q")]
    pair_rows = []
    lookup = merged.set_index("compound_id")
    for a, b in isotope_pairs:
        for score in ["Cumulative-2D-QSAR", "Cumulative-Chiral-QSAR", "BE_avg", "BE_min"]:
            pa, pb = float(lookup.loc[a, score]), float(lookup.loc[b, score])
            if score.startswith("BE_"):
                pa, pb = -pa, -pb
            true_delta = float(lookup.loc[b, "pEC50"] - lookup.loc[a, "pEC50"])
            pair_rows.append({"pair": f"{a}->{b}", "method": score,
                              "true_delta_pEC50": true_delta, "predicted_delta": pb - pa,
                              "direction_correct": bool(np.sign(true_delta) == np.sign(pb - pa))})
    pd.DataFrame(pair_rows).to_csv(OUT / "isotope_pair_audit.csv", index=False)

    audit = {
        "status": "generation4_stereochemical_mechanism_stress_test",
        "n_primary_stereoisomers": len(primary),
        "docking_complete": bool(primary.complete_ensemble.all()),
        "best_method_by_spearman": str(metrics.iloc[0].method),
        "primary_method_GaMD_BEavg_spearman": float(metrics.set_index("method").loc["GaMD-BEavg", "Spearman"]),
        "claim_boundary": "N=4 mechanism stress test; no SOTA or functional PAM claim.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = [
        "# Generation-4 stereochemical stress test", "",
        "Four non-deuterated 3-fluoropiperidine stereoisomers share one 2D graph.", "",
        "| Method | Spearman | exact permutation p | score range |", "|---|---:|---:|---:|",
    ]
    for row in metrics.itertuples():
        report.append(f"| {row.method} | {row.Spearman:.3f} | {row.exact_permutation_p_two_sided:.3f} | {row.prediction_range:.3f} |")
    report += ["", "Boundary: N=4; docking scores measure pose/affinity evidence, not PAM efficacy."]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
