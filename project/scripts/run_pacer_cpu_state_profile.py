# -*- coding: utf-8 -*-
"""CPU-only M4 PAM benchmark using receptor-state preference profiles.

This script deliberately treats ensemble docking as a binding/pose feature source,
not as direct evidence of PAM efficacy.  Validation is leave-one-matched-pair-out:
the positive and physicochemically matched negative of each pair are never seen in
training.  Chemical, state-profile, and late-fusion models use the same splits.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "results" / "m4_gamd_ensemble"
STATEWISE = ROOT / "validation_primary_docking_protocol_v127_ph7_rank1" / "statewise_scores.csv"
LABELS = ROOT / "validation_set_primary.csv"
ENSEMBLE = ROOT / "validation_ensemble_scores.csv"
OUT = PROJECT / "results" / "pacer_cpu_state_profile_v01"
N_CLUSTERS = 10
C_GRID = (0.01, 0.1, 1.0, 10.0)


def parse_contacts(value: object) -> set[int]:
    if pd.isna(value) or str(value).strip() == "":
        return set()
    return {int(x) for x in str(value).split(";") if x.strip()}


def softmax(x: np.ndarray) -> np.ndarray:
    z = x - np.max(x)
    e = np.exp(z)
    return e / e.sum()


def state_features(state: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    rows: list[dict[str, float | str]] = []
    feature_names: list[str] | None = None
    for molecule_id, part in state.groupby("molecule_id", sort=True):
        part = part.sort_values("cluster")
        if len(part) != N_CLUSTERS or part.cluster.nunique() != N_CLUSTERS:
            continue
        aff = part.vina_affinity.to_numpy(float)
        cov = part.native_pocket_coverage.to_numpy(float)
        ncon = part.n_contacted_residues.to_numpy(float)
        contacts = [parse_contacts(x) for x in part.contacted_residues]
        base = contacts[0]
        jaccard = np.asarray([
            len(base & x) / len(base | x) if (base | x) else 1.0
            for x in contacts[1:]
        ])
        # Relative profiles suppress the global ligand-size/affinity offset.
        values: dict[str, float | str] = {"molecule_id": molecule_id}
        values.update({f"aff_delta_c{i}": float(aff[i] - aff[0]) for i in range(1, N_CLUSTERS)})
        values.update({f"coverage_delta_c{i}": float(cov[i] - cov[0]) for i in range(1, N_CLUSTERS)})
        values.update({f"contact_delta_c{i}": float(ncon[i] - ncon[0]) for i in range(1, N_CLUSTERS)})
        values.update({f"contact_jaccard_c{i}": float(jaccard[i - 1]) for i in range(1, N_CLUSTERS)})
        p = softmax(-aff)  # fixed 1 kcal/mol scale; descriptive, not fitted.
        contact_union = set().union(*contacts)
        contact_core = set.intersection(*contacts) if contacts else set()
        values.update({
            "affinity_sd": float(np.std(aff)),
            "affinity_range": float(np.ptp(aff)),
            "state_entropy": float(-(p * np.log(p + 1e-12)).sum()),
            "state_max_probability": float(p.max()),
            "cluster0_probability": float(p[0]),
            "coverage_sd": float(np.std(cov)),
            "contact_count_sd": float(np.std(ncon)),
            "contact_union_size": float(len(contact_union)),
            "contact_core_size": float(len(contact_core)),
            "contact_jaccard_mean": float(jaccard.mean()),
        })
        if feature_names is None:
            feature_names = [k for k in values if k != "molecule_id"]
        rows.append(values)
    return pd.DataFrame(rows), feature_names or []


def ecfp(smiles: pd.Series, n_bits: int = 1024) -> np.ndarray:
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
    out = np.zeros((len(smiles), n_bits), dtype=np.float32)
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smi}")
        out[i] = np.asarray(gen.GetFingerprintAsNumPy(mol), dtype=np.float32)
    return out


def balanced_pair_folds(groups: np.ndarray, held_out: int) -> tuple[np.ndarray, np.ndarray]:
    test = groups == held_out
    return np.flatnonzero(~test), np.flatnonzero(test)


def choose_c(x: np.ndarray, y: np.ndarray, groups: np.ndarray) -> float:
    """Inner leave-one-pair-out selection using pairwise ranking accuracy."""
    unique = np.unique(groups)
    best = (-np.inf, C_GRID[0])
    for c in C_GRID:
        correct: list[float] = []
        for g in unique:
            tr, te = balanced_pair_folds(groups, g)
            model = make_pipeline(
                StandardScaler(),
                LogisticRegression(C=c, penalty="l2", solver="liblinear", max_iter=3000),
            )
            model.fit(x[tr], y[tr])
            pred = model.predict_proba(x[te])[:, 1]
            pos = pred[y[te] == 1]
            neg = pred[y[te] == 0]
            if len(pos) == 1 and len(neg) == 1:
                correct.append(float(pos[0] > neg[0]) + 0.5 * float(pos[0] == neg[0]))
        score = float(np.mean(correct)) if correct else -np.inf
        # Prefer stronger regularization when tied.
        if score > best[0] + 1e-12:
            best = (score, c)
    return float(best[1])


def nested_pair_oof(x: np.ndarray, y: np.ndarray, groups: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    pred = np.full(len(y), np.nan, float)
    folds: list[dict] = []
    for g in np.unique(groups):
        tr, te = balanced_pair_folds(groups, g)
        c = choose_c(x[tr], y[tr], groups[tr])
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(C=c, penalty="l2", solver="liblinear", max_iter=3000),
        )
        model.fit(x[tr], y[tr])
        pred[te] = model.predict_proba(x[te])[:, 1]
        folds.append({"held_out_pair": int(g), "selected_C": c, "n_train": int(len(tr)), "n_test": int(len(te))})
    if np.isnan(pred).any():
        raise RuntimeError("OOF prediction incomplete")
    return pred, folds


def ef10(y: np.ndarray, score: np.ndarray) -> float:
    k = max(1, int(math.ceil(0.1 * len(y))))
    hits = y[np.argsort(-score)[:k]].sum()
    return float((hits / k) / y.mean())


def binomial_two_sided(k: int, n: int) -> float:
    # Exact two-sided sign-test p-value under p=0.5.
    probs = [math.comb(n, i) / (2**n) for i in range(n + 1)]
    pk = probs[k]
    return float(min(1.0, sum(p for p in probs if p <= pk + 1e-15)))


def evaluate(y: np.ndarray, score: np.ndarray, groups: np.ndarray) -> dict:
    margins: list[float] = []
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        pos = score[idx][y[idx] == 1]
        neg = score[idx][y[idx] == 0]
        if len(pos) == 1 and len(neg) == 1:
            margins.append(float(pos[0] - neg[0]))
    arr = np.asarray(margins)
    wins = int((arr > 0).sum())
    ties = int((arr == 0).sum())
    return {
        "roc_auc_oof": float(roc_auc_score(y, score)),
        "pr_auc_oof": float(average_precision_score(y, score)),
        "ef10_oof": ef10(y, score),
        "matched_pair_accuracy": float(((arr > 0).sum() + 0.5 * (arr == 0).sum()) / len(arr)),
        "matched_pair_wins": wins,
        "matched_pair_ties": ties,
        "matched_pair_n": int(len(arr)),
        "sign_test_p_two_sided": binomial_two_sided(wins, len(arr) - ties) if len(arr) > ties else 1.0,
        "median_pair_margin": float(np.median(arr)),
    }


def paired_bootstrap_delta(
    y: np.ndarray,
    groups: np.ndarray,
    score_a: np.ndarray,
    score_b: np.ndarray,
    n_boot: int = 10000,
    seed: int = 20260920,
) -> dict:
    """Bootstrap matched pairs; compare within-pair margin B minus A."""
    margins_a, margins_b = [], []
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        margins_a.append(float(score_a[idx][y[idx] == 1][0] - score_a[idx][y[idx] == 0][0]))
        margins_b.append(float(score_b[idx][y[idx] == 1][0] - score_b[idx][y[idx] == 0][0]))
    a, b = np.asarray(margins_a), np.asarray(margins_b)
    rng = np.random.default_rng(seed)
    deltas = np.empty(n_boot)
    for i in range(n_boot):
        take = rng.integers(0, len(a), len(a))
        deltas[i] = np.mean(b[take] - a[take])
    return {
        "mean_margin_delta": float(np.mean(b - a)),
        "ci95": [float(x) for x in np.quantile(deltas, [0.025, 0.975])],
        "p_delta_gt_0": float(np.mean(deltas > 0)),
        "n_pairs": int(len(a)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    state = pd.read_csv(STATEWISE)
    state = state[(state.error.isna()) | state.error.eq("")].dropna(subset=["vina_affinity"])
    feat, state_cols = state_features(state)
    labels = pd.read_csv(LABELS)
    ensemble = pd.read_csv(ENSEMBLE)[["molecule_id", "cluster00_vina", "BE_min", "BE_avg"]]
    data = labels.merge(feat, left_on="canonical_molecule_id", right_on="molecule_id", validate="one_to_one")
    data = data.merge(ensemble, on="molecule_id", validate="one_to_one")
    data = data.sort_values(["match_pair", "target"], ascending=[True, False]).reset_index(drop=True)

    y = data.target.to_numpy(int)
    groups = data.match_pair.to_numpy(int)
    x_state = data[state_cols].to_numpy(float)
    x_chem = ecfp(data.canonical_smiles)
    x_fusion = np.concatenate([x_chem, x_state], axis=1)

    predictions: dict[str, np.ndarray] = {
        "single_cluster_vina": -data.cluster00_vina.to_numpy(float),
        "gamd_BE_min": -data.BE_min.to_numpy(float),
        "gamd_BE_avg": -data.BE_avg.to_numpy(float),
    }
    fold_details: dict[str, list[dict]] = {}
    for name, x in [("ecfp4", x_chem), ("state_profile", x_state), ("late_fusion", x_fusion)]:
        predictions[name], fold_details[name] = nested_pair_oof(x, y, groups)

    metrics = {name: evaluate(y, score, groups) for name, score in predictions.items()}
    comparisons = {
        "state_vs_single_vina": paired_bootstrap_delta(y, groups, predictions["single_cluster_vina"], predictions["state_profile"]),
        "fusion_vs_ecfp4": paired_bootstrap_delta(y, groups, predictions["ecfp4"], predictions["late_fusion"]),
        "fusion_vs_state": paired_bootstrap_delta(y, groups, predictions["state_profile"], predictions["late_fusion"]),
    }
    pred_df = data[["molecule_id", "canonical_smiles", "target", "match_pair", "source_component", "murcko_scaffold"]].copy()
    for name, score in predictions.items():
        pred_df[name] = score
    pred_df.to_csv(OUT / "oof_predictions.csv", index=False)
    pd.DataFrame(fold_details["state_profile"]).to_csv(OUT / "state_profile_fold_audit.csv", index=False)
    pd.DataFrame(fold_details["ecfp4"]).to_csv(OUT / "ecfp4_fold_audit.csv", index=False)
    pd.DataFrame(fold_details["late_fusion"]).to_csv(OUT / "late_fusion_fold_audit.csv", index=False)

    audit = {
        "protocol": "leave-one-matched-pair-out with inner pair-wise C selection",
        "n_molecules": int(len(data)),
        "n_pairs": int(len(np.unique(groups))),
        "n_positive": int(y.sum()),
        "n_negative": int((1 - y).sum()),
        "n_state_features": int(len(state_cols)),
        "state_feature_names": state_cols,
        "metrics": metrics,
        "comparisons": comparisons,
        "claim_boundary": (
            "Retrospective matched-pair benchmark. State profiles come from one known-PAM GaMD ensemble; "
            "results do not establish prospective PAM efficacy or replace functional experiments."
        ),
    }
    (OUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# PACER CPU state-profile benchmark", "",
        f"- Data: {len(data)} molecules, {int(y.sum())} PAM / {int((1-y).sum())} inactive, {len(np.unique(groups))} matched pairs.",
        "- Validation: leave-one-matched-pair-out; each positive and its matched negative are held out together.",
        "- Ensemble docking is treated as a binding/state feature source, not a PAM-efficacy oracle.", "",
        "| Method | OOF ROC-AUC | OOF PR-AUC | EF10% | Pair accuracy | Sign-test p |", "|---|---:|---:|---:|---:|---:|",
    ]
    labels_map = {
        "single_cluster_vina": "Single-state Vina",
        "gamd_BE_min": "GaMD BEmin",
        "gamd_BE_avg": "GaMD BEavg",
        "ecfp4": "ECFP4 logistic",
        "state_profile": "CPU state-profile",
        "late_fusion": "ECFP4 + state-profile",
    }
    for key in labels_map:
        m = metrics[key]
        lines.append(
            f"| {labels_map[key]} | {m['roc_auc_oof']:.3f} | {m['pr_auc_oof']:.3f} | {m['ef10_oof']:.2f} | "
            f"{m['matched_pair_accuracy']:.3f} ({m['matched_pair_wins']}/{m['matched_pair_n']}) | {m['sign_test_p_two_sided']:.3g} |"
        )
    lines += ["", "## Predeclared comparisons", ""]
    for key, value in comparisons.items():
        lines.append(
            f"- {key}: mean matched-pair margin delta {value['mean_margin_delta']:+.4f}; "
            f"95% CI [{value['ci95'][0]:+.4f}, {value['ci95'][1]:+.4f}]; P(delta>0)={value['p_delta_gt_0']:.3f}."
        )
    lines += ["", "## Claim boundary", "", audit["claim_boundary"], ""]
    (OUT / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
