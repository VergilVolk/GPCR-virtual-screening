# -*- coding: utf-8 -*-
"""PACER-XR: cross-target, M4-blind rank fusion of the official benchmark.

Glide/Vina and PDB/BEmin/BEavg score ranks are fused by logistic regression.
The regularisation strength is selected only by leave-one-target-out validation
on M2R, CCR2 and B2AR.  M4R labels are never used until the final evaluation.
Exact M4 SMILES are removed from development training as an extra leakage guard.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score


P = Path(__file__).resolve().parents[1]
ROOT = P / "tools" / "gpcr-am-ensemble-docking" / "docking_scores"
OUT = P / "results" / "pacer_xr_v01"
OUT.mkdir(parents=True, exist_ok=True)
DEV_TARGETS = ("M2R", "CCR2", "B2AR")
TEST_TARGET = "M4R"
BASE_FEATURES = (
    "Glide_PDB", "Glide_BEmin", "Glide_BEavg",
    "Vina_PDB", "Vina_BEmin", "Vina_BEavg",
)


def load_target(target):
    files = {
        "Glide_PDB": (ROOT / target / "PDB" / f"{target}_PDB_Glide_scores.csv", "glide_gscore"),
        "Glide_BEmin": (ROOT / target / "Ensemble" / f"{target}_Ensemble_Glide_BEmin_ranked.csv", "BE_min"),
        "Glide_BEavg": (ROOT / target / "Ensemble" / f"{target}_Ensemble_Glide_BEavg_ranked.csv", "BE_avg"),
        "Vina_PDB": (ROOT / target / "PDB" / f"{target}_PDB_Vina_scores.csv", "vina_score"),
        "Vina_BEmin": (ROOT / target / "Ensemble" / f"{target}_Ensemble_Vina_BEmin_ranked.csv", "BE_min"),
        "Vina_BEavg": (ROOT / target / "Ensemble" / f"{target}_Ensemble_Vina_BEavg_ranked.csv", "BE_avg"),
    }
    merged = None
    for name, (path, score_col) in files.items():
        d = pd.read_csv(path, usecols=["ligand_id", "SMILES", score_col]).drop_duplicates("ligand_id")
        d = d.rename(columns={score_col: name, "SMILES": f"SMILES_{name}"})
        merged = d if merged is None else merged.merge(d, on="ligand_id", how="inner")
    merged["SMILES"] = merged["SMILES_Glide_PDB"].astype(str)
    merged["target_name"] = target
    merged["target"] = merged.ligand_id.astype(str).str.startswith("ASD").astype(int)
    # Rank percentile: one is best, zero is worst.  This removes engine scale.
    for col in BASE_FEATURES:
        merged[col] = 1.0 - merged[col].rank(method="average", pct=True, ascending=True)
    merged["consensus_mean"] = merged[list(BASE_FEATURES)].mean(axis=1)
    merged["consensus_min"] = merged[list(BASE_FEATURES)].min(axis=1)
    merged["consensus_max"] = merged[list(BASE_FEATURES)].max(axis=1)
    merged["consensus_sd"] = merged[list(BASE_FEATURES)].std(axis=1)
    merged["engine_gap"] = (
        merged[[x for x in BASE_FEATURES if x.startswith("Glide")]].mean(axis=1)
        - merged[[x for x in BASE_FEATURES if x.startswith("Vina")]].mean(axis=1)
    )
    return merged


FEATURES = BASE_FEATURES + (
    "consensus_mean", "consensus_min", "consensus_max", "consensus_sd", "engine_gap",
)


def balanced(df, seed):
    pos = df[df.target == 1]
    neg = df[df.target == 0].sample(n=min(len(pos), (df.target == 0).sum()), random_state=seed)
    return pd.concat([pos, neg], ignore_index=True)


def fit_model(frames, c, seed):
    train = pd.concat([balanced(x, seed + i) for i, x in enumerate(frames)], ignore_index=True)
    model = LogisticRegression(C=c, max_iter=2000, class_weight=None, random_state=seed)
    model.fit(train[list(FEATURES)], train.target)
    return model, len(train)


def ef(y, score, fraction):
    n = max(1, int(np.ceil(len(y) * fraction)))
    order = np.argsort(-score)[:n]
    return float(y[order].mean() / y.mean())


def ef_prime(y, score, fraction):
    n = max(1, int(np.ceil(len(y) * fraction)))
    order = np.argsort(-score)
    top = order[:n]
    hit_positions = np.where(y[top] == 1)[0] + 1
    if not len(hit_positions):
        return 0.0
    apr = float(np.mean(hit_positions / len(y) * 100.0))
    return float((50.0 / apr) * (len(hit_positions) / y.sum()))


def official_logauc(y, score):
    """Match the released auc_calculator.py definition (reported as percent)."""
    order = np.argsort(-score)
    yy = y[order]
    pct_a = np.r_[0, np.cumsum(yy)] / yy.sum() * 100
    pct_d = np.r_[0, np.cumsum(1 - yy)] / (len(yy) - yy.sum()) * 100
    points = [[d / 100, a / 100] for d, a in zip(pct_d, pct_a) if .001 <= d / 100 <= 1.0]
    i = next((i for i, p in enumerate(points) if p[0] >= .001), None)
    if i:
        p, q = points[i], points[i - 1]
        slope = (p[1] - q[1]) / (p[0] - q[0])
        points.insert(i, [.00100001, slope * .00100001 + p[1] - slope * p[0]])
    area = sum(
        (p[1] - q[1]) / np.log(10)
        + (p[1] - (p[1] - q[1]) / (p[0] - q[0]) * p[0]) * (np.log10(p[0]) - np.log10(q[0]))
        for p, q in zip(points[1:], points[:-1]) if p[0] - q[0] > 1e-6
    )
    random = (1.0 - .001) / np.log(10) / np.log10(1000)
    return float((area / np.log10(1000) - random) * 100)


def auc_diff_ci(y, new, base):
    """Paired asymptotic AUC-difference CI using influence values."""
    def influence(score):
        # Keep the original positive/negative order so influence vectors from
        # two methods remain paired by molecule.  Sorting both vectors here
        # would destroy covariance and severely understate the delta CI.
        pos = score[y == 1]; neg = score[y == 0]
        neg_sorted = np.sort(neg); pos_sorted = np.sort(pos)
        vlo = np.searchsorted(neg_sorted, pos, side="left")
        vhi = np.searchsorted(neg_sorted, pos, side="right")
        v10 = (vlo + .5 * (vhi - vlo)) / len(neg)
        plo = np.searchsorted(pos_sorted, neg, side="left")
        phi = np.searchsorted(pos_sorted, neg, side="right")
        v01 = (len(pos) - phi + .5 * (phi - plo)) / len(pos)
        return v10, v01
    a10, a01 = influence(new); b10, b01 = influence(base)
    delta = float(a10.mean() - b10.mean())
    se = float(np.sqrt(np.var(a10 - b10, ddof=1) / len(a10) + np.var(a01 - b01, ddof=1) / len(a01)))
    return delta, [delta - 1.96 * se, delta + 1.96 * se]


def metrics(y, score):
    return {
        "ROC_AUC": float(roc_auc_score(y, score)),
        "PR_AUC": float(average_precision_score(y, score)),
        "logAUC_pct": official_logauc(y, score),
        "EF_0.5pct": ef(y, score, .005),
        "EFprime_0.5pct": ef_prime(y, score, .005),
        "EF_1pct": ef(y, score, .01),
        "EFprime_1pct": ef_prime(y, score, .01),
        "EF_2pct": ef(y, score, .02),
        "EF_5pct": ef(y, score, .05),
    }


def cascade_score(early_score, global_score, fraction=.01):
    """Preserve a fixed early-ranking list, then use global consensus.

    The 1% boundary is fixed from the paper's primary early-enrichment regime,
    not selected with M4 labels.
    """
    n = len(global_score)
    k = int(np.ceil(n * fraction))
    early_order = np.argsort(-early_score)
    head = early_order[:k]
    head_set = np.zeros(n, dtype=bool); head_set[head] = True
    tail = np.where(~head_set)[0]
    tail = tail[np.argsort(-global_score[tail])]
    order = np.r_[head, tail]
    score = np.empty(n, float)
    score[order] = np.arange(n, 0, -1, dtype=float)
    return score / n


def main():
    frames = {t: load_target(t) for t in DEV_TARGETS + (TEST_TARGET,)}
    m4_smiles = set(frames[TEST_TARGET].SMILES)
    overlap = {t: int(frames[t].SMILES.isin(m4_smiles).sum()) for t in DEV_TARGETS}
    dev_clean = {t: frames[t][~frames[t].SMILES.isin(m4_smiles)].copy() for t in DEV_TARGETS}

    c_rows = []
    for c in (.01, .03, .1, .3, 1.0, 3.0, 10.0):
        aucs = []
        for vi, valid_target in enumerate(DEV_TARGETS):
            train_targets = [t for t in DEV_TARGETS if t != valid_target]
            model, _ = fit_model([dev_clean[t] for t in train_targets], c, 100 + vi)
            valid = dev_clean[valid_target]
            score = model.predict_proba(valid[list(FEATURES)])[:, 1]
            aucs.append(float(roc_auc_score(valid.target, score)))
        c_rows.append({"C": c, "dev_LOTO_macro_AUC": float(np.mean(aucs)), **{
            f"AUC_{t}": aucs[i] for i, t in enumerate(DEV_TARGETS)
        }})
    cv = pd.DataFrame(c_rows).sort_values(["dev_LOTO_macro_AUC", "C"], ascending=[False, True])
    selected_c = float(cv.iloc[0].C)
    cv.to_csv(OUT / "development_target_cv.csv", index=False)

    model, n_train = fit_model([dev_clean[t] for t in DEV_TARGETS], selected_c, 829)
    test = frames[TEST_TARGET].copy()
    y = test.target.to_numpy(int)
    test["CrossTargetLogistic"] = model.predict_proba(test[list(FEATURES)])[:, 1]
    test["PACER_XR"] = test[list(BASE_FEATURES)].mean(axis=1)
    test["PACER_XR_Cascade1pct"] = cascade_score(
        test.Glide_BEmin.to_numpy(float), test.PACER_XR.to_numpy(float), .01
    )

    methods = list(BASE_FEATURES) + ["CrossTargetLogistic", "PACER_XR", "PACER_XR_Cascade1pct"]
    result = []
    for method in methods:
        z = metrics(y, test[method].to_numpy(float))
        result.append({"method": method, **z})
    result = pd.DataFrame(result).sort_values("ROC_AUC", ascending=False)
    result.to_csv(OUT / "m4_blind_test_metrics.csv", index=False)
    test[["ligand_id", "SMILES", "target"] + methods].to_csv(
        OUT / "m4_blind_test_predictions.csv", index=False
    )

    single = result[result.method.isin(BASE_FEATURES)].iloc[0]
    delta, ci = auc_diff_ci(
        y, test.PACER_XR.to_numpy(float), test[single.method].to_numpy(float)
    )
    cascade_delta, cascade_ci = auc_diff_ci(
        y, test.PACER_XR_Cascade1pct.to_numpy(float), test[single.method].to_numpy(float)
    )
    cascade_row = result[result.method == "PACER_XR_Cascade1pct"].iloc[0]
    coef = dict(zip(FEATURES, model.coef_[0].astype(float)))
    audit = {
        "protocol": "M4-blind cross-target rank fusion",
        "official_repository_commit": "44798c841ee77230b1f89fe41074e41b458c7070",
        "development_targets": list(DEV_TARGETS),
        "test_target": TEST_TARGET,
        "m4_labels_used_for_training_or_selection": False,
        "exact_m4_smiles_removed_from_development": overlap,
        "selected_C": selected_c,
        "balanced_development_training_n": n_train,
        "m4_intersection_n": len(test),
        "m4_actives": int(y.sum()),
        "best_single_method_same_intersection": str(single.method),
        "best_single_AUC": float(single.ROC_AUC),
        "PACER_XR_AUC": float(result.loc[result.method == "PACER_XR", "ROC_AUC"].iloc[0]),
        "delta_AUC_vs_best_single": delta,
        "delta_AUC_CI95": ci,
        "PACER_XR_Cascade1pct": {
            "ROC_AUC": float(cascade_row.ROC_AUC),
            "PR_AUC": float(cascade_row.PR_AUC),
            "logAUC_pct": float(cascade_row.logAUC_pct),
            "EF_0.5pct": float(cascade_row["EF_0.5pct"]),
            "EF_1pct": float(cascade_row["EF_1pct"]),
            "EF_2pct": float(cascade_row["EF_2pct"]),
            "EF_5pct": float(cascade_row["EF_5pct"]),
            "delta_AUC_vs_best_single": cascade_delta,
            "delta_AUC_CI95": cascade_ci,
            "selection_rule": "top 1% Glide-BEmin, remaining 99% PACER-XR consensus",
        },
        "coefficients": coef,
        "claim_boundary": (
            "This tests broad ASD allosteric-modulator retrieval, not functional M4 PAM identity or efficacy."
        ),
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(cv.to_string(index=False))
    print(result.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
