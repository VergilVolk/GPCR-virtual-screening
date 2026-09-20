"""Frozen ordinal stress test on the named US20260055116A1 subset."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import average_precision_score, balanced_accuracy_score
from sklearn.preprocessing import StandardScaler

from run_pacer_fs_baselines import adapted, choose_diverse, features, sims, train_delta


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "benchmarks" / "m4_pam_v1"
TRAIN = ROOT / "potency_molecules.csv"
SOURCE = ROOT / "external_us20260055116" / "named_structure_functional_subset.csv"
OUT = PROJECT / "results" / "pacer_external_us20260055116_v01"
SEED = 20260830
N_BOOT = 20_000
LATENT = {"A": 7.30, "B": 6.65, "C": 6.00, "D": 5.30}
ORDER = {"A": 3, "B": 2, "C": 1, "D": 0}


def collapse(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for smiles, part in frame[frame.high_confidence_named_structure.astype(bool)].groupby(
            "canonical_smiles", sort=False):
        first = part.iloc[0].copy()
        first["compound_id"] = "|".join(part.compound_id)
        first["source_row_count"] = len(part)
        for endpoint in ["human_m4_perk", "rat_m4_perk", "human_m4_gtpgs"]:
            values = [x for x in part[endpoint].astype(str) if x in ORDER]
            if len(set(values)) > 1:
                raise RuntimeError(f"Conflicting duplicate labels for {smiles}: {endpoint} {values}")
            first[endpoint] = values[0] if values else "N/A"
        rows.append(first)
    return pd.DataFrame(rows).reset_index(drop=True)


def lgbm(seed: int):
    return LGBMRegressor(n_estimators=350, learning_rate=.03, num_leaves=15,
                         max_depth=5, min_child_samples=12, reg_alpha=.1,
                         reg_lambda=1, random_state=seed, verbosity=-1, n_jobs=6)


def concordance(y, p):
    y, p = np.asarray(y), np.asarray(p, float)
    good, total = 0.0, 0
    for i in range(len(y)):
        for j in range(i + 1, len(y)):
            if y[i] == y[j]:
                continue
            total += 1
            direction = np.sign(y[i] - y[j]) * np.sign(p[i] - p[j])
            good += 1.0 if direction > 0 else .5 if direction == 0 else 0.0
    return float(good / total) if total else np.nan


def to_class(p):
    p = np.asarray(p, float)
    return np.where(p > 7.0, 3, np.where(p >= 6.30103, 2,
                    np.where(p >= 5.69897, 1, 0)))


def metrics(y, p):
    y = np.asarray(y, int); p = np.asarray(p, float)
    rho = spearmanr(y, p).statistic
    return {"ordinal_concordance": concordance(y, p),
            "Spearman": float(rho) if np.isfinite(rho) else 0.0,
            "balanced_accuracy": float(balanced_accuracy_score(y, to_class(p))),
            "A_vs_rest_PR_AUC": float(average_precision_score((y == 3).astype(int), p))}


def paired_bootstrap(y, candidate, reference):
    y = np.asarray(y, int)
    candidate = np.asarray(candidate, float)
    reference = np.asarray(reference, float)
    n = len(y)

    # A nonparametric sample bootstrap is equivalent to assigning each original
    # compound a multinomial count and weighting every unordered pair by the
    # product of its two counts.  Precomputing the pair outcomes preserves the
    # exact bootstrap estimand while avoiding millions of Python-level loops.
    left, right = np.triu_indices(n, k=1)
    comparable = y[left] != y[right]
    left, right = left[comparable], right[comparable]

    def pair_credit(prediction):
        direction = np.sign(y[left] - y[right]) * np.sign(
            prediction[left] - prediction[right])
        return np.where(direction > 0, 1.0, np.where(direction == 0, 0.5, 0.0))

    candidate_credit = pair_credit(candidate)
    reference_credit = pair_credit(reference)
    rng = np.random.default_rng(SEED)
    values = []
    for _ in range(N_BOOT):
        counts = np.bincount(rng.integers(0, n, n), minlength=n)
        pair_weights = counts[left] * counts[right]
        denominator = pair_weights.sum()
        if denominator:
            values.append(float(np.dot(pair_weights, candidate_credit - reference_credit)
                                / denominator))
    values = np.asarray(values)
    return {"estimate": float(concordance(y, candidate) - concordance(y, reference)),
            "ci95": [float(x) for x in np.quantile(values, [.025, .975])],
            "probability_gt_zero": float(np.mean(values > 0)),
            "valid_replicates": int(len(values))}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    train = pd.read_csv(TRAIN).reset_index(drop=True)
    ext = collapse(pd.read_csv(SOURCE))
    if ext.exact_train_overlap.any() or ext.canonical_smiles.nunique() != len(ext):
        raise RuntimeError("External overlap/deduplication contract violated")
    ext.to_csv(OUT / "external_unique_active_moieties.csv", index=False)

    combined = pd.concat([train.canonical_smiles, ext.canonical_smiles], ignore_index=True)
    X, D, bv = features(combined)
    n_train = len(train); tr = np.arange(n_train); target = np.arange(n_train, len(combined))
    y_train = train.pEC50.to_numpy(float)
    y_class = ext.human_m4_perk.map(ORDER).to_numpy(int)
    y_latent = ext.human_m4_perk.map(LATENT).to_numpy(float)

    selected = np.argsort(X[tr].var(axis=0))[-1024:]
    matrix = np.c_[X[:, selected], D]
    scaler = StandardScaler().fit(matrix[tr])
    scaled = scaler.transform(matrix)
    ridge = Ridge(alpha=10).fit(scaled[tr], y_train).predict(scaled[target])
    forest = ExtraTreesRegressor(n_estimators=500, min_samples_leaf=3,
                                 max_features="sqrt", random_state=SEED,
                                 n_jobs=6).fit(matrix[tr], y_train).predict(matrix[target])
    absolute = lgbm(11).fit(matrix[tr], y_train).predict(matrix[target])
    groups = train.source_component.astype(str).to_numpy()
    centered_y = np.zeros(n_train)
    for group in sorted(set(groups)):
        pos = np.where(groups == group)[0]
        centered_y[pos] = y_train[pos] - y_train[pos].mean()
    centered = lgbm(12).fit(matrix[tr], centered_y).predict(matrix[target]) + y_train.mean()
    similarity = sims(target, tr, bv)
    weights = np.maximum(similarity, 1e-6) ** 3
    knn = weights @ y_train / weights.sum(axis=1)
    nearest = y_train[np.argmax(similarity, axis=1)]
    zero = {"Constant": np.repeat(np.median(y_train), len(ext)), "ECFP4-1NN": nearest,
            "Similarity-kNN": knn, "Ridge": ridge, "ExtraTrees": forest,
            "Absolute-LGBM": absolute, "PACER-Centered-LGBM": centered}

    support_local = choose_diverse(np.arange(len(ext)), 7, bv[n_train:])
    # choose_diverse indexes the supplied local fingerprint list.
    support = target[support_local]
    query_local = np.asarray([i for i in range(len(ext)) if i not in set(support_local)], int)
    query = target[query_local]
    y_all = np.r_[y_train, y_latent]
    group_all = np.r_[groups, np.repeat("US20260055116", len(ext))]
    delta_pack = train_delta(tr, X, D, y_all, group_all, bv, seed=SEED, cliff_weight=0)
    base_full = np.full(len(combined), np.nan); base_full[target] = centered
    seven = {}
    for name, prediction in zero.items():
        offset = float(np.mean(y_latent[support_local] - prediction[support_local]))
        seven[name + "+7AnchorOffset"] = prediction[query_local] + offset
    for method in ["TanimotoGP", "DeltaSAR", "DeltaSARHybrid"]:
        seven[f"PACER-FS-{method}-7Anchor"] = adapted(
            method, base_full, support, query, y_all, X, D, bv, delta_pack)

    rows = []
    for protocol, truth, methods in [("zero_shot", y_class, zero),
                                      ("seven_anchor", y_class[query_local], seven)]:
        for name, prediction in methods.items():
            rows.append({"protocol": protocol, "method": name,
                         "n": len(truth), **metrics(truth, prediction)})
    metric_frame = pd.DataFrame(rows)
    metric_frame.to_csv(OUT / "metrics.csv", index=False)

    prediction = ext[["compound_id", "canonical_smiles", "human_m4_perk",
                      "rat_m4_perk", "human_m4_gtpgs", "max_train_tanimoto_ecfp4"]].copy()
    prediction["anchor_role"] = ["anchor" if i in set(support_local) else "query"
                                  for i in range(len(ext))]
    for name, value in zero.items(): prediction[name] = value
    for name, value in seven.items():
        prediction[name] = np.nan
        prediction.loc[query_local, name] = value
    prediction.to_csv(OUT / "predictions.csv", index=False)

    primary_name = "PACER-FS-DeltaSAR-7Anchor"
    primary = seven[primary_name]
    comparisons = {}
    for name, reference in seven.items():
        if name != primary_name:
            comparisons[name] = paired_bootstrap(y_class[query_local], primary, reference)
    all_lower_positive = all(value["ci95"][0] > 0 for value in comparisons.values())
    cross = {}
    for endpoint in ["rat_m4_perk", "human_m4_gtpgs"]:
        mask = ext[endpoint].isin(ORDER)
        cross[endpoint] = {
            "n": int(mask.sum()),
            "human_vs_endpoint_spearman": float(spearmanr(
                ext.loc[mask, "human_m4_perk"].map(ORDER),
                ext.loc[mask, endpoint].map(ORDER)).statistic),
            "exact_class_agreement": float(np.mean(
                ext.loc[mask, "human_m4_perk"] == ext.loc[mask, endpoint])),
        }
    audit = {
        "status": "frozen_after_source_qualification_ordinal_external_stress",
        "n_unique_active_moieties": len(ext), "n_query_after_7_anchors": len(query_local),
        "anchor_ids": ext.iloc[support_local].compound_id.tolist(),
        "median_max_train_tanimoto": float(ext.max_train_tanimoto_ecfp4.median()),
        "labels_visible_during_source_qualification": True,
        "model_or_feature_changed_after_label_view": False,
        "primary_method": primary_name, "paired_concordance_comparisons": comparisons,
        "primary_support_against_every_baseline": bool(all_lower_positive),
        "cross_readout_audit": cross,
        "claim_boundary": "Ordinal named-subset stress test; not pristine blind SOTA, exact potency, or prospective PAM confirmation.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# US20260055116 named-subset external stress test", "",
              f"Unique named active moieties: {len(ext)}; seven anchors; query N={len(query_local)}; exact historical overlap=0.", "",
              metric_frame.sort_values(["protocol", "ordinal_concordance"], ascending=[True, False]).to_markdown(index=False, floatfmt=".3f"), "",
              f"Primary support against every baseline: **{all_lower_positive}**.", "",
              "Labels were visible during source qualification, so even a positive result would not be pristine blind SOTA evidence."]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(metric_frame.sort_values(["protocol", "ordinal_concordance"], ascending=[True, False]).to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
