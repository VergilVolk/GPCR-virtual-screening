"""Post-hoc anchor robustness for the frozen Acadia 2025 campaign.

This is a stress test, not a second independent external validation.  Anchor
triples are sampled without looking at labels, and every method receives the
same three functional measurements.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

from pacer_acm_v2 import select_three_anchors
from run_pacer_fs_baselines import adapted, features, sims, train_delta


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "benchmarks" / "m4_pam_v1"
TRAIN = ROOT / "potency_molecules.csv"
EXTERNAL = ROOT / "external_acadia_2025" / "external_molecule_benchmark.csv"
OUT = PROJECT / "results" / "pacer_external_acadia_anchor_robustness_v01"
SEED = 20260830
N_RANDOM = 1000


def rho(y: np.ndarray, p: np.ndarray) -> float:
    value = spearmanr(y, p).statistic
    return float(value) if np.isfinite(value) else 0.0


def top_recall(y: np.ndarray, p: np.ndarray, fraction: float = .10) -> float:
    k = max(1, int(np.ceil(len(y) * fraction)))
    return float(len(set(np.argsort(y)[-k:]) & set(np.argsort(p)[-k:])) / k)


def centered_lgbm(matrix, train_index, test_index, y_train, groups):
    target = np.zeros(len(train_index), float)
    for group in sorted(set(groups)):
        pos = np.where(groups == group)[0]
        target[pos] = y_train[pos] - y_train[pos].mean()
    model = LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1,
        random_state=12, verbosity=-1, n_jobs=6,
    ).fit(matrix[train_index], target)
    return model.predict(matrix[test_index]) + y_train.mean()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    train = pd.read_csv(TRAIN).reset_index(drop=True)
    ext_all = pd.read_csv(EXTERNAL).reset_index(drop=True)
    ext = ext_all[ext_all.eligible_exact_potency.astype(bool)].reset_index(drop=True)

    combined = pd.concat([train.canonical_smiles, ext.canonical_smiles], ignore_index=True)
    X, D, bv = features(combined)
    n_train = len(train)
    tr = np.arange(n_train)
    target = np.arange(n_train, n_train + len(ext))
    y_train = train.pEC50.to_numpy(float)
    y = np.r_[y_train, ext.pEC50.to_numpy(float)]
    groups = train.source_component.astype(str).to_numpy()
    group_all = np.r_[groups, np.repeat("ACADIA_2025", len(ext))]

    selected = np.argsort(X[tr].var(axis=0))[-1024:]
    matrix = np.c_[X[:, selected], D]
    base = np.full(len(combined), np.nan)
    base[target] = centered_lgbm(matrix, tr, target, y_train, groups)
    similarity = sims(target, tr, bv)
    weights = np.maximum(similarity, 1e-6) ** 3
    knn = np.full(len(combined), np.nan)
    knn[target] = weights @ y_train / weights.sum(axis=1)
    delta_pack = train_delta(tr, X, D, y, group_all, bv, seed=SEED, cliff_weight=0)

    deterministic_local = np.asarray(select_three_anchors(ext.canonical_smiles.tolist()), int)
    rng = np.random.default_rng(SEED)
    triples = [("deterministic_diverse", 0, deterministic_local)]
    triples += [("random", i + 1, np.sort(rng.choice(len(ext), 3, replace=False)))
                for i in range(N_RANDOM)]

    rows = []
    for policy, repeat, support_local in triples:
        support = target[support_local]
        query_local = np.asarray([i for i in range(len(ext)) if i not in set(support_local)], int)
        query = target[query_local]
        truth = y[query]
        offset = float(np.mean(y[support] - knn[support]))
        predictions = {
            "Similarity-kNN+3AnchorOffset": knn[query] + offset,
            "PACER-FS-TanimotoGP-Centered": adapted(
                "TanimotoGP", base, support, query, y, X, D, bv, delta_pack),
            "PACER-FS-DeltaSAR-Centered": adapted(
                "DeltaSAR", base, support, query, y, X, D, bv, delta_pack),
            "PACER-FS-DeltaSARHybrid-Centered": adapted(
                "DeltaSARHybrid", base, support, query, y, X, D, bv, delta_pack),
        }
        baseline_rho = rho(truth, predictions["Similarity-kNN+3AnchorOffset"])
        baseline_mae = mean_absolute_error(truth, predictions["Similarity-kNN+3AnchorOffset"])
        anchor_ids = ";".join(ext.iloc[support_local].compound_id)
        anchor_span = float(np.ptp(y[support]))
        query_support_sims = sims(query, support, bv)
        coverage = float(np.mean(query_support_sims.max(axis=1) >= .30))
        anchor_similarity = sims(support, support, bv)
        off_diagonal = anchor_similarity[np.triu_indices(3, 1)]
        for method, pred in predictions.items():
            rows.append({
                "policy": policy, "repeat": repeat, "method": method,
                "anchor_ids": anchor_ids, "anchor_pEC50_span": anchor_span,
                "anchor_base_prediction_span": float(np.ptp(base[support])),
                "anchor_knn_prediction_span": float(np.ptp(knn[support])),
                "anchor_residual_sd": float(np.std(y[support] - base[support], ddof=1)),
                "anchor_pair_tanimoto_mean": float(off_diagonal.mean()),
                "anchor_pair_tanimoto_max": float(off_diagonal.max()),
                "mean_query_max_anchor_tanimoto": float(query_support_sims.max(axis=1).mean()),
                "query_anchor_coverage_ge_0.30": coverage, "n_query": len(query),
                "Spearman": rho(truth, pred),
                "MAE": float(mean_absolute_error(truth, pred)),
                "top10_recall": top_recall(truth, pred),
                "delta_Spearman_vs_knn": rho(truth, pred) - baseline_rho,
                "MAE_improvement_vs_knn": float(baseline_mae - mean_absolute_error(truth, pred)),
            })

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "repeats.csv", index=False)
    random = frame[frame.policy == "random"]
    summary_rows = []
    for method, part in random.groupby("method"):
        summary_rows.append({
            "method": method, "n_anchor_triples": len(part),
            "median_Spearman": part.Spearman.median(),
            "Spearman_q025": part.Spearman.quantile(.025),
            "Spearman_q975": part.Spearman.quantile(.975),
            "median_MAE": part.MAE.median(),
            "MAE_q025": part.MAE.quantile(.025),
            "MAE_q975": part.MAE.quantile(.975),
            "fraction_delta_Spearman_gt_0": (part.delta_Spearman_vs_knn > 0).mean(),
            "fraction_MAE_improvement_gt_0": (part.MAE_improvement_vs_knn > 0).mean(),
            "median_top10_recall": part.top10_recall.median(),
        })
    summary = pd.DataFrame(summary_rows).sort_values("median_Spearman", ascending=False)
    summary.to_csv(OUT / "summary.csv", index=False)

    primary = random[random.method == "PACER-FS-DeltaSARHybrid-Centered"]
    audit = {
        "status": "post_hoc_anchor_robustness_not_independent_validation",
        "seed": SEED, "n_random_label_blind_triples": N_RANDOM,
        "n_exact_molecules": len(ext),
        "deterministic_anchor_ids": ext.iloc[deterministic_local].compound_id.tolist(),
        "primary_fraction_rank_better_than_knn": float((primary.delta_Spearman_vs_knn > 0).mean()),
        "primary_fraction_mae_better_than_knn": float((primary.MAE_improvement_vs_knn > 0).mean()),
        "claim_boundary": "Same Acadia campaign and labels; robustness only, not new external evidence or SOTA proof.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    lines = ["# Acadia three-anchor robustness", "",
             "This is a post-hoc stress test on the same frozen campaign, not an independent external validation.", "",
             "| Method | Median rho | 95% anchor range | Median MAE | P(rho > kNN) | P(MAE < kNN) |",
             "|---|---:|---:|---:|---:|---:|"]
    for row in summary.itertuples():
        lines.append(f"| {row.method} | {row.median_Spearman:.3f} | [{row.Spearman_q025:.3f}, {row.Spearman_q975:.3f}] | {row.median_MAE:.3f} | {row.fraction_delta_Spearman_gt_0:.3f} | {row.fraction_MAE_improvement_gt_0:.3f} |")
    lines += ["", "Anchor sensitivity is part of the deployable uncertainty. A favorable median cannot replace an untouched campaign."]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
