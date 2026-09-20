# -*- coding: utf-8 -*-
"""PACER-BG: Bayesian few-shot gating of series-centred SAR experts.

The complete evaluation is nested.  For every held-out medicinal-chemistry
series, expert reliability priors are estimated only from the remaining
development series.  Three label-free predicted-span anchors then update a
Beta-Binomial concordance model.  The posterior weights combine rank-normalised
expert predictions; low-evidence episodes shrink back to the centred LightGBM
reference.

This script evaluates a fixed family of prior strengths.  The reported primary
configuration is selected inside each outer fold, never on the outer labels.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.special import expit, logit
from scipy.stats import rankdata
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from run_pacer_fs_baselines import (
    DATA, MIN_GROUP, base_model, choose_predicted_span, features, safe_rho,
)


P = Path(__file__).resolve().parents[1]
OUT = P / "results" / "pacer_bayes_gate_v01"
OUT.mkdir(parents=True, exist_ok=True)
EXPERTS = ("LightGBM", "RandomForest", "ExtraTrees", "Ridge")
PRIOR_STRENGTHS = (2.0, 4.0, 8.0, 16.0)
TEMPERATURES = (1.0, 2.0, 4.0)


def make_models(seed: int):
    return {
        "LightGBM": LGBMRegressor(
            n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
            min_child_samples=12, reg_alpha=.1, reg_lambda=1,
            random_state=seed, verbosity=-1, n_jobs=6,
        ),
        "RandomForest": RandomForestRegressor(
            n_estimators=400, min_samples_leaf=2, max_features=.35,
            random_state=seed, n_jobs=6,
        ),
        "ExtraTrees": ExtraTreesRegressor(
            n_estimators=400, min_samples_leaf=2, max_features=.35,
            random_state=seed, n_jobs=6,
        ),
        "Ridge": make_pipeline(StandardScaler(), Ridge(alpha=20.0)),
    }


def centered_target(y, groups, train):
    z = np.zeros(len(y), float)
    for group in set(groups[train]):
        idx = train[groups[train] == group]
        z[idx] = y[idx] - y[idx].mean()
    return z


def design_matrix(X, D, train):
    sel = np.argsort(X[train].var(0))[-1024:]
    return np.c_[X[:, sel], D]


def expert_predictions(X, D, y, groups, train, test, seed):
    xx = design_matrix(X, D, train)
    z = centered_target(y, groups, train)
    return {
        name: model.fit(xx[train], z[train]).predict(xx[test])
        for name, model in make_models(seed).items()
    }


def pair_concordance(y, score, indices):
    correct = total = 0.0
    for a, b in itertools.combinations(indices, 2):
        dy = y[b] - y[a]
        ds = score[b] - score[a]
        if dy == 0 or ds == 0:
            continue
        total += 1.0
        correct += float(np.sign(dy) == np.sign(ds))
    return correct, total


def balanced_prior_accuracy(y, groups, predictions, valid_groups):
    """Macro-average pair concordance so large series cannot dominate."""
    out = {}
    for name in EXPERTS:
        vals = []
        for group in valid_groups:
            idx = np.where(groups == group)[0]
            c, n = pair_concordance(y, predictions[name], idx)
            if n:
                vals.append(c / n)
        out[name] = float(np.mean(vals)) if vals else .5
    return out


def posterior_weights(prior_acc, support_counts, strength, temperature):
    means = []
    for name in EXPERTS:
        p = float(np.clip(prior_acc[name], .501, .999))
        c, n = support_counts[name]
        post = (p * strength + c) / (strength + n)
        means.append(post)
    # Evidence is expressed as posterior log-odds above chance.  Negative
    # evidence does not reverse a model; it suppresses it.  A small LightGBM
    # floor implements conservative fallback to the validated PACER-FS core.
    evidence = np.maximum(logit(np.clip(means, .501, .999)), 0.0)
    weights = np.exp(temperature * (evidence - evidence.max()))
    weights /= weights.sum()
    floor = .35
    weights *= 1.0 - floor
    weights[0] += floor
    return dict(zip(EXPERTS, weights)), dict(zip(EXPERTS, means))


def combine(expert_scores, target, weights):
    n = len(target)
    ranks = np.vstack([
        rankdata(expert_scores[name], method="average") / n for name in EXPERTS
    ])
    return np.asarray([weights[name] for name in EXPERTS]) @ ranks


def main():
    d = pd.read_csv(DATA).reset_index(drop=True)
    X, D, _ = features(d.canonical_smiles)
    y = d.pEC50.to_numpy(float)
    groups = d.source_component.astype(str).to_numpy()
    eligible = [x for x, n in d.source_component.value_counts().items() if n >= MIN_GROUP]

    outer_rows = []
    selection_rows = []
    weight_rows = []

    for oi, outer_group in enumerate(eligible):
        outer = np.where(groups == outer_group)[0]
        dev_groups = [x for x in eligible if x != outer_group]

        # Build true nested OOF predictions for prior estimation and inner
        # configuration selection.  Outer molecules are absent throughout.
        dev_pred = {name: np.full(len(d), np.nan) for name in EXPERTS}
        dev_abs = np.full(len(d), np.nan)
        for ii, inner_group in enumerate(dev_groups):
            inner = np.where(groups == inner_group)[0]
            train = np.where((groups != outer_group) & (groups != inner_group))[0]
            pred = expert_predictions(X, D, y, groups, train, inner, 1000 + oi * 100 + ii)
            for name in EXPERTS:
                dev_pred[name][inner] = pred[name]
            dev_abs[inner] = base_model(X, D, y, train, inner)

        prior_acc = balanced_prior_accuracy(y, groups, dev_pred, dev_groups)

        configs = [(s, t) for s in PRIOR_STRENGTHS for t in TEMPERATURES]
        inner_scores = {cfg: [] for cfg in configs}
        for inner_group in dev_groups:
            target = np.where(groups == inner_group)[0]
            support = choose_predicted_span(target, 3, dev_abs)
            query = np.asarray([i for i in target if i not in set(support)], int)
            counts = {name: pair_concordance(y, dev_pred[name], support) for name in EXPERTS}
            local = {idx: j for j, idx in enumerate(target)}
            ql = np.asarray([local[i] for i in query])
            scores = {name: dev_pred[name][target] for name in EXPERTS}
            for cfg in configs:
                weights, _ = posterior_weights(prior_acc, counts, *cfg)
                pred = combine(scores, target, weights)[ql]
                inner_scores[cfg].append(safe_rho(y[query], pred))

        ranked = sorted(
            [(float(np.mean(v)), -float(np.std(v)), cfg) for cfg, v in inner_scores.items()],
            reverse=True,
        )
        selected = ranked[0][2]
        for cfg, vals in inner_scores.items():
            selection_rows.append({
                "outer_group": outer_group, "prior_strength": cfg[0],
                "temperature": cfg[1], "inner_macro_Spearman": float(np.mean(vals)),
            })

        train = np.where(groups != outer_group)[0]
        outer_pred = expert_predictions(X, D, y, groups, train, outer, 9000 + oi)
        outer_abs = np.full(len(d), np.nan)
        outer_abs[outer] = base_model(X, D, y, train, outer)
        support = choose_predicted_span(outer, 3, outer_abs)
        query = np.asarray([i for i in outer if i not in set(support)], int)
        local = {idx: j for j, idx in enumerate(outer)}
        ql = np.asarray([local[i] for i in query])
        counts = {
            name: pair_concordance(y, outer_pred[name], [local[i] for i in support])
            for name in EXPERTS
        }
        weights, posterior = posterior_weights(prior_acc, counts, *selected)
        gated = combine(outer_pred, outer, weights)[ql]
        reference = outer_pred["LightGBM"][ql]
        outer_rows.append({
            "group": outer_group, "n_group": len(outer), "n_query": len(query),
            "prior_strength": selected[0], "temperature": selected[1],
            "Spearman": safe_rho(y[query], gated),
            "Reference_Spearman": safe_rho(y[query], reference),
            "delta_Spearman": safe_rho(y[query], gated) - safe_rho(y[query], reference),
            "MAE_rankscale": float(mean_absolute_error(rankdata(y[query]) / len(query), gated)),
            "support_ids": "|".join(d.loc[support, "canonical_molecule_id"].astype(str)),
        })
        for name in EXPERTS:
            weight_rows.append({
                "group": outer_group, "expert": name, "prior_accuracy": prior_acc[name],
                "support_correct": counts[name][0], "support_pairs": counts[name][1],
                "posterior_accuracy": posterior[name], "weight": weights[name],
            })
        print(outer_group, selected, f"rho={outer_rows[-1]['Spearman']:.3f}", flush=True)

    out = pd.DataFrame(outer_rows)
    out.to_csv(OUT / "outer_episodes.csv", index=False)
    pd.DataFrame(selection_rows).to_csv(OUT / "inner_selection.csv", index=False)
    pd.DataFrame(weight_rows).to_csv(OUT / "expert_posteriors.csv", index=False)
    metrics = {
        "protocol": "nested_leave_one_series_out_3shot_bayesian_expert_gating",
        "n_outer_series": len(out),
        "macro_Spearman": float(out.Spearman.mean()),
        "reference_macro_Spearman_same_queries": float(out.Reference_Spearman.mean()),
        "delta_macro_Spearman": float(out.delta_Spearman.mean()),
        "median_Spearman": float(out.Spearman.median()),
        "worst_Spearman": float(out.Spearman.min()),
        "positive_series": int((out.Spearman > 0).sum()),
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
