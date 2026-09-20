# -*- coding: utf-8 -*-
"""PACER-DeltaR: leakage-controlled residual delta learning for M4 PAM SAR."""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

RDLogger.DisableLog("rdApp.*")
P = Path(__file__).resolve().parents[1]
DATA = P / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
OUT = P / "results" / "pacer_deltar_v01"
OUT.mkdir(parents=True, exist_ok=True)
MIN_GROUP = 12
N_ANCHOR = 3
PAIR_LIMIT = 2500
SHRINKAGE = 1.5
DESC = [
    Descriptors.MolWt, Descriptors.MolLogP, Descriptors.NumHDonors,
    Descriptors.NumHAcceptors, Descriptors.TPSA, Descriptors.NumRotatableBonds,
    Descriptors.NumAromaticRings, Descriptors.FractionCSP3, Descriptors.HeavyAtomCount,
]


def features(smiles):
    fps, desc, bvs = [], [], []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        bv = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros(2048, np.float32)
        DataStructs.ConvertToNumpyArray(bv, arr)
        fps.append(arr); desc.append([f(mol) for f in DESC]); bvs.append(bv)
    return np.asarray(fps, np.float32), np.asarray(desc, np.float32), bvs


def fit_base(X, D, y, train, seed):
    sel = np.argsort(X[train].var(0))[-1024:]
    model = LGBMRegressor(
        n_estimators=350, learning_rate=.03, num_leaves=15, max_depth=5,
        min_child_samples=12, reg_alpha=.1, reg_lambda=1,
        random_state=seed, verbosity=-1, n_jobs=6,
    ).fit(np.c_[X[train][:, sel], D[train]], y[train])
    return model, sel


def predict_base(model, sel, X, D, idx):
    return model.predict(np.c_[X[idx][:, sel], D[idx]])


def pair_features(a, b, X, D, bvs, sel, dmean, dstd):
    da = X[b][:, sel] - X[a][:, sel]
    common = X[b][:, sel] * X[a][:, sel]
    dd = (D[b] - D[a]) / dstd
    mid = ((D[b] + D[a]) / 2 - dmean) / dstd
    sim = np.asarray([
        DataStructs.TanimotoSimilarity(bvs[i], bvs[j]) for i, j in zip(a, b)
    ])[:, None]
    return np.c_[da, common, dd, mid, sim].astype(np.float32)


def inner_oos_base(X, D, y, groups, outer_train, outer_group_index):
    """Base predictions for training molecules from models not seeing their series."""
    pred = np.full(len(y), np.nan)
    train_groups = sorted(set(groups[outer_train]))
    for inner_i, held in enumerate(train_groups):
        te = outer_train[groups[outer_train] == held]
        tr = outer_train[groups[outer_train] != held]
        model, sel = fit_base(X, D, y, tr, 1000 + outer_group_index * 100 + inner_i)
        pred[te] = predict_base(model, sel, X, D, te)
    if np.isnan(pred[outer_train]).any():
        raise RuntimeError("inner out-of-series base prediction is incomplete")
    return pred


def train_delta_residual(X, D, y, groups, bvs, outer_train, base_oos, seed):
    sel = np.argsort(X[outer_train].var(0))[-512:]
    dmean = D[outer_train].mean(0)
    dstd = D[outer_train].std(0); dstd[dstd < 1e-8] = 1
    residual = y - base_oos
    rng = np.random.default_rng(seed)
    aa, bb = [], []
    for group in sorted(set(groups[outer_train])):
        idx = outer_train[groups[outer_train] == group]
        pairs = np.asarray(list(itertools.combinations(idx, 2)), int)
        if len(pairs) > PAIR_LIMIT:
            pairs = pairs[rng.choice(len(pairs), PAIR_LIMIT, replace=False)]
        if len(pairs):
            aa.extend(pairs[:, 0]); bb.extend(pairs[:, 1])
            aa.extend(pairs[:, 1]); bb.extend(pairs[:, 0])
    aa = np.asarray(aa, int); bb = np.asarray(bb, int)
    z = pair_features(aa, bb, X, D, bvs, sel, dmean, dstd)
    target = residual[bb] - residual[aa]
    # High-similarity, large-residual changes carry the activity-cliff signal.
    weights = 1.0 + ((z[:, -1] >= .5) & (np.abs(target) >= .75)).astype(float)
    model = LGBMRegressor(
        n_estimators=400, learning_rate=.025, num_leaves=15, max_depth=6,
        min_child_samples=20, reg_alpha=.2, reg_lambda=2,
        random_state=seed, verbosity=-1, n_jobs=6,
    ).fit(z, target, sample_weight=weights)
    return model, sel, dmean, dstd, len(target)


def predicted_span(target, k, base):
    order = target[np.argsort(base[target])]
    positions = np.linspace(0, 1, k + 2)[1:-1]
    return np.asarray([order[int(round(q * (len(order) - 1)))] for q in positions], int)


def similarities(query, support, bvs):
    return np.asarray([
        [DataStructs.TanimotoSimilarity(bvs[q], bvs[a]) for a in support]
        for q in query
    ], float)


def predict_deltar(base, support, query, y, X, D, bvs, pack):
    model, sel, dmean, dstd, _ = pack
    aa = np.tile(support, len(query))
    qq = np.repeat(query, len(support))
    zaq = pair_features(aa, qq, X, D, bvs, sel, dmean, dstd)
    zqa = pair_features(qq, aa, X, D, bvs, sel, dmean, dstd)
    correction = .5 * (model.predict(zaq) - model.predict(zqa))
    correction = correction.reshape(len(query), len(support))
    anchor_residual = y[support] - base[support]
    sim = similarities(query, support, bvs)
    weights = np.maximum(sim, 1e-4) ** 3
    update = (weights * (anchor_residual[None, :] + correction)).sum(1)
    update /= weights.sum(1) + SHRINKAGE
    return base[query] + update


def safe_rho(y, pred):
    value = spearmanr(y, pred).statistic
    return float(value) if np.isfinite(value) else 0.0


def bootstrap_delta(rows, method, baseline="Absolute", n=50000, seed=829):
    pivot = rows.pivot(index="group", columns="method", values="Spearman")
    values = (pivot[method] - pivot[baseline]).to_numpy(float)
    rng = np.random.default_rng(seed)
    draws = values[rng.integers(0, len(values), size=(n, len(values)))].mean(1)
    return {
        "estimate": float(values.mean()),
        "ci95": [float(x) for x in np.quantile(draws, [.025, .975])],
        "probability_gt_zero": float((draws > 0).mean()),
    }


def main():
    data = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bvs = features(data.canonical_smiles)
    y = data.pEC50.to_numpy(float)
    groups = data.source_component.astype(str).to_numpy()
    eval_groups = [g for g, n in data.source_component.value_counts().items() if n >= MIN_GROUP]
    rows = []
    for outer_i, held in enumerate(eval_groups):
        target = np.flatnonzero(groups == held)
        train = np.flatnonzero(groups != held)
        base_model, base_sel = fit_base(X, D, y, train, 42 + outer_i)
        base = np.full(len(y), np.nan)
        base[target] = predict_base(base_model, base_sel, X, D, target)
        base_oos = inner_oos_base(X, D, y, groups, train, outer_i)
        pack = train_delta_residual(X, D, y, groups, bvs, train, base_oos, 2000 + outer_i)
        support = predicted_span(target, N_ANCHOR, base)
        query = np.asarray([x for x in target if x not in set(support)], int)
        sim = similarities(query, support, bvs)
        w = np.maximum(sim, 1e-4) ** 3
        offset = (w * (y[support] - base[support])[None, :]).sum(1) / (w.sum(1) + SHRINKAGE)
        predictions = {
            "Absolute": base[query],
            "AnchorResidual": base[query] + offset,
            "PACER-DeltaR": predict_deltar(base, support, query, y, X, D, bvs, pack),
        }
        for method, pred in predictions.items():
            rows.append({
                "group": held, "n_group": len(target), "n_query": len(query),
                "method": method, "Spearman": safe_rho(y[query], pred),
                "MAE": float(mean_absolute_error(y[query], pred)),
                "support_ids": "|".join(data.loc[support, "canonical_molecule_id"]),
                "delta_training_pairs": pack[-1],
            })
        print(held, "query", len(query), "pairs", pack[-1], flush=True)

    result = pd.DataFrame(rows)
    result.to_csv(OUT / "series_results.csv", index=False)
    summary = result.groupby("method").agg(
        macro_Spearman=("Spearman", "mean"), median_Spearman=("Spearman", "median"),
        worst_Spearman=("Spearman", "min"),
        positive_series=("Spearman", lambda x: int((x > 0).sum())),
        macro_MAE=("MAE", "mean"), n_series=("Spearman", "size"),
    ).reset_index().sort_values("macro_Spearman", ascending=False)
    summary.to_csv(OUT / "summary.csv", index=False)
    deltar = bootstrap_delta(result, "PACER-DeltaR")
    anchor = bootstrap_delta(result, "AnchorResidual")
    audit = {
        "protocol": "11-series LOSO; 3 predicted-span anchors; inner-LOSO residual targets",
        "outer_labels_used_for_model_or_anchor_selection": False,
        "n_molecules": int(len(data)), "n_series": int(len(eval_groups)),
        "PACER_DeltaR_vs_Absolute": deltar,
        "AnchorResidual_vs_Absolute": anchor,
        "preregistered_threshold": {
            "delta_ci_lower_gt_zero": bool(deltar["ci95"][0] > 0),
            "macro_spearman_not_below_frozen_PACER_FS_0.263": bool(
                summary.set_index("method").loc["PACER-DeltaR", "macro_Spearman"] >= .263
            ),
        },
        "claim_boundary": "internal series-held-out functional potency benchmark; not external SOTA or PAM confirmation",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
