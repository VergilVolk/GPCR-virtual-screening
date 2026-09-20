"""Strict leave-series-out PACER functional hard-negative metric learning."""

from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error
from torch import nn
from torch.nn import functional as F

from run_pacer_fs_baselines import base_model, choose_predicted_span, features


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
OUT = PROJECT / "results" / "pacer_functional_triplet_v01"
SEEDS = [17, 29, 43]
EPOCHS = 160
MARGIN = .20
METRIC_WEIGHT = .50


class MetricRegressor(nn.Module):
    def __init__(self, n_features: int):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_features, 256), nn.ReLU(), nn.Dropout(.10),
            nn.Linear(256, 64), nn.ReLU(),
        )
        self.metric = nn.Linear(64, 32)
        self.regressor = nn.Linear(64, 1)

    def forward(self, x):
        hidden = self.encoder(x)
        embedding = F.normalize(self.metric(hidden), dim=1)
        prediction = self.regressor(hidden).squeeze(1)
        return embedding, prediction


def safe_rho(y, p):
    value = spearmanr(y, p).statistic
    return float(value) if np.isfinite(value) else 0.0


def make_input(X, D, train):
    mean = D[train].mean(axis=0)
    std = D[train].std(axis=0)
    std[std < 1e-8] = 1
    return np.c_[X, (D - mean) / std].astype(np.float32)


def pair_lists(train, y, groups, bv):
    positives = {}
    hard = {}
    broad = {}
    train_set = set(map(int, train))
    for group in sorted(set(groups[train])):
        idx = [int(i) for i in np.where(groups == group)[0] if int(i) in train_set]
        for i in idx:
            positives.setdefault(i, []); hard.setdefault(i, []); broad.setdefault(i, [])
        for a_pos, i in enumerate(idx):
            for j in idx[a_pos + 1:]:
                from rdkit import DataStructs
                similarity = DataStructs.TanimotoSimilarity(bv[i], bv[j])
                delta = abs(y[i] - y[j])
                if similarity >= .50 and delta <= .30:
                    positives[i].append(j); positives[j].append(i)
                if similarity >= .55 and delta >= 1.0:
                    hard[i].append(j); hard[j].append(i)
                if delta > .30:
                    broad[i].append(j); broad[j].append(i)
    return positives, hard, broad


def build_triplets(train, y, groups, bv, kind, seed):
    positives, hard, broad = pair_lists(train, y, groups, bv)
    rng = np.random.default_rng(seed)
    rows = []
    for anchor in sorted(positives):
        if not positives[anchor]:
            continue
        candidates = hard[anchor] if kind == "hard" else broad[anchor]
        if not candidates:
            continue
        for positive in positives[anchor]:
            chosen = candidates if kind == "hard" else [int(rng.choice(candidates))]
            for negative in chosen:
                rows.append((anchor, positive, negative))
    if not rows:
        return np.empty((0, 3), int)
    rows = np.asarray(rows, int)
    if len(rows) > 6000:
        rows = rows[rng.choice(len(rows), 6000, replace=False)]
    return rows


def train_neural(inputs, centered_y, train, groups, bv, method, seed):
    torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
    torch.set_num_threads(6)
    model = MetricRegressor(inputs.shape[1])
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    x = torch.from_numpy(inputs)
    target = torch.from_numpy(centered_y.astype(np.float32))
    train_tensor = torch.as_tensor(train, dtype=torch.long)
    triplets = None
    if method != "MLP-Reg":
        kind = "hard" if method == "PACER-FM-HardTriplet" else "random"
        triplets = build_triplets(train, centered_y, groups, bv, kind, seed)
        triplets = torch.as_tensor(triplets, dtype=torch.long)
    for _ in range(EPOCHS):
        model.train(); optimizer.zero_grad()
        embedding, prediction = model(x)
        regression = F.mse_loss(prediction[train_tensor], target[train_tensor])
        metric = torch.tensor(0.0)
        if triplets is not None and len(triplets):
            a, p, n = triplets.T
            dap = 1 - (embedding[a] * embedding[p]).sum(dim=1)
            dan = 1 - (embedding[a] * embedding[n]).sum(dim=1)
            metric = F.relu(MARGIN + dap - dan).mean()
        loss = regression + METRIC_WEIGHT * metric
        loss.backward(); optimizer.step()
    model.eval()
    with torch.no_grad():
        _, prediction = model(x)
    return prediction.numpy(), 0 if triplets is None else len(triplets)


def cliff_accuracy(query, y, bv):
    pairs = []
    query = list(map(int, query))
    from rdkit import DataStructs
    for pos, i in enumerate(query):
        for j in query[pos + 1:]:
            if DataStructs.TanimotoSimilarity(bv[i], bv[j]) >= .55 and abs(y[i] - y[j]) >= 1.0:
                pairs.append((i, j))
    return pairs


def group_bootstrap(frame, candidate, reference, n_boot=20_000):
    pivot = frame.pivot(index="group", columns="method", values="Spearman")
    delta = pivot[candidate] - pivot[reference]
    rng = np.random.default_rng(20260830)
    values = np.empty(n_boot)
    raw = delta.to_numpy(float)
    for i in range(n_boot):
        values[i] = rng.choice(raw, len(raw), replace=True).mean()
    return {"estimate": float(raw.mean()),
            "ci95": [float(x) for x in np.quantile(values, [.025, .975])],
            "probability_gt_zero": float(np.mean(values > 0))}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bv = features(data.canonical_smiles)
    y = data.pEC50.to_numpy(float)
    groups = data.source_component.astype(str).to_numpy()
    outer_groups = [g for g, n in data.source_component.value_counts().items() if n >= 12]
    records, predictions, cliff_records = [], [], []
    neural_methods = ["MLP-Reg", "MLP-RandomTriplet", "PACER-FM-HardTriplet"]
    for fold, group in enumerate(outer_groups):
        train = np.where(groups != group)[0]
        target = np.where(groups == group)[0]
        absolute = np.full(len(data), np.nan)
        absolute[target] = base_model(X, D, y, train, target)
        support = choose_predicted_span(target, 3, absolute)
        query = np.asarray([i for i in target if i not in set(support)], int)

        centered_y = np.zeros(len(data), float)
        for train_group in sorted(set(groups[train])):
            idx = np.where((groups == train_group) & np.isin(np.arange(len(data)), train))[0]
            centered_y[idx] = y[idx] - y[idx].mean()
        centered = np.full(len(data), np.nan)
        centered[target] = base_model(X, D, centered_y, train, target)
        classical = {"Absolute-LightGBM": absolute.copy(), "PACER-FS": centered.copy()}
        inputs = make_input(X, D, train)
        for method in neural_methods:
            seed_predictions, counts = [], []
            for seed in SEEDS:
                pred, n_triplets = train_neural(inputs, centered_y, train, groups, bv, method, seed)
                seed_predictions.append(pred); counts.append(n_triplets)
            classical[method] = np.mean(seed_predictions, axis=0)
            print(group, method, "triplets", counts, flush=True)

        cliffs = cliff_accuracy(query, y, bv)
        for method, full_prediction in classical.items():
            offset = float(np.mean(y[support] - full_prediction[support]))
            pred = full_prediction[query] + offset
            rho = safe_rho(y[query], pred)
            mae = float(mean_absolute_error(y[query], pred))
            records.append({"group": group, "method": method, "n_query": len(query),
                            "Spearman": rho, "MAE": mae,
                            "support_ids": "|".join(data.iloc[support].canonical_molecule_id)})
            for idx, value in zip(query, pred):
                predictions.append({"group": group, "method": method,
                                    "canonical_molecule_id": data.iloc[idx].canonical_molecule_id,
                                    "observed": y[idx], "predicted": float(value)})
            correct = [np.sign(y[i] - y[j]) == np.sign(full_prediction[i] - full_prediction[j])
                       for i, j in cliffs]
            cliff_records.append({"group": group, "method": method,
                                  "n_cliff_pairs": len(correct),
                                  "cliff_direction_accuracy": float(np.mean(correct)) if correct else np.nan})

    episodes = pd.DataFrame(records)
    episodes.to_csv(OUT / "outer_episodes.csv", index=False)
    pd.DataFrame(predictions).to_csv(OUT / "outer_predictions.csv", index=False)
    cliffs = pd.DataFrame(cliff_records)
    cliffs.to_csv(OUT / "outer_cliff_metrics.csv", index=False)
    summary = episodes.groupby("method").agg(
        macro_Spearman=("Spearman", "mean"), worst_Spearman=("Spearman", "min"),
        positive_series=("Spearman", lambda x: int((x > 0).sum())),
        macro_MAE=("MAE", "mean")).reset_index()
    cliff_summary = cliffs.groupby("method").agg(
        evaluable_series=("cliff_direction_accuracy", "count"),
        macro_cliff_accuracy=("cliff_direction_accuracy", "mean"),
        total_query_cliff_pairs=("n_cliff_pairs", "sum")).reset_index()
    summary = summary.merge(cliff_summary, on="method", how="left")
    summary.to_csv(OUT / "summary.csv", index=False)
    comparisons = {reference: group_bootstrap(episodes,
                    "PACER-FM-HardTriplet", reference)
                   for reference in ["PACER-FS", "MLP-Reg", "MLP-RandomTriplet"]}
    indexed = summary.set_index("method")
    hard = indexed.loc["PACER-FM-HardTriplet"]
    fs = indexed.loc["PACER-FS"]
    internal_go = bool(comparisons["PACER-FS"]["ci95"][0] > 0 and
                       hard.worst_Spearman >= fs.worst_Spearman and
                       hard.macro_cliff_accuracy > indexed.loc["MLP-Reg"].macro_cliff_accuracy)
    audit = {"outer_groups": len(outer_groups), "three_anchor_labels_used_for_rank": False,
             "triplets_constructed_after_outer_split": True,
             "hyperparameters_selected_from_outer_results": False,
             "comparisons": comparisons, "internal_go": internal_go,
             "external_promotion": False,
             "claim_boundary": "Internal series-held-out development; not external SOTA or PAM confirmation."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# PACER Functional Triplet v0.1", "", summary.to_markdown(index=False, floatfmt=".3f"),
              "", "```json", json.dumps(audit, indent=2), "```"]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(summary.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
