"""Preregistered DeepRLI frozen-encoder functional adapter benchmark."""

from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from torch import nn
from torch.nn import functional as F

from run_pacer_fs_baselines import base_model, choose_predicted_span, features


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
EMBED = PROJECT / "results" / "pacer_deeprli_m4_v01"
OUT = PROJECT / "results" / "pacer_deeprli_adapter_v01"
STATES = ["7TRQ", "7TRP", "7TRS"]
SEEDS = [17, 29, 43]
EPOCHS = 160


class BagAdapter(nn.Module):
    def __init__(self):
        super().__init__()
        self.adapter = nn.Sequential(nn.Linear(64, 32), nn.ReLU(), nn.Linear(32, 64))
        nn.init.zeros_(self.adapter[-1].weight); nn.init.zeros_(self.adapter[-1].bias)
        self.attention = nn.Linear(64, 1)
        self.metric = nn.Linear(64, 32)
        self.regressor = nn.Linear(64, 1)

    def forward(self, x):
        residual = self.adapter(x)
        h = x + residual
        weights = torch.softmax(self.attention(h).squeeze(-1), dim=1)
        pooled = (h * weights[..., None]).sum(1)
        metric = F.normalize(self.metric(pooled), dim=1)
        return metric, self.regressor(pooled).squeeze(1), residual, weights


def safe_rho(y, p):
    value = spearmanr(y, p).statistic
    return float(value) if np.isfinite(value) else 0.0


def load_bags(data):
    scores = pd.read_csv(EMBED / "deeprli_scores.csv")
    zframe = pd.read_csv(EMBED / "deeprli_z64.csv")
    zcols = [f"z{i:02d}" for i in range(64)]
    score_map = scores.set_index("complex_path")
    zmap = zframe.set_index("complex_path")
    bags, score_tensor = [], []
    for mol_id in data.canonical_molecule_id.astype(str):
        keys = [f"{state}/{mol_id}" for state in STATES]
        bags.append(zmap.loc[keys, zcols].to_numpy(np.float32))
        score_tensor.append(score_map.loc[keys, ["scoring_score", "docking_score", "screening_score"]].to_numpy(float))
    return np.asarray(bags), np.asarray(score_tensor)


def partner_lists(train, y, groups, bv, frozen_mean):
    train_set = set(map(int, train)); positive, hard, broad = {}, {}, {}
    norm = frozen_mean / np.maximum(np.linalg.norm(frozen_mean, axis=1, keepdims=True), 1e-8)
    for group in sorted(set(groups[train])):
        idx = [int(i) for i in np.where(groups == group)[0] if int(i) in train_set]
        for i in idx:
            positive[i], hard[i], broad[i] = [], [], []
        for pos, i in enumerate(idx):
            for j in idx[pos + 1:]:
                sim = DataStructs.TanimotoSimilarity(bv[i], bv[j]); delta = abs(y[i] - y[j])
                if sim >= .50 and delta <= .30:
                    positive[i].append(j); positive[j].append(i)
                if sim >= .55 and delta >= 1.0:
                    hard[i].append(j); hard[j].append(i)
                if delta > .30:
                    broad[i].append(j); broad[j].append(i)
    # Freeze the nearest binding-space cliff for each anchor.
    nearest = {}
    for i, candidates in hard.items():
        if candidates:
            nearest[i] = min(candidates, key=lambda j: 1 - float(norm[i] @ norm[j]))
    return positive, nearest, broad


def train_adapter(bags, centered_y, train, groups, bv, method, seed):
    torch.manual_seed(seed); np.random.seed(seed); random.seed(seed); torch.set_num_threads(6)
    mean = bags[train].reshape(-1, 64).mean(0); std = bags[train].reshape(-1, 64).std(0); std[std < 1e-6] = 1
    x = torch.from_numpy(((bags - mean) / std).astype(np.float32))
    target = torch.from_numpy(centered_y.astype(np.float32)); tr = torch.as_tensor(train, dtype=torch.long)
    model = BagAdapter(); opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    positive, nearest, broad = partner_lists(train, centered_y, groups, bv, bags.mean(1))
    rng = np.random.default_rng(seed)
    usable = [i for i in train if positive.get(int(i)) and
              ((method == "HardTriplet-Adapter" and int(i) in nearest) or
               (method == "RandomTriplet-Adapter" and broad.get(int(i))))]
    for _ in range(EPOCHS):
        model.train(); opt.zero_grad(); metric_z, pred, residual, _ = model(x)
        loss = F.mse_loss(pred[tr], target[tr]) + .10 * residual[tr].pow(2).mean()
        if method != "Adapter-Reg" and usable:
            by_group = []
            for group in sorted(set(groups[usable])):
                anchors = [int(i) for i in usable if groups[int(i)] == group]
                triplets = []
                for a in anchors:
                    p = int(rng.choice(positive[a]))
                    n = int(nearest[a]) if method == "HardTriplet-Adapter" else int(rng.choice(broad[a]))
                    triplets.append((a, p, n))
                t = torch.as_tensor(triplets, dtype=torch.long)
                a, p, n = t.T
                dap = 1 - (metric_z[a] * metric_z[p]).sum(1)
                dan = 1 - (metric_z[a] * metric_z[n]).sum(1)
                by_group.append(F.relu(.2 + dap - dan).mean())
            loss = loss + .50 * torch.stack(by_group).mean()
        loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        _, pred, _, attention = model(x)
    return pred.numpy(), attention.numpy(), len(usable)


def group_bootstrap(frame, candidate, reference, n_boot=20000):
    pivot = frame.pivot(index="group", columns="method", values="Spearman")
    raw = (pivot[candidate] - pivot[reference]).to_numpy(float)
    rng = np.random.default_rng(20260830)
    samples = np.asarray([rng.choice(raw, len(raw), replace=True).mean() for _ in range(n_boot)])
    return {"estimate": float(raw.mean()), "ci95": np.quantile(samples, [.025, .975]).tolist(),
            "probability_gt_zero": float(np.mean(samples > 0))}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(DATA).reset_index(drop=True)
    X, D, bv = features(data.canonical_smiles); y = data.pEC50.to_numpy(float)
    groups = data.source_component.astype(str).to_numpy(); bags, scores = load_bags(data)
    outer_groups = [g for g, n in data.source_component.value_counts().items() if n >= 12]
    records, predictions, attention_rows = [], [], []
    methods_neural = ["Adapter-Reg", "RandomTriplet-Adapter", "HardTriplet-Adapter"]
    for group in outer_groups:
        train = np.where(groups != group)[0]; target_idx = np.where(groups == group)[0]
        absolute = np.full(len(data), np.nan); absolute[target_idx] = base_model(X, D, y, train, target_idx)
        support = choose_predicted_span(target_idx, 3, absolute)
        query = np.asarray([i for i in target_idx if i not in set(support)], int)
        centered_y = np.zeros(len(data))
        for g in sorted(set(groups[train])):
            idx = np.where((groups == g) & np.isin(np.arange(len(data)), train))[0]
            centered_y[idx] = y[idx] - y[idx].mean()
        fs = np.full(len(data), np.nan); fs[target_idx] = base_model(X, D, centered_y, train, target_idx)
        mean_z = bags.mean(1).reshape(len(data), -1)
        ridge = Ridge(alpha=10).fit(mean_z[train], centered_y[train]).predict(mean_z)
        models = {"PACER-FS": fs, "FrozenZ-Ridge": ridge,
                  "DeepRLI-Scoring-Mean": -scores[:, :, 0].mean(1),
                  "DeepRLI-Screening-Mean": -scores[:, :, 2].mean(1)}
        for method in methods_neural:
            ps, ats, counts = [], [], []
            for seed in SEEDS:
                pred, attention, count = train_adapter(bags, centered_y, train, groups, bv, method, seed)
                ps.append(pred); ats.append(attention); counts.append(count)
            models[method] = np.mean(ps, axis=0)
            mean_attention = np.mean(ats, axis=0)
            for i in query:
                attention_rows.append({"group": group, "method": method,
                                       "canonical_molecule_id": data.iloc[i].canonical_molecule_id,
                                       **{f"attention_{state}": mean_attention[i, s] for s, state in enumerate(STATES)}})
            print(group, method, "usable_anchors", counts, flush=True)
        for method, full in models.items():
            offset = float(np.mean(y[support] - full[support])); pred = full[query] + offset
            rho = safe_rho(y[query], pred); mae = float(mean_absolute_error(y[query], pred))
            records.append({"group": group, "method": method, "n_query": len(query),
                            "Spearman": rho, "MAE": mae})
            for i, value in zip(query, pred):
                predictions.append({"group": group, "method": method,
                                    "canonical_molecule_id": data.iloc[i].canonical_molecule_id,
                                    "observed": y[i], "predicted": float(value)})
    episodes = pd.DataFrame(records); episodes.to_csv(OUT / "outer_episodes.csv", index=False)
    pd.DataFrame(predictions).to_csv(OUT / "outer_predictions.csv", index=False)
    pd.DataFrame(attention_rows).to_csv(OUT / "state_attention.csv", index=False)
    summary = episodes.groupby("method").agg(macro_Spearman=("Spearman", "mean"),
        worst_Spearman=("Spearman", "min"), positive_series=("Spearman", lambda x: int((x > 0).sum())),
        macro_MAE=("MAE", "mean")).reset_index()
    summary.to_csv(OUT / "summary.csv", index=False)
    comparisons = {r: group_bootstrap(episodes, "HardTriplet-Adapter", r)
                   for r in ["PACER-FS", "Adapter-Reg", "RandomTriplet-Adapter", "FrozenZ-Ridge"]}
    hard = summary.set_index("method").loc["HardTriplet-Adapter"]
    fsrow = summary.set_index("method").loc["PACER-FS"]
    audit = {"outer_groups": len(outer_groups), "encoder_frozen": True,
             "triplets_built_after_outer_split": True, "state_bag_not_independent_samples": True,
             "comparisons": comparisons,
             "internal_go": bool(comparisons["PACER-FS"]["ci95"][0] > 0 and hard.worst_Spearman >= fsrow.worst_Spearman),
             "external_promotion": False,
             "claim_boundary": "Internal frozen-encoder transfer; no external SOTA or confirmed PAM."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    (OUT / "REPORT.md").write_text("# PACER-DeepRLI Adapter v0.1\n\n" +
        summary.to_markdown(index=False, floatfmt=".3f") + "\n\n```json\n" + json.dumps(audit, indent=2) + "\n```\n", encoding="utf-8")
    print(summary.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
