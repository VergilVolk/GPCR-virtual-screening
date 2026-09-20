"""Leave-one-replica-out replica-invariant state representation pilot."""

from __future__ import annotations

import json
import random
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.linalg import eigh
from scipy.spatial.distance import jensenshannon
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.preprocessing import StandardScaler
from torch import nn
from torch.nn import functional as F


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "results" / "m4_gamd_public_recluster_v01" / "aligned_pocket_features.npy"
OUT = PROJECT / "results" / "pacer_replica_invariant_state_v02"
N_REP, N_FRAME, K, LAG = 6, 500, 10, 5
SEEDS = [17, 29, 43]


class StateEncoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(32, 64), nn.ReLU(), nn.Linear(64, 16))
        self.decoder = nn.Sequential(nn.Linear(16, 64), nn.ReLU(), nn.Linear(64, 32))

    def forward(self, x):
        z = F.normalize(self.encoder(x), dim=1)
        return z, self.decoder(z)


def make_partners(x, reps, train, seed):
    rng = np.random.default_rng(seed); train = np.asarray(train, int)
    temporal, cross, negative = [], [], []
    for a in train:
        rep = reps[a]
        near = train[(reps[train] == rep) & (np.abs(train - a) >= 1) & (np.abs(train - a) <= LAG)]
        temporal.append(int(rng.choice(near)))
        other = train[reps[train] != rep]
        cross.append(int(other[np.argmin(np.sum((x[other] - x[a]) ** 2, axis=1))]))
        far = train[(reps[train] == rep) & (np.abs(train - a) >= 50)]
        negative.append(int(far[np.argmax(np.sum((x[far] - x[a]) ** 2, axis=1))]))
    return train, np.asarray(temporal), np.asarray(cross), np.asarray(negative)


def train_encoder(x, reps, train, seed):
    torch.manual_seed(seed); np.random.seed(seed); random.seed(seed); torch.set_num_threads(6)
    a, pt, pc, n = make_partners(x, reps, train, seed)
    a, pt, pc, n = [torch.as_tensor(v, dtype=torch.long) for v in (a, pt, pc, n)]
    tx = torch.from_numpy(x.astype(np.float32)); model = StateEncoder()
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    train_reps = sorted(set(reps[train]))
    for _ in range(120):
        model.train(); opt.zero_grad(); z, recon = model(tx)
        dap_t = 1 - (z[a] * z[pt]).sum(1); dap_c = 1 - (z[a] * z[pc]).sum(1)
        dan = 1 - (z[a] * z[n]).sum(1)
        metric = .5 * F.relu(.2 + dap_t - dan).mean() + .5 * F.relu(.2 + dap_c - dan).mean()
        reconstruction = F.mse_loss(recon[a], tx[a])
        means = torch.stack([z[torch.as_tensor(np.where(reps == r)[0])].mean(0) for r in train_reps])
        invariance = means.var(0).mean()
        variance = F.relu(.10 - torch.sqrt(z[a].var(0) + 1e-4)).mean()
        loss = metric + .10 * reconstruction + .10 * invariance + .10 * variance
        loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        return model(tx)[0].numpy()


def tica(train_x, all_x, reps, train, lag=LAG):
    x0, xt = [], []
    for rep in sorted(set(reps[train])):
        idx = train[reps[train] == rep]
        x0.append(train_x[np.searchsorted(train, idx[:-lag])])
        xt.append(train_x[np.searchsorted(train, idx[lag:])])
    x0, xt = np.vstack(x0), np.vstack(xt); mean = np.vstack([x0, xt]).mean(0)
    x0, xt = x0 - mean, xt - mean
    c00 = (x0.T @ x0 + xt.T @ xt) / (2 * len(x0)) + np.eye(x0.shape[1]) * 1e-5
    c01 = (x0.T @ xt + xt.T @ x0) / (2 * len(x0))
    values, vectors = eigh(c01, c00); order = np.argsort(np.abs(values))[::-1][:5]
    return (all_x - mean) @ vectors[:, order]


def evaluate(method, z, reps, train, test, seed):
    km = KMeans(n_clusters=K, random_state=seed, n_init=20).fit(z[train])
    labels = km.predict(z)
    train_occ = np.bincount(labels[train], minlength=K) / len(train)
    test_occ = np.bincount(labels[test], minlength=K) / len(test)
    per_train = []
    for rep in sorted(set(reps[train])):
        idx = train[reps[train] == rep]
        per_train.append(np.bincount(labels[idx], minlength=K) / len(idx))
    entropy = -np.sum(test_occ[test_occ > 0] * np.log(test_occ[test_occ > 0]))
    return {
        "method": method, "heldout_replica": int(reps[test[0]] + 1), "seed": seed,
        "heldout_occupancy_JSD": float(jensenshannon(test_occ, train_occ, base=2) ** 2),
        "heldout_temporal_persistence": float(np.mean(labels[test][1:] == labels[test][:-1])),
        "common_state_coverage": int(np.sum((test_occ >= .01) & np.all(np.asarray(per_train) >= .01, axis=0))),
        "heldout_effective_states": float(np.exp(entropy)),
        "heldout_max_state_fraction": float(test_occ.max()),
    }, labels


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    raw = np.load(SOURCE)
    if raw.shape != (3000, 717): raise RuntimeError(raw.shape)
    reps = np.repeat(np.arange(N_REP), N_FRAME)
    rows, assignments, seed_amis = [], [], []
    for heldout in range(N_REP):
        test = np.where(reps == heldout)[0]; train = np.where(reps != heldout)[0]
        scaler = StandardScaler().fit(raw[train]); scaled = scaler.transform(raw)
        pca = PCA(n_components=32, random_state=17).fit(scaled[train]); pca32 = pca.transform(scaled).astype(np.float32)
        row, labels = evaluate("PCA-10", pca32[:, :10], reps, train, test, 17)
        rows.append(row)
        tica5 = tica(pca32[train], pca32, reps, train)
        row, labels = evaluate("tICA-5", tica5, reps, train, test, 17); rows.append(row)
        fold_labels = []
        for seed in SEEDS:
            z = train_encoder(pca32, reps, train, seed)
            row, labels = evaluate("ReplicaInvariantTriplet-16", z, reps, train, test, seed)
            rows.append(row); fold_labels.append(labels[test])
            assignments.extend({"heldout_replica": heldout + 1, "seed": seed,
                                "global_frame": int(i), "state": int(labels[i])} for i in test)
        amis = [adjusted_mutual_info_score(fold_labels[i], fold_labels[j])
                for i, j in combinations(range(len(SEEDS)), 2)]
        seed_amis.extend({"heldout_replica": heldout + 1, "AMI": float(x)} for x in amis)
        print(f"heldout_replica={heldout+1} mean_seed_AMI={np.mean(amis):.3f}", flush=True)
    metrics = pd.DataFrame(rows); metrics.to_csv(OUT / "fold_metrics.csv", index=False)
    pd.DataFrame(assignments).to_csv(OUT / "heldout_assignments.csv", index=False)
    pd.DataFrame(seed_amis).to_csv(OUT / "seed_stability.csv", index=False)
    summary = metrics.groupby("method").agg(
        mean_heldout_JSD=("heldout_occupancy_JSD", "mean"),
        mean_temporal_persistence=("heldout_temporal_persistence", "mean"),
        mean_common_coverage=("common_state_coverage", "mean"),
        mean_effective_states=("heldout_effective_states", "mean"),
        mean_max_state_fraction=("heldout_max_state_fraction", "mean")).reset_index()
    summary.to_csv(OUT / "summary.csv", index=False)
    triplet = summary.set_index("method").loc["ReplicaInvariantTriplet-16"]
    pca_row = summary.set_index("method").loc["PCA-10"]
    mean_ami = float(pd.DataFrame(seed_amis).AMI.mean())
    audit = {"outer_folds": N_REP, "functional_labels_used": False,
             "single_pam_chemotype": True, "mean_seed_AMI": mean_ami,
             "technical_go": bool(mean_ami >= .60 and triplet.mean_effective_states >= 3 and
                 triplet.mean_max_state_fraction <= .90 and triplet.mean_common_coverage >= 3 and
                 triplet.mean_heldout_JSD <= pca_row.mean_heldout_JSD and
                 triplet.mean_temporal_persistence >= pca_row.mean_temporal_persistence),
             "functional_statemil_go": False,
             "claim_boundary": "Leave-one-replica-out single-system state pilot; no PAM efficacy inference."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    (OUT / "REPORT.md").write_text("# PACER Replica-Invariant StateMetric v0.2\n\n" +
        summary.to_markdown(index=False, floatfmt=".4f") + "\n\n```json\n" + json.dumps(audit, indent=2) + "\n```\n", encoding="utf-8")
    print(summary.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
