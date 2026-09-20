"""Frozen unsupervised state-representation pilot on six public M4 GaMD replicas."""

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
OUT = PROJECT / "results" / "pacer_unsupervised_state_pilot_v01"
N_REPLICAS = 6
FRAMES_PER_REPLICA = 500
K = 10
LAG = 5
SEEDS = [17, 29, 43]


class TemporalEncoder(nn.Module):
    def __init__(self, n_features):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(n_features, 64), nn.ReLU(),
                                 nn.Linear(64, 32), nn.ReLU(), nn.Linear(32, 16))

    def forward(self, x):
        return F.normalize(self.net(x), dim=1)


def replica_ids(n):
    return np.repeat(np.arange(N_REPLICAS), FRAMES_PER_REPLICA)[:n]


def temporal_pairs(n, seed):
    rng = np.random.default_rng(seed)
    anchors, positives, negatives = [], [], []
    for i in range(n):
        rep = i // FRAMES_PER_REPLICA
        lo = rep * FRAMES_PER_REPLICA
        hi = min((rep + 1) * FRAMES_PER_REPLICA, n)
        offsets = [d for d in range(-LAG, LAG + 1) if d and lo <= i + d < hi]
        if not offsets:
            continue
        p = i + int(rng.choice(offsets))
        other = int(rng.choice([r for r in range(N_REPLICAS) if r != rep]))
        nidx = other * FRAMES_PER_REPLICA + int(rng.integers(FRAMES_PER_REPLICA))
        anchors.append(i); positives.append(p); negatives.append(nidx)
    return np.asarray(anchors), np.asarray(positives), np.asarray(negatives)


def train_temporal(x, seed):
    torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
    torch.set_num_threads(6)
    model = TemporalEncoder(x.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    tx = torch.from_numpy(x.astype(np.float32))
    a, p, n = [torch.as_tensor(v, dtype=torch.long) for v in temporal_pairs(len(x), seed)]
    for _ in range(120):
        model.train(); opt.zero_grad()
        z = model(tx)
        dap = 1 - (z[a] * z[p]).sum(1)
        dan = 1 - (z[a] * z[n]).sum(1)
        loss = F.relu(.2 + dap - dan).mean()
        # Variance floor prevents all frames from collapsing onto one point.
        std = torch.sqrt(z.var(0) + 1e-4)
        loss = loss + .10 * F.relu(.10 - std).mean()
        loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        return model(tx).numpy()


def tica(x, reps, lag=LAG, n_components=5):
    x0, xt = [], []
    for rep in range(N_REPLICAS):
        idx = np.where(reps == rep)[0]
        x0.append(x[idx[:-lag]]); xt.append(x[idx[lag:]])
    x0, xt = np.vstack(x0), np.vstack(xt)
    mean = np.vstack([x0, xt]).mean(0)
    x0, xt = x0 - mean, xt - mean
    c00 = (x0.T @ x0 + xt.T @ xt) / (2 * len(x0))
    c01 = (x0.T @ xt + xt.T @ x0) / (2 * len(x0))
    c00 += np.eye(c00.shape[0]) * 1e-5
    values, vectors = eigh(c01, c00)
    order = np.argsort(np.abs(values))[::-1][:n_components]
    return (x - mean) @ vectors[:, order], values[order]


def evaluate(name, z, seed, reps):
    labels = KMeans(n_clusters=K, random_state=seed, n_init=20).fit_predict(z)
    occupancy = np.zeros((N_REPLICAS, K), float)
    for rep in range(N_REPLICAS):
        occupancy[rep] = np.bincount(labels[reps == rep], minlength=K)
        occupancy[rep] /= occupancy[rep].sum()
    jsd = [jensenshannon(occupancy[a], occupancy[b], base=2) ** 2
           for a, b in combinations(range(N_REPLICAS), 2)]
    global_occ = np.bincount(labels, minlength=K) / len(labels)
    entropy = -np.sum(global_occ[global_occ > 0] * np.log(global_occ[global_occ > 0]))
    transitions = []
    for rep in range(N_REPLICAS):
        seq = labels[reps == rep]
        transitions.extend(seq[1:] == seq[:-1])
    row = {
        "method": name, "seed": seed,
        "mean_replica_JSD": float(np.mean(jsd)),
        "max_replica_JSD": float(np.max(jsd)),
        "all_replica_state_coverage": int(np.sum(np.all(occupancy >= .01, axis=0))),
        "temporal_persistence": float(np.mean(transitions)),
        "effective_state_count": float(np.exp(entropy)),
        "max_state_fraction": float(global_occ.max()),
    }
    return row, labels, occupancy


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    raw = np.load(SOURCE)
    if raw.shape != (3000, 717):
        raise RuntimeError(f"unexpected feature shape: {raw.shape}")
    reps = replica_ids(len(raw))
    scaled = StandardScaler().fit_transform(raw)
    pca32_model = PCA(n_components=32, random_state=17).fit(scaled)
    pca32 = pca32_model.transform(scaled).astype(np.float32)
    representations = {"PCA-10": [(17, pca32[:, :10])]}
    tica5, tica_eigenvalues = tica(pca32, reps)
    representations["tICA-5"] = [(17, tica5)]
    representations["TemporalTriplet-16"] = [(seed, train_temporal(pca32, seed)) for seed in SEEDS]

    rows, label_map, occupancy_rows = [], {}, []
    for method, variants in representations.items():
        for seed, z in variants:
            row, labels, occupancy = evaluate(method, z, seed, reps)
            rows.append(row); label_map[(method, seed)] = labels
            for rep in range(N_REPLICAS):
                for state in range(K):
                    occupancy_rows.append({"method": method, "seed": seed, "replica": rep + 1,
                                           "state": state, "fraction": occupancy[rep, state]})

    temporal_labels = [label_map[("TemporalTriplet-16", s)] for s in SEEDS]
    amis = [adjusted_mutual_info_score(temporal_labels[i], temporal_labels[j])
            for i, j in combinations(range(len(SEEDS)), 2)]
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "method_metrics.csv", index=False)
    pd.DataFrame(occupancy_rows).to_csv(OUT / "replica_state_occupancy.csv", index=False)
    for (method, seed), labels in label_map.items():
        safe = method.lower().replace("-", "_")
        pd.DataFrame({"global_frame": np.arange(len(labels)), "replica": reps + 1,
                      "state": labels}).to_csv(OUT / f"assignments_{safe}_seed{seed}.csv", index=False)

    pca = result[result.method == "PCA-10"].iloc[0]
    temporal = result[result.method == "TemporalTriplet-16"]
    audit = {
        "public_frames": len(raw), "replicas": N_REPLICAS,
        "single_pam_chemotype": True, "functional_labels_used": False,
        "tica_lag_ns": LAG,
        "tica_eigenvalues": [float(x) for x in tica_eigenvalues],
        "temporal_triplet_seed_AMI": [float(x) for x in amis],
        "temporal_triplet_mean_AMI": float(np.mean(amis)),
        "technical_go": bool(np.mean(amis) >= .60 and
                             temporal.effective_state_count.mean() >= 3 and
                             temporal.max_state_fraction.mean() <= .90 and
                             temporal.all_replica_state_coverage.mean() >= 3 and
                             temporal.mean_replica_JSD.mean() <= pca.mean_replica_JSD and
                             temporal.temporal_persistence.mean() >= pca.temporal_persistence),
        "functional_statemil_go": False,
        "claim_boundary": "Unsupervised single-system replica-consistency pilot; no PAM efficacy inference.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# PACER Unsupervised State Pilot v0.1", "", result.to_markdown(index=False, floatfmt=".4f"),
              "", "```json", json.dumps(audit, indent=2), "```", "",
              "The pilot does not contain matched functional controls and cannot identify a PAM state."]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(result.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
