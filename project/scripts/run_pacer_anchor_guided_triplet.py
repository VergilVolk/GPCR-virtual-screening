"""LOPO anchor-guided triplet representation for M4 GaMD snapshots."""

from __future__ import annotations

import json
import random
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.spatial.distance import jensenshannon
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.preprocessing import StandardScaler
from torch import nn
from torch.nn import functional as F


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "pacer_anchor_guided_triplet_v03"
FRAME = PROJECT / "results" / "pacer_consensus_pam_signature_v01" / "frame_signature_scores.csv"
FEATURES = PROJECT / "results" / "pacer_consensus_pam_signature_v01" / "signature_features.csv"
PAIR_DISTANCE = PROJECT / "results" / "pacer_consensus_pam_signature_v01" / "ca_pair_distance_features_210.npy"
SEEDS = [17, 29, 43]
K, N_REP, N_FRAME = 6, 6, 500


class Encoder(nn.Module):
    def __init__(self):
        super().__init__()
        self.enc = nn.Sequential(nn.Linear(32, 64), nn.ReLU(), nn.Linear(64, 16))
        self.dec = nn.Sequential(nn.Linear(16, 64), nn.ReLU(), nn.Linear(64, 32))

    def forward(self, x):
        z = F.normalize(self.enc(x), dim=1)
        return z, self.dec(z)


def eta_squared(v, labels):
    total = np.sum((v - v.mean()) ** 2)
    return float(sum(np.sum(labels == c) * (v[labels == c].mean() - v.mean()) ** 2
                     for c in np.unique(labels)) / total) if total > 0 else 0.0


def select_pair_distances() -> np.ndarray:
    distances = np.load(PAIR_DISTANCE)
    if distances.shape != (3000, 210):
        raise RuntimeError(distances.shape)
    selected = pd.read_csv(FEATURES)
    static_ids = [89, 92, 93, 96, 97, 98, 108, 183, 184, 185, 186,
                  187, 190, 416, 420, 423, 432, 433, 435, 436, 439]
    pair_to_i = {(static_ids[a], static_ids[b]): i
                 for i, (a, b) in enumerate(combinations(range(21), 2))}
    cols = [pair_to_i[(int(r.residue_i), int(r.residue_j))] for r in selected.itertuples()]
    if len(cols) != 141:
        raise RuntimeError(f"Expected 141 frozen features, got {len(cols)}")
    return distances[:, cols].astype(np.float32)


def partners(x, score, reps, train):
    train = np.asarray(train, int)
    positives, negatives, temporal = [], [], []
    for a in train:
        other = train[reps[train] != reps[a]]
        order = other[np.argsort(np.abs(score[other] - score[a]))[:20]]
        positives.append(int(order[np.argmin(np.sum((x[order] - x[a]) ** 2, axis=1))]))
        hard = train[np.abs(score[train] - score[a]) >= .50]
        if len(hard) == 0:
            hard = train[np.argsort(np.abs(score[train] - score[a]))[-100:]]
        negatives.append(int(hard[np.argmin(np.sum((x[hard] - x[a]) ** 2, axis=1))]))
        same = train[(reps[train] == reps[a]) & (np.abs(train - a) >= 1) & (np.abs(train - a) <= 5)]
        temporal.append(int(same[np.argmin(np.abs(same - a))]))
    return tuple(torch.as_tensor(v, dtype=torch.long) for v in
                 (train, np.asarray(positives), np.asarray(negatives), np.asarray(temporal)))


def train_encoder(x, score, reps, train, seed):
    torch.manual_seed(seed); np.random.seed(seed); random.seed(seed); torch.set_num_threads(6)
    a, p, n, t = partners(x, score, reps, train)
    tx = torch.from_numpy(x.astype(np.float32)); model = Encoder()
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    for _ in range(120):
        model.train(); opt.zero_grad(); z, recon = model(tx)
        dap = 1 - (z[a] * z[p]).sum(1); dan = 1 - (z[a] * z[n]).sum(1)
        dat = 1 - (z[a] * z[t]).sum(1)
        metric = F.relu(.20 + dap - dan).mean() + .25 * F.relu(.10 + dat - dan).mean()
        reconstruction = F.mse_loss(recon[a], tx[a])
        variance = F.relu(.10 - torch.sqrt(z[a].var(0) + 1e-4)).mean()
        loss = metric + .10 * reconstruction + .10 * variance
        loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        return model(tx)[0].numpy()


def metrics(method, z, score, reps, train, test, seed):
    km = KMeans(n_clusters=K, random_state=seed, n_init=30).fit(z[train])
    labels = km.predict(z)
    train_occ = np.bincount(labels[train], minlength=K) / len(train)
    test_occ = np.bincount(labels[test], minlength=K) / len(test)
    entropy = -np.sum(test_occ[test_occ > 0] * np.log(test_occ[test_occ > 0]))
    return {"method": method, "heldout_replica": int(reps[test[0]]), "seed": seed,
            "eta_squared": eta_squared(score[test], labels[test]),
            "temporal_persistence": float(np.mean(labels[test][1:] == labels[test][:-1])),
            "occupancy_JSD": float(jensenshannon(test_occ, train_occ, base=2) ** 2),
            "effective_states": float(np.exp(entropy)), "max_state_fraction": float(test_occ.max())}, labels[test]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    xraw = select_pair_distances()
    frame = pd.read_csv(FRAME)
    score = frame.signature_projection.to_numpy(float)
    reps = frame.replica.to_numpy(int)
    rows, stability, assignments = [], [], []
    for heldout in range(1, N_REP + 1):
        train, test = np.where(reps != heldout)[0], np.where(reps == heldout)[0]
        scaler = StandardScaler().fit(xraw[train]); scaled = scaler.transform(xraw)
        pca = PCA(n_components=32, random_state=17).fit(scaled[train])
        x = pca.transform(scaled).astype(np.float32)
        row, lab = metrics("PCA-10", x[:, :10], score, reps, train, test, 17); rows.append(row)
        assignments.extend({"method":"PCA-10","heldout_replica":heldout,"seed":17,
                            "global_frame":int(g),"state":int(y)} for g,y in zip(test,lab))
        seed_labels = []
        for seed in SEEDS:
            z = train_encoder(x, score, reps, train, seed)
            row, lab = metrics("AnchorGuidedTriplet-16", z, score, reps, train, test, seed)
            rows.append(row); seed_labels.append(lab)
            assignments.extend({"method":"AnchorGuidedTriplet-16","heldout_replica":heldout,"seed":seed,
                                "global_frame":int(g),"state":int(y)} for g,y in zip(test,lab))
        for a, b in combinations(range(3), 2):
            stability.append({"heldout_replica": heldout,
                              "AMI": adjusted_mutual_info_score(seed_labels[a], seed_labels[b])})
        print(f"heldout={heldout} seed_AMI={np.mean([r['AMI'] for r in stability if r['heldout_replica']==heldout]):.3f}", flush=True)
    folds = pd.DataFrame(rows); folds.to_csv(OUT / "fold_metrics.csv", index=False)
    pd.DataFrame(assignments).to_csv(OUT / "heldout_assignments.csv", index=False)
    stab = pd.DataFrame(stability); stab.to_csv(OUT / "seed_stability.csv", index=False)
    summary = folds.groupby("method").agg(mean_eta_squared=("eta_squared", "mean"),
        mean_temporal_persistence=("temporal_persistence", "mean"), mean_JSD=("occupancy_JSD", "mean"),
        mean_effective_states=("effective_states", "mean"), mean_max_state_fraction=("max_state_fraction", "mean")).reset_index()
    summary.to_csv(OUT / "summary.csv", index=False)
    pivot = folds.groupby(["method", "heldout_replica"]).eta_squared.mean().unstack(0)
    delta = pivot["AnchorGuidedTriplet-16"] - pivot["PCA-10"]
    tri = summary.set_index("method").loc["AnchorGuidedTriplet-16"]
    pca_row = summary.set_index("method").loc["PCA-10"]
    mean_ami = float(stab.AMI.mean())
    audit = {"outer_folds": 6, "heldout_signature_used_in_training": False,
             "pseudo_label_is_static_structural_signature": True, "mean_seed_AMI": mean_ami,
             "triplet_minus_pca_eta_squared": float(delta.mean()),
             "triplet_better_folds": int((delta > 0).sum()),
             "technical_go": bool(delta.mean() >= .05 and (delta > 0).sum() >= 4 and mean_ami >= .60 and
                 tri.mean_temporal_persistence >= pca_row.mean_temporal_persistence - .05 and
                 tri.mean_effective_states >= 3 and tri.mean_max_state_fraction <= .90),
             "functional_state_go": False,
             "claim_boundary": "Weakly supervised structural-state metric; no PAM efficacy inference."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = "# PACER Anchor-Guided Trajectory Triplet v0.3\n\n" + summary.to_markdown(index=False, floatfmt=".4f")
    report += "\n\n## Held-out eta-squared\n\n" + pivot.reset_index().to_markdown(index=False, floatfmt=".4f")
    report += "\n\n```json\n" + json.dumps(audit, indent=2) + "\n```\n"
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(summary.to_string(index=False)); print(pivot.to_string()); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
