"""Annotate frozen trajectory clusters with an independently constructed PAM signature."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import run_pacer_replica_invariant_state as state_v02


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "pacer_triplet_cluster_anchor_alignment_v01"
RAW = PROJECT / "results" / "m4_gamd_public_recluster_v01" / "aligned_pocket_features.npy"
SIG = PROJECT / "results" / "pacer_consensus_pam_signature_v01" / "frame_signature_scores.csv"
TRIPLET = PROJECT / "results" / "pacer_replica_invariant_state_v02" / "heldout_assignments.csv"
N_PER_REP = 500
N_PERM = 1000


def eta_squared(values: np.ndarray, labels: np.ndarray) -> float:
    grand = float(values.mean())
    total = float(np.sum((values - grand) ** 2))
    if total <= 0:
        return 0.0
    between = 0.0
    for label in np.unique(labels):
        group = values[labels == label]
        between += len(group) * float((group.mean() - grand) ** 2)
    return between / total


def frozen_baseline_labels(raw: np.ndarray, reps0: np.ndarray) -> dict[str, dict[int, np.ndarray]]:
    output = {"PCA-10": {}, "tICA-5": {}}
    for heldout in range(6):
        test = np.where(reps0 == heldout)[0]
        train = np.where(reps0 != heldout)[0]
        scaler = StandardScaler().fit(raw[train])
        scaled = scaler.transform(raw)
        pca = PCA(n_components=32, random_state=17).fit(scaled[train])
        pca32 = pca.transform(scaled).astype(np.float32)
        km = KMeans(n_clusters=10, random_state=17, n_init=20).fit(pca32[train, :10])
        output["PCA-10"][heldout + 1] = km.predict(pca32[test, :10])
        tica5 = state_v02.tica(pca32[train], pca32, reps0, train)
        km = KMeans(n_clusters=10, random_state=17, n_init=20).fit(tica5[train])
        output["tICA-5"][heldout + 1] = km.predict(tica5[test])
    return output


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = np.load(RAW)
    sig = pd.read_csv(SIG)
    if raw.shape[0] != len(sig) or len(sig) != 3000:
        raise RuntimeError((raw.shape, len(sig)))
    values = sig["signature_projection"].to_numpy(float)
    reps0 = np.repeat(np.arange(6), N_PER_REP)
    base = frozen_baseline_labels(raw, reps0)
    trip = pd.read_csv(TRIPLET)

    assignments: dict[str, dict[int, list[np.ndarray]]] = {
        "PCA-10": {r: [base["PCA-10"][r]] for r in range(1, 7)},
        "tICA-5": {r: [base["tICA-5"][r]] for r in range(1, 7)},
        "ReplicaInvariantTriplet-16": {},
    }
    for r in range(1, 7):
        assignments["ReplicaInvariantTriplet-16"][r] = []
        subset = trip[trip.heldout_replica == r]
        for seed in sorted(subset.seed.unique()):
            x = subset[subset.seed == seed].sort_values("global_frame")
            if len(x) != N_PER_REP:
                raise RuntimeError(f"Replica {r}, seed {seed}: {len(x)} assignments")
            assignments["ReplicaInvariantTriplet-16"][r].append(x.state.to_numpy(int))

    rows = []
    observed = {}
    for method, by_rep in assignments.items():
        fold_values = []
        for r in range(1, 7):
            v = values[(r - 1) * N_PER_REP:r * N_PER_REP]
            seed_etas = [eta_squared(v, labels) for labels in by_rep[r]]
            fold_eta = float(np.mean(seed_etas))
            fold_values.append(fold_eta)
            rows.append({"method": method, "heldout_replica": r,
                         "eta_squared": fold_eta, "seed_eta_sd": float(np.std(seed_etas))})
        observed[method] = float(np.mean(fold_values))

    rng = np.random.default_rng(20260901)
    nulls = {m: [] for m in assignments}
    for _ in range(N_PERM):
        shifted = {}
        for r in range(1, 7):
            v = values[(r - 1) * N_PER_REP:r * N_PER_REP]
            shifted[r] = np.roll(v, int(rng.integers(25, N_PER_REP - 25)))
        for method, by_rep in assignments.items():
            fold_etas = []
            for r in range(1, 7):
                fold_etas.append(float(np.mean([eta_squared(shifted[r], lab) for lab in by_rep[r]])))
            nulls[method].append(float(np.mean(fold_etas)))

    summary = []
    for method in assignments:
        null = np.asarray(nulls[method])
        summary.append({"method": method, "mean_eta_squared": observed[method],
                        "circular_shift_null_mean": float(null.mean()),
                        "empirical_p": float((1 + np.sum(null >= observed[method])) / (N_PERM + 1))})
    summary_df = pd.DataFrame(summary)
    fold_df = pd.DataFrame(rows)
    fold_df.to_csv(OUT / "fold_alignment.csv", index=False)
    summary_df.to_csv(OUT / "summary.csv", index=False)
    pivot = fold_df.pivot(index="heldout_replica", columns="method", values="eta_squared")
    delta = pivot["ReplicaInvariantTriplet-16"] - pivot["PCA-10"]
    audit = {
        "triplet_trained_before_anchor_signature": True,
        "retrospective_external_annotation": True,
        "triplet_retrained": False,
        "triplet_minus_pca_mean_eta_squared": float(delta.mean()),
        "triplet_better_than_pca_folds": int((delta > 0).sum()),
        "triplet_increment_supported": bool(delta.mean() > 0 and (delta > 0).sum() >= 4),
        "claim_boundary": "Cluster alignment with a static structural signature; not PAM efficacy validation.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = "# PACER Triplet-Cluster External Annotation v0.1\n\n"
    report += summary_df.to_markdown(index=False, floatfmt=".4f")
    report += "\n\n## Held-out replicas\n\n" + pivot.reset_index().to_markdown(index=False, floatfmt=".4f")
    report += "\n\n```json\n" + json.dumps(audit, indent=2) + "\n```\n"
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(summary_df.to_string(index=False))
    print(pivot.to_string())
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
