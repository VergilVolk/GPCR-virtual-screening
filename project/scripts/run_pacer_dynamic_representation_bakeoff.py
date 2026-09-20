"""Leakage-controlled dynamic representation benchmark on six public M4 GaMD replicas.

This benchmark is deliberately label-free. It asks whether a representation learned on
five independent trajectories transfers to the sixth trajectory and preserves slow,
predictive dynamics. It does not evaluate PAM efficacy.
"""

from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.spatial.distance import jensenshannon
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "results" / "m4_gamd_public_recluster_v01" / "aligned_pocket_features.npy"
OUT = PROJECT / "results" / "pacer_dynamic_representation_bakeoff_v01"
N_REP, N_FRAME, LAG, K = 6, 500, 5, 8
RIDGE = 1e-5
SEED = 17


def lagged_pairs(x: np.ndarray, reps: np.ndarray, allowed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x0, x1 = [], []
    allowed_set = set(int(i) for i in allowed)
    for rep in sorted(set(reps[allowed])):
        idx = np.where(reps == rep)[0]
        keep = [(int(a), int(b)) for a, b in zip(idx[:-LAG], idx[LAG:])
                if int(a) in allowed_set and int(b) in allowed_set]
        if keep:
            a, b = zip(*keep)
            x0.append(x[np.asarray(a)]); x1.append(x[np.asarray(b)])
    return np.vstack(x0), np.vstack(x1)


def invsqrt_psd(c: np.ndarray, eps: float = RIDGE) -> np.ndarray:
    values, vectors = eigh((c + c.T) / 2)
    keep = values > eps
    if not np.any(keep):
        raise RuntimeError("covariance has no retained dimensions")
    return (vectors[:, keep] / np.sqrt(values[keep])) @ vectors[:, keep].T


def fit_tica(x: np.ndarray, reps: np.ndarray, train: np.ndarray, n_components: int = 5):
    x0, x1 = lagged_pairs(x, reps, train)
    mean = np.vstack([x0, x1]).mean(0)
    a, b = x0 - mean, x1 - mean
    c00 = (a.T @ a + b.T @ b) / (2 * len(a)) + np.eye(a.shape[1]) * RIDGE
    c01 = (a.T @ b + b.T @ a) / (2 * len(a))
    values, vectors = eigh(c01, c00)
    order = np.argsort(np.abs(values))[::-1][:n_components]
    return mean, vectors[:, order], values[order]


def fit_vamp(x: np.ndarray, reps: np.ndarray, train: np.ndarray, n_components: int = 5):
    x0, x1 = lagged_pairs(x, reps, train)
    mean0, mean1 = x0.mean(0), x1.mean(0)
    a, b = x0 - mean0, x1 - mean1
    c00 = a.T @ a / len(a) + np.eye(a.shape[1]) * RIDGE
    c11 = b.T @ b / len(b) + np.eye(b.shape[1]) * RIDGE
    c01 = a.T @ b / len(a)
    w0, w1 = invsqrt_psd(c00), invsqrt_psd(c11)
    u, s, vt = np.linalg.svd(w0 @ c01 @ w1, full_matrices=False)
    left = w0 @ u[:, :n_components]
    right = w1 @ vt.T[:, :n_components]
    return mean0, mean1, left, right, s[:n_components]


def predictive_r2(z: np.ndarray, reps: np.ndarray, train: np.ndarray, test: np.ndarray) -> float:
    tr0, tr1 = lagged_pairs(z, reps, train)
    te0, te1 = lagged_pairs(z, reps, test)
    model = Ridge(alpha=1e-3).fit(tr0, tr1)
    pred = model.predict(te0)
    numerator = float(np.sum((te1 - pred) ** 2))
    denominator = float(np.sum((te1 - tr1.mean(0)) ** 2))
    return 1.0 - numerator / denominator if denominator > 0 else float("nan")


def state_transfer(z: np.ndarray, reps: np.ndarray, train: np.ndarray, test: np.ndarray):
    km = KMeans(n_clusters=K, random_state=SEED, n_init=20).fit(z[train])
    labels = km.predict(z)
    train_occ = np.bincount(labels[train], minlength=K) / len(train)
    test_occ = np.bincount(labels[test], minlength=K) / len(test)
    seq = labels[test]
    return {
        "heldout_state_JSD": float(jensenshannon(test_occ, train_occ, base=2) ** 2),
        "heldout_temporal_persistence": float(np.mean(seq[1:] == seq[:-1])),
        "heldout_effective_states": float(np.exp(-np.sum(test_occ[test_occ > 0] * np.log(test_occ[test_occ > 0])))),
        "heldout_max_state_fraction": float(test_occ.max()),
    }


def exact_sign_flip_p(differences: np.ndarray) -> float:
    observed = abs(float(np.mean(differences)))
    null = [abs(float(np.mean(differences * np.asarray(signs))))
            for signs in product((-1.0, 1.0), repeat=len(differences))]
    return float(np.mean(np.asarray(null) >= observed - 1e-12))


def bootstrap_ci(differences: np.ndarray, seed: int = 1701):
    rng = np.random.default_rng(seed)
    means = [np.mean(rng.choice(differences, size=len(differences), replace=True)) for _ in range(20000)]
    return [float(x) for x in np.quantile(means, [0.025, 0.975])]


def markdown_table(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |",
             "| " + " | ".join("---" for _ in columns) + " |"]
    for _, row in frame.iterrows():
        values = []
        for column in columns:
            value = row[column]
            values.append(f"{value:.4f}" if isinstance(value, (float, np.floating)) else str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    raw = np.load(SOURCE)
    if raw.shape != (N_REP * N_FRAME, 717):
        raise RuntimeError(f"unexpected source shape {raw.shape}")
    reps = np.repeat(np.arange(N_REP), N_FRAME)
    rows, spectra = [], []
    for heldout in range(N_REP):
        test = np.where(reps == heldout)[0]
        train = np.where(reps != heldout)[0]
        scaler = StandardScaler().fit(raw[train])
        scaled = scaler.transform(raw)
        pca = PCA(n_components=32, random_state=SEED).fit(scaled[train])
        x = pca.transform(scaled)

        representations = {"PCA-10": x[:, :10]}
        mean, vectors, eig = fit_tica(x, reps, train)
        representations["tICA-5"] = (x - mean) @ vectors
        mean0, mean1, left, right, sing = fit_vamp(x, reps, train)
        # State coordinates use the forward (left) singular functions.
        representations["VAMP-5"] = (x - mean0) @ left
        for i, value in enumerate(eig):
            spectra.append({"heldout_replica": heldout + 1, "method": "tICA-5", "component": i + 1,
                            "value": float(value)})
        for i, value in enumerate(sing):
            spectra.append({"heldout_replica": heldout + 1, "method": "VAMP-5", "component": i + 1,
                            "value": float(value)})

        for method, z in representations.items():
            row = {"heldout_replica": heldout + 1, "method": method,
                   "heldout_lagged_R2": predictive_r2(z, reps, train, test)}
            row.update(state_transfer(z, reps, train, test))
            rows.append(row)

    folds = pd.DataFrame(rows)
    folds.to_csv(OUT / "fold_metrics.csv", index=False)
    pd.DataFrame(spectra).to_csv(OUT / "spectra.csv", index=False)
    summary = folds.groupby("method").agg(
        mean_lagged_R2=("heldout_lagged_R2", "mean"),
        sd_lagged_R2=("heldout_lagged_R2", "std"),
        mean_state_JSD=("heldout_state_JSD", "mean"),
        mean_temporal_persistence=("heldout_temporal_persistence", "mean"),
        mean_effective_states=("heldout_effective_states", "mean"),
        mean_max_state_fraction=("heldout_max_state_fraction", "mean"),
    ).reset_index()
    summary.to_csv(OUT / "summary.csv", index=False)

    pivot = folds.pivot(index="heldout_replica", columns="method", values="heldout_lagged_R2")
    comparisons = {}
    for baseline in ("PCA-10", "tICA-5"):
        d = (pivot["VAMP-5"] - pivot[baseline]).to_numpy()
        comparisons[f"VAMP-5_minus_{baseline}"] = {
            "mean_delta_lagged_R2": float(np.mean(d)),
            "bootstrap_CI95": bootstrap_ci(d),
            "exact_two_sided_sign_flip_p": exact_sign_flip_p(d),
            "better_folds": int(np.sum(d > 0)),
            "n_folds": int(len(d)),
        }
    best = summary.sort_values("mean_lagged_R2", ascending=False).iloc[0]
    audit = {
        "source": str(SOURCE.relative_to(PROJECT)),
        "frames": int(len(raw)),
        "independent_replicas": N_REP,
        "single_system": "M4 + iperoxo + MK-97",
        "functional_labels_used": False,
        "outer_split": "leave-one-replica-out",
        "lag_ns": LAG,
        "best_mean_predictive_method": str(best.method),
        "comparisons": comparisons,
        "functional_claim_allowed": False,
        "claim_boundary": "Single-system kinetic representation benchmark; cannot identify or rank functional PAMs.",
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = [
        "# PACER Dynamic Representation Bake-off v0.1", "",
        markdown_table(summary), "",
        "```json", json.dumps(audit, indent=2), "```", "",
        "This is a leakage-controlled representation benchmark on six independent replicas of one published PAM system.",
        "It evaluates transfer of slow/predictive dynamics, not PAM efficacy or cross-chemotype generalization.",
    ]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
