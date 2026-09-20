"""Frozen internal-coordinate test of a cross-chemotype M4 PAM signature."""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.preprocessing import StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "pacer_consensus_pam_signature_v01"
STATIC = PROJECT / "results" / "structure" / "ensemble"
PUBLIC = PROJECT / "data" / "m4_gamd_figshare_33283491"
PUBLIC_RESIDS = [687, 690, 691, 694, 695, 696, 706, 781, 782, 783, 784,
                 785, 788, 852, 856, 859, 868, 869, 871, 872, 875]
STATIC_RESIDS = [89, 92, 93, 96, 97, 98, 108, 183, 184, 185, 186,
                 187, 190, 416, 420, 423, 432, 433, 435, 436, 439]
THRESHOLDS = [0.0, 0.05, 0.10, 0.20]
PRIMARY = 0.10
K = 6


def static_ca(path: Path) -> dict[int, tuple[str, np.ndarray]]:
    atoms = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[21:22] != "R":
            continue
        if line[12:16].strip() != "CA":
            continue
        resid = int(line[22:26])
        atoms[resid] = (line[17:20].strip(), np.array([
            float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float))
    return atoms


def pair_features(coords: np.ndarray) -> np.ndarray:
    pairs = np.asarray(list(combinations(range(coords.shape[-2]), 2)), dtype=int)
    delta = coords[..., pairs[:, 0], :] - coords[..., pairs[:, 1], :]
    return np.sqrt(np.sum(delta * delta, axis=-1))


def load_static() -> tuple[np.ndarray, np.ndarray, np.ndarray, list[tuple[int, int]]]:
    maps = [static_ca(STATIC / f"{p}_R.pdb") for p in ("7TRQ", "7TRP", "7TRS")]
    static_ids = STATIC_RESIDS
    missing = [[r for r in static_ids if r not in m] for m in maps]
    if any(missing):
        raise RuntimeError(f"Missing static CA residues: {missing}")
    coords = [np.stack([m[r][1] for r in static_ids]) for m in maps]
    residue_pairs = [(static_ids[i], static_ids[j])
                     for i, j in combinations(range(len(static_ids)), 2)]
    return pair_features(coords[0]), pair_features(coords[1]), pair_features(coords[2]), residue_pairs


def load_trajectories() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    all_features, replicas, times = [], [], []
    for rep in range(1, 7):
        u = mda.Universe(str(PUBLIC / "sys-protein.pdb"),
                         str(PUBLIC / f"run{rep}-GaMD-stride1000ps.nc"))
        ag = u.select_atoms("name CA and resid " + " ".join(map(str, PUBLIC_RESIDS)))
        by_resid = {int(a.resid): a for a in ag}
        if set(by_resid) != set(PUBLIC_RESIDS):
            raise RuntimeError(f"Replica {rep}: CA selection mismatch {sorted(by_resid)}")
        frames = []
        for ts in u.trajectory:
            frames.append(np.stack([by_resid[r].position.copy() for r in PUBLIC_RESIDS]))
        feats = pair_features(np.asarray(frames, dtype=float))
        all_features.append(feats)
        replicas.extend([rep] * len(feats))
        times.extend(range(len(feats)))
    return np.vstack(all_features), np.asarray(replicas), np.asarray(times)


def signature(q: np.ndarray, p: np.ndarray, s: np.ndarray,
              x: np.ndarray, threshold: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dq, dp = q - s, p - s
    mask = (dq * dp > 0) & (np.abs(dq) >= threshold) & (np.abs(dp) >= threshold)
    d = (dq + dp) / 2
    agreement = np.minimum(np.abs(dq), np.abs(dp)) / np.maximum(np.abs(dq), np.abs(dp))
    v = d[mask] * np.sqrt(agreement[mask])
    denom = float(np.dot(v, v))
    if denom <= 0:
        raise RuntimeError(f"No signature at threshold {threshold}")
    score = ((x[:, mask] - s[mask]) * np.sqrt(agreement[mask])) @ v / denom
    anchor_scores = np.array([
        np.dot((q[mask] - s[mask]) * np.sqrt(agreement[mask]), v) / denom,
        np.dot((p[mask] - s[mask]) * np.sqrt(agreement[mask]), v) / denom,
        0.0,
    ])
    return score, mask, anchor_scores


def loro_clusters(x: np.ndarray, replicas: np.ndarray) -> np.ndarray:
    labels = np.full(len(x), -1, dtype=int)
    for heldout in range(1, 7):
        train, test = replicas != heldout, replicas == heldout
        scaler = StandardScaler().fit(x[train])
        train_x, all_x = scaler.transform(x[train]), scaler.transform(x)
        full = PCA(random_state=17).fit(train_x)
        ncomp = min(10, int(np.searchsorted(np.cumsum(full.explained_variance_ratio_), .90) + 1))
        pca = PCA(n_components=ncomp, random_state=17).fit(train_x)
        z = pca.transform(all_x)
        km = KMeans(n_clusters=K, random_state=17, n_init=30).fit(z[train])
        # Cluster ids are fold-local and are used only for held-out enrichment summaries.
        labels[test] = km.predict(z[test])
    if np.any(labels < 0):
        raise RuntimeError("Unassigned held-out frames")
    return labels


def global_clusters(x: np.ndarray) -> np.ndarray:
    scaled = StandardScaler().fit_transform(x)
    full = PCA(random_state=17).fit(scaled)
    ncomp = min(10, int(np.searchsorted(np.cumsum(full.explained_variance_ratio_), .90) + 1))
    z = PCA(n_components=ncomp, random_state=17).fit_transform(scaled)
    return KMeans(n_clusters=K, random_state=17, n_init=30).fit_predict(z)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    q, p, s, pairs = load_static()
    x, replicas, times = load_trajectories()
    np.save(OUT / "ca_pair_distance_features_210.npy", x.astype(np.float32))
    sensitivity, primary_score, primary_mask = [], None, None
    for threshold in THRESHOLDS:
        score, mask, anchors = signature(q, p, s, x, threshold)
        medians = [float(np.median(score[replicas == r])) for r in range(1, 7)]
        sensitivity.append({
            "threshold_A": threshold, "features": int(mask.sum()),
            "7TRQ_anchor_score": float(anchors[0]), "7TRP_anchor_score": float(anchors[1]),
            "positive_replica_medians": int(sum(v > 0 for v in medians)),
            **{f"replica_{r}_median": medians[r - 1] for r in range(1, 7)},
        })
        if threshold == PRIMARY:
            primary_score, primary_mask = score, mask
    pd.DataFrame(sensitivity).to_csv(OUT / "threshold_sensitivity.csv", index=False)

    selected_x = x[:, primary_mask]
    labels = loro_clusters(selected_x, replicas)
    global_labels = global_clusters(selected_x)
    frame = pd.DataFrame({"replica": replicas, "frame_ns": times,
                          "signature_projection": primary_score,
                          "fold_local_cluster": labels,
                          "global_cluster": global_labels})
    frame.to_csv(OUT / "frame_signature_scores.csv", index=False)
    rep_rows = []
    for r in range(1, 7):
        vals = primary_score[replicas == r]
        rep_rows.append({"replica": r, "mean": float(vals.mean()),
                         "median": float(np.median(vals)),
                         "fraction_positive": float(np.mean(vals > 0))})
    rep_df = pd.DataFrame(rep_rows)
    rep_df.to_csv(OUT / "replica_summary.csv", index=False)

    # Fold-local cluster ids are not comparable across held-out folds. AMI is permutation invariant.
    cluster_rows = []
    for r in range(1, 7):
        for c in range(K):
            idx = (replicas == r) & (labels == c)
            if not idx.any():
                continue
            cluster_rows.append({"heldout_replica": r, "cluster": c, "frames": int(idx.sum()),
                                 "fraction_of_replica": float(idx.mean() * 6),
                                 "median_signature": float(np.median(primary_score[idx])),
                                 "fraction_positive": float(np.mean(primary_score[idx] > 0))})
    cluster_df = pd.DataFrame(cluster_rows)
    cluster_df.to_csv(OUT / "heldout_cluster_summary.csv", index=False)

    global_rows = []
    shared_positive = []
    for c in range(K):
        idx = global_labels == c
        occupancy = [float(np.mean(global_labels[replicas == r] == c)) for r in range(1, 7)]
        row = {"cluster": c, "frames": int(idx.sum()), "replicas_ge_1pct": int(sum(v >= .01 for v in occupancy)),
               "median_signature": float(np.median(primary_score[idx])),
               "fraction_positive": float(np.mean(primary_score[idx] > 0)),
               **{f"replica_{r}_occupancy": occupancy[r-1] for r in range(1, 7)}}
        global_rows.append(row)
        if all(v >= .01 for v in occupancy) and row["median_signature"] > 0:
            shared_positive.append(c)
    global_df = pd.DataFrame(global_rows)
    global_df.to_csv(OUT / "global_cluster_summary.csv", index=False)
    loro_ami = float(adjusted_mutual_info_score(global_labels, labels))

    sens = pd.DataFrame(sensitivity).set_index("threshold_A")
    sensitivity_consistent = bool(all(sens.loc[t, "positive_replica_medians"] == 6
                                      for t in (0.05, 0.20)))
    common_cluster_supported = bool(shared_positive)
    audit = {
        "frames": int(len(x)), "distance_features_total": int(x.shape[1]),
        "primary_signature_features": int(primary_mask.sum()),
        "all_six_replica_median_positive": bool((rep_df["median"] > 0).all()),
        "sensitivity_direction_consistent": sensitivity_consistent,
        "loro_vs_global_assignment_AMI": loro_ami,
        "shared_positive_global_clusters": [int(v) for v in shared_positive],
        "shared_cluster_claim_evaluable": True,
        "common_cluster_supported": common_cluster_supported,
        "premise_supported": bool(primary_mask.sum() >= 10 and
                                   (rep_df["median"] > 0).all() and
                                   sensitivity_consistent and common_cluster_supported),
        "claim_boundary": "Static-anchor structural direction only; no PAM efficacy or functional-state inference.",
        "note": "Global clusters are unsupervised descriptors; biological validation still needs matched multi-system trajectories.",
    }
    # The preregistered cluster Go condition is deliberately not relaxed post hoc.
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    selected_pairs = pd.DataFrame([{"residue_i": pairs[i][0], "residue_j": pairs[i][1],
                                    "7TRQ_minus_7TRS_A": float((q-s)[i]),
                                    "7TRP_minus_7TRS_A": float((p-s)[i])}
                                   for i in np.where(primary_mask)[0]])
    selected_pairs.to_csv(OUT / "signature_features.csv", index=False)
    report = "# PACER Consensus-PAM Trajectory Signature v0.1\n\n## Replica summary\n\n"
    report += rep_df.to_markdown(index=False, floatfmt=".4f")
    report += "\n\n## Threshold sensitivity\n\n" + pd.DataFrame(sensitivity).to_markdown(index=False, floatfmt=".4f")
    report += "\n\n## Global clusters\n\n" + global_df.to_markdown(index=False, floatfmt=".4f")
    report += "\n\n```json\n" + json.dumps(audit, indent=2) + "\n```\n"
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(rep_df.to_string(index=False))
    print(pd.DataFrame(sensitivity).to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
