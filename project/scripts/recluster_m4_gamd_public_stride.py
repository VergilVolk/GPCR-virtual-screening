# -*- coding: utf-8 -*-
"""Independent clustering of the public 1-ns-stride M4 GaMD trajectories.

This mirrors the paper/repository pocket selection and average-linkage target
of ten clusters, but operates on the 3000 public protein-only frames rather
than the unpublished/full-resolution ~3 million-frame stream used upstream.
"""
from __future__ import annotations

import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from sklearn.cluster import AgglomerativeClustering


P = Path(__file__).resolve().parents[1]
ROOT = P / "data" / "m4_gamd_figshare_33283491"
OUT = P / "results" / "m4_gamd_public_recluster_v01"
REP = OUT / "representatives"
REP.mkdir(parents=True, exist_ok=True)
RESIDS = [687, 690, 691, 694, 695, 696, 706, 781, 782, 783, 784, 785, 788,
          852, 856, 859, 868, 869, 871, 872, 875, 907, 908]
OFFICIAL_FRACTIONS = np.array([.795, .144, .032, .019, .003, .002, .002, .001, .001, .001])


def aligned(mobile, reference, weights):
    w = weights / weights.sum()
    mc = np.sum(mobile * w[:, None], axis=0)
    rc = np.sum(reference * w[:, None], axis=0)
    rot, _ = Rotation.align_vectors(reference - rc, mobile - mc, weights=weights)
    return rot.apply(mobile - mc) + rc, rot, mc, rc


def main():
    traj = sorted(ROOT.glob("run*-GaMD-stride1000ps.nc"))
    u = mda.Universe(str(ROOT / "sys-protein.pdb"), [str(x) for x in traj])
    selection = "resid " + " ".join(map(str, RESIDS)) + " and not name H*"
    pocket = u.select_atoms(selection)
    if pocket.n_atoms != 239 or len(u.trajectory) != 3000:
        raise RuntimeError(f"unexpected public data dimensions: atoms={pocket.n_atoms}, frames={len(u.trajectory)}")
    u.trajectory[0]
    reference = pocket.positions.astype(float).copy()
    weights = pocket.masses.astype(float)
    features = np.empty((len(u.trajectory), pocket.n_atoms * 3), np.float32)
    for i, _ in enumerate(u.trajectory):
        fit, _, _, _ = aligned(pocket.positions.astype(float), reference, weights)
        features[i] = fit.ravel()
        if (i + 1) % 500 == 0:
            print("aligned", i + 1, flush=True)
    np.save(OUT / "aligned_pocket_features.npy", features)

    model = AgglomerativeClustering(n_clusters=10, metric="euclidean", linkage="average")
    raw = model.fit_predict(features)
    counts = pd.Series(raw).value_counts().sort_values(ascending=False)
    mapping = {old: new for new, old in enumerate(counts.index)}
    labels = np.asarray([mapping[x] for x in raw], int)

    rows = []
    for cluster in range(10):
        idx = np.where(labels == cluster)[0]
        centroid = features[idx].mean(axis=0)
        rep_idx = int(idx[np.argmin(np.linalg.norm(features[idx] - centroid, axis=1))])
        rep_rmsd = float(np.linalg.norm(features[rep_idx] - features[0]) / np.sqrt(pocket.n_atoms))
        rows.append({
            "cluster": cluster, "frames": len(idx), "fraction": len(idx) / len(labels),
            "representative_global_frame": rep_idx,
            "representative_time_ps_within_run": int((rep_idx % 500 + 1) * 1000),
            "representative_run": int(rep_idx // 500 + 1),
            "pocket_RMSD_to_first_A": rep_rmsd,
            "official_fullstream_fraction_rankmatched": OFFICIAL_FRACTIONS[cluster],
        })
        u.trajectory[rep_idx]
        mobile = pocket.positions.astype(float)
        _, rot, mc, rc = aligned(mobile, reference, weights)
        original = u.atoms.positions.copy()
        u.atoms.positions = rot.apply(original - mc) + rc
        u.atoms.write(str(REP / f"public_cluster{cluster:02d}_representative.pdb"))
        u.atoms.positions = original

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "cluster_summary.csv", index=False)
    pd.DataFrame({"global_frame": np.arange(len(labels)), "cluster": labels}).to_csv(
        OUT / "frame_assignments.csv", index=False
    )
    audit = {
        "figshare_article": 33283491,
        "public_trajectory_count": len(traj),
        "public_frames": len(labels),
        "stride_ps": 1000,
        "aggregate_public_sampling_ns": 3000,
        "pocket_heavy_atoms": pocket.n_atoms,
        "resids": RESIDS,
        "algorithm": "mass-weighted reference alignment plus average-linkage agglomerative clustering, k=10",
        "dominant_fraction_public": float(table.iloc[0].fraction),
        "dominant_fraction_official_fullstream": float(OFFICIAL_FRACTIONS[0]),
        "fraction_L1_rankmatched": float(np.abs(table.fraction.to_numpy() - OFFICIAL_FRACTIONS).sum()),
        "scope_boundary": (
            "Independent robustness reproduction on public 1-ns-stride coordinates; not an exact replay "
            "of CPPTRAJ sieve-500 clustering over the full-resolution trajectory."
        ),
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(table.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
