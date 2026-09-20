"""Test whether a third-PAM GaMD trajectory approaches two independent PAM-bound anchors."""

from __future__ import annotations

import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "data" / "m4_gamd_figshare_33283491"
STRUCT = PROJECT / "results" / "structure" / "ensemble"
OUT = PROJECT / "results" / "pacer_static_anchored_trajectory_v01"
PUBLIC_RESIDS = [687, 690, 691, 694, 695, 696, 706, 781, 782, 783, 784, 785, 788,
                 852, 856, 859, 868, 869, 871, 872, 875]
STATIC = {"7TRQ_VU0467154": STRUCT / "7TRQ_R.pdb",
          "7TRP_LY2033298": STRUCT / "7TRP_R.pdb",
          "7TRS_ACh_only": STRUCT / "7TRS_R.pdb"}


def static_atoms(path):
    atoms = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[21] != "R": continue
        resid = int(line[22:26]); name = line[12:16].strip()
        element = line[76:78].strip().upper()
        if resid not in {x - 598 for x in PUBLIC_RESIDS} or element == "H" or name.startswith("H"): continue
        atoms[(resid, name)] = np.asarray([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return atoms


def align(mobile, reference):
    mc, rc = mobile.mean(0), reference.mean(0)
    rot, _ = Rotation.align_vectors(reference - rc, mobile - mc)
    return rot.apply(mobile - mc) + rc


def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    trajectories = sorted(ROOT.glob("run*-GaMD-stride1000ps.nc"))
    u = mda.Universe(str(ROOT / "sys-protein.pdb"), [str(x) for x in trajectories])
    public = {}
    for atom in u.atoms:
        if atom.resid in PUBLIC_RESIDS and not atom.name.startswith("H"):
            public[(int(atom.resid - 598), atom.name)] = atom.index
    anchors_raw = {name: static_atoms(path) for name, path in STATIC.items()}
    common = sorted(set(public).intersection(*(set(x) for x in anchors_raw.values())))
    if len(common) < 100: raise RuntimeError(f"too few common pocket atoms: {len(common)}")
    indices = np.asarray([public[k] for k in common], int)
    u.trajectory[0]; reference = u.atoms[indices].positions.astype(float).copy()
    anchors = {name: align(np.asarray([atoms[k] for k in common]), reference)
               for name, atoms in anchors_raw.items()}
    q, p, s = anchors["7TRQ_VU0467154"], anchors["7TRP_LY2033298"], anchors["7TRS_ACh_only"]
    static_distances = {"7TRQ_vs_7TRP": rmsd(q, p), "7TRQ_vs_7TRS": rmsd(q, s),
                        "7TRP_vs_7TRS": rmsd(p, s)}
    rows = []
    for global_frame, _ in enumerate(u.trajectory):
        frame = align(u.atoms[indices].positions.astype(float), reference)
        dq, dp, ds = rmsd(frame, q), rmsd(frame, p), rmsd(frame, s)
        rows.append({"global_frame": global_frame, "replica": global_frame // 500 + 1,
                     "frame_within_replica": global_frame % 500 + 1,
                     "RMSD_7TRQ": dq, "RMSD_7TRP": dp, "RMSD_7TRS": ds,
                     "PAM_likeness": ds - .5 * (dq + dp)})
    frame = pd.DataFrame(rows)
    clusters = pd.read_csv(PROJECT / "results" / "m4_gamd_public_recluster_v01" / "frame_assignments.csv")
    frame = frame.merge(clusters, on="global_frame", validate="one_to_one")
    frame.to_csv(OUT / "frame_anchor_scores.csv", index=False)
    replica = frame.groupby("replica").agg(mean_PAM_likeness=("PAM_likeness", "mean"),
        median_PAM_likeness=("PAM_likeness", "median"), fraction_PAM_like=("PAM_likeness", lambda x: float((x > 0).mean()))).reset_index()
    replica.to_csv(OUT / "replica_summary.csv", index=False)
    cluster = frame.groupby("cluster").agg(frames=("global_frame", "count"),
        replicas=("replica", "nunique"), mean_PAM_likeness=("PAM_likeness", "mean"),
        median_PAM_likeness=("PAM_likeness", "median"), fraction_PAM_like=("PAM_likeness", lambda x: float((x > 0).mean()))).reset_index()
    cluster.to_csv(OUT / "cluster_anchor_summary.csv", index=False)
    static_coherent = static_distances["7TRQ_vs_7TRP"] < min(static_distances["7TRQ_vs_7TRS"], static_distances["7TRP_vs_7TRS"])
    shared_cluster = cluster[(cluster.replicas == 6) & (cluster.median_PAM_likeness > 0)]
    audit = {"common_protein_pocket_atoms": len(common), "static_RMSD_A": static_distances,
             "pam_anchor_geometry_coherent": bool(static_coherent),
             "all_six_replica_median_positive": bool((replica.median_PAM_likeness > 0).all()),
             "shared_positive_clusters": shared_cluster.cluster.astype(int).tolist(),
             "premise_supported": bool(static_coherent and (replica.median_PAM_likeness > 0).all() and len(shared_cluster) > 0),
             "claim_boundary": "Three-static-anchor structural hypothesis test; not functional-state or efficacy validation."}
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    report = ["# PACER Static-Anchored Trajectory Test v0.1", "", "## Replica summary", "",
              replica.to_markdown(index=False, floatfmt=".4f"), "", "## Cluster summary", "",
              cluster.to_markdown(index=False, floatfmt=".4f"), "", "```json", json.dumps(audit, indent=2), "```"]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(replica.to_string(index=False)); print(cluster.to_string(index=False)); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
