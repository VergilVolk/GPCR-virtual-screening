# -*- coding: utf-8 -*-
"""Match public-stride recluster representatives to released full-stream clusters."""
from __future__ import annotations

import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation


P = Path(__file__).resolve().parents[1]
OFFICIAL = P / "tools" / "gpcr-am-ensemble-docking" / "ensembles" / "M4R" / "M4R_ensemble_clusters.pdb"
ROOT = P / "results" / "m4_gamd_public_recluster_v01"
RESIDS = [687, 690, 691, 694, 695, 696, 706, 781, 782, 783, 784, 785, 788,
          852, 856, 859, 868, 869, 871, 872, 875, 907, 908]
SEL = "resid " + " ".join(map(str, RESIDS)) + " and not name H*"


def keyed(group):
    return {(int(a.resid), str(a.name)): a.position.astype(float).copy() for a in group}


def rmsd(a, b):
    ac = a - a.mean(0); bc = b - b.mean(0)
    rot, _ = Rotation.align_vectors(ac, bc)
    fit = rot.apply(bc)
    return float(np.sqrt(np.mean(np.sum((ac - fit) ** 2, axis=1))))


def main():
    off = mda.Universe(str(OFFICIAL))
    official = []
    keys = None
    for ts in off.trajectory:
        k = keyed(off.select_atoms(SEL))
        keys = sorted(k) if keys is None else keys
        official.append(np.asarray([k[x] for x in keys]))
    public = []
    files = sorted((ROOT / "representatives").glob("public_cluster*_representative.pdb"))
    for path in files:
        u = mda.Universe(str(path))
        k = keyed(u.select_atoms(SEL))
        if sorted(k) != keys:
            raise RuntimeError(f"atom-key mismatch: {path.name}")
        public.append(np.asarray([k[x] for x in keys]))
    matrix = np.asarray([[rmsd(a, b) for b in official] for a in public])
    rr, cc = linear_sum_assignment(matrix)
    rows = [{
        "public_cluster": int(i), "official_cluster": int(j), "pocket_RMSD_A": float(matrix[i, j])
    } for i, j in zip(rr, cc)]
    pd.DataFrame(matrix, columns=[f"official_{i:02d}" for i in range(10)]).to_csv(
        ROOT / "representative_pairwise_RMSD.csv", index_label="public_cluster"
    )
    matched = pd.DataFrame(rows).sort_values("public_cluster")
    matched.to_csv(ROOT / "representative_hungarian_match.csv", index=False)
    audit = {
        "matched_mean_RMSD_A": float(matched.pocket_RMSD_A.mean()),
        "matched_median_RMSD_A": float(matched.pocket_RMSD_A.median()),
        "matched_max_RMSD_A": float(matched.pocket_RMSD_A.max()),
        "matches_below_2A": int((matched.pocket_RMSD_A < 2.0).sum()),
        "n_clusters": 10,
        "interpretation": "Lower values indicate public-stride clusters recover released full-stream pocket states.",
    }
    (ROOT / "representative_match_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(matched.to_string(index=False))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
