#!/usr/bin/env python3
"""Stage 3D - C: structural plausibility gate from real Stage 3A docking evidence.

Docking here is a STRUCTURAL PLAUSIBILITY GATE only. Vina values are used as
structural support evidence and are never interpreted as PAM efficacy/potency.

Pre-declared thresholds (fixed before inspecting per-candidate outcomes):
  G1 cluster_coverage_fraction      >= 0.80   (>= 8/10 M4R GaMD conformations gave a valid pose)
  G2 median native_pocket_coverage  >= 0.50   (>= 7/14 declared allosteric-pocket residues contacted)
  G3 median n_contacted_residues    >= 8      (key-contact evidence present)
  G4 max centroid outlier distance  <= 20.0 A (no implausible pose far outside the ensemble region)
"""
from __future__ import annotations
import hashlib, itertools, json
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DOCK = ROOT / "project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01"
from identity_contract import validate_inputs, graph
ROUTER_PATH = ROOT / "project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv"
validate_inputs(ROOT, DOCK, ROUTER_PATH)
OUT  = ROOT / "project/results/pacer_stage3d_corrected_20261002_v01"
TH = dict(g1_cluster_coverage=0.80, g2_median_pocket_coverage=0.50,
          g3_median_contacts=8, g4_max_outlier_distance_A=20.0)

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

led = pd.read_json(DOCK / "ledger.jsonl", lines=True)
ok = led[led.error.isna() | led.error.astype(str).str.strip().eq("")].copy()
port = pd.read_csv(ROOT / "project/results/pacer_candidates_v01/predock_portfolio.csv")

xyz = ok[["pose_centroid_x", "pose_centroid_y", "pose_centroid_z"]].to_numpy(float)
ok["_res"] = ok.contacted_residues.fillna("").astype(str)

rows = []
for cid, g in ok.groupby("molecule_id"):
    g = g.sort_values("cluster")
    pts = g[["pose_centroid_x", "pose_centroid_y", "pose_centroid_z"]].to_numpy(float)
    med = np.median(pts, axis=0)
    dist = np.linalg.norm(pts - med, axis=1)
    sets = [set(s.split(";")) - {""} for s in g._res]
    jac = [len(a & b) / len(a | b) for a, b in itertools.combinations(sets, 2) if (a | b)] or [np.nan]
    v = g.vina_affinity.to_numpy(float)
    rows.append(dict(
        candidate_id=cid,
        n_valid_states=int(g.state_id.nunique()),
        n_clusters_with_valid_pose=int(g.cluster.nunique()),
        cluster_coverage_fraction=g.cluster.nunique() / 10.0,
        best_vina=float(v.min()), median_vina=float(np.median(v)),
        mean_vina=float(v.mean()), vina_std=float(v.std(ddof=0)),
        median_native_pocket_coverage=float(np.median(g.native_pocket_coverage)),
        min_native_pocket_coverage=float(g.native_pocket_coverage.min()),
        median_n_contacted_residues=float(np.median(g.n_contacted_residues)),
        max_n_contacted_residues=int(g.n_contacted_residues.max()),
        centroid_dispersion_mean_A=float(dist.mean()),
        centroid_dispersion_max_A=float(dist.max()),
        pocket_contact_jaccard_mean=float(np.nanmean(jac)),
        microstate_sensitivity="not_available (single pH7 microstate per candidate)",
    ))
st = pd.DataFrame(rows).merge(port[["candidate_id", "canonical_smiles"]], on="candidate_id", how="left")
st = st.sort_values("candidate_id").reset_index(drop=True)

st["gate_g1_cluster_coverage"]  = st.cluster_coverage_fraction >= TH["g1_cluster_coverage"]
st["gate_g2_pocket_occupancy"]  = st.median_native_pocket_coverage >= TH["g2_median_pocket_coverage"]
st["gate_g3_key_contacts"]      = st.median_n_contacted_residues >= TH["g3_median_contacts"]
st["gate_g4_no_catastrophic_outlier"] = st.centroid_dispersion_max_A <= TH["g4_max_outlier_distance_A"]
st["structural_gate_pass"] = st[["gate_g1_cluster_coverage","gate_g2_pocket_occupancy",
                                "gate_g3_key_contacts","gate_g4_no_catastrophic_outlier"]].all(axis=1)
def reason(r):
    fails = []
    if not r.gate_g1_cluster_coverage: fails.append("G1_low_cluster_coverage")
    if not r.gate_g2_pocket_occupancy: fails.append("G2_pocket_occupancy_below_0.5")
    if not r.gate_g3_key_contacts: fails.append("G3_median_contacts_below_8")
    if not r.gate_g4_no_catastrophic_outlier: fails.append("G4_implausible_outlier_pose")
    return "PASS" if not fails else ";".join(fails)
st["structural_gate_reason"] = st.apply(reason, axis=1)

cols = ["candidate_id","canonical_smiles","n_valid_states","n_clusters_with_valid_pose",
        "cluster_coverage_fraction","best_vina","median_vina","mean_vina","vina_std",
        "median_native_pocket_coverage","min_native_pocket_coverage","median_n_contacted_residues",
        "max_n_contacted_residues","centroid_dispersion_mean_A","centroid_dispersion_max_A",
        "pocket_contact_jaccard_mean","microstate_sensitivity",
        "gate_g1_cluster_coverage","gate_g2_pocket_occupancy","gate_g3_key_contacts",
        "gate_g4_no_catastrophic_outlier","structural_gate_pass","structural_gate_reason"]
st[cols].to_csv(OUT / "pacer200_structural_gate.csv", index=False)
print(json.dumps({
  "thresholds": TH,
  "n_candidates": len(st),
  "structural_gate_pass": int(st.structural_gate_pass.sum()),
  "structural_gate_fail": int((~st.structural_gate_pass).sum()),
  "fail_by_reason": st.loc[~st.structural_gate_pass, "structural_gate_reason"].value_counts().to_dict(),
  "gate_component_pass": {c: int(st[c].sum()) for c in
      ["gate_g1_cluster_coverage","gate_g2_pocket_occupancy","gate_g3_key_contacts","gate_g4_no_catastrophic_outlier"]},
  "output_sha256": sha(OUT / "pacer200_structural_gate.csv"),
}, indent=2))
print()
print("=== observed ranges (context for the thresholds) ===")
for c in ["cluster_coverage_fraction","median_native_pocket_coverage","median_n_contacted_residues",
          "centroid_dispersion_max_A","pocket_contact_jaccard_mean","vina_std"]:
    s = st[c]
    print(f"  {c:32s} min={s.min():.4f} p25={s.quantile(.25):.4f} med={s.median():.4f} p75={s.quantile(.75):.4f} max={s.max():.4f}")
