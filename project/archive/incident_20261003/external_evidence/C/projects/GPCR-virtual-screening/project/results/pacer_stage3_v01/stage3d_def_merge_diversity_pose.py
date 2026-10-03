#!/usr/bin/env python3
"""Stage 3D - D/E/F: merge with the frozen M4-safe router, scaffold-diverse selection,
and freezing of concrete docking poses for the recommended four-context MD set.

Rules enforced:
  * final_rank is taken verbatim from the frozen router output (2023 M4-LOTO only).
  * structural gate only REMOVES; it never reorders binding priority.
  * score_2026_13t_famaug / rankpct_2026_13t_famaug are supporting-only fields.
  * no weighted fusion, no new score, no refitting.
"""
from __future__ import annotations
import hashlib, json, shutil
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[3]
OUT  = ROOT / "project/results/pacer_stage3_v01"
DOCK = ROOT / "project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01"
MD_N = 3          # mainline doc: "输出 1 至 3 个高优先级候选"

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def scaff(s):
    m = Chem.MolFromSmiles(str(s))
    return MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else ""

router = pd.read_csv(OUT / "pacer200_m4_safe_routed_ranking.csv")
gate   = pd.read_csv(OUT / "pacer200_structural_gate.csv")
led    = pd.read_json(DOCK / "ledger.jsonl", lines=True)

# ---------- D. merge ----------
router["candidate_id"] = router.pair_id.astype(str)   # bundle uses pair_id == PACER candidate id
m = router.merge(gate.drop(columns=["canonical_smiles"]), on="candidate_id", how="left", validate="1:1")
assert len(m) == 200 and m.structural_gate_pass.notna().all()
m["stage3d_status"] = np.where(m.structural_gate_pass, "ELIGIBLE", "EXCLUDED_STRUCTURAL")
m["pacer_binding_score_unchanged_from_router"] = np.isclose(m.pacer_binding_score, m.rankpct_2023_gpcr_loto)
m["final_rank_unchanged_from_router"] = True
prio = m.loc[m.stage3d_status == "ELIGIBLE"].sort_values("final_rank", kind="mergesort")
m["stage3d_priority"] = pd.Series({cid: i for i, cid in enumerate(prio.pair_id, 1)})
m["stage3d_priority"] = m.stage3d_priority.astype("Int64")
m = m.sort_values(["final_rank"], kind="mergesort").reset_index(drop=True)
m.to_csv(OUT / "pacer200_stage3d_merged.csv", index=False)

# verify frozen order untouched
chk = m.sort_values("pair_id").reset_index(drop=True)
r2  = router.sort_values("pair_id").reset_index(drop=True)
order_ok = bool((chk.final_rank.to_numpy() == r2.final_rank.to_numpy()).all()
                and np.allclose(chk.pacer_binding_score, r2.pacer_binding_score)
                and (chk.route == "m4_2023_gpcr_loto").all())

# ---------- E. scaffold diversity (deterministic) ----------
elig = m[m.stage3d_status == "ELIGIBLE"].sort_values(["final_rank", "pair_id"], kind="mergesort").copy()
elig["murcko_scaffold"] = elig.canonical_smiles.map(scaff)
seen, sel = set(), []
for row in elig.itertuples():
    if row.murcko_scaffold in seen:
        continue
    seen.add(row.murcko_scaffold)
    sel.append(row.Index)
div = elig.loc[sel].copy().sort_values("final_rank", kind="mergesort").reset_index(drop=True)
div.insert(0, "diverse_rank", range(1, len(div) + 1))
div["selection_rule"] = "one representative per Murcko scaffold, ordered by frozen M4-safe final_rank"
div["md_set_member"] = div.diverse_rank <= MD_N
div["md_set_reason"] = np.where(div.md_set_member,
    "top-%d scaffold-diverse candidate by frozen M4-safe final_rank" % MD_N, "")
div.to_csv(OUT / "pacer_stage3d_diverse_shortlist.csv", index=False)

# ---------- F. freeze poses for the recommended MD set ----------
mdset = div[div.md_set_member].sort_values("diverse_rank")
posedir = OUT / "md_shortlist_poses"; posedir.mkdir(parents=True, exist_ok=True)
manifest = []
for row in mdset.itertuples():
    g = led[(led.molecule_id == row.candidate_id) &
            (led.error.isna() | led.error.astype(str).str.strip().eq(""))].copy()
    pts = g[["pose_centroid_x", "pose_centroid_y", "pose_centroid_z"]].to_numpy(float)
    med = np.median(pts, axis=0)
    g["dist_to_median_A"] = np.linalg.norm(pts - med, axis=1)
    # ensemble-consistent / plausible pose; deliberately NOT lowest-Vina-first
    g = g.sort_values(["dist_to_median_A", "native_pocket_coverage", "vina_affinity", "cluster"],
                      ascending=[True, False, True, True], kind="mergesort")
    pick = g.iloc[0]
    key  = hashlib.sha1(str(row.candidate_id).encode()).hexdigest()[:12]
    src  = DOCK / "poses" / f"{key}_s{int(pick.state_id):02d}_c{int(pick.cluster):02d}.pdbqt"
    assert src.exists(), src
    dstname = f"{row.candidate_id}_c{int(pick.cluster):02d}_s{int(pick.state_id):02d}.pdbqt"
    dst = posedir / dstname
    shutil.copyfile(src, dst)
    assert sha(src) == sha(dst)
    manifest.append(dict(
        candidate_id=row.candidate_id, canonical_smiles=row.canonical_smiles,
        source_cluster=int(pick.cluster), source_state=int(pick.state_id),
        source_state_smiles=pick.state_smiles,
        pose_file=str(dst.relative_to(ROOT)).replace("\\", "/"),
        pose_source_file=str(src.relative_to(ROOT)).replace("\\", "/"),
        pose_file_sha256=sha(dst), docking_run="pacer_stage3_200_docking_stage3_v01",
        selection_reason=("structural gate PASS; pose nearest the candidate's ensemble-median centroid "
                          "(ensemble-consistent), then highest native pocket coverage, then Vina; "
                          "not selected by lowest Vina alone"),
        router_final_rank=int(row.final_rank), pacer_binding_score=float(row.pacer_binding_score),
        route=row.route, diverse_rank=int(row.diverse_rank), murcko_scaffold=row.murcko_scaffold,
        structural_gate_pass=bool(row.structural_gate_pass),
        structural_evidence_summary=dict(
            cluster_coverage_fraction=float(row.cluster_coverage_fraction),
            median_native_pocket_coverage=float(row.median_native_pocket_coverage),
            median_n_contacted_residues=float(row.median_n_contacted_residues),
            pocket_contact_jaccard_mean=float(row.pocket_contact_jaccard_mean),
            centroid_dispersion_max_A=float(row.centroid_dispersion_max_A),
            vina_affinity=float(pick.vina_affinity),
            native_pocket_coverage=float(pick.native_pocket_coverage),
            n_contacted_residues=int(pick.n_contacted_residues),
            contacted_residues=str(pick.contacted_residues)),
        claim_boundary="Structural plausibility evidence only; not PAM function, efficacy or potency.",
    ))
pd.DataFrame(manifest).to_csv(OUT / "md_shortlist_pose_manifest.csv", index=False)

print(json.dumps({
  "merged_rows": len(m),
  "final_rank_untouched_vs_router": order_ok,
  "eligible": int((m.stage3d_status == "ELIGIBLE").sum()),
  "excluded_structural": int((m.stage3d_status == "EXCLUDED_STRUCTURAL").sum()),
  "router_top50_excluded_by_structure": int(((m.final_rank <= 50) & (~m.structural_gate_pass)).sum()),
  "router_top25_excluded_by_structure": int(((m.final_rank <= 25) & (~m.structural_gate_pass)).sum()),
  "n_distinct_scaffolds_eligible": int(div.murcko_scaffold.nunique()),
  "diverse_shortlist_rows": len(div),
  "recommended_md_set": mdset.candidate_id.tolist(),
  "md_set_poses": [f"{r['candidate_id']} -> cluster {r['source_cluster']} / {Path(r['pose_file']).name}" for r in manifest],
  "outputs": {n: sha(OUT / n) for n in
      ["pacer200_stage3d_merged.csv","pacer_stage3d_diverse_shortlist.csv","md_shortlist_pose_manifest.csv"]},
}, indent=2))
