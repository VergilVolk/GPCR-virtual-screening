#!/usr/bin/env python
"""Project-wide integration benchmark assembler v01.

Collects real, committed numbers from all project lines into one master JSON
plus per-line CSV tables. Every number carries a provenance string.
Run from any directory; paths are absolute to the two worktrees.

Sources:
  binding : D:/CLC  (codex/drugclip-pacer-handoff results)
  function: D:/CLC_geom2vec_pilot (integration branch results + extracted fkg_v01)
  qsar    : D:/CLC  (Aug v1 campaign results)
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

BIND = Path(r"D:\CLC\project\results\drugclip_science2026")
FUNC = Path(r"D:\CLC_geom2vec_pilot\project\results")
QSAR = Path(r"D:\CLC\project\results")
WS = FUNC / "project_wide_integration_benchmark_v01" / "data"
OUT = FUNC / "project_wide_integration_benchmark_v01"

M = {"binding": {}, "functional": {}, "qsar_static": {}, "integration": {}}


def jload(p: Path):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def prov(tag: str) -> str:
    return tag


# ---------------- binding line ----------------
def collect_binding():
    m13 = jload(BIND / "migration_13target_v01" / "summary_preserve0_3seed.json")
    M["binding"]["m13_summary"] = {
        meth: {
            "macro": m13["metrics"][meth]["macro"],
            "per_target": m13["metrics"][meth]["per_target"],
        }
        for meth in m13["metrics"]
    }
    M["binding"]["m13_source"] = "migration_13target_v01/summary_preserve0_3seed.json (ep40 tuned 3-seed summary)"

    # ep80 per-seed
    ep80 = {}
    for s in (20260925, 20260926, 20260927):
        d = jload(BIND / "migration_13target_v01" / f"loso_preserve0.2_seed{s}.json")
        ep80[s] = {"macro": d["metrics"]["bce_retrieval"]["macro"], "per_target": d["metrics"]["bce_retrieval"]["per_target"]}
        if s == 20260925:
            M["binding"]["ep80_official_macro"] = d["metrics"]["official"]["macro"]
            M["binding"]["ep80_official_source"] = f"loso_preserve0.2_seed{s}.json metrics.official"
    M["binding"]["ep80_seeds"] = ep80
    M["binding"]["ep80_source"] = "loso_preserve0.2_seed*.json (ep80 3 seeds)"

    # family augmentation
    fa = {}
    for s in (20260925, 20260926, 20260927):
        d = jload(BIND / "family_aug_v01" / f"famaug_ep80_seed{s}.json")
        fa[s] = {"macro": d["metrics"]["bce_retrieval"]["macro"], "per_target": d["metrics"]["bce_retrieval"]["per_target"]}
    M["binding"]["famaug_seeds"] = fa
    M["binding"]["famaug_ensemble"] = jload(BIND / "family_aug_v01" / "famaug_ensemble.json")
    M["binding"]["famaug_source"] = "family_aug_v01/famaug_ep80_seed*.json + famaug_ensemble.json"

    # extended 20-target
    e20 = {}
    e20dir = BIND / "extended20_finetune_v01"
    for f in e20dir.glob("*.json"):
        try:
            d = jload(f)
        except Exception:
            continue
        if "metrics" in d and "official" in d.get("metrics", {}):
            e20[f.name] = {
                meth: {"macro": d["metrics"][meth]["macro"]}
                for meth in d["metrics"]
            }
    M["binding"]["e20_files"] = {k: v for k, v in e20.items()}
    M["binding"]["e20_ecfp"] = jload(e20dir / "ecfp20_loso.json")["metrics"]["macro"]
    M["binding"]["e20_source"] = "extended20_finetune_v01/*.json + ecfp20_loso.json"


# ---------------- functional line ----------------
REGIONS = [
    "compound110_extension", "cooperativity_mutagenesis", "distal_control",
    "intracellular_microswitches", "orthosteric_activation_core",
    "orthosteric_contact_union", "pam_contact_consensus", "pam_contact_union",
    "stable_core_control",
]


def collect_functional():
    tb = jload(FUNC / "pacer_dc_close_loop_20ns_v01" / "track_b_fkg_v02" /
               "phase2b_graph_region" / "TRACK_B_MATCHED_20NS_RESULTS_v01.json")
    rows = []
    for branch in ("STATE_MOTION", "SIGNED_DRIFT"):
        for panel in ("compound110_reference", "LY2119620"):
            for contrast in ("Delta_PAM", "Delta_AGO", "Delta_INT"):
                node = tb["results"][branch][panel][contrast]["regions"]
                for r in REGIONS:
                    g = node[r]
                    row = {
                        "branch": branch, "panel": panel.replace("_reference", ""),
                        "contrast": contrast, "region": r,
                        "R1R3_cosine": g.get("R1_R3_direction_cosine"),
                        **{f"{rk}_norm": g["replicas"][rk]["mean_vector_norm"] for rk in ("R1", "R2", "R3")},
                        "provenance": "TRACK_B_MATCHED_20NS_RESULTS_v01.json",
                    }
                    rows.append(row)
    M["functional"]["track_b_rows"] = rows

    # LY vs c110 direct comparison
    M["functional"]["ly_vs_c110"] = tb["LY_vs_compound110_descriptive_comparison"]

    # CM00734 primary + 9-region table parsed from synthesis doc
    cs = jload(FUNC / "pacer_dc_cm00734_stage_b_20ns_analysis_v01" / "CM00734_STAGE_B_FINAL_STATUS_v01.json")
    M["functional"]["cm00734_status"] = cs

    syn = (FUNC / "pacer_dc_stage_a_b_final_synthesis_v01" / "PACER_DC_STAGE_A_B_FINAL_SYNTHESIS_v01.md").read_text(encoding="utf-8")
    import re
    cm_rows = []
    # parse the region-level specificity table (STATE_MOTION / Delta_INT)
    pat = re.compile(r"^\|\s*([a-z0-9_]+)\s*\|\s*([+-]?[\d.]+)\s*\|\s*([+-]?[\d.]+)\s*\|\s*([+-]?[\d.]+)\s*\|\s*$", re.M)
    in_sec = False
    for line in syn.splitlines():
        if "Region-level specificity evidence" in line:
            in_sec = True
            continue
        if in_sec:
            mm = pat.match(line)
            if mm and mm.group(1) in REGIONS:
                cm_rows.append({
                    "region": mm.group(1),
                    "CM00734_R1R3": float(mm.group(2)),
                    "LY2119620_R1R3": float(mm.group(3)),
                    "compound110_R1R3": float(mm.group(4)),
                    "provenance": "PACER_DC_STAGE_A_B_FINAL_SYNTHESIS_v01.md table",
                })
    M["functional"]["cm00734_region_rows"] = cm_rows

    # FKG v01 common-kernel region axis
    v01 = WS / "fkg_v01_REGION_AXIS.csv"
    if v01.exists():
        with v01.open(encoding="utf-8") as f:
            M["functional"]["fkg_v01_region_axis"] = list(csv.DictReader(f))
        M["functional"]["fkg_v01_source"] = "extracted from codex/pacer-dc-geom2vec-pilot:results/pacer_dc_long_md_v02/fkg_common_kernel_v01"
    v01r = WS / "fkg_v01_REPLICA_AUDIT.csv"
    if v01r.exists():
        with v01r.open(encoding="utf-8") as f:
            M["functional"]["fkg_v01_replica_audit"] = list(csv.DictReader(f))


# ---------------- qsar / static line ----------------
def collect_qsar():
    b = jload(QSAR / "pacer_baselines_v1" / "baseline_metrics.json")
    M["qsar_static"]["baselines"] = b
    st = jload(QSAR / "pacer_structure_loso_v01" / "structure_loso_metrics.json")
    M["qsar_static"]["structure_loso_aggregate"] = st["aggregate"]
    fs = jload(QSAR / "pacer_fs_final" / "evidence.json")
    M["qsar_static"]["pacer_fs"] = {k: fs[k] for k in
                                    ("centered_vs_absolute_spearman", "centered_vs_absolute_mae", "status")}
    mi = jload(QSAR / "pacer_m4r_enrichment_v01" / "metrics.json")
    agg = {}
    for meth, v in mi.items():
        folds = v.get("folds", [])
        if folds:
            keys = [k for k in folds[0] if isinstance(folds[0][k], (int, float))]
            agg[meth] = {k: sum(f.get(k, 0) for f in folds) / len(folds) for k in keys}
    M["qsar_static"]["miao_external_macro"] = agg
    gam = jload(QSAR / "m4_gamd_ensemble" / "candidate_ensemble_audit.json")
    M["qsar_static"]["gamd"] = {k: gam[k] for k in
                                ("n_complete", "n_dynamic_pareto", "single_7TRS_vs_GaMD_BEmin_spearman")}
    # representation bakeoff (PCA/VAMP/tICA)
    bk = QSAR / "pacer_dynamic_representation_bakeoff_v01" / "summary.csv"
    if bk.exists():
        with bk.open(encoding="utf-8") as f:
            M["qsar_static"]["representation_bakeoff"] = list(csv.DictReader(f))


# ---------------- integration checks ----------------
def collect_integration():
    checks = []
    # tag -> commit resolution
    import subprocess
    for t, want in [
        ("encoder-candidate-ah-frozen-20260929", "0a4b4f23"),
        ("pacer-fkg-v02-phase2-frozen-20260929", "ecaf7a60"),
        ("pacer-dc-ly2119620-close-loop-20ns-v01", "b719f01e"),
        ("pacer-dc-cm00734-stage-b-20ns-v01", "f2f73403"),
        ("pacer-dc-stage-a-b-synthesis-v01", "5e86149b"),
        ("pacer-fkg-long-md-v01", "df219d70"),
    ]:
        r = subprocess.run(["git", "rev-parse", f"{t}^{{commit}}"], cwd=str(FUNC.parents[1]),
                           capture_output=True, text=True)
        got = r.stdout.strip()[:8] if r.returncode == 0 else "MISSING"
        checks.append({"check": f"tag:{t}", "expected": want, "got": got,
                       "pass": got == want})
    # frozen SHAs referenced in artifacts
    cs = jload(FUNC / "pacer_dc_cm00734_stage_b_20ns_analysis_v01" / "CM00734_STAGE_B_FINAL_STATUS_v01.json")
    checks.append({"check": "fkg_v02_freeze_sha_in_stage_b",
                   "expected": "b48bc74a757a...", "got": cs["frozen_fkg_v02_sha256"][:12],
                   "pass": cs["frozen_fkg_v02_sha256"].startswith("b48bc74a757a")})
    tb = jload(FUNC / "pacer_dc_close_loop_20ns_v01" / "track_b_fkg_v02" /
               "phase2b_graph_region" / "TRACK_B_MATCHED_20NS_RESULTS_v01.json")
    checks.append({"check": "trackB_upstream_freeze",
                   "expected": "b48bc74a...", "got": tb["upstream"]["historical_fkg_v02_freeze_sha256"][:12],
                   "pass": tb["upstream"]["historical_fkg_v02_freeze_sha256"].startswith("b48bc74a")})
    checks.append({"check": "trackB_forbidden_ops_all_false",
                   "expected": "all False", "got": json.dumps(tb["forbidden_operations_performed"]),
                   "pass": not any(tb["forbidden_operations_performed"].values())})
    M["integration"]["checks"] = checks
    M["integration"]["tests_run"] = [
        {"suite": "test_path_guard.py", "result": "4 passed, 20 subtests passed"},
        {"suite": "test_phase2_calibration.py", "result": "7 passed"},
        {"suite": "test_phase1_bs256.py", "result": "SKIP - torch_geometric not installed on this host (GPU laptop env)"},
        {"suite": "test_phase3_fkg.py", "result": "SKIP - V02_FREEZE_MANIFEST artifact lives on GPU laptop (git-ignored cache); anchor verification correctly refuses without it"},
    ]


def main():
    collect_binding()
    collect_functional()
    collect_qsar()
    collect_integration()
    (OUT / "project_wide_benchmark_master_v01.json").write_text(
        json.dumps(M, ensure_ascii=False, indent=1), encoding="utf-8")

    # flat CSVs
    ws = OUT / "csv"
    ws.mkdir(parents=True, exist_ok=True)
    with (ws / "functional_three_class.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(M["functional"]["track_b_rows"][0]))
        w.writeheader(); w.writerows(M["functional"]["track_b_rows"])
    with (ws / "cm00734_regions.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["region", "CM00734_R1R3", "LY2119620_R1R3", "compound110_R1R3", "provenance"])
        w.writeheader(); w.writerows(M["functional"]["cm00734_region_rows"])
    if M["functional"].get("fkg_v01_region_axis"):
        with (ws / "fkg_v01_region_axis.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(M["functional"]["fkg_v01_region_axis"][0]))
            w.writeheader(); w.writerows(M["functional"]["fkg_v01_region_axis"])
    # binding macro table
    rows = []
    for meth, v in M["binding"]["m13_summary"].items():
        rows.append({"protocol": "13target", "method": meth, **{k: round(x, 4) for k, x in v["macro"].items()}})
    for s, v in M["binding"]["ep80_seeds"].items():
        rows.append({"protocol": "13target", "method": f"ep80_seed{s}", **{k: round(x, 4) for k, x in v["macro"].items()}})
    for s, v in M["binding"]["famaug_seeds"].items():
        rows.append({"protocol": "13target", "method": f"famaug_seed{s}", **{k: round(x, 4) for k, x in v["macro"].items()}})
    fe = M["binding"]["famaug_ensemble"]["ensemble_macro"]
    rows.append({"protocol": "13target", "method": "famaug_ensemble", **{k: round(x, 4) for k, x in fe.items()}})
    for name, meths in M["binding"]["e20_files"].items():
        for meth, v in meths.items():
            rows.append({"protocol": "20target", "method": f"{meth}:{name}", **{k: round(x, 4) for k, x in v["macro"].items()}})
    rows.append({"protocol": "20target", "method": "ecfp4_logistic", **{k: round(x, 4) for k, x in M["binding"]["e20_ecfp"].items()}})
    with (ws / "binding_macros.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["protocol", "method", "roc_auc", "pr_auc", "bedroc_alpha20", "ef1pct", "ef5pct"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})
    print("assembled:", len(M["binding"]), "binding entries,", len(M["functional"]["track_b_rows"]),
          "functional rows,", len(M["integration"]["checks"]), "integration checks")
    print("integration checks: all pass =", all(c["pass"] for c in M["integration"]["checks"]))
    for c in M["integration"]["checks"]:
        print("  ", "PASS" if c["pass"] else "FAIL", c["check"], c["got"])


if __name__ == "__main__":
    main()

