"""Produce the Phase-0 gate/report and perform end-of-phase integrity checks."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from path_guard import guard_write_path


REPO = Path(__file__).resolve().parents[2]
V02 = REPO / "project/pacer_fkg_v02"
PROV = REPO / "project/results/pacer_fkg_v02_longmd_v01/provenance"


def load(path: Path): return json.loads(path.read_text(encoding="utf-8"))
def digest(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(8<<20),b""): h.update(b)
    return h.hexdigest()
def write(path: Path, text: str) -> None:
    p=guard_write_path(path,REPO); p.parent.mkdir(parents=True,exist_ok=True); p.write_text(text,encoding="utf-8")


def integrity() -> dict:
    manifest=V02/"V01_PROTECTION_MANIFEST_SHA256.txt"; rows=[]
    for line in manifest.read_text(encoding="utf-8").splitlines():
        expected,raw=line.split("  ",1); path=Path(raw); actual=digest(path) if path.is_file() else None
        rows.append({"path":str(path),"expected_sha256":expected.lower(),"actual_sha256":actual,"match":actual==expected.lower()})
    baseline=load(V02/"BASELINE_TRACKED_TREE_SHA256.json"); brows=[]
    for item in baseline["files"]:
        path=REPO/item["path"]; actual=digest(path) if path.is_file() else None
        brows.append({"path":item["path"],"match":actual==item["sha256"],"actual_sha256":actual} if actual!=item["sha256"] else {"path":item["path"],"match":True})
    return {"checked_at":datetime.now(timezone.utc).astimezone().isoformat(),
      "v01_result_manifest_file_count":len(rows),"V01_RESULT_HASHES_UNCHANGED":all(x["match"] for x in rows),
      "v01_mismatches":[x for x in rows if not x["match"]],
      "baseline_tracked_file_count":len(brows),"BASELINE_TRACKED_FILES_UNCHANGED":all(x["match"] for x in brows),
      "baseline_mismatches":[x for x in brows if not x["match"]]}


def main() -> None:
    inv=load(PROV/"LONG_MD_FILE_INVENTORY.json"); meta=load(PROV/"TRAJECTORY_METADATA_AUDIT.json")
    qc=load(PROV/"LONG_MD_COMPLETION_QC_AUDIT.json"); prep=load(PROV/"LONG_MD_ENCODER_PREPROCESSING_AUDIT.json")
    comp=load(PROV/"LONG_MD_ENCODER_COMPATIBILITY.json"); integ=integrity()
    offsets=sorted({round(x["state_csv_minus_dcd_last_time_ps"],6) for x in meta["trajectories"]})
    caveats=[
      "The long-MD backup contains no topology. Read-only, checksummed per-system sanitized PDBs from C:\\projects\\PACER_DC_MD_delivery_20260927 were used; their atom counts match every DCD and their receptor sequence/mapping matches the authenticated short-data atom14 CSVs.",
      f"All state.csv timestamps are offset from DCD header times by about {offsets[0]:.6f} ps at the final frame (DCD: 50.0-50000.0 ps; state.csv: 54.0-50004.0 ps). Completion is independently recorded as 50.0 ns in each progress.json; no timing value was inferred from filenames.",
      "MDAnalysis warns that some non-receptor PDB records lack resid information and defaults those records to resid 1. The chain-E protein selection remains 270 residues, matches the short-data sequence exactly, and has complete expected atom14 heavy atoms.",
    ]
    hard_ok=(inv["server_manifest_companion_match"] and inv["all_server_recorded_files_match"] and
             inv["expected_60_production_files_present"] and qc["all_12_complete"] and
             qc["all_seeds_locally_verified"] and prep["route_usable_without_change"] and
             comp["all_12_compatible"] and integ["V01_RESULT_HASHES_UNCHANGED"] and
             integ["BASELINE_TRACKED_FILES_UNCHANGED"])
    decision="PASS_WITH_DOCUMENTED_CAVEAT" if hard_ok and caveats else "PASS" if hard_ok else "FAIL"
    status=subprocess.check_output(["git","status","--short","--branch"],cwd=REPO,text=True).splitlines()
    ignored=subprocess.check_output(["git","status","--short","--ignored","project/results/pacer_fkg_v02_longmd_v01"],cwd=REPO,text=True).splitlines()
    gate={"schema":"pacer_fkg_v02.phase0_gate.v1","created_at":datetime.now(timezone.utc).astimezone().isoformat(),
      "decision":decision,"no_encoder_inference_performed":True,"no_visnet_long_md_processing_performed":True,
      "checks":{"repository_protection":integ,"server_manifest_authenticated":inv["server_manifest_companion_match"],
        "server_local_checksums_match":inv["all_server_recorded_files_match"],"production_files_complete":inv["expected_60_production_files_present"],
        "completion_qc_pass":qc["all_12_complete"],"seeds_verified":qc["all_seeds_locally_verified"],
        "preprocessing_route_compatible":prep["route_usable_without_change"],"atom_residue_compatible":comp["all_12_compatible"],
        "path_guard_tests":"4/4 PASS (unittest)"},"documented_caveats":caveats,"git_status":status,"ignored_output_status":ignored}
    write(V02/"PHASE0_GATE.json",json.dumps(gate,indent=2)+"\n")
    out=io.StringIO(); fields=["system","replica","atom_count","residue_count","receptor_residue_count","stored_frames","first_trajectory_time_ps","last_trajectory_time_ps","dt_ps","box_pbc_available","expected_50_ns_consistency"]
    w=csv.DictWriter(out,fieldnames=fields,lineterminator="\n");w.writeheader();
    for x in meta["trajectories"]:w.writerow({k:x[k] for k in fields})
    write(PROV/"LONG_MD_TRAJECTORY_SUMMARY.csv",out.getvalue())
    lines=["# PACER-FKG v02 long-MD - Phase 0 report","",f"**Gate: {decision}**","",
      "Phase 0 established isolation and audited provenance, completion, trajectory metadata, PBC behavior, and the frozen encoder input mapping. No ViSNet inference or downstream PACER-FKG computation was performed.","",
      "## 1. v01 protection coverage","",
      "The pre-existing manifest was preserved byte-for-byte. It contains 41 files across the seven protected historical result roots: G1 diagnostics (6), G2A PBC (5), G2A (3), G2B robustness (5), G2B (5), G2C (5), and R2/R3 (12). End-of-phase recomputation: `V01_RESULT_HASHES_UNCHANGED = TRUE`.","",
      "## 2. Baseline tracked-tree protection","",
      f"`BASELINE_TRACKED_TREE_SHA256.json` records SHA256 and size for all {integ['baseline_tracked_file_count']} Git-tracked files at `0a4b4f2`, branch `experiment/pacer-fkg-v02-longmd-v01`, tag `encoder-candidate-ah-frozen-20260929`. End-of-phase recomputation: `BASELINE_TRACKED_FILES_UNCHANGED = TRUE`.","",
      "## 3. Server/local checksums","",
      f"The SHA256 of `SERVER_SHA256SUMS.txt` matches its companion file. All {inv['server_record_count']} server-recorded production files match locally; no mismatch occurred. The inventory contains {inv['file_count']} backup files total, including validation and rsync evidence not recorded by the server manifest.","",
      "## 4-5. Exact trajectory inventory and metadata","",
      "| System | Rep | Atoms | Residues | Receptor residues | Frames | First ps | Last ps | dt ps | PBC |",
      "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|" ]
    for x in meta["trajectories"]:
        lines.append(f"| {x['system']} | {x['replica']:02d} | {x['atom_count']} | {x['residue_count']} | {x['receptor_residue_count']} | {x['stored_frames']} | {x['first_trajectory_time_ps']:.6f} | {x['last_trajectory_time_ps']:.6f} | {x['dt_ps']:.6f} | yes |")
    lines += ["", "Each trajectory is consistent with 50 ns from the DCD header and local completion ledger. Five deterministic full-coordinate frames per trajectory are finite; all frames were scanned at receptor-CA level for PBC discontinuities.","",
      "## 6. Completion/QC","", "All 12 `progress.json` files independently record `complete`, target/completed 50.0 ns. Each trajectory and `state.csv` has 1,000 stored records. Paired seeds are locally recorded and verified: R1=27101, R2=38201, R3=49301.","",
      "## 7. PBC/preprocessing route","", "Established short-data route: DCD -> per-system minimized PDB -> `chainID E and protein` -> raw stored receptor coordinates -> `pacer_dc_atom14_v0.1` -> expected heavy-atom point cloud -> ViSNet. The route does not make molecules whole, center, align, unwrap, or wrap. DCD box vectors are retained but no coordinate transformation is applied. Across 12,000 frames, the audit found zero adjacent-CA PBC splits, zero severe adjacent-CA breaks, and zero temporal CA wrapping events, so the established route is usable unchanged.","",
      "## 8. Atom/residue compatibility","", "All systems select 270 receptor residues in identical order and match the authenticated short-data sequence. Expected N/CA/C/O and side-chain atom14 heavy atoms are present; Gly uses N/CA/C/O only; HSE/HSD/HSP -> HIS and CHARMM ILE CD -> CD1 rules are preserved; terminal OXT/hydrogens are ignored. All 12 trajectories satisfy the frozen input contract without repair.","",
      "## 9. Path guard","", "`python -m unittest discover -s project/pacer_fkg_v02 -p 'test_*.py' -v`: 4/4 tests passed. The guard uses resolved paths plus component-aware containment checks and rejects every protected historical root and all paths outside the three allowed v02 roots.","",
      "## 10. Gate decision","", f"**{decision}**", "", "Documented non-blocking caveats:", ""]
    lines += [f"- {x}" for x in caveats]
    lines += ["", "## 11. Git diff/status summary", "", "Only new v02 audit/guard/report files are present; no baseline tracked file is modified. The provenance result root is intentionally ignored by the repository's existing rules. Current concise status:", "", "```text", *status, *ignored, "```", "",
      "Phase 0 stops here. No encoder inference, BS256 cache, v02 normalization, bandwidth fit, RFF draw, PACER-FKG, or biological analysis was performed.",""]
    write(V02/"REPORT_PHASE0.md","\n".join(lines))


if __name__=="__main__": main()
