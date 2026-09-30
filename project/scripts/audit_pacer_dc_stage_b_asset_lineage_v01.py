#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import socket
import subprocess
from pathlib import Path

CANDIDATES = ("LY2119620", "CM00734", "compound110")
KEY_FILENAMES = {
    "audit.json", "preflight.json", "progress.json", "state.csv",
    "minimized.pdb", "system.xml", "state.xml", "production_start_state.xml",
    "latest_state.xml", "checkpoint.chk", "trajectory.dcd",
}
SEARCH_NAMES = (
    "7TRS.pdb", "7TRS_OPM.pdb", "7V68.pdb", "7V6A.pdb",
    "2CU_ideal.sdf", "ACH_ideal.sdf", "7TRS_R.pdb",
    "7TRS_R_meeko.pdbqt", "CM00734_7TRS_redocked.sdf",
    "compound110_7TRS_redocked.sdf",
)

def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def git(repo_root: Path, args: list[str]) -> str | None:
    try:
        p = subprocess.run(
            ["git", *args], cwd=repo_root, capture_output=True, text=True, check=False
        )
        return p.stdout.strip() if p.returncode == 0 else None
    except Exception:
        return None

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default=socket.gethostname())
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument("--search-root", type=Path, action="append", default=[])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    project = repo_root / "project"

    asset_roots = {
        "reference_complexes": project / "results" / "pacer_dc_reference_complexes_v01",
        "openmm_reference_qc": project / "results" / "pacer_dc_openmm_reference_qc_v01",
        "membrane_reference": project / "results" / "pacer_dc_membrane_reference_v01",
        "membrane_smoke": project / "results" / "pacer_dc_membrane_smoke_v01",
        "short_equilibration": project / "results" / "pacer_dc_short_equilibration_v01",
        "restraint_release": project / "results" / "pacer_dc_restraint_release_v01",
        "production": project / "results" / "pacer_dc_production_v01",
        "close_loop_20ns": project / "results" / "pacer_dc_close_loop_20ns_v01" / "production",
    }
    source_files = {
        "7TRS": project / "data" / "pdb" / "7TRS.pdb",
        "7TRS_OPM": project / "data" / "pdb" / "7TRS_OPM.pdb",
        "LY_source_7V68": project / "data" / "pdb" / "m4_external_structures" / "7V68.pdb",
        "compound110_source_7V6A": project / "data" / "pdb" / "m4_external_structures" / "7V6A.pdb",
        "LY_template_2CU": project / "data" / "pdb" / "m4_ligands" / "2CU_ideal.sdf",
        "ACH_template": project / "data" / "pdb" / "m4_ligands" / "ACH_ideal.sdf",
        "CM_pose_sdf": project / "results" / "pacer_dc_reference_complexes_v01" / "CM00734_7TRS_redocked.sdf",
        "compound110_pose_sdf": project / "results" / "pacer_dc_reference_complexes_v01" / "compound110_7TRS_redocked.sdf",
        "docking_receptor_pdb": project / "results" / "structure" / "ensemble" / "7TRS_R.pdb",
        "docking_receptor_pdbqt": project / "results" / "structure" / "ensemble" / "7TRS_R_meeko.pdbqt",
        "vina_windows": project / "tools" / "vina.exe",
    }
    pipeline_scripts = (
        project / "scripts" / "build_pacer_dc_reference_complexes.py",
        project / "scripts" / "build_pacer_dc_openmm_reference_systems.py",
        project / "scripts" / "build_pacer_dc_membrane_reference.py",
        project / "scripts" / "run_pacer_dc_membrane_smoke.py",
        project / "scripts" / "run_pacer_dc_short_equilibration.py",
        project / "scripts" / "run_pacer_dc_restraint_release.py",
        project / "scripts" / "run_pacer_dc_production_md.py",
        project / "scripts" / "prepare_cm00734_stage_b_pose_v01.py",
    )

    def rel(path: Path) -> str:
        try:
            return str(path.relative_to(project))
        except ValueError:
            return str(path)

    def file_record(path: Path, *, hash_file: bool = True) -> dict:
        rec = {
            "path": rel(path),
            "exists": path.exists(),
            "is_file": path.is_file(),
            "is_dir": path.is_dir(),
        }
        if path.is_file():
            rec["size_bytes"] = path.stat().st_size
            if hash_file:
                rec["sha256"] = sha256(path)
        return rec

    def stage_summary(candidate: str) -> dict:
        names = (f"{candidate}__candidate_no_probe", f"{candidate}__candidate_probe")
        out = {}
        for label, root in asset_roots.items():
            rows = []
            for name in names:
                p = root / name
                files = []
                if p.exists():
                    for x in sorted((q for q in p.rglob("*") if q.is_file()), key=lambda q: str(q)):
                        if x.name in KEY_FILENAMES or x.suffix.lower() in {".json", ".csv"}:
                            files.append(file_record(x, hash_file=False))
                rows.append({
                    "context": name,
                    "exists": p.exists(),
                    "key_files": files,
                    "file_count": sum(1 for x in p.rglob("*") if x.is_file()) if p.exists() else 0,
                })
            out[label] = rows
        return out

    def named_hits(candidate: str) -> list[dict]:
        hits = []
        needle = candidate.lower()
        for label, root in asset_roots.items():
            if not root.exists():
                continue
            for p in root.rglob("*"):
                if p.is_file() and needle in str(p).lower():
                    hits.append({"stage": label, **file_record(p, hash_file=False)})
        return hits

    def script_contract(path: Path) -> dict:
        rec = file_record(path, hash_file=True)
        if not path.is_file():
            return rec
        text = path.read_text(encoding="utf-8", errors="ignore")
        rec["mentions"] = {c: text.count(c) for c in CANDIDATES}
        rec["contains_7V68"] = "7V68" in text
        rec["contains_7V6A"] = "7V6A" in text
        rec["contains_7TRS_common_pocket_pose_proposal"] = (
            "7TRS_common_pocket_pose_proposal" in text
        )
        return rec

    manifest = project / "config" / "pacer_dc_pilot_md_manifest_v2.csv"
    manifest_rows = []
    if manifest.exists():
        for line in manifest.read_text(encoding="utf-8-sig", errors="ignore").splitlines():
            if any(c in line for c in CANDIDATES):
                manifest_rows.append(line)

    external_hits = {}
    for root in [p.resolve() for p in args.search_root]:
        root_hits = []
        if root.exists():
            for target in SEARCH_NAMES:
                for p in root.rglob(target):
                    if p.is_file():
                        root_hits.append(file_record(p, hash_file=False))
            for candidate in CANDIDATES:
                for p in root.rglob(f"*{candidate}*"):
                    if p.exists():
                        root_hits.append(file_record(p, hash_file=False))
        external_hits[str(root)] = root_hits

    report = {
        "schema": "pacer_dc.stage_b_asset_lineage_audit.v3",
        "label": args.label,
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "cwd": os.getcwd(),
        "repo_root": str(repo_root),
        "project_root": str(project),
        "git": {
            "branch": git(repo_root, ["branch", "--show-current"]),
            "head": git(repo_root, ["rev-parse", "HEAD"]),
            "status_short": git(repo_root, ["status", "--short"]),
        },
        "source_files": {k: file_record(v) for k, v in source_files.items()},
        "manifest": {**file_record(manifest), "candidate_rows": manifest_rows},
        "pipeline_scripts": [script_contract(p) for p in pipeline_scripts],
        "candidates": {},
        "external_search_hits": external_hits,
    }

    for candidate in CANDIDATES:
        report["candidates"][candidate] = {
            "stage_summary": stage_summary(candidate),
            "named_result_file_hits": named_hits(candidate),
        }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"AUDIT_LABEL={args.label}")
    print(f"REPO_ROOT={repo_root}")
    print(f"GIT_BRANCH={report['git']['branch']}")
    print(f"GIT_HEAD={report['git']['head']}")
    for candidate in CANDIDATES:
        stages = report["candidates"][candidate]["stage_summary"]
        present = []
        for stage, rows in stages.items():
            for row in rows:
                if row["exists"]:
                    present.append(f"{stage}:{row['context']}")
        print(f"{candidate}_PRESENT_STAGES=" + (";".join(present) if present else "NONE"))
    for key, rec in report["source_files"].items():
        print(f"SOURCE_{key}={'PRESENT' if rec['exists'] else 'MISSING'}")
    for root, hits in external_hits.items():
        print(f"SEARCH_ROOT={root} HITS={len(hits)}")
    print(f"AUDIT_JSON={args.out}")

if __name__ == "__main__":
    main()
