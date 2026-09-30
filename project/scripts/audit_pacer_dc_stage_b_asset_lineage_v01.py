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

ROOT = Path(__file__).resolve().parents[1]

CANDIDATES = ("LY2119620", "CM00734", "compound110")

ASSET_ROOTS = {
    "reference_complexes": ROOT / "results" / "pacer_dc_reference_complexes_v01",
    "openmm_reference_qc": ROOT / "results" / "pacer_dc_openmm_reference_qc_v01",
    "membrane_reference": ROOT / "results" / "pacer_dc_membrane_reference_v01",
    "membrane_smoke": ROOT / "results" / "pacer_dc_membrane_smoke_v01",
    "short_equilibration": ROOT / "results" / "pacer_dc_short_equilibration_v01",
    "restraint_release": ROOT / "results" / "pacer_dc_restraint_release_v01",
    "production": ROOT / "results" / "pacer_dc_production_v01",
    "close_loop_20ns": ROOT / "results" / "pacer_dc_close_loop_20ns_v01",
}

SOURCE_FILES = {
    "7TRS": ROOT / "data" / "pdb" / "7TRS.pdb",
    "7TRS_OPM": ROOT / "data" / "pdb" / "7TRS_OPM.pdb",
    "LY_source_7V68": ROOT / "data" / "pdb" / "m4_external_structures" / "7V68.pdb",
    "compound110_source_7V6A": ROOT / "data" / "pdb" / "m4_external_structures" / "7V6A.pdb",
    "LY_template_2CU": ROOT / "data" / "pdb" / "m4_ligands" / "2CU_ideal.sdf",
    "ACH_template": ROOT / "data" / "pdb" / "m4_ligands" / "ACH_ideal.sdf",
    "CM_pose_sdf": ROOT / "results" / "pacer_dc_reference_complexes_v01" / "CM00734_7TRS_redocked.sdf",
    "compound110_pose_sdf": ROOT / "results" / "pacer_dc_reference_complexes_v01" / "compound110_7TRS_redocked.sdf",
    "docking_receptor_pdb": ROOT / "results" / "structure" / "ensemble" / "7TRS_R.pdb",
    "docking_receptor_pdbqt": ROOT / "results" / "structure" / "ensemble" / "7TRS_R_meeko.pdbqt",
    "vina_windows": ROOT / "tools" / "vina.exe",
}

PIPELINE_SCRIPTS = (
    ROOT / "scripts" / "build_pacer_dc_reference_complexes.py",
    ROOT / "scripts" / "build_pacer_dc_openmm_reference_systems.py",
    ROOT / "scripts" / "build_pacer_dc_membrane_reference.py",
    ROOT / "scripts" / "run_pacer_dc_membrane_smoke.py",
    ROOT / "scripts" / "run_pacer_dc_short_equilibration.py",
    ROOT / "scripts" / "run_pacer_dc_restraint_release.py",
    ROOT / "scripts" / "run_pacer_dc_production_md.py",
    ROOT / "scripts" / "prepare_cm00734_stage_b_pose_v01.py",
)

def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def file_record(path: Path) -> dict:
    rec = {
        "path": str(path.relative_to(ROOT)) if path.is_absolute() else str(path),
        "exists": path.exists(),
        "is_file": path.is_file(),
        "is_dir": path.is_dir(),
    }
    if path.is_file():
        rec["size_bytes"] = path.stat().st_size
        rec["sha256"] = sha256(path)
    return rec

def git(args: list[str]) -> str | None:
    try:
        p = subprocess.run(
            ["git", *args], cwd=ROOT.parent, capture_output=True, text=True, check=False
        )
        return p.stdout.strip() if p.returncode == 0 else None
    except Exception:
        return None

def candidate_paths(candidate: str) -> list[dict]:
    hits = []
    needle = candidate.lower()
    for label, root in ASSET_ROOTS.items():
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if needle in p.name.lower() or needle in str(p.parent).lower():
                if p.is_file():
                    hits.append({
                        "stage": label,
                        **file_record(p),
                    })
    return hits

def stage_summary(candidate: str) -> dict:
    names = (
        f"{candidate}__candidate_no_probe",
        f"{candidate}__candidate_probe",
    )
    out = {}
    for label, root in ASSET_ROOTS.items():
        rows = []
        for name in names:
            p = root / name
            rows.append({
                "context": name,
                "exists": p.exists(),
                "files": sorted(
                    str(x.relative_to(ROOT))
                    for x in p.rglob("*")
                    if x.is_file()
                ) if p.exists() else [],
            })
        out[label] = rows
    return out

def script_contract(path: Path) -> dict:
    rec = file_record(path)
    if not path.is_file():
        return rec
    text = path.read_text(encoding="utf-8", errors="ignore")
    rec["mentions"] = {
        c: text.count(c) for c in CANDIDATES
    }
    rec["contains_7V68"] = "7V68" in text
    rec["contains_7V6A"] = "7V6A" in text
    rec["contains_7TRS_common_pocket_pose_proposal"] = "7TRS_common_pocket_pose_proposal" in text
    return rec

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default=socket.gethostname())
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    manifest = ROOT / "config" / "pacer_dc_pilot_md_manifest_v2.csv"
    manifest_rows = []
    if manifest.exists():
        for line in manifest.read_text(encoding="utf-8-sig", errors="ignore").splitlines():
            if any(c in line for c in ("LY2119620", "CM00734", "compound110")):
                manifest_rows.append(line)

    report = {
        "schema": "pacer_dc.stage_b_asset_lineage_audit.v1",
        "label": args.label,
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "cwd": os.getcwd(),
        "repo_root": str(ROOT.parent),
        "git": {
            "branch": git(["branch", "--show-current"]),
            "head": git(["rev-parse", "HEAD"]),
            "status_short": git(["status", "--short"]),
        },
        "source_files": {k: file_record(v) for k, v in SOURCE_FILES.items()},
        "manifest": {
            **file_record(manifest),
            "candidate_rows": manifest_rows,
        },
        "pipeline_scripts": [script_contract(p) for p in PIPELINE_SCRIPTS],
        "candidates": {},
    }

    for candidate in CANDIDATES:
        report["candidates"][candidate] = {
            "stage_summary": stage_summary(candidate),
            "all_result_file_hits": candidate_paths(candidate),
        }

    text = json.dumps(report, indent=2, ensure_ascii=False)
    print(text)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")

if __name__ == "__main__":
    main()
