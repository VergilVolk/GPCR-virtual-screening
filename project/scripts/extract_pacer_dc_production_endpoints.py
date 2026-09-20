#!/usr/bin/env python
"""Extract matched endpoints from PACER-DC GPU production jobs after the frozen QC gate."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "results" / "pacer_dc_membrane_reference_v01"
COMMON = ROOT / "results" / "pacer_dc_common_protein_v01" / "7TRS_common_protein_pH74.pdb"
MODEL = ROOT / "config" / "pacer_dc_coupling_coordinate_v01.json"
DEFAULT_MANIFEST = ROOT / "config" / "pacer_dc_gpu_phase1_manifest.csv"
DEFAULT_RUNS = ROOT / "results" / "pacer_dc_production_v01"
DEFAULT_QC = ROOT / "results" / "pacer_dc_production_phase1_audit_v01" / "job_qc.csv"


def ca_residues(path: Path, chain: str) -> list[tuple[int, str]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if (
            line.startswith(("ATOM", "HETATM")) and len(line) >= 26
            and line[21] == chain and line[12:16].strip() == "CA"
        ):
            rows.append((int(line[22:26]), line[17:20].strip()))
    return rows


def system_spec(system: str) -> dict[str, str]:
    if system == "apo":
        return {"candidate": "SHARED_CONTROL", "context": "apo"}
    if system == "probe_only":
        return {
            "candidate": "SHARED_CONTROL", "context": "probe_only",
            "probe": "chainID J and resname UNK and not name H*",
        }
    candidate, context = system.split("__", 1)
    if context == "candidate_no_probe":
        return {
            "candidate": candidate, "context": context,
            "ligand": "chainID J and resname UNK and not name H*",
        }
    if context == "candidate_probe":
        return {
            "candidate": candidate, "context": context,
            "probe": "chainID J and resname UNK and not name H*",
            "ligand": "chainID K and resname UNK and not name H*",
        }
    raise ValueError(system)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--runs", type=Path, default=DEFAULT_RUNS)
    parser.add_argument("--qc", type=Path, default=DEFAULT_QC)
    parser.add_argument("--evidence-level", choices=("production_pilot", "trajectory"), default="production_pilot")
    parser.add_argument("--allow-formal-trajectory", action="store_true")
    args = parser.parse_args()
    if args.evidence_level == "trajectory" and not args.allow_formal_trajectory:
        raise SystemExit("Formal trajectory evidence requires explicit --allow-formal-trajectory after convergence audit.")
    manifest = pd.read_csv(args.manifest)
    qc = pd.read_csv(args.qc)
    merged = manifest.merge(qc[["array_id", "qc_pass"]], on="array_id", how="left")
    if len(merged) != len(manifest) or not merged.qc_pass.fillna(False).all():
        raise RuntimeError("All manifest jobs must pass the frozen production QC before endpoint extraction.")

    model = json.loads(MODEL.read_text(encoding="utf-8"))
    model_resids = sorted({p["resid_i"] for p in model["pairs"]} | {p["resid_j"] for p in model["pairs"]})
    original = ca_residues(COMMON, "R")
    renumbered = ca_residues(REFERENCE / "apo" / "minimized.pdb", "E")
    if len(original) != len(renumbered) or [x[1] for x in original] != [x[1] for x in renumbered]:
        raise RuntimeError("Audited receptor residue mapping failed.")
    mapping = {old[0]: new[0] for old, new in zip(original, renumbered)}
    receptor_resids = ",".join(str(mapping[x]) for x in model_resids)
    extractor = ROOT / "scripts" / "extract_pacer_dc_trajectory_endpoints.py"
    evidence_paths = []

    for job in merged.itertuples(index=False):
        spec = system_spec(job.system)
        run = args.runs / job.system / f"replica_{int(job.replica):02d}"
        out = run / "endpoints"
        out.mkdir(parents=True, exist_ok=True)
        clean = out / "topology_mdanalysis.pdb"
        clean.write_text(
            "\n".join(
                line for line in (REFERENCE / job.system / "minimized.pdb").read_text(encoding="utf-8").splitlines()
                if not line.startswith("CONECT")
            ) + "\n", encoding="utf-8",
        )
        command = [
            sys.executable, str(extractor),
            "--topology", str(clean),
            "--trajectory", str(run / "trajectory.dcd"),
            "--candidate-id", spec["candidate"],
            "--context", spec["context"],
            "--replicate-id", f"paired_seed_{int(job.replica)}",
            "--source-id", f"pacer_dc_production_v01_{float(job.target_ns):g}ns",
            "--evidence-level", args.evidence_level,
            "--receptor-resids", receptor_resids,
            "--receptor-chain", "E",
            "--reference-topology", str(clean),
            "--outdir", str(out),
        ]
        if spec.get("probe"):
            command.extend(["--probe-selection", spec["probe"]])
        if spec.get("ligand"):
            command.extend(["--allosteric-selection", spec["ligand"]])
        subprocess.run(command, check=True)
        evidence_paths.append(out / "replica_evidence.csv")

    evidence = pd.concat([pd.read_csv(path) for path in evidence_paths], ignore_index=True)
    output = args.runs / f"{args.evidence_level}_evidence.csv"
    evidence.to_csv(output, index=False)
    audit = {
        "jobs": len(merged),
        "evidence_rows": len(evidence),
        "evidence_level": args.evidence_level,
        "eligible_for_functional_scoring": args.evidence_level == "trajectory",
        "claim_boundary": (
            "Five-nanosecond production pilots are QC/signal-discovery data, not converged PAM evidence."
            if args.evidence_level == "production_pilot" else
            "Formal trajectory label records protocol eligibility only; biological claims still require the frozen statistics."
        ),
    }
    (args.runs / f"{args.evidence_level}_endpoint_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
