#!/usr/bin/env python
"""Extract QC-only endpoints from the six matched PACER-DC equilibration pilots."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "results" / "pacer_dc_membrane_reference_v01"
RUNS = ROOT / "results" / "pacer_dc_short_equilibration_v01"
MODEL = ROOT / "config" / "pacer_dc_coupling_coordinate_v01.json"
COMMON_PROTEIN = ROOT / "results" / "pacer_dc_common_protein_v01" / "7TRS_common_protein_pH74.pdb"

SYSTEMS = {
    "apo": {"candidate": "SHARED_CONTROL", "context": "apo"},
    "probe_only": {
        "candidate": "SHARED_CONTROL",
        "context": "probe_only",
        "probe": "chainID J and resname UNK and not name H*",
    },
    "LY2119620__candidate_no_probe": {
        "candidate": "LY2119620",
        "context": "candidate_no_probe",
        "ligand": "chainID J and resname UNK and not name H*",
    },
    "LY2119620__candidate_probe": {
        "candidate": "LY2119620",
        "context": "candidate_probe",
        "probe": "chainID J and resname UNK and not name H*",
        "ligand": "chainID K and resname UNK and not name H*",
    },
    "compound110__candidate_no_probe": {
        "candidate": "compound110",
        "context": "candidate_no_probe",
        "ligand": "chainID J and resname UNK and not name H*",
    },
    "compound110__candidate_probe": {
        "candidate": "compound110",
        "context": "candidate_probe",
        "probe": "chainID J and resname UNK and not name H*",
        "ligand": "chainID K and resname UNK and not name H*",
    },
}


def main() -> None:
    model = json.loads(MODEL.read_text(encoding="utf-8"))
    resids = sorted({p["resid_i"] for p in model["pairs"]} | {p["resid_j"] for p in model["pairs"]})
    def ca_residues(path: Path, chain: str) -> list[tuple[int, str]]:
        rows = []
        for line in path.read_text(encoding="utf-8").splitlines():
            if (
                line.startswith(("ATOM", "HETATM"))
                and len(line) >= 26
                and line[21] == chain
                and line[12:16].strip() == "CA"
            ):
                rows.append((int(line[22:26]), line[17:20].strip()))
        return rows

    original = ca_residues(COMMON_PROTEIN, "R")
    renumbered = ca_residues(REFERENCE / "apo" / "minimized.pdb", "E")
    if len(original) != len(renumbered) or [x[1] for x in original] != [x[1] for x in renumbered]:
        raise RuntimeError("Cannot establish the audited original-R to membrane-E receptor residue map.")
    residue_map = {old[0]: new[0] for old, new in zip(original, renumbered)}
    missing = [resid for resid in resids if resid not in residue_map]
    if missing:
        raise RuntimeError(f"Model residues absent from membrane receptor: {missing}")
    receptor_resids = ",".join(str(residue_map[x]) for x in resids)
    audits = []
    evidence_paths = []
    extractor = ROOT / "scripts" / "extract_pacer_dc_trajectory_endpoints.py"
    for system_name, spec in SYSTEMS.items():
        source = REFERENCE / system_name
        run = RUNS / system_name / "replica_01"
        audit_path = run / "audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if not audit.get("passed", False):
            raise RuntimeError(f"Equilibration gate not passed: {system_name}")
        out = run / "endpoints"
        out.mkdir(parents=True, exist_ok=True)
        # OpenMM writes hexadecimal atom serials above 99,999.  MDAnalysis can
        # preserve the atom order but cannot parse hexadecimal CONECT records.
        # Removing CONECT is safe here because all endpoints are coordinate/contact based.
        clean_topology = out / "topology_mdanalysis.pdb"
        clean_topology.write_text(
            "\n".join(
                line for line in (source / "minimized.pdb").read_text(encoding="utf-8").splitlines()
                if not line.startswith("CONECT")
            ) + "\n",
            encoding="utf-8",
        )
        command = [
            sys.executable,
            str(extractor),
            "--topology", str(clean_topology),
            "--trajectory", str(run / "npt.dcd"),
            "--candidate-id", spec["candidate"],
            "--context", spec["context"],
            "--replicate-id", "paired_seed_1",
            "--source-id", "pacer_dc_short_equilibration_v01",
            "--evidence-level", "equilibration_pilot",
            "--receptor-resids", receptor_resids,
            "--receptor-chain", "E",
            "--reference-topology", str(clean_topology),
            "--outdir", str(out),
        ]
        if spec.get("probe"):
            command.extend(["--probe-selection", spec["probe"]])
        if spec.get("ligand"):
            command.extend(["--allosteric-selection", spec["ligand"]])
        subprocess.run(command, check=True)
        evidence_paths.append(out / "replica_evidence.csv")
        audits.append(json.loads((out / "audit.json").read_text(encoding="utf-8")))

    import pandas as pd

    evidence = pd.concat([pd.read_csv(path) for path in evidence_paths], ignore_index=True)
    evidence.to_csv(RUNS / "equilibration_pilot_evidence.csv", index=False)
    summary = {
        "n_systems": len(audits),
        "n_evidence_rows": len(evidence),
        "evidence_level": "equilibration_pilot",
        "eligible_for_functional_scoring": False,
        "claim_boundary": (
            "Endpoints are plumbing/QC outputs from 1 ps NPT segments. They must not be used to "
            "estimate PAM activity, model performance, convergence, or biological mechanism."
        ),
    }
    (RUNS / "endpoint_extraction_audit.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
