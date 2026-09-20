#!/usr/bin/env python
"""Build auditable PACER-DC reference starting complexes on a common 7TRS base.

Experimental allosteric ligand coordinates are transferred by a Kabsch fit of
common M4 receptor (chain R) C-alpha atoms.  The script only constructs starting
coordinates; membrane building, force-field assignment and MD are separate QC
stages.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "LY2119620": ("7V68.pdb", "2CU"),
    "compound110": ("7V6A.pdb", "5XI"),
}
POCKET_RESIDS = {89, 92, 93, 96, 184, 186, 190, 423, 432, 433, 435, 436}


def xyz(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def receptor_ca(lines: list[str], chain: str = "R") -> dict[tuple[int, str], np.ndarray]:
    result = {}
    for line in lines:
        if not line.startswith("ATOM  ") or line[21:22] != chain:
            continue
        if line[12:16].strip() == "CA":
            result[(int(line[22:26]), line[26:27])] = xyz(line)
    return result


def fit_transform(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    h = (source - source_center).T @ (target - target_center)
    u, _, vt = np.linalg.svd(h)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = vt.T @ u.T
    translation = target_center - source_center @ rotation.T
    fitted = source @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def replace_xyz(line: str, coordinates: np.ndarray, serial: int, chain: str = "L") -> str:
    return (
        f"HETATM{serial:5d}" + line[11:21] + chain + f"{1:4d}" + line[26:30]
        + f"{coordinates[0]:8.3f}{coordinates[1]:8.3f}{coordinates[2]:8.3f}"
        + line[54:]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=ROOT / "data" / "pdb" / "7TRS.pdb")
    parser.add_argument(
        "--source-dir", type=Path,
        default=ROOT / "data" / "pdb" / "m4_external_structures",
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=ROOT / "results" / "pacer_dc_reference_complexes_v01",
    )
    args = parser.parse_args()

    base_lines = args.base.read_text(encoding="utf-8", errors="ignore").splitlines()
    base_ca = receptor_ca(base_lines)
    protein = [line for line in base_lines if line.startswith("ATOM  ")]
    protein_heavy_xyz = np.stack([
        xyz(line) for line in protein if not line[76:78].strip().upper().startswith("H")
    ])
    ach = [
        line for line in base_lines
        if line.startswith("HETATM") and line[17:20].strip() == "ACH"
    ]
    args.outdir.mkdir(parents=True, exist_ok=True)
    audits = []

    # Shared controls contain no transferred allosteric ligand.
    for context, include_probe in (("probe_only", True), ("apo", False)):
        output = args.outdir / f"SHARED_CONTROL__{context}.pdb"
        lines = protein + (ach if include_probe else []) + ["END"]
        output.write_text("\n".join(lines) + "\n", encoding="utf-8")
        audits.append({
            "candidate_id": "SHARED_CONTROL", "context": context,
            "output": str(output), "protein_atoms": len(protein),
            "probe_atoms": len(ach) if include_probe else 0,
            "allosteric_atoms": 0,
        })

    for candidate, (filename, residue_name) in SOURCES.items():
        source_lines = (args.source_dir / filename).read_text(
            encoding="utf-8", errors="ignore"
        ).splitlines()
        source_ca = receptor_ca(source_lines)
        common = sorted(set(base_ca) & set(source_ca))
        if len(common) < 100:
            raise ValueError(f"{candidate}: only {len(common)} common receptor CA atoms")
        source_xyz = np.stack([source_ca[key] for key in common])
        target_xyz = np.stack([base_ca[key] for key in common])
        _, _, global_fit_rmsd = fit_transform(source_xyz, target_xyz)
        local_common = [key for key in common if key[0] in POCKET_RESIDS]
        if len(local_common) < 6:
            raise ValueError(f"{candidate}: only {len(local_common)} common pocket CA atoms")
        local_source = np.stack([source_ca[key] for key in local_common])
        local_target = np.stack([base_ca[key] for key in local_common])
        rotation, translation, fit_rmsd = fit_transform(local_source, local_target)
        ligand = [
            line for line in source_lines
            if line.startswith("HETATM") and line[17:20].strip() == residue_name
        ]
        if not ligand:
            raise ValueError(f"{candidate}: residue {residue_name} not found in {filename}")
        transferred = []
        transferred_xyz = []
        serial = len(protein) + len(ach) + 1
        for line in ligand:
            transformed = xyz(line) @ rotation.T + translation
            transferred_xyz.append(transformed)
            transferred.append(replace_xyz(line, transformed, serial))
            serial += 1
        minimum_protein_distance = float(np.min(np.linalg.norm(
            np.stack(transferred_xyz)[:, None, :] - protein_heavy_xyz[None, :, :], axis=2
        )))

        for context, include_probe in (("candidate_probe", True), ("candidate_no_probe", False)):
            output = args.outdir / f"{candidate}__{context}.pdb"
            lines = protein + (ach if include_probe else []) + transferred + ["END"]
            output.write_text("\n".join(lines) + "\n", encoding="utf-8")
            audits.append({
                "candidate_id": candidate, "context": context,
                "source_structure": filename, "source_ligand_resname": residue_name,
                "common_receptor_ca": len(common),
                "global_alignment_rmsd_A": global_fit_rmsd,
                "pocket_alignment_ca": len(local_common),
                "pocket_alignment_rmsd_A": fit_rmsd,
                "minimum_ligand_protein_distance_A": minimum_protein_distance,
                "severe_clash": minimum_protein_distance < 1.0,
                "output": str(output), "protein_atoms": len(protein),
                "probe_atoms": len(ach) if include_probe else 0,
                "allosteric_atoms": len(transferred),
            })

    audit = {
        "method": "PACER-DC common-base reference complex construction",
        "base_structure": str(args.base),
        "base_probe": "ACH",
        "complex_count": len(audits),
        "records": audits,
        "claim_boundary": "Starting-coordinate construction only; no force-field or biological claim.",
    }
    (args.outdir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
