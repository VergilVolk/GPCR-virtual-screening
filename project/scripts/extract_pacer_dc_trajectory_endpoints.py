#!/usr/bin/env python
"""Extract PACER-DC frame traces and one replica-level evidence record.

Selections are explicit command-line inputs because residue/ligand naming differs across
trajectory engines.  The output is compatible with ``pacer_dc_score.py``.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation


ROOT = Path(__file__).resolve().parents[1]


def align_coordinates(mobile: np.ndarray, reference: np.ndarray, coords: np.ndarray) -> np.ndarray:
    mc, rc = mobile.mean(axis=0), reference.mean(axis=0)
    rotation, _ = Rotation.align_vectors(reference - rc, mobile - mc)
    return rotation.apply(coords - mc) + rc


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def nearest_periodic_image(coords: np.ndarray, reference_center: np.ndarray, box_lengths: np.ndarray) -> np.ndarray:
    """Translate one molecule to the periodic image nearest a reference center."""
    delta = coords.mean(axis=0) - reference_center
    shift = box_lengths * np.round(delta / box_lengths)
    return coords - shift


def pbc_distances(a: np.ndarray, b: np.ndarray, box_lengths: np.ndarray) -> np.ndarray:
    delta = a[:, None, :] - b[None, :, :]
    delta -= box_lengths * np.round(delta / box_lengths)
    return np.linalg.norm(delta, axis=2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument("--context", choices=["candidate_probe", "candidate_no_probe", "probe_only", "apo"], required=True)
    parser.add_argument("--replicate-id", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument(
        "--evidence-level",
        choices=["trajectory", "production_pilot", "equilibration_pilot"],
        default="trajectory",
        help="Only trajectory satisfies scoring gates; pilot levels are protocol/QC evidence.",
    )
    parser.add_argument("--receptor-resids", required=True,
                        help="21 topology residue ids ordered like model receptor residues")
    parser.add_argument("--receptor-chain", default="",
                        help="Topology chain ID for the receptor; required when residue ids repeat across chains")
    parser.add_argument("--probe-selection", default="")
    parser.add_argument("--allosteric-selection", default="")
    parser.add_argument("--pocket-resids", default="89,92,93,96,184,186,190,423,432,433,435,436")
    parser.add_argument("--reference-topology", type=Path)
    parser.add_argument("--model", type=Path,
                        default=ROOT / "config" / "pacer_dc_coupling_coordinate_v01.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    model = json.loads(args.model.read_text(encoding="utf-8"))
    model_resids = sorted({p["resid_i"] for p in model["pairs"]} | {p["resid_j"] for p in model["pairs"]})
    topology_resids = [int(x) for x in args.receptor_resids.split(",")]
    if len(topology_resids) != len(model_resids):
        raise ValueError(f"Expected {len(model_resids)} receptor ids, got {len(topology_resids)}")
    resid_map = dict(zip(model_resids, topology_resids))
    pair_indices = [(topology_resids.index(resid_map[p["resid_i"]]),
                     topology_resids.index(resid_map[p["resid_j"]])) for p in model["pairs"]]
    baseline = np.array([p["baseline_distance_A"] for p in model["pairs"]])
    weights = np.array([p["weight"] for p in model["pairs"]])

    universe = mda.Universe(str(args.topology), str(args.trajectory))
    chain_clause = f"chainID {args.receptor_chain} and " if args.receptor_chain else ""
    ca = universe.select_atoms(chain_clause + "name CA and resid " + " ".join(map(str, topology_resids)))
    ca_by_resid = {int(atom.resid): atom for atom in ca}
    if set(ca_by_resid) != set(topology_resids):
        raise ValueError(f"CA residue mismatch: got {sorted(ca_by_resid)}")
    ordered_ca = [ca_by_resid[r] for r in topology_resids]

    reference_path = args.reference_topology or args.topology
    reference = mda.Universe(str(reference_path))
    ref_ca_group = reference.select_atoms(
        chain_clause + "name CA and resid " + " ".join(map(str, topology_resids))
    )
    ref_ca_by_resid = {int(atom.resid): atom for atom in ref_ca_group}
    ref_ca = np.stack([ref_ca_by_resid[r].position.copy() for r in topology_resids])

    probe = universe.select_atoms(args.probe_selection) if args.probe_selection else None
    ref_probe = reference.select_atoms(args.probe_selection) if args.probe_selection else None
    if probe is not None and (not len(probe) or len(probe) != len(ref_probe)):
        raise ValueError(f"Probe selection mismatch: trajectory={len(probe)}, reference={len(ref_probe)}")
    ref_probe_xyz = ref_probe.positions.copy() if probe is not None else None

    ligand = universe.select_atoms(args.allosteric_selection) if args.allosteric_selection else None
    pocket_model = [int(x) for x in args.pocket_resids.split(",") if x]
    pocket_topology = [resid_map[r] for r in pocket_model if r in resid_map]
    pocket = universe.select_atoms(
        chain_clause + "resid " + " ".join(map(str, pocket_topology)) + " and not name H*"
    ) if ligand is not None else None
    if ligand is not None and (not len(ligand) or not len(pocket)):
        raise ValueError("Allosteric ligand or pocket selection is empty")

    frame_rows = []
    for frame_index, ts in enumerate(universe.trajectory):
        if ts.dimensions is None or np.any(np.asarray(ts.dimensions[:3]) <= 0):
            raise ValueError("Periodic box dimensions are required for ligand endpoint extraction")
        box_lengths = np.asarray(ts.dimensions[:3], dtype=float)
        coords = np.stack([atom.position.copy() for atom in ordered_ca])
        distances = np.array([np.linalg.norm(coords[i] - coords[j]) for i, j in pair_indices])
        coupling = float((distances - baseline) @ weights)
        row = {"frame": frame_index, "time_ps": float(ts.time), "coupling_coordinate": coupling}
        if probe is not None:
            transformed = align_coordinates(coords, ref_ca, probe.positions.copy())
            transformed = nearest_periodic_image(transformed, ref_probe_xyz.mean(axis=0), box_lengths)
            row["orthosteric_pose_rmsd_A"] = rmsd(transformed, ref_probe_xyz)
        if ligand is not None:
            ligand_xyz = ligand.positions.astype(float)
            # Per-residue contact coverage, averaged later as a 0..1 compatibility gate.
            contacted = 0
            for residue in pocket.residues:
                heavy = residue.atoms.select_atoms("not name H*").positions.astype(float)
                if np.min(pbc_distances(ligand_xyz, heavy, box_lengths)) <= 4.5:
                    contacted += 1
            row["binding_compatibility"] = contacted / max(len(pocket.residues), 1)
        frame_rows.append(row)

    frames = pd.DataFrame(frame_rows)
    evidence = []
    for metric in ("coupling_coordinate", "orthosteric_pose_rmsd_A", "binding_compatibility"):
        if metric in frames:
            evidence.append({
                "candidate_id": args.candidate_id,
                "context": args.context,
                "replicate_id": args.replicate_id,
                "metric": metric,
                "value": float(frames[metric].mean()),
                "evidence_level": args.evidence_level,
                "source_id": args.source_id,
            })
    args.outdir.mkdir(parents=True, exist_ok=True)
    frames.to_csv(args.outdir / "frame_endpoints.csv", index=False)
    pd.DataFrame(evidence).to_csv(args.outdir / "replica_evidence.csv", index=False)
    audit = {
        "frames": int(len(frames)),
        "candidate_id": args.candidate_id,
        "context": args.context,
        "replicate_id": args.replicate_id,
        "metrics": [x["metric"] for x in evidence],
        "evidence_level": args.evidence_level,
        "inference_unit": "one trajectory replica",
        "eligible_for_functional_scoring": args.evidence_level == "trajectory",
    }
    (args.outdir / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
