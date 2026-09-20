#!/usr/bin/env python
"""Freeze the probe-matched M4 internal-distance coupling coordinate as JSON."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESIDS = [89, 92, 93, 96, 97, 98, 108, 183, 184, 185, 186,
          187, 190, 416, 420, 423, 432, 433, 435, 436, 439]


def ca_coordinates(path: Path) -> dict[int, np.ndarray]:
    atoms: dict[int, np.ndarray] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[21:22] != "R":
            continue
        if line[12:16].strip() != "CA":
            continue
        atoms[int(line[22:26])] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
    missing = sorted(set(RESIDS) - set(atoms))
    if missing:
        raise ValueError(f"{path.name}: missing chain-R CA residues {missing}")
    return atoms


def distances(atom_map: dict[int, np.ndarray]) -> tuple[list[tuple[int, int]], np.ndarray]:
    pairs = list(combinations(RESIDS, 2))
    values = np.array([np.linalg.norm(atom_map[a] - atom_map[b]) for a, b in pairs])
    return pairs, values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threshold-A", type=float, default=0.1)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "config" / "pacer_dc_coupling_coordinate_v01.json",
    )
    args = parser.parse_args()
    static = ROOT / "results" / "structure" / "ensemble"
    pairs, q = distances(ca_coordinates(static / "7TRQ_R.pdb"))
    _, p = distances(ca_coordinates(static / "7TRP_R.pdb"))
    _, s = distances(ca_coordinates(static / "7TRS_R.pdb"))
    dq, dp = q - s, p - s
    mask = (dq * dp > 0) & (np.abs(dq) >= args.threshold_A) & (np.abs(dp) >= args.threshold_A)
    direction = (dq + dp) / 2.0
    agreement = np.minimum(np.abs(dq), np.abs(dp)) / np.maximum(np.abs(dq), np.abs(dp))
    v = direction[mask] * np.sqrt(agreement[mask])
    denom = float(v @ v)
    weights = direction[mask] * agreement[mask] / denom
    selected_pairs = [pairs[i] for i in np.where(mask)[0]]
    baseline = s[mask]

    def score(x: np.ndarray) -> float:
        return float((x[mask] - baseline) @ weights)

    artifact = {
        "model_id": "PACER-DC-coupling-coordinate-v0.1",
        "threshold_A": args.threshold_A,
        "source_contrasts": ["7TRQ-7TRS", "7TRP-7TRS"],
        "source_probe": "iperoxo",
        "representation": "receptor_CA_pair_distances",
        "residue_numbering": "M4 receptor PDB numbering",
        "feature_count": int(mask.sum()),
        "pairs": [
            {"resid_i": int(a), "resid_j": int(b),
             "baseline_distance_A": float(base), "weight": float(weight)}
            for (a, b), base, weight in zip(selected_pairs, baseline, weights)
        ],
        "anchor_scores": {
            "7TRS_probe_only": score(s),
            "7TRQ_VU0467154_PAM": score(q),
            "7TRP_LY2033298_PAM": score(p),
        },
        "claim_boundary": (
            "Static probe-matched allosteric coupling direction; not PAM-specific and not efficacy."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(args.output), **artifact["anchor_scores"]}, indent=2))


if __name__ == "__main__":
    main()
