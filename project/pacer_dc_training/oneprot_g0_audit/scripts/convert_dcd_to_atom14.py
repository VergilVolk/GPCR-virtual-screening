#!/usr/bin/env python
"""Convert the public M4 xanomeline trajectory (PSF + DCD) into the atom14
[n_frames, L, 14, 3] representation MDGen consumes, in Angstrom.

One-time data prep for the G0 embedding smoke test.  Produces:
  artifacts/m4_atom14/m4xan_i1.npy          atom14 coords, float32, Angstrom
  artifacts/m4_atom14/m4xan_test.csv        name + one-letter seqres

CHARMM naming quirks handled here:
  * HSE/HSD/HSP -> HIS
  * ILE delta carbon is "CD" in CHARMM but "CD1" in the atom14 encoding
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import MDAnalysis as mda

ROOT = Path("project/tools/oneprot-embeddings").resolve()
sys.path.insert(0, str(ROOT / "external" / "mdgen"))
from mdgen import residue_constants as rc  # noqa: E402

TRAJ = ROOT / "artifacts" / "m4_traj"
OUT = ROOT / "artifacts" / "m4_atom14"
OUT.mkdir(parents=True, exist_ok=True)

PSF = TRAJ / "M4_xan_ortho_and_allo.psf"
DCD = TRAJ / "M4_xan_ortho_and_allo_1.dcd"
NUM_FRAMES = 100  # >= num_frames of the encoder config

RESNAME_MAP = {"HSE": "HIS", "HSD": "HIS", "HSP": "HIS", "HSD": "HIS"}
# charmm atom name -> atom14 atom name (residue-specific quirk only)
CHARMM_ATOM_RENAMES = {"ILE": {"CD": "CD1"}}


def main() -> None:
    u = mda.Universe(str(PSF), str(DCD))
    # protein residues only: segments A+B, drop ACE cap and XAN ligands
    protein = u.select_atoms("segid A B and not resname ACE XAN")
    residues = list(protein.residues)
    L = len(residues)
    print(f"protein residues : {L}")

    # Precompute per-residue atom14 slot -> global atom index.
    atom14_idx = np.full((L, 14), -1, dtype=np.int64)
    one_letters = []
    for j, res in enumerate(residues):
        std = RESNAME_MAP.get(res.resname, res.resname)
        names = rc.restype_name_to_atom14_names.get(std)
        if names is None:
            raise SystemExit(f"non-standard residue {res.resname} in protein selection")
        ol = rc.restype_3to1.get(std)
        if ol is None:
            raise SystemExit(f"no one-letter code for {std}")
        one_letters.append(ol)
        rename = CHARMM_ATOM_RENAMES.get(std, {})
        for slot, nm in enumerate(names):
            if not nm:
                continue
            for atom in res.atoms:
                an = rename.get(atom.name, atom.name)
                if an == nm:
                    atom14_idx[j, slot] = atom.index
                    break

    n_frames = min(NUM_FRAMES, u.trajectory.n_frames)
    atom14 = np.zeros((n_frames, L, 14, 3), dtype=np.float32)
    for fi, ts in enumerate(u.trajectory[:n_frames]):
        pos = ts.positions  # [N_atoms, 3]
        valid = atom14_idx >= 0
        idx = np.where(valid, atom14_idx, 0)
        atom14[fi] = pos[idx] * valid[..., None]

    missing = (atom14_idx == -1).sum(axis=1)
    print(f"frames          : {n_frames}")
    print(f"atoms missing/残基 (应为 GLY 的 CB 等) : {missing.min()}-{missing.max()}")

    npy = OUT / "m4xan_i1.npy"
    np.save(npy, atom14)
    print(f"wrote {npy}  shape={atom14.shape}  dtype={atom14.dtype}")

    seqres = "".join(one_letters)
    csv = OUT / "m4xan_test.csv"
    with csv.open("w", encoding="utf-8") as fh:
        fh.write("name,seqres\n")
        fh.write(f"m4xan,{seqres}\n")
    print(f"wrote {csv}  seqres length={len(seqres)}")
    print("first 20:", seqres[:20])
    print("CONVERT_DONE")


if __name__ == "__main__":
    main()
