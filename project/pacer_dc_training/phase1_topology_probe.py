#!/usr/bin/env python3
"""READ-ONLY: identify the receptor selection inside a PACER-DC membrane topology."""
import sys
from pathlib import Path

import MDAnalysis as mda

W = Path("/root/work/GPCR-virtual-screening")
SYS = sys.argv[1] if len(sys.argv) > 1 else "apo"
PDB = W / "project/results/pacer_dc_membrane_reference_v01" / SYS / "minimized.pdb"

u = mda.Universe(str(PDB))
print(f"topology      : {PDB}")
print(f"total atoms   : {len(u.atoms)}")
print(f"segments      : {[ (s.segid, len(s.residues), len(s.atoms)) for s in u.segments ]}")
print()
print("=== per-chain residue / atom counts ===")
from collections import OrderedDict
chains = OrderedDict()
for a in u.atoms:
    chains.setdefault(a.chainID, {"res": set(), "atoms": 0, "names": set()})
    chains[a.chainID]["res"].add(a.resid)
    chains[a.chainID]["atoms"] += 1
for cid, d in chains.items():
    print(f"  chain {cid!r}: {len(d['res']):>5} residues  {d['atoms']:>7} atoms")

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "HSE": "H", "HSD": "H", "HSP": "H", "ILE": "I", "LEU": "L",
    "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
print()
print("=== protein-like chains (standard residues only) + sequence ===")
for cid, d in chains.items():
    sel = u.select_atoms(f"chainID {cid}")
    resnames = sorted({r.resname for r in sel.residues})
    if not all(rn in AA3 for rn in resnames):
        print(f"  chain {cid!r}: not all standard -> {resnames[:8]}{'...' if len(resnames) > 8 else ''}")
        continue
    seq = "".join(AA3[r.resname] for r in sel.residues)
    print(f"  chain {cid!r}: {len(sel.residues)} residues, {len(sel.atoms)} atoms")
    print(f"      seq[:60] = {seq[:60]}")
    print(f"      seq[-30:] = {seq[-30:]}")
    if len(seq) > 250:
        (Path("/tmp") / f"receptor_{SYS}.seq").write_text(seq + "\n")
        print(f"      full sequence written to /tmp/receptor_{SYS}.seq")
