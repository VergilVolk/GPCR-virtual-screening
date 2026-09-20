"""Build auditable DeepRLI complex inputs from existing M4 Vina poses."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
OUT = PROJECT / "data" / "deeprli_m4_v01"
STATES = {
    "7TRQ": (PROJECT / "results" / "pacer_structure_loso_v01" / "poses",
             PROJECT / "results" / "structure" / "7trq_R_receptor.pdb"),
    "7TRP": (PROJECT / "results" / "pacer_structure_7trp_v01" / "poses",
             PROJECT / "results" / "structure" / "ensemble" / "7TRP_R.pdb"),
    "7TRS": (PROJECT / "results" / "pacer_structure_7trs_v01" / "poses",
             PROJECT / "results" / "structure" / "ensemble" / "7TRS_R.pdb"),
}


def pose_path(pose_dir: Path, mol_id: str, seed: int = 42):
    key = hashlib.sha1(f"{mol_id}|{seed}".encode()).hexdigest()[:14]
    return pose_dir / f"{key}.pose.pdbqt"


def pdbqt_heavy(path: Path):
    serial_to_atom, smiles_to_serial = {}, {}
    amap = {"A": "C", "C": "C", "N": "N", "NA": "N", "O": "O", "OA": "O",
            "S": "S", "SA": "S", "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I"}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("REMARK SMILES IDX"):
            values = [int(x) for x in line.split()[3:]]
            for smiles_idx, serial in zip(values[::2], values[1::2]):
                smiles_to_serial[smiles_idx] = serial
        if not line.startswith(("ATOM", "HETATM")):
            continue
        serial = int(line[6:11])
        atom_type = line.split()[-1].upper()
        if atom_type in {"H", "HD", "HS"}:
            continue
        serial_to_atom[serial] = (amap.get(atom_type, atom_type.title()),
                                  [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not smiles_to_serial:
        raise ValueError("Meeko SMILES-IDX map missing")
    return smiles_to_serial, serial_to_atom


def ligand_from_smiles_pose(smiles: str, pose: Path):
    # Isotopic hydrogens are absent from the DeepRLI heavy-atom graph and from
    # Meeko's docked coordinate list; normalize D to H while preserving indices.
    params = Chem.SmilesParserParams(); params.removeHs = False
    original = Chem.MolFromSmiles(smiles.replace("[2H]", "[H]"), params)
    if original is None:
        raise ValueError("SMILES parse failed")
    smiles_to_serial, serial_to_atom = pdbqt_heavy(pose)
    heavy_original = [i for i, atom in enumerate(original.GetAtoms()) if atom.GetAtomicNum() != 1]
    ordered = [serial_to_atom[smiles_to_serial[i + 1]] for i in heavy_original]
    xyz = np.asarray([x[1] for x in ordered], np.float64)
    elements = [x[0] for x in ordered]
    mol = Chem.RemoveHs(original)
    expected = [a.GetSymbol() for a in mol.GetAtoms()]
    if len(xyz) != mol.GetNumAtoms() or elements != expected:
        raise ValueError(f"pose/SMILES atom-order mismatch: pose={elements} smiles={expected}")
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, point in enumerate(xyz):
        conf.SetAtomPosition(i, point)
    mol.RemoveAllConformers(); mol.AddConformer(conf)
    return mol, xyz


def receptor_residue_blocks(path: Path):
    blocks, current_key, current_lines, current_xyz = [], None, [], []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines(True):
        if not line.startswith(("ATOM", "HETATM")) or line[17:20] == "HOH":
            continue
        element = line[76:78].strip().upper()
        if element == "H" or line[12:16].strip().startswith("H"):
            continue
        key = (line[21], line[22:26], line[26])
        if current_key is not None and key != current_key:
            blocks.append((current_lines, np.asarray(current_xyz, float)))
            current_lines, current_xyz = [], []
        current_key = key
        current_lines.append(line)
        current_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if current_lines:
        blocks.append((current_lines, np.asarray(current_xyz, float)))
    return blocks


def pocket_text(blocks, ligand_xyz, cutoff=6.5):
    selected = []
    for lines, xyz in blocks:
        if np.linalg.norm(xyz[:, None, :] - ligand_xyz[None, :, :], axis=2).min() < cutoff:
            selected.extend(lines)
    return "".join(selected) + "END\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--ids", default="")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--clean", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    if args.clean and out.exists():
        shutil.rmtree(out)
    (out / "raw").mkdir(parents=True, exist_ok=True)
    (out / "index").mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(DATA)
    if args.ids:
        requested = {x.strip() for x in args.ids.split(",") if x.strip()}
        data = data[data.canonical_molecule_id.astype(str).isin(requested)]
    if args.limit:
        data = data.head(args.limit)
    rows, audit_rows = [], []
    blocks = {state: receptor_residue_blocks(receptor) for state, (_, receptor) in STATES.items()}
    for state, (pose_dir, _) in STATES.items():
        for row in data.itertuples():
            mol_id = str(row.canonical_molecule_id)
            complex_path = f"{state}/{mol_id}"
            destination = out / "raw" / complex_path
            destination.mkdir(parents=True, exist_ok=True)
            pose = pose_path(pose_dir, mol_id)
            try:
                mol, xyz = ligand_from_smiles_pose(row.canonical_smiles, pose)
                writer = Chem.SDWriter(str(destination / "ligand.sdf")); writer.write(mol); writer.close()
                pocket = pocket_text(blocks[state], xyz)
                (destination / "pocket_6.50A.pdb").write_text(pocket, encoding="utf-8")
                n_res = len({(x[21], x[22:26], x[26]) for x in pocket.splitlines() if x.startswith(("ATOM", "HETATM"))})
                rows.append({"complex_path": complex_path, "pKd": row.pEC50, "DeltaG": np.nan})
                audit_rows.append({"state": state, "canonical_molecule_id": mol_id, "status": "ok",
                                   "ligand_atoms": mol.GetNumAtoms(), "pocket_residues": n_res})
            except Exception as exc:
                audit_rows.append({"state": state, "canonical_molecule_id": mol_id,
                                   "status": "failed", "error": repr(exc)[:500]})
    pd.DataFrame(rows).to_csv(out / "index" / "complexes.csv", index=False)
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(out / "preparation_audit.csv", index=False)
    summary = {"requested_molecules": len(data), "requested_complexes": len(data) * len(STATES),
               "prepared_complexes": len(rows), "failed_complexes": int((audit.status == "failed").sum()),
               "states": list(STATES), "seed": 42, "pocket_cutoff_A": 6.5,
               "bond_order_source": "canonical_smiles", "coordinates_source": "Vina pose"}
    (out / "audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if summary["failed_complexes"]:
        print(audit[audit.status == "failed"].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
