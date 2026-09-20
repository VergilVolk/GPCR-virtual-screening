# -*- coding: utf-8 -*-
"""Generate label-free 7TRQ structural features for the M4 PAM potency set.

Vina is used only to propose poses.  Features are explicit and auditable:
atom-pair contacts, residue contact fingerprints, crystal-IFP Jaccard,
crystal-pose centroid displacement, and TYR439 distance.  TYR439 is treated as
a coupling *hypothesis*, never as a positive label.

The run is resumable.  Usage:
  python scripts/dock_potency_coupling_features.py --workers 4 --exhaustiveness 4
  python scripts/dock_potency_coupling_features.py --limit 24 --seeds 42,43,44
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "benchmarks" / "m4_pam_v1" / "potency_molecules.csv"
STRUCT = PROJECT / "results" / "structure"
VINA = PROJECT / "tools" / "vina.exe"
STATE = os.environ.get("M4_STATE", "7TRQ").upper()
STATE_CONFIG = {
    "7TRQ": {"receptor": STRUCT/"7trq_R_meeko.pdbqt", "pdb": STRUCT/"7trq_R_receptor.pdb",
              "ligand": STRUCT/"7trq_IUI_ligand.pdb", "center": (107.948,85.737,70.412)},
    "7TRP": {"receptor": STRUCT/"ensemble"/"7TRP_R_meeko.pdbqt", "pdb": STRUCT/"ensemble"/"7TRP_R.pdb",
              "ligand": STRUCT/"ensemble"/"7TRP_IUE.pdb", "center": (113.7571,93.2032,74.5355)},
    "7TRS": {"receptor": STRUCT/"ensemble"/"7TRS_R_meeko.pdbqt", "pdb": STRUCT/"ensemble"/"7TRS_R.pdb",
              "ligand": None, "center": (110.2121366,107.8369051,68.1181557)},
}
if STATE not in STATE_CONFIG: raise ValueError(f"unknown M4_STATE={STATE}")
CFG=STATE_CONFIG[STATE]; RECEPTOR=CFG["receptor"]; RECEPTOR_PDB=CFG["pdb"]; CRYSTAL=CFG["ligand"]
OUTDIR = PROJECT / "results" / ("pacer_structure_loso_v01" if STATE=="7TRQ" else f"pacer_structure_{STATE.lower()}_v01")
OUTCSV = OUTDIR / "docking_features.csv"
OUTJSONL = OUTDIR / "docking_features.jsonl"
TMP = OUTDIR / "poses"
POCKET = (89, 92, 93, 96, 184, 186, 190, 423, 432, 433, 435, 436, 439)
HUB = 439
CONTACT_CUT = 4.5
BOX_CENTER = CFG["center"]
BOX_SIZE = (22.0, 22.0, 22.0)


def pdb_atoms(path: Path, first_model=True):
    atoms, in_model, saw_model = [], False, False
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("MODEL"):
            saw_model = True; in_model = True; continue
        if line.startswith("ENDMDL") and first_model:
            break
        if not line.startswith(("ATOM", "HETATM")) or (saw_model and not in_model):
            continue
        atype = line[76:79].strip().upper()
        name = line[12:16].strip()
        if atype in {"H", "HD"} or name.startswith("H"):
            continue
        atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54]),
                      line[17:20].strip(), int(line[22:26])))
    return atoms


def receptor_arrays():
    all_atoms = pdb_atoms(RECEPTOR_PDB)
    xyz = np.asarray([a[:3] for a in all_atoms], dtype=np.float32)
    nums = np.asarray([a[4] for a in all_atoms], dtype=np.int32)
    return xyz, nums


REC_XYZ, REC_NUMS = receptor_arrays()
CRYSTAL_XYZ = np.asarray([a[:3] for a in pdb_atoms(CRYSTAL)], dtype=np.float32) if CRYSTAL else np.empty((0,3))
CRYSTAL_CENTROID = CRYSTAL_XYZ.mean(0) if len(CRYSTAL_XYZ) else np.asarray(BOX_CENTER,dtype=np.float32)


def pose_features(pose_path: Path):
    lig = np.asarray([a[:3] for a in pdb_atoms(pose_path)], dtype=np.float32)
    if not len(lig):
        raise ValueError("empty pose")
    pocket_mask = np.isin(REC_NUMS, POCKET)
    pxyz, pnums = REC_XYZ[pocket_mask], REC_NUMS[pocket_mask]
    dist = np.linalg.norm(lig[:, None, :] - pxyz[None, :, :], axis=2)
    contacts = dist < CONTACT_CUT
    contacted = sorted({int(x) for x in pnums[np.any(contacts, axis=0)]})
    ref = set(POCKET); obs = set(contacted)
    per_res_min = {r: float(dist[:, pnums == r].min()) for r in POCKET}
    per_res_pairs = {r: int(contacts[:, pnums == r].sum()) for r in POCKET}
    feat = {
        "n_heavy_atoms": int(len(lig)),
        "pocket_atom_pair_contacts": int(contacts.sum()),
        "ligand_atoms_contacting_pocket": int(np.any(contacts, axis=1).sum()),
        "unique_pocket_residue_contacts": int(len(contacted)),
        "crystal_ifp_jaccard": float(len(obs & ref) / len(obs | ref)) if obs | ref else 0.0,
        "crystal_centroid_distance_A": float(np.linalg.norm(lig.mean(0) - CRYSTAL_CENTROID)),
        "tyr439_min_distance_A": per_res_min[HUB],
        "tyr439_atom_pair_contacts": per_res_pairs[HUB],
        "contacted_residues": ";".join(map(str, contacted)),
    }
    for r in POCKET:
        feat[f"res{r}_min_A"] = per_res_min[r]
        feat[f"res{r}_pairs"] = per_res_pairs[r]
    return feat


def make_ligand(smiles: str, path: Path, seed: int):
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("invalid SMILES")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3(); params.randomSeed = int(seed)
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError("3D embedding failed")
    try: AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    except Exception: pass
    setup = MoleculePreparation().prepare(mol)[0]
    text, ok, error = PDBQTWriterLegacy.write_string(setup)
    if not ok: raise ValueError(str(error))
    path.write_text(text, encoding="utf-8")


def dock_one(task):
    mol_id, smiles, seed, exhaustiveness = task
    key = hashlib.sha1(f"{mol_id}|{seed}".encode()).hexdigest()[:14]
    lig, pose = TMP / f"{key}.lig.pdbqt", TMP / f"{key}.pose.pdbqt"
    try:
        if not pose.exists():
            make_ligand(smiles, lig, seed)
            cmd = [str(VINA), "--receptor", str(RECEPTOR), "--ligand", str(lig),
                   "--center_x", str(BOX_CENTER[0]), "--center_y", str(BOX_CENTER[1]),
                   "--center_z", str(BOX_CENTER[2]), "--size_x", str(BOX_SIZE[0]),
                   "--size_y", str(BOX_SIZE[1]), "--size_z", str(BOX_SIZE[2]),
                   "--exhaustiveness", str(exhaustiveness), "--num_modes", "1",
                   "--cpu", "1", "--seed", str(seed), "--out", str(pose)]
            run = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if run.returncode != 0: raise RuntimeError(run.stderr[-500:])
        text = pose.read_text(encoding="utf-8", errors="ignore")
        affinity = np.nan
        for line in text.splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                affinity = float(line.split()[3]); break
        return {"canonical_molecule_id": mol_id, "seed": seed,
                "vina_affinity": affinity, **pose_features(pose), "error": ""}
    except Exception as exc:
        return {"canonical_molecule_id": mol_id, "seed": seed, "error": repr(exc)[:500]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--seeds", default="42")
    args = ap.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True); TMP.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA)
    if args.limit: df = df.head(args.limit)
    seeds = [int(x) for x in args.seeds.split(",")]
    done = set()
    if OUTJSONL.exists():
        old_rows = [json.loads(line) for line in OUTJSONL.read_text(encoding="utf-8").splitlines() if line]
        done = {(str(r["canonical_molecule_id"]), int(r["seed"])) for r in old_rows}
    tasks = [(str(r.canonical_molecule_id), r.canonical_smiles, seed, args.exhaustiveness)
             for r in df.itertuples() for seed in seeds
             if (str(r.canonical_molecule_id), seed) not in done]
    print(f"scheduled={len(tasks)} already_done={len(done)} receptor={RECEPTOR.name}", flush=True)
    t0=time.time(); completed=0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(dock_one, t) for t in tasks]
        for future in as_completed(futures):
            row=future.result(); completed += 1
            with OUTJSONL.open("a", encoding="utf-8") as fh:
                fh.write(json.dumps(row, ensure_ascii=False) + "\n")
            if completed % 20 == 0 or completed == len(tasks):
                print(f"completed={completed}/{len(tasks)} elapsed_s={time.time()-t0:.0f}", flush=True)
    meta={"state":STATE,"receptor":str(RECEPTOR),"n_molecules":int(len(df)),"seeds":seeds,
          "exhaustiveness":args.exhaustiveness,"contact_cut_A":CONTACT_CUT,
          "pocket_residues":POCKET,
          "caution":"TYR439 features test a coupling hypothesis; they are not PAM labels."}
    (OUTDIR/"run_metadata.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
    rows=[json.loads(line) for line in OUTJSONL.read_text(encoding="utf-8").splitlines() if line]
    pd.DataFrame(rows).to_csv(OUTCSV,index=False)


if __name__ == "__main__": main()
