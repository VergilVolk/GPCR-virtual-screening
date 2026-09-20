# -*- coding: utf-8 -*-
"""CPU ensemble docking against the published 10-cluster M4R GaMD ensemble."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
P = Path(__file__).resolve().parents[1]
ROOT = P / "results" / "m4_gamd_ensemble"
AUDIT = ROOT / "ensemble_preparation_audit.json"
VINA = P / "tools" / "vina_1.2.7.exe"
AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HID","HIE","HIP","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}


def pdb_atoms(path, protein_only=False, first_model_only=False):
    atoms = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if first_model_only and line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        if protein_only and resname not in AA:
            continue
        element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
        if element == "H":
            continue
        atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54]), int(line[22:26]), resname))
    return atoms


def native_pocket(cluster):
    rec = pdb_atoms(ROOT / "receptors" / f"cluster_{cluster:02d}_receptor.pdb", protein_only=True)
    lig = pdb_atoms(ROOT / "receptors" / f"cluster_{cluster:02d}_MK97.pdb")
    rx = np.asarray([x[:3] for x in rec]); lx = np.asarray([x[:3] for x in lig])
    close = np.linalg.norm(rx[:, None, :] - lx[None, :, :], axis=2).min(1) <= 5.0
    return sorted({rec[i][3] for i in np.where(close)[0]})


def prepare_ligand(smiles, path, seed=42):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3(); params.randomSeed = seed
    if AllChem.EmbedMolecule(m, params) != 0:
        raise ValueError("3D embedding failed")
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=300)
    except Exception:
        pass
    setup = MoleculePreparation().prepare(m)[0]
    text, ok, error = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise ValueError(str(error))
    path.write_text(text, encoding="utf-8")


def enumerate_ph7_states(smiles, max_states=16):
    """Molscrub 0.1.1 protonation/tautomer enumeration used by the reference protocol."""
    from scrubber import Scrub
    scrub = Scrub(ph_low=7.0, ph_high=7.0, skip_gen3d=True)
    source = Chem.MolFromSmiles(smiles)
    states = scrub(source) if source is not None else []
    unique = {}
    for mol in states:
        try:
            smi = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)
            if Chem.MolFromSmiles(smi) is not None:
                unique.setdefault(smi, None)
        except Exception:
            continue
    if not unique:
        unique[Chem.MolToSmiles(source, canonical=True, isomericSmiles=True)] = None
    return list(unique)[:max_states]


def pose_features(path, cluster, pocket):
    # Vina writes several MODEL blocks; geometry must describe rank-1 only.
    lig = np.asarray([x[:3] for x in pdb_atoms(path, first_model_only=True)], float)
    rec = pdb_atoms(ROOT / "receptors" / f"cluster_{cluster:02d}_receptor.pdb", protein_only=True)
    rx = np.asarray([x[:3] for x in rec]); rn = np.asarray([x[3] for x in rec])
    d = np.linalg.norm(lig[:, None, :] - rx[None, :, :], axis=2)
    contacted = sorted({int(x) for x in rn[np.any(d <= 4.5, axis=0)]})
    overlap = set(contacted) & set(pocket)
    return {
        "pose_centroid_x": float(lig[:, 0].mean()), "pose_centroid_y": float(lig[:, 1].mean()),
        "pose_centroid_z": float(lig[:, 2].mean()), "n_contacted_residues": len(contacted),
        "native_pocket_coverage": len(overlap) / len(pocket) if pocket else 0.0,
        "contacted_residues": ";".join(map(str, contacted)),
    }


def one(task):
    molecule_id, state_id, state_smiles, cluster, receptor, ligand, center, size, pose, exhaustiveness, num_modes, pocket = task
    try:
        cmd = [str(VINA), "--receptor", receptor, "--ligand", ligand,
               "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
               "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
               "--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes), "--cpu", "1",
               "--energy_range", "3.0", "--seed", "42", "--out", pose]
        run = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        if run.returncode:
            raise RuntimeError((run.stdout + run.stderr)[-800:])
        affinity = np.nan
        for line in Path(pose).read_text(errors="ignore").splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                affinity = float(line.split()[3]); break
        return {"molecule_id": molecule_id, "state_id": state_id, "state_smiles": state_smiles,
                "cluster": cluster, "vina_affinity": affinity,
                **pose_features(pose, cluster, pocket), "error": ""}
    except Exception as exc:
        return {"molecule_id": molecule_id, "state_id": state_id, "state_smiles": state_smiles,
                "cluster": cluster, "error": repr(exc)[:800]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["validation", "candidates", "custom"], required=True)
    ap.add_argument("--source", type=Path, help="CSV source for --set custom")
    ap.add_argument("--id-col", default="compound_id", help="molecule id column for --set custom")
    ap.add_argument("--custom-name", default="custom", help="output prefix for --set custom")
    ap.add_argument("--validation-scope", choices=["primary", "full"], default="primary")
    ap.add_argument("--clusters", default="0-9")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--run-name", default="protocol_v127_ph7_rank1")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    if args.set == "validation":
        source = ROOT / ("validation_set_primary.csv" if args.validation_scope == "primary" else "validation_set.csv")
        idcol = "canonical_molecule_id"
    elif args.set == "candidates":
        source = P / "results" / "pacer_candidates_v01" / "final" / "final_candidate_hypotheses.csv"
        idcol = "candidate_id"
    else:
        if args.source is None:
            ap.error("--source is required with --set custom")
        source = args.source
        idcol = args.id_col
    df = pd.read_csv(source)
    if args.limit:
        df = df.head(args.limit)
    if "-" in args.clusters:
        a, b = map(int, args.clusters.split("-")); clusters = list(range(a, b + 1))
    else:
        clusters = [int(x) for x in args.clusters.split(",")]
    meta = json.loads(AUDIT.read_text())
    cmeta = {x["cluster"]: x for x in meta["clusters"]}
    if args.set == "validation":
        set_name = f"validation_{args.validation_scope}"
    elif args.set == "custom":
        set_name = args.custom_name
    else:
        set_name = args.set
    out = ROOT / f"{set_name}_docking_{args.run_name}"; ligdir = out / "ligands"; posedir = out / "poses"
    ligdir.mkdir(parents=True, exist_ok=True); posedir.mkdir(parents=True, exist_ok=True)
    ledger = out / "ledger.jsonl"
    old = [json.loads(x) for x in ledger.read_text().splitlines() if x] if ledger.exists() else []
    done = {(str(x["molecule_id"]), int(x.get("state_id", 0)), int(x["cluster"]))
            for x in old if not x.get("error")}
    ligand_paths = {}
    for row in df.itertuples():
        mid = str(getattr(row, idcol)); key = hashlib.sha1(mid.encode()).hexdigest()[:12]
        ligand_paths[mid] = []
        for sid, state_smiles in enumerate(enumerate_ph7_states(row.canonical_smiles)):
            lp = ligdir / f"{key}_s{sid:02d}.pdbqt"
            if not lp.exists():
                prepare_ligand(state_smiles, lp, seed=42 + sid)
            ligand_paths[mid].append((sid, state_smiles, lp))
    pockets = {c: native_pocket(c) for c in clusters}
    tasks = []
    for row in df.itertuples():
        mid = str(getattr(row, idcol))
        for sid, state_smiles, ligand_path in ligand_paths[mid]:
            for c in clusters:
                if (mid, sid, c) in done:
                    continue
                cm = cmeta[c]
                pose = posedir / f"{hashlib.sha1(mid.encode()).hexdigest()[:12]}_s{sid:02d}_c{c:02d}.pdbqt"
                tasks.append((mid, sid, state_smiles, c, cm["pdbqt"], str(ligand_path),
                              cm["box_center"], cm["box_size"], str(pose), args.exhaustiveness,
                              args.num_modes, pockets[c]))
    print(f"set={args.set} molecules={len(df)} clusters={clusters} scheduled={len(tasks)}", flush=True)
    start = time.time(); completed = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(one, task) for task in tasks]):
            result = future.result(); completed += 1
            with ledger.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(result) + "\n")
            if completed % 20 == 0 or completed == len(tasks):
                print(f"completed={completed}/{len(tasks)} elapsed_s={time.time()-start:.0f}", flush=True)
    rows = [json.loads(x) for x in ledger.read_text().splitlines() if x]
    states = pd.DataFrame(rows).drop_duplicates(["molecule_id", "state_id", "cluster"], keep="last")
    states.to_csv(out / "statewise_scores.csv", index=False)
    valid = states[(states.error.isna() | states.error.eq("")) & states.vina_affinity.notna()].copy()
    best_idx = valid.groupby(["molecule_id", "cluster"]).vina_affinity.idxmin()
    valid.loc[best_idx].sort_values(["molecule_id", "cluster"]).to_csv(out / "framewise_scores.csv", index=False)
    run_meta = {"set": args.set, "source": str(source),
                "validation_scope": args.validation_scope if args.set == "validation" else None,
                "n_molecules": len(df), "clusters": clusters, "exhaustiveness": args.exhaustiveness,
                "num_modes": args.num_modes, "vina_version": "1.2.7", "workers": args.workers,
                "ligand_states": "molscrub 0.1.1, pH 7.0, protonation and tautomer enumeration; best state per cluster",
                "box": "per-cluster co-simulated MK-97 centroid, 30 A cube",
                "interpretation": "ensemble pose/affinity evidence only; not PAM functional efficacy"}
    (out / "metadata.json").write_text(json.dumps(run_meta, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
