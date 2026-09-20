# -*- coding: utf-8 -*-
"""DEPRECATED exploratory candidate docking script.

Do not use its fixed 0.8/0.2 fusion or historical AUC comments as validation.
The audited implementation is dock_potency_coupling_features.py and must first
pass leave-one-series-out ablation before any structural score is promoted.

对 PEAM-VS 候选 (算法线 QSAR/生成) 执行:
  embed -> meeko PDBQT -> vina pose (seed=42, 仅 pose 生成)
  -> 几何接触证据打分 (n_contacts / hub_contacts / buried_ratio / pocket_density)
  -> 与 QSAR 排名融合为 PEAM-VS v2 (结构证据入列)

用法: python scripts/dock_candidates_geometry.py [N] [start]
输出: results/qsar/peamvs_candidates_geometry.csv (+ .json 汇总)
"""
from __future__ import annotations
import sys, json, subprocess, hashlib, time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = Path("results/qsar")
STRUCT = Path("results/structure")
BOX = STRUCT / "allosteric_box.config"
VINA = Path("tools/vina.exe")
SEED = 42

# 7TRQ 别构口袋 13 残基 (validate_allosteric_pocket.py 验证, 接触 13/13)
POCKET_RESIDUES = {89, 92, 93, 96, 184, 186, 190, 423, 432, 433, 435, 436, 439}
HUB_RESIDUES = {439}  # TYR439 contact is a hypothesis, not a PAM label.
CONTACT_CUT = 4.5
BURIED_CUT = 6.0


def load_receptor(pdbqt_path):
    atoms = []
    with open(pdbqt_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append({
                    "name": line[12:16].strip(), "res": int(line[22:26]),
                    "resname": line[17:20].strip(),
                    "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
                })
    return atoms


def parse_pose(pdbqt_path):
    """MODEL 1 重原子."""
    atoms = []
    with open(pdbqt_path, encoding="utf-8") as f:
        in_model = False
        for line in f:
            if line.startswith("MODEL"):
                in_model = True
                continue
            if line.startswith("ENDMDL"):
                break
            if in_model and line.startswith(("ATOM", "HETATM")):
                elem = line[76:78].strip().upper()
                if elem == "H":
                    continue
                atoms.append({"x": float(line[30:38]), "y": float(line[38:46]),
                              "z": float(line[46:54])})
    return atoms


def geometric_features(lig_atoms, rec_atoms):
    n_contacts = 0
    hub_contacts = 0
    buried = 0
    for la in lig_atoms:
        min_d = 1e9
        hub_d = 1e9
        for ra in rec_atoms:
            d = ((la["x"] - ra["x"]) ** 2 + (la["y"] - ra["y"]) ** 2 + (la["z"] - ra["z"]) ** 2) ** 0.5
            if d < min_d:
                min_d = d
            if ra["res"] in HUB_RESIDUES and d < hub_d:
                hub_d = d
        if min_d < CONTACT_CUT:
            n_contacts += 1
            if hub_d < CONTACT_CUT:
                hub_contacts += 1
        if min_d < BURIED_CUT:
            buried += 1
    n = max(len(lig_atoms), 1)
    return {"n_contacts": n_contacts, "hub_contacts": hub_contacts,
            "buried_ratio": round(buried / n, 3), "pocket_density": round(n_contacts / n, 3),
            "n_lig_atoms": len(lig_atoms)}


def prepare_and_dock(smiles, tmp):
    """meeko 制备 + vina 对接(仅 pose), 返回 (vina_affinity, 特征)."""
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None:
        return None, None
    if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
        try:
            if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(randomSeed=SEED)) != 0:
                return None, None
        except Exception:
            return None, None
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        return None, None
    pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return None, None
    key = hashlib.md5(smiles.encode()).hexdigest()[:10]
    lig = tmp / f"{key}.pdbqt"
    lig.write_text(pdbqt_str, encoding="utf-8")
    cfg_lines = [l for l in BOX.read_text(encoding="utf-8").splitlines()
                 if not l.strip().startswith(("exhaustiveness", "num_modes"))]
    cfg_lines += [f"exhaustiveness = 8", "num_modes = 1", f"seed = {SEED}"]
    cfg = tmp / f"{key}.cfg"
    cfg.write_text("\n".join(cfg_lines) + "\n", encoding="utf-8")
    out_pdbqt = lig.with_suffix(".out.pdbqt")
    r = subprocess.run([str(VINA), "--receptor", str(STRUCT / "7trq_R_receptor.pdbqt"),
                        "--ligand", str(lig), "--config", str(cfg), "--out", str(out_pdbqt)],
                       capture_output=True, text=True)
    if r.returncode != 0 or not out_pdbqt.exists():
        return None, None
    affinity = None
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            affinity = float(parts[1])
            break
    lig_atoms = parse_pose(out_pdbqt)
    if not lig_atoms:
        return affinity, None
    return affinity, geometric_features(lig_atoms, REC_ATOMS)


def zscore(s):
    return (s - s.mean()) / (s.std() + 1e-9)


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    START = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    SOURCE = sys.argv[3] if len(sys.argv) > 3 else None  # fragment / lstm
    df = pd.read_csv(OUT / "peamvs_candidates_ranked.csv")
    if SOURCE:
        df = df[df["source"] == SOURCE].reset_index(drop=True)
    sub = df.iloc[START:START + N].reset_index(drop=True)
    tmp = STRUCT / "_tmp_geomdock"
    tmp.mkdir(parents=True, exist_ok=True)

    rows = []
    t0 = time.time()
    for i, r in enumerate(sub.itertuples()):
        aff, feats = prepare_and_dock(r.smiles, tmp)
        if feats is None:
            continue
        row = {"smiles": r.smiles, "source": r.source,
               "pEC50_pred": round(r.pEC50_pred, 3), "PAM_prob": round(r.PAM_prob, 3),
               "PEAMVS_score": round(r.PEAMVS_score, 3),
               "vina_affinity_ref": round(aff, 2) if aff else None}
        row.update(feats)
        rows.append(row)
        if (i + 1) % 20 == 0:
            print(f"  [{START+i+1}/{START+N}] {time.time()-t0:.0f}s 成功 {len(rows)}", flush=True)

    if not rows:
        print("无成功对接")
        return
    g = pd.DataFrame(rows)
    # Historical exploratory fusion only; weights were not independently validated.
    g["z_n_contacts"] = zscore(g["n_contacts"])
    g["z_hub"] = zscore(g["hub_contacts"])
    g["geom_score"] = 0.6 * g["z_n_contacts"] + 0.4 * g["z_hub"]
    g["z_peamvs"] = zscore(g["PEAMVS_score"])
    g["PEAMVS_v2"] = 0.8 * g["z_peamvs"] + 0.2 * g["geom_score"]
    g = g.sort_values("PEAMVS_v2", ascending=False).reset_index(drop=True)
    g["rank_v2"] = g.index + 1

    g.to_csv(OUT / "peamvs_candidates_geometry.csv", index=False, encoding="utf-8-sig")
    summary = {
        "n_docked": len(g), "n_total_requested": min(N, len(sub)),
        "mean_vina_ref": round(g["vina_affinity_ref"].mean(), 2) if g["vina_affinity_ref"].notna().any() else None,
        "mean_n_contacts": round(g["n_contacts"].mean(), 1),
        "mean_hub_contacts": round(g["hub_contacts"].mean(), 2),
        "hub_contact_rate": round(float((g["hub_contacts"] > 0).mean()), 3),
        "top10": g.head(10)[["rank_v2", "source", "PEAMVS_score", "n_contacts", "hub_contacts", "PEAMVS_v2"]].to_dict("records"),
        "note": "vina_affinity_ref 仅参考(pose 生成副产品), 排名完全由几何证据+QSAR 决定",
    }
    (OUT / "peamvs_geometry_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    REC_ATOMS = load_receptor(STRUCT / "7trq_R_receptor.pdbqt")
    print(f"受体原子: {len(REC_ATOMS)}")
    main()
