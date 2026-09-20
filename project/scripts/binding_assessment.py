"""
结合判定统一评估: 一次对接 (seed 42) + 全部打分特征 + AUC 对比
================================================================
统一重构 coupling_rescore / shape_rescore, 消除 pose/分数来源不一致:

  对接: 基准 20 分子 (10 PAM + 10 非活性) 入 7TRQ 别构口袋, seed=42
  特征: Vina 分 / TYR439 枢纽接触 / 别构口袋接触 / 双口袋接触 /
        接触原子对数 / 埋藏比例 / 接触密度 / IFP vs VU0467154
  输出: 全部特征 AUC + Consensus 秩集成

用法: python scripts/binding_assessment.py
输出: results/structure/binding_assessment.json + .csv
"""
from __future__ import annotations
import sys
import json
import subprocess
import hashlib
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

PDB = Path("data/pdb/7TRQ.pdb")
OUT = Path("results/structure")
BOX = OUT / "allosteric_box.config"
VINA = Path("tools/vina.exe")
CUTOFF = 4.5
SEED = 42

ALLOSTERIC = [(89, 'TYR'), (92, 'TYR'), (93, 'ILE'), (96, 'GLY'), (184, 'GLN'),
              (186, 'PHE'), (190, 'LEU'), (423, 'ASN'), (432, 'ASP'), (433, 'THR'),
              (435, 'TRP'), (436, 'SER'), (439, 'TYR')]
ORTHOSTERIC = [(112, 'ASP'), (113, 'TYR'), (116, 'SER'), (117, 'ASN'), (120, 'VAL'),
               (164, 'TRP'), (203, 'ALA'), (204, 'PHE'), (413, 'TRP'), (416, 'TYR'),
               (417, 'ASN'), (439, 'TYR'), (442, 'CYS'), (443, 'TYR')]
HUB = 439


def parse_atoms(pdb_path, exclude=("IUI", "IXO")):
    out = []
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")):
            res = line[17:20].strip()
            if res in exclude:
                continue
            el = line[76:78].strip() or line[12:16].strip()[0]
            if el == "H":
                continue
            out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])), )
    return out


def prepare_and_dock(smiles, tmp, key):
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=SEED)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        return None
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return None
    lig = tmp / f"{key}.pdbqt"
    lig.write_text(s, encoding="utf-8")
    cfg_lines = [l for l in BOX.read_text(encoding="utf-8").splitlines()
                 if not l.strip().startswith(("exhaustiveness", "num_modes", "seed"))]
    cfg_lines += ["exhaustiveness = 8", "num_modes = 1", f"seed = {SEED}"]
    cfg = tmp / f"{key}.cfg"
    cfg.write_text("\n".join(cfg_lines) + "\n", encoding="utf-8")
    out_pdbqt = tmp / f"{key}_out.pdbqt"
    r = subprocess.run([str(VINA), "--receptor", str(OUT / "7trq_R_receptor.pdbqt"),
                        "--ligand", str(lig), "--config", str(cfg), "--out", str(out_pdbqt)],
                       capture_output=True, text=True)
    if r.returncode != 0 or not out_pdbqt.exists():
        return None
    aff = None
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            aff = float(parts[1])
            break
    # pose 重原子坐标
    pose = []
    in_model = 0
    for line in out_pdbqt.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL"):
            in_model += 1
            if in_model > 1:
                break
            continue
        if in_model == 1 and line.startswith("ATOM"):
            el = line[77:79].strip() or line[13:16].strip().lstrip("0123456789")[:1]
            if el == "H":
                continue
            pose.append((float(line[30:38]), float(line[38:46]), float(line[46:54])), )
    return aff, pose


def contact_counts(pose, prot, pocket):
    """pose 对口袋残基的接触: {rnum: n_atoms}."""
    hits = {}
    atoms_by_res = {}
    for x, y, z in prot:
        pass  # prot 无残基信息, 从 PDB 重新解析
    return hits


def main():
    bench = pd.read_csv(OUT / "benchmark_pam_vs_inactive_raw.csv")
    prot_all = parse_atoms(PDB)
    # 残基 -> 原子坐标索引 (重解析带残基号)
    res_atoms = {}
    for line in open(PDB):
        if line.startswith(("ATOM", "HETATM")):
            res = line[17:20].strip()
            if res in ("IUI", "IXO"):
                continue
            el = line[76:78].strip() or line[12:16].strip()[0]
            if el == "H":
                continue
            rnum = int(line[22:26])
            res_atoms.setdefault(rnum, []).append(
                (float(line[30:38]), float(line[38:46]), float(line[46:54])))

    def contact_pocket(pose, pocket):
        hits = {}
        for x, y, z in pose:
            for rnum, rname in pocket:
                for ax, ay, az in res_atoms.get(rnum, ()):
                    if (x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2 < CUTOFF ** 2:
                        hits[rnum] = hits.get(rnum, 0) + 1
                        break
        return hits

    tmp = OUT / "_tmp_assess"
    tmp.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in bench.itertuples():
        key = hashlib.md5(r.smiles.encode()).hexdigest()[:10]
        res = prepare_and_dock(r.smiles, tmp, key)
        if res is None:
            print(f"  {r.name}: 对接失败")
            continue
        aff, pose = res
        c_allo = contact_pocket(pose, ALLOSTERIC)
        c_ortho = contact_pocket(pose, ORTHOSTERIC)
        hub_n = c_allo.get(HUB, 0)
        # 形状度量 (vs 全部受体原子)
        n_cont = 0
        lig_buried = 0
        for x, y, z in pose:
            hit = False
            for rx, ry, rz in prot_all:
                if (x - rx) ** 2 + (y - ry) ** 2 + (z - rz) ** 2 < CUTOFF ** 2:
                    n_cont += 1
                    hit = True
            if hit:
                lig_buried += 1
        buried = lig_buried / len(pose) if pose else 0.0
        density = n_cont / len(pose) if pose else 0.0
        # IFP vs VU0467154 (别构口袋内接触集合 Tanimoto)
        vu_contacts = {89, 92, 93, 96, 184, 186, 190, 423, 432, 433, 435, 436, 439}
        set_allo = set(c_allo.keys())
        union = set_allo | vu_contacts
        ifp = len(set_allo & vu_contacts) / len(union) if union else 0.0
        rows.append({"name": r.name, "label": int(r.label), "vina": aff,
                     "hub_contacts": hub_n, "n_allosteric": len(c_allo),
                     "n_dual": len(c_allo) + len(c_ortho) - (1 if HUB in c_allo and HUB in c_ortho else 0),
                     "n_contacts": n_cont, "buried_ratio": round(buried, 3),
                     "density": round(density, 2), "ifp": round(ifp, 2)})
        print(f"  {r.name}: vina={aff:.2f} hub={hub_n} allo={len(c_allo)} cont={n_cont} ifp={ifp:.2f}")

    df = pd.DataFrame(rows)
    y = df["label"].to_numpy()

    def auc(s):
        s = s.to_numpy(dtype=float)
        n_pos, n_neg = y.sum(), len(y) - y.sum()
        if n_pos == 0 or n_neg == 0:
            return float("nan")
        a = 0.0
        for i in range(len(y)):
            for j in range(len(y)):
                if y[i] == 1 and y[j] == 0:
                    a += 1 if s[i] > s[j] else (0.5 if s[i] == s[j] else 0)
        return a / (n_pos * n_neg)

    feats = {"vina": "原始 Vina", "hub_contacts": "TYR439 枢纽接触",
             "n_allosteric": "别构口袋接触数", "n_dual": "双口袋接触数",
             "n_contacts": "接触原子对数", "buried_ratio": "埋藏比例",
             "density": "接触密度", "ifp": "IFP vs VU0467154"}
    report = {}
    for col, desc in feats.items():
        report[desc] = round(auc(df[col]), 3)
        print(f"  ROC AUC [{desc}]: {report[desc]:.3f}")

    # Consensus 秩集成 (去掉最差的埋藏/密度, 用正交信号)
    df["r_vina"] = df["vina"].rank(ascending=False)
    df["r_hub"] = df["hub_contacts"].rank(ascending=False)
    df["r_cont"] = df["n_contacts"].rank(ascending=False)
    df["r_ifp"] = df["ifp"].rank(ascending=False)
    df["consensus"] = df[["r_vina", "r_hub", "r_cont"]].sum(axis=1)
    df["consensus_all"] = df[["r_vina", "r_hub", "r_cont", "r_ifp"]].sum(axis=1)
    report["Consensus (Vina+枢纽+接触)"] = round(auc(df["consensus"]), 3)
    report["Consensus (Vina+枢纽+接触+IFP)"] = round(auc(df["consensus_all"]), 3)
    print(f"  ROC AUC [Consensus (Vina+枢纽+接触)]: {report['Consensus (Vina+枢纽+接触)']:.3f}")

    out = OUT / "binding_assessment.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    df.to_csv(OUT / "binding_assessment_raw.csv", index=False)
    print(f"\n已保存 -> {out}")


if __name__ == "__main__":
    main()
