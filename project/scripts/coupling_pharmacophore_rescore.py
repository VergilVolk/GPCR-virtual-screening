"""
变构枢纽药效团重打分 — 方向 1 的方法落地
==========================================
机制发现: VU0467154 与 iperoxo 相距 8.03 A 不直接接触, TYR439 是唯一共享残基
(变构枢纽)。因此"PAM 效力"在结构上要求分子占据别构口袋并接触传导枢纽。

方法: 对基准分子(10 PAM + 10 非活性)重对接, 计算每个 pose 的:
  - 接触 TYR439 的原子数 (枢纽接触)
  - 接触别构口袋残基数 (13 残基命中)
  - 接触双口袋残基数 (26 残基命中)
  - 与 VU0467154 的接触模式 Tanimoto (IFP)
用这些特征重打分, 对比原始 Vina 的 ROC AUC (0.68) 是否提升
=> "从机制到方法": 打分失效的修复方向不是改力场, 而是编码变构耦合的先验

用法: python scripts/coupling_pharmacophore_rescore.py
输出: results/structure/benchmark_coupling_rescore.json
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

# 从 coupling_analysis.json 读口袋
COUPLING = OUT / "coupling_analysis.json"

# 别构口袋残基 (VU0467154 接触) — 从 coupling 报告
ALLOSTERIC = [(89, 'TYR'), (92, 'TYR'), (93, 'ILE'), (96, 'GLY'), (184, 'GLN'),
              (186, 'PHE'), (190, 'LEU'), (423, 'ASN'), (432, 'ASP'), (433, 'THR'),
              (435, 'TRP'), (436, 'SER'), (439, 'TYR')]
ORTHOSTERIC = [(112, 'ASP'), (113, 'TYR'), (116, 'SER'), (117, 'ASN'), (120, 'VAL'),
               (164, 'TRP'), (203, 'ALA'), (204, 'PHE'), (413, 'TRP'), (416, 'TYR'),
               (417, 'ASN'), (439, 'TYR'), (442, 'CYS'), (443, 'TYR')]
HUB = (439, 'TYR')  # 变构枢纽

CUTOFF = 4.5


def parse_atoms(pdb_path):
    out = []
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")):
            res = line[17:20].strip()
            el = line[76:78].strip() or line[12:16].strip()[0]
            if el == "H":
                continue
            out.append((float(line[30:38]), float(line[38:46]), float(line[46:54]),
                        res, int(line[22:26])))
    return out


def pocket_contacts(lig_atoms, pocket):
    """配体对口袋残基的接触: {残基号: 接触原子数}."""
    hits = {}
    for x, y, z, _, _ in lig_atoms:
        for rnum, rname in pocket:
            for line in open(PDB):
                pass  # 直接下面从 PDB 取残基原子 (避免重复读, 改为传入)
    return hits


def residue_atoms(prot_all, rnum):
    return [(x, y, z) for x, y, z, rname, rn in prot_all if rn == rnum]


def contacts_for(lig_atoms, prot_all, pocket):
    hits = {}
    for x, y, z, _, _ in lig_atoms:
        for rnum, rname in pocket:
            for ax, ay, az in residue_atoms(prot_all, rnum):
                if (x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2 < CUTOFF ** 2:
                    hits[rnum] = hits.get(rnum, 0) + 1
                    break
    return hits


def prepare_and_dock(smiles, tmp, key):
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=42)
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
                 if not l.strip().startswith(("exhaustiveness", "num_modes"))]
    cfg_lines += ["exhaustiveness = 8", "num_modes = 1", "seed = 42"]
    cfg = tmp / f"{key}.cfg"
    cfg.write_text("\n".join(cfg_lines) + "\n", encoding="utf-8")
    out_pdbqt = tmp / f"{key}_out.pdbqt"
    r = subprocess.run([str(VINA), "--receptor", str(OUT / "7trq_R_receptor.pdbqt"),
                        "--ligand", str(lig), "--config", str(cfg), "--out", str(out_pdbqt)],
                       capture_output=True, text=True)
    if r.returncode != 0 or not out_pdbqt.exists():
        return None
    # 解析 pose 坐标 (pdbqt MODEL 1, 取 ATOM 行)
    pose_atoms = []
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
            pose_atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54]), "LIG", 0))
    # 亲和力
    aff = None
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            aff = float(parts[1])
            break
    return aff, pose_atoms


def roc_auc(y, s):
    n_pos = y.sum()
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    auc = 0.0
    for i in range(len(y)):
        for j in range(len(y)):
            if y[i] == 1 and y[j] == 0:
                if s[i] > s[j]:
                    auc += 1
                elif s[i] == s[j]:
                    auc += 0.5
    return auc / (n_pos * n_neg)


def main():
    bench = pd.read_csv(OUT / "benchmark_pam_vs_inactive_raw.csv")
    prot_all = parse_atoms(PDB)
    # VU0467154 接触模式 (IFP 参考)
    vu_atoms = [(x, y, z, "LIG", 0) for x, y, z, rname, rn in parse_atoms(PDB) if rname == "IUI"]
    vu_contacts = set(contacts_for(vu_atoms, prot_all, ALLOSTERIC).keys())

    tmp = OUT / "_tmp_rescore"
    tmp.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in bench.itertuples():
        key = hashlib.md5(r.smiles.encode()).hexdigest()[:10]
        res = prepare_and_dock(r.smiles, tmp, key)
        if res is None:
            print(f"  {r.name}: 对接失败")
            continue
        aff, pose = res
        c_allo = contacts_for(pose, prot_all, ALLOSTERIC)
        c_ortho = contacts_for(pose, prot_all, ORTHOSTERIC)
        hub_n = c_allo.get(HUB[0], 0)
        n_allo = len(c_allo)
        n_dual = len(c_allo) + len(c_ortho) - (1 if HUB[0] in c_allo and HUB[0] in c_ortho else 0)
        # IFP vs VU0467154 (别构口袋内残基命中集合的 Tanimoto)
        set_allo = set(c_allo.keys())
        union = set_allo | vu_contacts
        inter = set_allo & vu_contacts
        ifp = len(inter) / len(union) if union else 0.0
        rows.append({"name": r.name, "label": int(r.label), "vina": aff,
                     "hub_contacts": hub_n, "n_allosteric": n_allo,
                     "n_dual_pocket": n_dual, "ifp_vs_VU": ifp})
        print(f"  {r.name}: vina={aff:.2f} hub={hub_n} allo={n_allo} dual={n_dual} ifp={ifp:.2f}")

    df = pd.DataFrame(rows)
    y = df["label"].to_numpy()
    report = {}
    for col, desc in [("vina", "原始 Vina"), ("hub_contacts", "TYR439 枢纽接触"),
                      ("n_allosteric", "别构口袋残基接触数"),
                      ("n_dual_pocket", "双口袋残基接触数"), ("ifp_vs_VU", "IFP vs VU0467154")]:
        s = df[col].to_numpy(dtype=float)
        report[desc] = round(roc_auc(y, s), 3)
        print(f"  ROC AUC [{desc}]: {report[desc]:.3f}")

    out = OUT / "benchmark_coupling_rescore.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    df.to_csv(OUT / "benchmark_coupling_rescore_raw.csv", index=False)
    print(f"\n已保存 -> {out}")


if __name__ == "__main__":
    main()
