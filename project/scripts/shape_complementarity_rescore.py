"""
形状/接触互补打分 — 结合判定的物理版本 (不经过力场)
=====================================================
对 docked pose 直接计算几何度量 (RDKit/坐标, 与 Vina 误差正交):

  1. 接触原子对数: 配体-受体原子对 < 4.5 A (埋藏程度)
  2. 埋藏比例: 配体原子中被接触的比例
  3. 平均接触密度: 每配体原子接触的受体原子数 (形状互补代理)
  4. 接触残基数 (已有, 补充)

对比原始 Vina 的 ROC AUC, 并测试 Consensus (Vina + 形状 + 枢纽) 集成。

用法: python scripts/shape_complementarity_rescore.py
输出: results/structure/benchmark_shape_rescore.json
"""
from __future__ import annotations
import sys
import json
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PDB = Path("data/pdb/7TRQ.pdb")
OUT = Path("results/structure")
POSEDIR = OUT / "_tmp_rescore"
CUTOFF = 4.5


def parse_atoms_heavy(pdb_path, exclude_res=("IUI", "IXO")):
    out = []
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")):
            res = line[17:20].strip()
            if res in exclude_res:
                continue
            el = line[76:78].strip() or line[12:16].strip()[0]
            if el == "H":
                continue
            out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])), )
    return out


def parse_pdbqt_pose(pdbqt_path):
    """vina 输出 pdbqt 的 MODEL 1 重原子坐标."""
    atoms = []
    in_model = 0
    for line in open(pdbqt_path):
        if line.startswith("MODEL"):
            in_model += 1
            if in_model > 1:
                break
            continue
        if in_model == 1 and line.startswith("ATOM"):
            el = line[77:79].strip() or line[13:16].strip().lstrip("0123456789")[:1]
            if el == "H":
                continue
            atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])), )
    return atoms


def shape_metrics(lig_atoms, rec_atoms, cutoff=CUTOFF):
    """几何度量: 接触对 / 埋藏比例 / 平均密度."""
    n_lig = len(lig_atoms)
    if n_lig == 0:
        return 0, 0.0, 0.0
    # 受体原子快查 (粗网格哈希)
    contacts_per_lig = [0] * n_lig
    rec_list = rec_atoms
    for i, (lx, ly, lz) in enumerate(lig_atoms):
        for rx, ry, rz in rec_list:
            if (lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2 < cutoff ** 2:
                contacts_per_lig[i] += 1
    n_contacts = sum(contacts_per_lig)
    buried_ratio = sum(1 for c in contacts_per_lig if c > 0) / n_lig
    density = n_contacts / n_lig
    return n_contacts, buried_ratio, density


def roc_auc(y, s):
    n_pos = y.sum()
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    auc = 0.0
    for i in range(len(y)):
        for j in range(len(y)):
            if y[i] == 1 and y[j] == 0:
                auc += 1 if s[i] > s[j] else (0.5 if s[i] == s[j] else 0)
    return auc / (n_pos * n_neg)


def main():
    bench = pd.read_csv(OUT / "benchmark_pam_vs_inactive_raw.csv")
    rec_atoms = parse_atoms_heavy(PDB)

    rows = []
    for r in bench.itertuples():
        key = hashlib.md5(r.smiles.encode()).hexdigest()[:10]
        pose_file = POSEDIR / f"{key}_out.pdbqt"
        if not pose_file.exists():
            print(f"  {r.name}: pose 缺失 {pose_file.name}")
            continue
        lig = parse_pdbqt_pose(pose_file)
        n_cont, buried, density = shape_metrics(lig, rec_atoms)
        rows.append({"name": r.name, "label": int(r.label), "vina": r.vina_affinity,
                     "n_contacts": n_cont, "buried_ratio": round(buried, 3),
                     "density": round(density, 3)})
        print(f"  {r.name}: vina={r.vina_affinity:.2f} contacts={n_cont} buried={buried:.2f} density={density:.1f}")

    df = pd.DataFrame(rows)
    y = df["label"].to_numpy()
    # 与枢纽特征合并 (从 coupling_rescore_raw)
    hub = pd.read_csv(OUT / "benchmark_coupling_rescore_raw.csv")[["name", "hub_contacts"]]
    df = df.merge(hub, on="name", how="left")

    report = {}
    for col, desc in [("vina", "原始 Vina"), ("n_contacts", "接触原子对数"),
                      ("buried_ratio", "埋藏比例"), ("density", "接触密度"),
                      ("hub_contacts", "TYR439 枢纽接触")]:
        s = df[col].to_numpy(dtype=float)
        report[desc] = round(roc_auc(y, s), 3)
        print(f"  ROC AUC [{desc}]: {report[desc]:.3f}")

    # Consensus: 秩平均 (Vina + 枢纽 + 埋藏)
    df["r_vina"] = df["vina"].rank(ascending=False)
    df["r_hub"] = df["hub_contacts"].rank(ascending=False)
    df["r_buried"] = df["buried_ratio"].rank(ascending=False)
    df["r_contacts"] = df["n_contacts"].rank(ascending=False)
    df["consensus3"] = df["r_vina"] + df["r_hub"] + df["r_buried"]
    df["consensus4"] = df["r_vina"] + df["r_hub"] + df["r_buried"] + df["r_contacts"]
    for col, desc in [("consensus3", "Consensus (Vina+枢纽+埋藏)"),
                      ("consensus4", "Consensus (Vina+枢纽+埋藏+接触)")]:
        s = df[col].to_numpy(dtype=float)
        report[desc] = round(roc_auc(y, -s), 3)  # 秩小=好
        print(f"  ROC AUC [{desc}]: {report[desc]:.3f}")

    out = OUT / "benchmark_shape_rescore.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    df.to_csv(OUT / "benchmark_shape_rescore_raw.csv", index=False)
    print(f"\n已保存 -> {out}")


if __name__ == "__main__":
    main()
