"""
CHRM4 别构口袋打分基准: PAM vs 严格非活性 (打分失效实证)
==========================================================
科学问题: 传统对接打分(Vina)能否区分「M4 PAM」与「M4 实验确认非活性分子」?
  - 若 ROC AUC ≈ 0.5: 打分对别构 PAM 失效 -> 支撑「PAM 效力 != 亲和力」核心论点
  - 这是 MRGPRX1 别构 (AUC=0.44) 之外的第 2 个靶点实证 (双靶点结论)

数据集:
  阳性: modeling_potency_exact_calcium.csv pEC50 最高的分子 (PAM)
  阴性: modeling_pam_vs_inactive.csv 的 molecule_strict_label != positive (严格非活性)
口袋: results/structure/allosteric_box.config (7TRQ VU0467154 位点)

用法: python scripts/dock_benchmark_chrm4.py
输出: results/structure/benchmark_pam_vs_inactive.json + 明细 CSV
"""
from __future__ import annotations
import sys
import os
import json
import subprocess
import tempfile
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

DATA = Path(r"D:\CLC\CHRM4_PAM_Modeling_Handoff_v1.0\data")
OUT = Path("results/structure")
BOX = OUT / "allosteric_box.config"
VINA = Path("tools/vina.exe")
N_POS, N_NEG = 10, 10
EXHAUSTIVENESS = 8


def load_positives(n=10):
    df = pd.read_csv(DATA / "modeling_potency_exact_calcium.csv")
    cmp = pd.read_csv(DATA / "compounds.csv")[["canonical_molecule_id", "canonical_smiles"]]
    df = df.merge(cmp, on="canonical_molecule_id", how="inner")
    df = df.drop_duplicates("measurement_group_id")
    agg = df.groupby("canonical_molecule_id").agg(
        pEC50=("aggregated_pEC50", "mean"), smiles=("canonical_smiles", "first")).reset_index()
    top = agg.nlargest(n, "pEC50")
    return [(r.smiles, f"PAM_pEC50_{r.pEC50:.2f}") for r in top.itertuples()]


def load_negatives(n=10):
    df = pd.read_csv(DATA / "modeling_pam_vs_inactive.csv")
    cmp = pd.read_csv(DATA / "compounds.csv")[["canonical_molecule_id", "canonical_smiles"]]
    df = df.merge(cmp, on="canonical_molecule_id", how="inner")
    mol = df.drop_duplicates("canonical_molecule_id")
    neg = mol[mol["molecule_strict_label"] != "positive"]
    neg = neg[neg["canonical_smiles"].apply(lambda s: Chem.MolFromSmiles(s) is not None)]
    picked = neg.head(n)
    return [(r.canonical_smiles, f"INACTIVE_{i+1}") for i, r in enumerate(picked.itertuples())]


def prepare_ligand(smiles, pdbqt_path):
    """meeko SMILES -> PDBQT (单分子)."""
    from meeko import MoleculePreparation
    from meeko.preparation import PDBQTWriterLegacy
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=42)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        return False
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return False
    pdbqt_path.write_text(s, encoding="utf-8")
    return True


def dock(pdbqt_path):
    """Vina 对接, 返回最佳亲和力 (kcal/mol)."""
    out = pdbqt_path.with_suffix(".out.pdbqt")
    cfg_lines = [l for l in BOX.read_text(encoding="utf-8").splitlines()
                 if not l.strip().startswith(("exhaustiveness", "num_modes"))]
    cfg_lines.append(f"exhaustiveness = {EXHAUSTIVENESS}")
    cfg_lines.append("num_modes = 1")
    cfg_lines.append("seed = 42")  # 固定随机种子: n=20 小样本下随机波动极大, 必须可复现
    cfg = "\n".join(cfg_lines) + "\n"
    cfg_path = pdbqt_path.with_suffix(".cfg")
    cfg_path.write_text(cfg, encoding="utf-8")
    r = subprocess.run([str(VINA), "--receptor", str(OUT / "7trq_R_receptor.pdbqt"),
                        "--ligand", str(pdbqt_path), "--config", str(cfg_path),
                        "--out", str(out)], capture_output=True, text=True)
    if r.returncode != 0:
        return None
    # 解析 "   1       -9.736          0          0" (mode 行无单位, 只有表头有 kcal/mol)
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            return float(parts[1])
    return None


def main():
    pos = load_positives(N_POS)
    neg = load_negatives(N_NEG)
    print(f"阳性 PAM: {len(pos)} | 阴性非活性: {len(neg)}")

    rows = []
    tmp = OUT / "_tmp_benchmark"
    tmp.mkdir(parents=True, exist_ok=True)
    for smiles, name in pos + neg:
        pdbqt = tmp / f"{name.replace(' ', '_').replace('.', '_')}.pdbqt"
        if not prepare_ligand(smiles, pdbqt):
            print(f"  {name}: 配体制备失败")
            continue
        aff = dock(pdbqt)
        if aff is None:
            print(f"  {name}: 对接失败")
            continue
        label = 1 if name.startswith("PAM") else 0
        rows.append({"name": name, "label": label, "vina_affinity": aff,
                     "smiles": smiles})
        print(f"  {name}: {aff:.2f} kcal/mol")

    df = pd.DataFrame(rows)
    if len(df) == 0:
        print("无结果")
        return
    df.to_csv(OUT / "benchmark_pam_vs_inactive_raw.csv", index=False)

    # ROC AUC (正类=1 PAM)
    y = df["label"].to_numpy()
    s = -df["vina_affinity"].to_numpy()  # 越负越好 -> 越大越可能是阳性
    order = np.argsort(-s)
    n_pos = y.sum()
    n_neg = len(y) - n_pos
    # 简单 AUC (Wilcoxon-Mann-Whitney)
    auc = 0.0
    for i in range(len(y)):
        for j in range(len(y)):
            if y[i] == 1 and y[j] == 0:
                if s[i] > s[j]:
                    auc += 1
                elif s[i] == s[j]:
                    auc += 0.5
    auc /= (n_pos * n_neg) if n_pos * n_neg else 1
    # EF(20%)
    k = max(1, len(y) // 5)
    top_k = y[order[:k]]
    ef = top_k.sum() / (k * n_pos / len(y)) if n_pos else 0

    result = {"n_pos": int(n_pos), "n_neg": int(n_neg), "ROC_AUC": float(auc),
              "EF_20%": float(ef)}
    (OUT / "benchmark_pam_vs_inactive.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    print("\n===== 结果 =====")
    print(json.dumps(result, indent=2))
    print("\n排序表 (分数从高到低):")
    for i in order:
        mark = "PAM" if y[i] else "非活性"
        print(f"  {i+1:>2}. [{mark}] {df.iloc[i]['name']:<20} {df.iloc[i]['vina_affinity']:.2f}")


if __name__ == "__main__":
    main()
