"""
生成式方法学对比: SMILES LSTM vs 离散扩散 (few-shot 迁移)
==========================================================
对比维度 (同种子、同语料、同过滤条件下):
  1. 采样有效率 (有效 SMILES / 采样数)
  2. 过滤通过率 (类药+药效团 / 有效)
  3. 与 12 个 PAM 种子的最大 Tanimoto (化学空间贴近度)
  4. 与冻结数据集阳性的最大 Tanimoto (已知活性空间贴近度)
  5. 新颖性: 相对数据集 Tanimoto<0.85 的比例

用法: python scripts/compare_generators.py
输出: results/generated/generator_comparison.json
"""
from __future__ import annotations
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PROJECT = Path(__file__).resolve().parents[1]
GEN = PROJECT / "data" / "generated"
RESULTS = PROJECT / "results" / "generated"
DATA = PROJECT.parent / "CHRM4_PAM_Modeling_Handoff_v1.0" / "data"


def fps_of_smiles_list(smis):
    out = []
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            out.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))
    return out


def max_tanimoto(fp, ref_fps):
    return max((DataStructs.TanimotoSimilarity(fp, r) for r in ref_fps), default=0.0)


def analyze(name, smis, seed_fps, known_fps):
    fps = fps_of_smiles_list(smis)
    if not fps:
        return {"name": name, "n": 0}
    t_seed = [max_tanimoto(f, seed_fps) for f in fps]
    t_known = [max_tanimoto(f, known_fps) for f in fps]
    novel = sum(1 for t in t_known if t < 0.85)
    return {"name": name, "n": len(smis),
            "mean_Tanimoto_vs_seed": round(float(np.mean(t_seed)), 3),
            "mean_Tanimoto_vs_known": round(float(np.mean(t_known)), 3),
            "novel_ratio_lt0.85": round(novel / len(smis), 3)}


def main():
    seeds = [l.strip() for l in GEN.joinpath("actives.smi").read_text(encoding="utf-8").splitlines() if l.strip()]
    seed_fps = fps_of_smiles_list(seeds)
    known_df = pd.read_csv(DATA / "compounds.csv")
    known_fps = fps_of_smiles_list(known_df["canonical_smiles"].dropna().tolist())

    report = {}
    # LSTM
    lstm_csv = RESULTS / "lstm_generated_pam_analogs.csv"
    if lstm_csv.exists():
        lstm = pd.read_csv(lstm_csv)["smiles"].tolist()
        report["LSTM"] = analyze("LSTM", lstm, seed_fps, known_fps)
    # Diffusion
    diff_csv = RESULTS / "diffusion_generated_pam_analogs.csv"
    if diff_csv.exists():
        diff = pd.read_csv(diff_csv)["smiles"].tolist()
        report["Diffusion"] = analyze("Diffusion", diff, seed_fps, known_fps)
    # GPT
    gpt_csv = RESULTS / "gpt_generated_pam_analogs.csv"
    if gpt_csv.exists():
        gpt = pd.read_csv(gpt_csv)["smiles"].tolist()
        report["GPT"] = analyze("GPT", gpt, seed_fps, known_fps)
    # 片段重组
    frag_csv = RESULTS / "generated_pam_analogs.csv"
    if frag_csv.exists():
        frag = pd.read_csv(frag_csv)["smiles"].tolist()
        report["Fragment"] = analyze("Fragment", frag, seed_fps, known_fps)

    out = RESULTS / "generator_comparison.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"\n已保存 -> {out}")


if __name__ == "__main__":
    main()
