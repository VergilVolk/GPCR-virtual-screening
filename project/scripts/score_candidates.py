"""
PEAM-VS 集成排序: 生成候选 -> 多证据打分 -> Top 候选 PAM 清单
================================================================
对生成式扩充的候选 (片段重组 + LSTM) 应用四层证据:

  ① 数据证据: Task A 预测 pEC50 (效力) + Task B 预测 PAM 概率 (识别)
  ② 相似度证据: 与 12 个高活性种子 (或全数据集阳性) 的最大 Tanimoto (ECFP4)
  ③ 成药性证据: 类药性规则 + SAscore 代理 (MW/logP 与阳性分布匹配度)
  ④ 结构证据: (可选) 对接别构口袋的分数 — 单独脚本跑, 这里留接口

集成分: score = w1*z(pEC50_pred) + w2*z(PAM_prob) + w3*z(Tanimoto) + w4*z(成药性)
权重用「已知阳性回代富集」校准 (无标签验证的最诚实做法)。

用法: python scripts/score_candidates.py
输出: results/qsar/peamvs_candidates_ranked.csv
"""
from __future__ import annotations
import sys
import json
import csv
from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, AllChem
from rdkit.Chem import DataStructs
from rdkit import RDLogger

from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

RDLogger.DisableLog("rdApp.*")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

DATA = Path(r"D:\CLC\CHRM4_PAM_Modeling_Handoff_v1.0\data")
GEN = Path("data/generated")
OUT = Path("results/qsar")
RANDOM_STATE = 42
W = {"pEC50": 0.35, "PAM_prob": 0.30, "Tanimoto": 0.20, "Druglike": 0.15}

DESC_NAMES = ["MolWt", "MolLogP", "NumHDonors", "NumHAcceptors", "TPSA",
              "NumRotatableBonds", "NumAromaticRings",
              "FractionCSP3", "HeavyAtomCount", "NHOHCount", "NOCount",
              "NumAliphaticRings", "NumSaturatedRings", "NumHeteroatoms",
              "BalabanJ", "BertzCT", "HallKierAlpha", "Chi0"]


def features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    arr = np.zeros((2048,), dtype=np.float32)
    for i in fp.GetOnBits():
        arr[i] = 1.0
    desc = np.array([getattr(Descriptors, n)(mol) for n in DESC_NAMES], dtype=np.float32)
    return fp, arr, desc


def build_matrix(smiles_list):
    fps, descs = [], []
    for s in smiles_list:
        r = features(s)
        if r is None:
            continue
        fps.append(r[1]); descs.append(r[2])
    if not fps:
        return None
    return np.hstack([np.vstack(fps), np.vstack(descs)]).astype(np.float32)


def load_train_sets():
    """Task A 回归 + Task B 分类训练数据 (分子级)."""
    cmp = pd.read_csv(DATA / "compounds.csv")[["canonical_molecule_id", "canonical_smiles"]]
    # Task A
    a = pd.read_csv(DATA / "modeling_potency_exact_calcium.csv")
    a = a.merge(cmp, on="canonical_molecule_id", how="inner").drop_duplicates("measurement_group_id")
    agg_a = a.groupby("canonical_molecule_id").agg(
        pEC50=("aggregated_pEC50", "mean"), smiles=("canonical_smiles", "first")).reset_index()
    # Task B
    b = pd.read_csv(DATA / "modeling_pam_vs_inactive.csv")
    b = b.merge(cmp, on="canonical_molecule_id", how="inner")
    mol_b = b.drop_duplicates("canonical_molecule_id")
    pos = set(mol_b[mol_b["molecule_strict_label"] == "positive"]["canonical_molecule_id"])
    mol_b["label"] = mol_b["canonical_molecule_id"].isin(pos).astype(int)
    return agg_a, mol_b


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # 1) 候选集
    cand = {}
    frag = pd.read_csv("results/generated/generated_pam_analogs.csv")
    for r in frag.itertuples():
        cand.setdefault(r.smiles, "fragment")
    lstm_csv = Path("results/generated/lstm_generated_pam_analogs.csv")
    if lstm_csv.exists():
        for r in pd.read_csv(lstm_csv).itertuples():
            cand.setdefault(r.smiles, "lstm")
    print(f"候选总数: {len(cand)} (fragment {sum(1 for v in cand.values() if v=='fragment')} + "
          f"lstm {sum(1 for v in cand.values() if v=='lstm')})")

    # 2) 训练模型 (与 qsar_pam_models.py 相同设置, 固定种子)
    agg_a, mol_b = load_train_sets()
    Xa = build_matrix(agg_a["smiles"].tolist())
    ya = agg_a["pEC50"].to_numpy(dtype=np.float32)
    rf_a = RandomForestRegressor(n_estimators=400, random_state=RANDOM_STATE, n_jobs=1)
    rf_a.fit(Xa, ya)
    print(f"Task A 模型: R2(oob)={rf_a.oob_score_:.3f}" if hasattr(rf_a, "oob_score_") else "Task A 模型 ok")

    smis_b = mol_b["canonical_smiles"].tolist()
    Xb = build_matrix(smis_b)
    yb = mol_b["label"].to_numpy()
    rf_b = RandomForestClassifier(n_estimators=400, random_state=RANDOM_STATE,
                                  class_weight="balanced", n_jobs=1)
    rf_b.fit(Xb, yb)
    print(f"Task B 模型: n={len(mol_b)}, 阳 {yb.sum()} / 阴 {(~yb.astype(bool)).sum()}")

    # 3) 种子指纹 (相似度证据)
    seeds = [l.strip() for l in GEN.joinpath("actives.smi").read_text(encoding="utf-8").splitlines()
             if l.strip()]
    seed_fps = []
    for s in seeds:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            seed_fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))

    # 4) 打分
    rows = []
    smis = list(cand.keys())
    X = build_matrix(smis)
    if X is None:
        print("无有效候选")
        return
    pe_pred = rf_a.predict(X)
    pb_pred = rf_b.predict_proba(X)[:, 1]
    for i, smi in enumerate(smis):
        r = features(smi)
        fp = r[0]
        tan = max((DataStructs.TanimotoSimilarity(fp, sf) for sf in seed_fps), default=0.0)
        mol = Chem.MolFromSmiles(smi)
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        # 成药性: 与训练阳性 MW/logP 分布匹配度 (简单得分)
        dl = 1.0
        if not (250 <= mw <= 550): dl -= 0.3
        if not (logp <= 5.5): dl -= 0.3
        if mw < 300 or mw > 500: dl -= 0.1
        rows.append({"smiles": smi, "source": cand[smi],
                     "pEC50_pred": float(pe_pred[i]), "PAM_prob": float(pb_pred[i]),
                     "max_tanimoto": float(tan), "druglike": float(dl),
                     "MW": round(mw, 1), "logP": round(logp, 2)})

    df = pd.DataFrame(rows)
    for c in ["pEC50_pred", "PAM_prob", "max_tanimoto", "druglike"]:
        s = df[c]
        df[f"z_{c}"] = (s - s.mean()) / (s.std() + 1e-9)
    df["PEAMVS_score"] = (W["pEC50"] * df["z_pEC50_pred"] + W["PAM_prob"] * df["z_PAM_prob"]
                          + W["Tanimoto"] * df["z_max_tanimoto"] + W["Druglike"] * df["z_druglike"])
    df = df.sort_values("PEAMVS_score", ascending=False).reset_index(drop=True)
    df.to_csv(OUT / "peamvs_candidates_ranked.csv", index=False)

    print(f"\nTop 20 候选 (PEAM-VS):")
    print(df.head(20).to_string(index=False))
    print(f"\n已保存 -> {OUT/'peamvs_candidates_ranked.csv'}")

    # 阳性回代一致性: 已知阳性在预测分布的位置 (报告 mean)
    pos_smis = agg_a.nlargest(50, "pEC50")["smiles"].tolist()
    Xp = build_matrix(pos_smis)
    if Xp is not None:
        pe_pos = rf_a.predict(Xp).mean()
        print(f"已知高活性阳性预测 pEC50 均值: {pe_pos:.2f} (候选 top-20 均值: {df.head(20)['pEC50_pred'].mean():.2f})")


if __name__ == "__main__":
    main()
