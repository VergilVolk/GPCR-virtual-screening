"""
Task A/B/C QSAR — CHRM4 PAM 效力/识别模型 (纯 CPU)
====================================================
解决核心问题「怎么判断一个分子是不是(多强的) M4 PAM」的三个可计算子问题:

  Task A 回归   : 结构 -> human M4 PAM pEC50 (效力数值)
  Task B 分类   : PAM vs 实验确认非活性 (严格识别)
  Task C 分类   : 10 uM 阈值内有效力 vs 未达阈值 (预筛)

关键设计:
  - 最小隔离单元 = canonical_molecule_id (禁止行级随机划分, 防同一分子跨区泄漏)
  - 两种划分策略对比: (1) 分子级随机 5-fold CV (2) Bemis-Murcko 骨架划分 (泛化压力)
  - 特征: ECFP4 (Morgan 2048) + 19 个理化描述符
  - 模型: RandomForest + GradientBoosting (基线), 输出特征重要性(可解释性)

用法: python scripts/qsar_pam_models.py
输出: results/qsar/{metrics, predictions, feature_importance}.json/csv
"""
from __future__ import annotations
import sys
import json
import csv
from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger

from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.metrics import (mean_absolute_error, mean_squared_error, r2_score,
                             roc_auc_score, average_precision_score,
                             matthews_corrcoef, balanced_accuracy_score)

RDLogger.DisableLog("rdApp.*")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

DATA = Path(r"D:\CLC\CHRM4_PAM_Modeling_Handoff_v1.0\data")
OUT = Path("results/qsar")
RANDOM_STATE = 42

# ---------------- feature builders ----------------

DESC_NAMES = ["MolWt", "MolLogP", "NumHDonors", "NumHAcceptors", "TPSA",
              "NumRotatableBonds", "NumAromaticRings",
              "FractionCSP3", "HeavyAtomCount", "NHOHCount", "NOCount",
              "NumAliphaticRings", "NumSaturatedRings", "NumHeteroatoms",
              "BalabanJ", "BertzCT", "HallKierAlpha", "Chi0"]


def mol_features(smiles: str) -> tuple[np.ndarray, np.ndarray]:
    """ECFP4 (2048) + 19 physchem descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    fp_arr = np.zeros((2048,), dtype=np.float32)
    for i in fp.GetOnBits():
        fp_arr[i] = 1.0
    desc = np.array([getattr(Descriptors, n)(mol) for n in DESC_NAMES], dtype=np.float32)
    return fp_arr, desc


def build_matrix(smiles_list):
    fps, descs = [], []
    for s in smiles_list:
        f, d = mol_features(s)
        fps.append(f)
        descs.append(d)
    X = np.hstack([np.vstack(fps), np.vstack(descs)])
    return X.astype(np.float32)


# ---------------- splits ----------------

def scaffold_groups(smiles_list):
    groups = []
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            groups.append(None)
            continue
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        groups.append(scaf)
    return groups


def group_kfold(groups, n_splits=5, random_state=RANDOM_STATE):
    """按骨架/分子组做 KFold: 同一 group 只出现在同一折."""
    unique = sorted({g for g in groups if g is not None})
    rng = np.random.RandomState(random_state)
    rng.shuffle(unique)
    folds = np.array_split(unique, n_splits)
    fold_of_group = {}
    for i, fs in enumerate(folds):
        for g in fs:
            fold_of_group[g] = i
    # 未分到组的(None)随机分配
    none_idxs = [i for i, g in enumerate(groups) if g is None]
    rng.shuffle(none_idxs)
    for k, i in enumerate(none_idxs):
        fold_of_group[(None, i)] = k % n_splits
    return np.array([fold_of_group.get((g, i) if g is None else g, 0)
                     for i, g in enumerate(groups)])


# ---------------- evaluation ----------------

def eval_reg(y_true, y_pred):
    return {
        "MAE": float(mean_absolute_error(y_true, y_pred)),
        "RMSE": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "R2": float(r2_score(y_true, y_pred)),
        "Spearman_rho": float(pd.Series(y_true).corr(pd.Series(y_pred), method="spearman")),
    }


def eval_clf(y_true, y_prob, y_pred):
    return {
        "ROC_AUC": float(roc_auc_score(y_true, y_prob)),
        "PR_AUC": float(average_precision_score(y_true, y_prob)),
        "MCC": float(matthews_corrcoef(y_true, y_pred)),
        "BalancedAcc": float(balanced_accuracy_score(y_true, y_pred)),
    }


# ---------------- data loading ----------------

def load_task_a():
    """Task A: 分子级 pEC50 (钙流 FLIPR, 探针=乙酰胆碱)."""
    df = pd.read_csv(DATA / "modeling_potency_exact_calcium.csv")
    cmp = pd.read_csv(DATA / "compounds.csv")[["canonical_molecule_id", "canonical_smiles", "identity_status"]]
    df = df.merge(cmp, on="canonical_molecule_id", how="inner")
    df = df[df["identity_status"] == "confirmed_structure"]
    # measurement_group_id 去重(表示级重复) -> 分子级聚合
    df = df.drop_duplicates("measurement_group_id")
    agg = df.groupby("canonical_molecule_id").agg(
        pEC50=("aggregated_pEC50", "mean"),
        smiles=("canonical_smiles", "first"),
        n_measurements=("aggregated_pEC50", "count"),
    ).reset_index()
    return agg


def load_task_bc(fname, label_col, conflict_col):
    df = pd.read_csv(DATA / fname)
    cmp = pd.read_csv(DATA / "compounds.csv")[["canonical_molecule_id", "canonical_smiles", "identity_status"]]
    df = df.merge(cmp, on="canonical_molecule_id", how="inner")
    df = df[df["identity_status"] == "confirmed_structure"]
    # 分子级标签
    mol = df.drop_duplicates("canonical_molecule_id")
    # 排除标签冲突的分子
    if conflict_col in mol.columns:
        mol = mol[mol[conflict_col] != True]  # noqa: E712
    pos = mol[mol[label_col] == "positive"]["canonical_molecule_id"]
    neg = mol[mol[label_col] != "positive"]["canonical_molecule_id"]
    return mol, pos, neg


# ---------------- main ----------------

def run_task(name, X, y, groups, clf=False):
    """5-fold 分子级 CV + 5-fold 骨架 CV, 两模型."""
    models_reg = {
        "RandomForest": RandomForestRegressor(n_estimators=400, random_state=RANDOM_STATE, n_jobs=1),
        "GradientBoosting": GradientBoostingRegressor(n_estimators=300, random_state=RANDOM_STATE),
    }
    models_clf = {
        "RandomForest": RandomForestClassifier(n_estimators=400, random_state=RANDOM_STATE,
                                               class_weight="balanced", n_jobs=1),
        "GradientBoosting": GradientBoostingClassifier(n_estimators=300, random_state=RANDOM_STATE),
    }
    models = models_clf if clf else models_reg
    results = {}
    # 两种划分: (1) 分子级随机 5-fold (2) 骨架级 5-fold (泛化压力)
    rng = np.random.RandomState(RANDOM_STATE)
    mol_folds = np.zeros(len(y), dtype=int)
    perm = rng.permutation(len(y))
    mol_folds[perm] = np.arange(len(y)) % 5
    scaf_folds = group_kfold(groups)
    for split_name, fold_ids in [("molecule_level", mol_folds), ("scaffold_level", scaf_folds)]:
        for mname, model in models.items():
            # 组级 KFold 折叠 (fold_ids 同组同折, 但 KFold 需 fold 编号)
            metrics, oof_true, oof_prob = [], [], []
            for fold in range(5):
                tr = fold_ids != fold
                te = fold_ids == fold
                if clf:
                    if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
                        continue
                    model.fit(X[tr], y[tr])
                    prob = model.predict_proba(X[te])[:, 1]
                    pred = (prob >= 0.5).astype(int)
                    oof_true.extend(y[te].tolist())
                    oof_prob.extend(prob.tolist())
                    metrics.append(eval_clf(y[te], prob, pred))
                else:
                    model.fit(X[tr], y[tr])
                    pred = model.predict(X[te])
                    oof_true.extend(y[te].tolist())
                    oof_prob.extend(pred.tolist())
                    metrics.append(eval_reg(y[te], pred))
            if not metrics:
                continue
            mean_m = {k: float(np.mean([m[k] for m in metrics])) for k in metrics[0]}
            std_m = {k: float(np.std([m[k] for m in metrics])) for k in metrics[0]}
            results[f"{split_name}|{mname}"] = {"mean": mean_m, "std": std_m}

    # 特征重要性 (RF, 分子级划分全量训练)
    rf = models["RandomForest"]
    rf.fit(X, y)
    imp = np.asarray(rf.feature_importances_)
    return results, imp


def save_predictions(task, X, y, groups, clf=False):
    """骨架划分 holdout: 训练集全量训练, 测试集(20% 骨架)预测, 存 CSV 供论文画图."""
    fold_ids = group_kfold(groups)
    te = fold_ids == 4
    tr = ~te
    m = (RandomForestClassifier(n_estimators=400, random_state=RANDOM_STATE,
                                class_weight="balanced", n_jobs=1)
         if clf else RandomForestRegressor(n_estimators=400, random_state=RANDOM_STATE, n_jobs=1))
    m.fit(X[tr], y[tr])
    prob = m.predict_proba(X[te])[:, 1] if clf else m.predict(X[te])
    return y[te], prob


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {}

    # ---- Task A ----
    print("=" * 60)
    print("Task A: PAM pEC50 回归 (calcium, 分子级)")
    a = load_task_a()
    smis = a["smiles"].tolist()
    Xa = build_matrix(smis)
    ya = a["pEC50"].to_numpy(dtype=np.float32)
    print(f"  分子数 {len(a)} | X {Xa.shape}")
    res_a, imp_a = run_task("A", Xa, ya, scaffold_groups(smis))
    report["TaskA_regression"] = res_a
    yt, yp = save_predictions("A", Xa, ya, scaffold_groups(smis))
    pd.DataFrame({"true_pEC50": yt, "pred_pEC50": yp}).to_csv(OUT / "taskA_predictions.csv", index=False)
    imp_df = pd.DataFrame({"feature": [f"FP{i}" for i in range(2048)] + DESC_NAMES,
                           "importance": imp_a}).sort_values("importance", ascending=False)
    imp_df.head(25).to_csv(OUT / "taskA_feature_importance_top25.csv", index=False)
    for k, v in res_a.items():
        print(f"  {k}: mean {v['mean']}")

    # ---- Task B ----
    print("=" * 60)
    print("Task B: PAM vs 实验非活性 (严格二分类)")
    mol_b, pos_b, neg_b = load_task_bc("modeling_pam_vs_inactive.csv", "molecule_strict_label", "molecule_label_conflict")
    pos_b = set(pos_b)
    mask = mol_b["canonical_molecule_id"].isin(pos_b)
    yb = mask.astype(int).to_numpy()
    smis_b = mol_b["canonical_smiles"].tolist()
    Xb = build_matrix(smis_b)
    print(f"  分子 {len(mol_b)} (阳 {yb.sum()} / 阴 {(~mask).sum()}) | X {Xb.shape}")
    res_b, imp_b = run_task("B", Xb, yb, scaffold_groups(smis_b), clf=True)
    report["TaskB_PAM_vs_inactive"] = res_b
    yt_b, yp_b = save_predictions("B", Xb, yb, scaffold_groups(smis_b), clf=True)
    pd.DataFrame({"true": yt_b, "prob": yp_b}).to_csv(OUT / "taskB_predictions.csv", index=False)
    for k, v in res_b.items():
        print(f"  {k}: mean {v['mean']}")

    # ---- Task C ----
    print("=" * 60)
    print("Task C: 10 uM 阈值分类")
    mol_c, pos_c, neg_c = load_task_bc("modeling_threshold_10uM.csv", "class_label", None)
    pos_c = set(mol_c[mol_c["class_label"] == "threshold_positive"]["canonical_molecule_id"])
    mask_c = mol_c["canonical_molecule_id"].isin(pos_c)
    yc = mask_c.astype(int).to_numpy()
    smis_c = mol_c["canonical_smiles"].tolist()
    Xc = build_matrix(smis_c)
    print(f"  分子 {len(mol_c)} (阳 {yc.sum()} / 阴 {(~mask_c).sum()}) | X {Xc.shape}")
    res_c, imp_c = run_task("C", Xc, yc, scaffold_groups(smis_c), clf=True)
    report["TaskC_threshold_10uM"] = res_c
    yt_c, yp_c = save_predictions("C", Xc, yc, scaffold_groups(smis_c), clf=True)
    pd.DataFrame({"true": yt_c, "prob": yp_c}).to_csv(OUT / "taskC_predictions.csv", index=False)
    for k, v in res_c.items():
        print(f"  {k}: mean {v['mean']}")

    with open(OUT / "metrics.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=2)
    print(f"\n已保存 -> {OUT}/metrics.json")


if __name__ == "__main__":
    main()
