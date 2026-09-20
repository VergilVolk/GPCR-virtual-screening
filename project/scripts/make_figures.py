"""
论文核心图 (matplotlib) — PEAM-VS 结果可视化
图1: QSAR 指标双划分对比 (Task A/B/C)
图2: 打分失效: Vina 分数 PAM vs 非活性分布
图3: Top 候选 2D 结构网格
用法: python scripts/make_figures.py
输出: results/figures/
"""
from __future__ import annotations
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

FIG = Path("results/figures")
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False


def fig1_qsar_metrics():
    """QSAR 指标: Task B ROC-AUC / PR-AUC 双划分对比."""
    m = json.loads(Path("results/qsar/metrics.json").read_text(encoding="utf-8"))
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    # Task A: Spearman + R2
    ta = m["TaskA_regression"]
    tasks = [("Task A\npEC50 回归", "Spearman_rho", "R2"),
             ("Task B\nPAM vs 非活性", "ROC_AUC", "PR_AUC"),
             ("Task C\n10µM 阈值", "ROC_AUC", "PR_AUC")]
    for ax, (title, m1, m2) in zip(axes, tasks):
        key = "TaskA_regression" if "Task A" in title else ("TaskB_PAM_vs_inactive" if "Task B" in title else "TaskC_threshold_10uM")
        data = m[key]
        mol = data["molecule_level|RandomForest"]["mean"]
        sca = data["scaffold_level|RandomForest"]["mean"]
        labels = [m1.replace("_", " "), m2.replace("_", " ")]
        x = np.arange(2)
        w = 0.35
        ax.bar(x - w / 2, [mol[m1], mol[m2]], w, label="分子级划分", color="#4C72B0")
        ax.bar(x + w / 2, [sca[m1], sca[m2]], w, label="骨架划分", color="#DD8452")
        ax.set_title(title, fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_ylim(0, 1.15)
        for xi, (a, b) in enumerate(zip([mol[m1], mol[m2]], [sca[m1], sca[m2]])):
            ax.text(xi - w / 2, a + 0.02, f"{a:.2f}", ha="center", fontsize=9)
            ax.text(xi + w / 2, b + 0.02, f"{b:.2f}", ha="center", fontsize=9)
    axes[0].legend(fontsize=9, loc="lower right")
    fig.suptitle("QSAR: 分子级 vs 骨架划分 (5-fold CV, RandomForest)", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_qsar_metrics.png", dpi=200)
    print("fig1 saved")


def fig2_scoring_failure():
    """打分失效: Vina 分数分布 (PAM vs 非活性) + 排名."""
    df = pd.read_csv("results/structure/benchmark_pam_vs_inactive_raw.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    # 左: 分布
    pam = df[df["label"] == 1]["vina_affinity"]
    neg = df[df["label"] == 0]["vina_affinity"]
    axes[0].hist(pam, bins=8, alpha=0.6, color="#4C72B0", label=f"PAM (n={len(pam)})")
    axes[0].hist(neg, bins=8, alpha=0.6, color="#C44E52", label=f"实验非活性 (n={len(neg)})")
    axes[0].set_xlabel("Vina 亲和力 (kcal/mol, 越负越强)")
    axes[0].set_ylabel("分子数")
    axes[0].legend()
    axes[0].set_title("别构口袋对接打分分布 (7TRQ VU0467154 位点)")
    # 右: 排名条
    df = df.sort_values("vina_affinity").reset_index(drop=True)
    colors = ["#4C72B0" if r.label == 1 else "#C44E52" for r in df.itertuples()]
    axes[1].barh(range(len(df)), df["vina_affinity"], color=colors)
    axes[1].set_yticks(range(len(df)))
    axes[1].set_yticklabels([r.name[:22] for r in df.itertuples()], fontsize=7)
    axes[1].set_xlabel("Vina 亲和力")
    axes[1].set_title("按 Vina 分数排序 (蓝=PAM, 红=非活性)")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_scoring_failure.png", dpi=200)
    print("fig2 saved")


def fig3_top_candidates():
    """Top 候选 2D 结构网格."""
    df = pd.read_csv("results/qsar/peamvs_candidates_ranked.csv")
    top = df.head(10)
    mols = [Chem.MolFromSmiles(s) for s in top["smiles"]]
    labels = [f"#{i+1} pEC50={r.pEC50_pred:.2f} PAM={r.PAM_prob:.2f}"
              for i, r in enumerate(top.itertuples())]
    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(350, 300),
                               legends=labels)
    img.save(str(FIG / "fig3_top10_candidates.png"))
    print("fig3 saved")


if __name__ == "__main__":
    fig1_qsar_metrics()
    fig2_scoring_failure()
    fig3_top_candidates()
