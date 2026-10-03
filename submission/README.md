# CHRM4 PAM 虚拟筛选 — 代码提交包（附件5 规范）

一键复现 M4 靶点虚拟筛选：DrugCLIP 检索 + 微调投影 + 标准化候选清单。

## 环境

- Python 3.11，依赖见 `requirements.txt`（`pip install -r requirements.txt`）
- CPU 即可运行 demo 与筛选（无需 GPU；全量新分子编码需 GPU，见下）
- 操作系统：Windows 10/11 或 Linux；无 CUDA 依赖（torch CPU 版）

## 快速开始

```bash
pip install -r requirements.txt

# 演示（自带 50 分子示例，约 5 秒）：
python screen.py --demo

# 全量 28,519 分子库（需主仓库在 ../project 下）：
python screen.py --full

# 复现训练（微调投影权重，需主仓库数据，约 3 分钟/种子 CPU）：
python train.py --full-data --seed 20260925

# 复现基准（13 靶 LOSO ×3 种子，约 40 分钟 CPU）：
python train.py --loso --seed 20260925
```

## 输入输出

- 输入：分子与口袋的 DrugCLIP 512 维表征（npz；demo 自带，全量在主仓库）
- 输出：`results/results.csv` — 标准候选清单（候选编号/赛道/SMILES/结合
  分数/模型版本/备注），UTF-8，按分数降序

## 目录

```
submission/
├── README.md            # 本文件
├── requirements.txt     # 依赖
├── screen.py            # 一键筛选（主运行入口）
├── train.py             # 一键训练/基准复现入口
├── src/                 # 投影头与筛选核心（含注释）
├── models/              # 3 种子最终权重 + MODEL_CARD.md
├── data/                # demo 数据 + 全量数据来源/许可/防泄漏说明
├── notebooks/           # demo_screening.ipynb 可执行演示
├── results/             # results.csv 候选清单
└── logs/                # 训练日志/复现权重输出位置
```

## 结果说明

排序逻辑：10 个 GaMD M4 口袋构象余弦分 mean 池化 → 3 随机种子等权平均。
不确定性：种子间分数标准差（notebook 中示例）。候选为**结合检索候选**，
PAM 功能结论须走四上下文 MD 复核（备注列已注明）。

## 本项目实际贡献（相对开源基座）

1. 毒蕈碱族 triplet 双侧投影微调协议（3 种子，ChEMBL 标签）；
2. 同测试集 head-to-head 证明 2023 基座 + mean 池化为 M4 最优配置
   （ROC 0.7136 vs Science-2026 基座 0.5095，详见 models/MODEL_CARD.md）；
3. 完整防泄漏评估链（LOSO + scaffold purge + frozen SHA）。

## 已知局限

top-1% 切片富集（EF1%）不稳定（0.87），建议保留 top-5% 候选池；
模型不预测 PAM 功能与效力。
