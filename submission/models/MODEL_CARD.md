# Model Card — M4 虚拟筛选模型 v01

## 模型概述
- **名称**：DrugCLIP(2023 基座) + muscarinic-triplet 双侧投影头
- **文件**：`models/external_result.seed2026092{4,5,6}.projection.pt`（3 随机种子集成，每个 2.5MB）
- **类型**：开源预训练分子-口袋检索模型（DrugCLIP，基座为 Uni-Mol 分子/口袋编码器）之上的轻量投影微调层（2 层线性 + ReLU + L2 归一化），非完整 backbone 微调
- **调用方式**：`python screen.py --demo` / `--full`

## 训练
- **微调数据**：ChEMBL 毒蕈碱族（M1-M5）triplet 对（内部靶留出协议）；**未使用** M4R GaMD 筛选评测对（该测试集用于 head-to-head 评估时模型从未见过）
- **训练配置**：三随机种子 20260924/25/26；关键超参与种子记录见主仓库 `project/results/muscarinic_triplet_expansion_v01/`
- **训练入口**：主仓库 `project/scripts/`（`train.py` 为一键包装）

## 评估（同测试集对决，2026-10-01）
| 配置 | ROC | BEDROC20 | PR-AUC |
|---|---:|---:|---:|
| **本模型 + mean 池化（采用）** | **0.7136** | **0.1802** | **0.1547** |
| 本模型 + max 池化 | 0.6673 | 0.1263 | 0.1295 |
| 对照：Science-2026 基座 famaug LOSO | 0.5095 | 0.1118 | 0.0926 |

测试集：GaMD M4R 筛选对 25,249 行（2,301 活性）。

## 输入输出
- **输入**：分子 512 维 DrugCLIP 表征 + 10 个 GaMD M4 口袋表征（npz）
- **输出**：每分子余弦结合分（10 口袋 mean 池化，3 种子平均），降序候选清单

## 适用范围与已知局限
- 适用于 CHRM4(M4R) 结合候选的大规模初筛排序（顶层切片建议 ≤5%）
- **不预测变构正性调节（PAM）功能**：候选需经四上下文 MD（PACER-DC）功能复核
- EF1%（top-1% 切片）在同测试集上低于随机基线（0.87），早期识别请以 BEDROC/PR 为参考并保留较大候选池
- 表征依赖 DrugCLIP/Uni-Mol 编码环境（完整新分子编码流程见主仓库）

## 第三方组件
| 组件 | 版本/来源 | 用途 | 许可 |
|---|---|---|---|
| DrugCLIP | 2023 官方 ckpt（主仓库 data/external/） | 分子-口袋检索基座 | Apache-2.0 |
| Uni-Mol / Uni-Core | 主仓库 tools/Uni-Core | 表征编码器 | BSD-3 |
| GaMD M4 构象 | figshare 33283491 + 内部聚类 | 口袋集合 | CC-BY(公开部分) |
| ChEMBL | 33 | 微调标签 | CC-BY-4.0 |
