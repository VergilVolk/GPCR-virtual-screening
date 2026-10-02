# DrugCLIP 代际与微调模型统一评测协议 v01

## 1. 目的

回答两个不同问题：

1. 在完全相同的分子、靶点、口袋构象和指标下，2023 与 2026 DrugCLIP 哪个更适合 GPCR 虚拟筛选？
2. 现有 GPCR 微调、13 靶点微调、family augmentation 和 20 靶点微调，是否真正优于原始模型及传统基线？

本评测只衡量**结合候选检索**，不把结果解释为 PAM 功能、效力或协同性预测。

## 2. 两层证据，禁止混用

### A. 开发/选择基准

- 4 个 GPCR：B2AR、CCR2、M2R、M4R；
- 28,618 个可评估靶点—分子对；
- 2023 与 2026 使用完全相同的 28,519 个分子和 40 个受体构象；
- 该数据已参与模型开发，只能用于模型诊断、消融和部署路线选择；
- 不能称为外部盲测，也不能据此宣称 SOTA。

### B. 一次性锁箱

- 分子及标签来自未用于训练、调参或选择模型的新来源；
- 标签由非建模队员保管；
- 建模侧只收到 `pair_id、target、smiles、pocket`；
- 在揭盲前冻结模型、聚合规则、预测文件和 SHA256；
- 只允许一次揭盲；失败结果同样保留。

## 3. 冻结模型清单

| 类型 | 模型 | 作用 |
|---|---|---|
| 原始模型 | DrugCLIP 2023 raw | 旧版通用结合表征 |
| 原始模型 | DrugCLIP Science 2026 raw | 新版 LIT-PCBA 表征 |
| 迁移模型 | 2023 GPCR LOTO | 三靶点训练、整靶点外推 |
| 迁移模型 | 2026 GPCR-only ep80 | 检验纯 GPCR 域适配 |
| 迁移模型 | 2026 13-target | 混合靶点域适配 |
| 迁移模型 | 2026 13-target family-aug | 加入家族上下文增强 |
| 迁移模型 | 2026 20-target | 更大范围多靶点训练 |
| 探索性集成 | 2023 GPCR LOTO + 2026 13T family-aug | 靶点内百分位秩固定 1:1 融合 |
| 冻结路由 | M4 使用 2023 GPCR LOTO，其余使用固定融合 | 延续既有 M4 退化审计结论，避免负迁移 |

探索性集成是在开发集上提出的，必须在锁箱前冻结，不能在锁箱上再选权重。

## 4. Baseline

- DrugCLIP 两代原始 checkpoint；
- ECFP4 logistic，整靶点留出并清除重叠 Murcko 骨架；
- Glide/Vina 单结构与 ensemble docking，在各自实际有结果的共同子集上配对比较；
- IID 随机分数；
- 随机标签训练控制。

注意：当前活性/decoy 集存在明显二维性质捷径，ECFP 的高分不能直接理解为真实前瞻筛选能力。

## 5. 指标与统计

主要指标：

- BEDROC α=20；
- EF1%；
- EF5%。

次要指标：

- PR-AUC；
- ROC-AUC。

所有总结果先逐靶点计算再做 macro average。置信区间按靶点分层，并以 Murcko scaffold 为重采样单位；禁止将每个分子行当作独立样本。

## 6. 决策规则

1. 虚拟筛选主排序优先看 BEDROC、EF1、EF5，不以 ROC-AUC 单独定胜负。
2. 若新模型只提高 ROC-AUC、却损害头部富集，不替换主筛模型。
3. 若两个模型互补，只允许使用锁箱前冻结的固定融合；不得按每个外部靶点事后挑最好模型。
4. M4 功能性 PAM 判别仍由 PACER-FKG 四上下文动态模块承担；DrugCLIP 只负责前级结合候选压缩。
5. 只有一次性锁箱对主要指标给出稳定增益，才升级为外部泛化结论。

## 7. 运行入口

统一入口：`project/scripts/run_drugclip_generation_bakeoff.py`。

输出：

- `macro_metrics.csv`；
- `per_target_metrics.csv`；
- `scaffold_bootstrap_ci.csv`；
- `matched_docking_comparison.csv`；
- `predictions.csv`；
- 三张标准图；
- 含全部输入 SHA256 的 `audit.json`。

## 8. 当前主张边界

该评测可以支持“不同 DrugCLIP 代际及适配策略在 GPCR 检索上存在可重复的能力差异和互补性”，但尚不能支持“新模型已达到 SOTA”“候选物是 PAM”或“候选物具有实验效力”。这些结论分别需要一次性外部锁箱和湿实验功能验证。
