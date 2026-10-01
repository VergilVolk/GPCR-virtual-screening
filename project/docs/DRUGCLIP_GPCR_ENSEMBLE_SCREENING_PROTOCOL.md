# GPCR 别构调节剂的动态集合 DrugCLIP 筛选协议

状态：开发集已冻结，DrugCLIP 计算进行中（2026-09-25）

## 1. 研究问题

本实验不预测 PAM 功能效力，而检验一个范围更清楚的问题：

> 将 GaMD 得到的多个 GPCR 别构口袋构象及其 cluster population 引入 DrugCLIP，能否比单构象 DrugCLIP、静态 docking 和传统 ensemble docking 更好地富集已知别构调节剂？

## 2. 数据

来源为公开的四靶点 GPCR allosteric-modulator benchmark：B2AR、CCR2、M2R、M4R。开发集保留全部无标签冲突的 active，并为每个靶点用预先固定的 SHA256 规则选择十倍 property-matched decoy。

| 靶点 | Active | Decoy |
|---|---:|---:|
| B2AR | 64 | 640 |
| CCR2 | 16 | 160 |
| M2R | 230 | 2,300 |
| M4R | 2,305 | 23,050 |

总计 28,765 个靶点—分子对、28,629 个唯一分子、19,448 个 Murcko scaffold。9 个同一靶点内 active/decoy 标签冲突分子已删除。

切分使用全局 Murcko scaffold：同一分子或骨架即使出现在多个靶点，也只能属于一个外层折。模型、ECFP 和所有消融共用同一五折切分。

## 3. 方法

### 3.1 单构象基线

- 官方 DrugCLIP cluster-0；
- 双侧 rank-4 projection LoRA；
- 目标和标签双平衡 BCE；
- 两个 Uni-Mol encoder 保持冻结。

### 3.2 动态集合方法

每个靶点使用公开 GaMD 聚类的 10 个代表构象及 population \(\pi_{tc}\)。DrugCLIP 给出分子 \(m\) 与构象 \(c\) 的余弦分数 \(s_{mtc}\)。集合分数定义为：

\[
S(m,t)=\tau_t\log\sum_c\pi_{tc}\exp(s_{mtc}/\tau_t).
\]

该形式在低温时接近最优构象，在高温时接近 population-weighted average。每个 \(\tau_t\) 只能在当前外层训练折学习，不得查看测试标签。训练仍只更新双侧 projection LoRA 和四个温度参数。

### 3.3 必须比较的基线

1. 官方 DrugCLIP：cluster-0、population mean、max、固定温度 population-LSE；
2. 单构象 LoRA；
3. 动态集合 LoRA；
4. 随机标签动态集合 LoRA；
5. ECFP4 logistic；
6. 公开静态 Glide/Vina；
7. 公开 GaMD ensemble Glide/Vina 的 BEmin/BEavg。

## 4. 指标

每个靶点独立计算 ROC-AUC、PR-AUC、BEDROC20、EF1% 和 EF5%，主结果为四靶点宏平均。任何 docking 缺失值均置于该方法最差分数之后，禁止只评价成功 docking 的分子。

开发集已经复算的 docking 基线：

| 方法 | 宏 AUC | 宏 PR-AUC | BEDROC20 | EF1% | EF5% |
|---|---:|---:|---:|---:|---:|
| 静态 Glide | 0.663 | 0.198 | 0.246 | 3.98 | 2.16 |
| 静态 Vina | 0.645 | 0.187 | 0.238 | 3.17 | 2.01 |
| Ensemble Glide BEavg | **0.690** | 0.277 | 0.332 | **6.43** | 3.26 |
| Ensemble Glide BEmin | 0.677 | **0.299** | **0.373** | 6.31 | **4.05** |
| Ensemble Vina BEavg | 0.654 | 0.182 | 0.246 | 2.27 | 2.38 |
| Ensemble Vina BEmin | 0.657 | 0.196 | 0.264 | 2.65 | 2.80 |

## 5. 晋级标准

只有同时满足以下条件才进入完整约 12 万 decoy 的最终实验：

1. 动态集合 LoRA 在至少三个靶点方向一致；
2. 宏 BEDROC20 或宏 EF1% 高于单构象 LoRA；
3. 优于随机标签 LoRA，并有 scaffold-cluster bootstrap 支持；
4. 不以降低多个靶点性能换取 M4R 样本量主导的总体提升；
5. 对 ECFP 的优势或互补性必须单独证明，不能只比较 docking。

## 6. 主张边界

通过该实验最多可以主张“GPCR 别构 active/decoy 的回顾性早期富集改善”。它不能证明候选是 PAM，不能预测 EC50、协同性、内在激动或临床效应。PETA 等 2026 年参数高效适配工作必须作为相关方法讨论；本工作的区别必须落在 GaMD population-aware conformational aggregation，而不是笼统的 LoRA 微调。
