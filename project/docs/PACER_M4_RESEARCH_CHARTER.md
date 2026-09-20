# PACER-M4 研究章程（冻结版 v1.0）

> 冻结日期：2026-08-28  
> 用途：约束科学问题、方法开发、验证口径与对外表述，防止项目重新退化为工具堆叠。  
> 变更规则：任何偏离本章程的主线变更，必须说明新增假设、验证数据和预期消融实验。

## 1. 项目目标

开发 **PACER-M4（Probe-Aware Allosteric Coupling and Efficacy Ranker）**：一个面向人源 CHRM4/M4 正性别构调节剂（PAM）的正构探针感知、受体状态感知、变构耦合感知的功能效力预测与候选生成框架。

项目追求的是 **M4 PAM 功能预测的 target-specific SOTA**，不是宣称建立适用于所有 GPCR 的通用模型。SOTA 只有在冻结 benchmark、公开 baseline、严格外推划分、置信区间和统计检验均完成后才能使用。

## 2. 唯一核心科学问题

> 能否显式建模“候选 PAM—正构探针—实验上下文—受体变构耦合”，从分子结构预测 M4 PAM 的功能效力，并在新骨架、新药化系列和时间外推测试上，稳定超过 ligand-only QSAR、传统 docking、相互作用指纹和通用分子深度学习模型？

### 2.1 PAM 的计算定义

PAM 不是“能进入别构口袋”或“具有较高结合亲和力”的同义词。PAM 是条件性功能表型，至少与以下因素有关：

- PAM 对受体的结合；
- 与正构激动剂的结合协同性；
- 对正构激动剂功能效力的协同性；
- PAM 自身的内在激动活性；
- 正构探针、探针浓度、受体状态、细胞系统、信号通路和 readout。

因此，任何结构分数只能作为结构证据，不能单独作为 PAM ground truth。

## 3. PACER-M4 方法架构

### 3.1 分子表示层

- 传统表示：ECFP/Morgan + 理化描述符；
- 图表示：D-MPNN/GIN 等低容量图模型；
- 预训练表示：MolFormer/ChemBERTa 类编码器冻结或参数高效微调；
- 430 分子规模下禁止从头训练高容量深度模型后直接声称提升。

### 3.2 Probe/assay-aware 条件层

模型输入应尽可能保留：orthosteric agonist、浓度或 ECxx、assay family、readout、cell system、G 蛋白/嵌合系统、species、endpoint type。主任务使用 assay-harmonized calcium/ACh 数据；异质 assay 作为辅助任务或条件任务，不与主标签无条件混合。

目标形式：

`y_hat = f(PAM molecule, orthosteric probe, assay context)`

**数据可识别性限制（2026-08-28 审计）**：当前扩展精确效力视图中，已标注正构探针只有 acetylcholine 一种；`agonist_conc_or_ecxx` 覆盖约 24.8%，且已知值只有 EC20。因此现阶段可以验证 assay-aware/source-invariant learning，但不能从现有数据验证跨探针 probe dependence。模型保留 probe 输入接口；只有补充至少两种探针且具有足够重叠化学空间的数据后，才允许声称 probe-aware 增益。

### 3.3 变构耦合结构图

以 7TRQ 等 M4 复合物和受体构象集合建立联合图：

`PAM 原子 -> 别构口袋残基 -> 残基传播网络 -> 正构口袋残基 -> 正构配体`

候选特征包括：

- PAM—残基相互作用；
- 残基间空间邻接与理化关系；
- 到正构口袋/候选枢纽残基的路径；
- 多受体构象中的接触稳定性；
- 正构配体存在条件下的联合相互作用图；
- 结构特征的不确定性与 pose 稳定性。

TYR439 等残基是待检验的机制假设，不是预先认定的真理。所有 hub/subpocket 结论必须在独立集上验证。

### 3.4 潜在亚位点混合专家

不手工用少数阳性定义亚位点后回测。先对已知阳性和严格阴性的多构象 pose/IFP 无监督聚类，再建立 latent-subpocket mixture-of-experts：

`y_hat = sum_k gate(k | molecule, pose) * expert_k(molecule, IFP_k)`

亚位点数量、稳定性和增益必须通过嵌套交叉验证与消融实验确定。

### 3.5 SAR 差分与 activity-cliff 学习

从 matched molecular pairs 学习局部结构变化导致的 `delta pEC50`。联合目标至少包含：

- 绝对 pEC50 回归；
- PAM vs 严格非活性；
- 10 uM 阈值任务；
- pairwise potency ranking；
- activity-cliff 判别；
- 可选的 coupling-consistency 正则。

首轮冻结基准显示，差分 SAR 在 scaffold-held-out 中改善部分 cliff 指标，但在 source-held-out 中增益有限或退化。因此 pairwise loss 是待消融的辅助目标，不预设为必然有效；只有同时改善来源外推与校准时才保留在最终模型中。

PACER-M4 v0.1 的三阴性系列留出实验进一步表明：分类 + pEC50 + delta-SAR 联合训练相对同架构分类器获得统计稳定的 ROC-AUC 增益，但尚未在 PR-AUC/Brier 上超过 Extra Trees。该结果将 pairwise 功能效力学习晋级为候选核心模块，但最终保留仍需 nested loss-weight selection、更多时间/系列外推和早期富集验证。

### 3.6 不确定性与适用域

最终输出必须包括预测值、区间、校准概率、适用域和模型分歧。候选生成奖励中必须惩罚高不确定性和超出适用域的分子，防止生成器 reward hacking。

## 4. 双线生成

### 路线 A：SAR/MMP 规则驱动生成

从已知药化系列学习有利局部变换，进行受控片段替换与重组。重点是可解释、可合成、低风险的邻域优化。

### 路线 B：条件生成模型

在大规模类药语料预训练，在 M4 PAM 空间微调，以 PACER-M4 的效力、耦合、选择性、成药性、不确定性和合成可及性形成多目标条件。模型类型可以是离散扩散、生成式 Transformer 或其他经 baseline 证明合适的方法；“使用 diffusion”本身不构成创新。

两条路线最终使用 Pareto front 和化学簇多样性选择，不采用未经校准的固定加权总分。

## 5. Benchmark v1.0：SOTA 的唯一坐标系

### 5.1 必备测试

1. molecule-level split：仅作宽松参考；
2. Bemis-Murcko scaffold split：测试新骨架；
3. leave-one-source/series-out：整篇论文或药化系列留出；
4. temporal split：早期文献训练，后期文献盲测；
5. activity-cliff/matched-pair test：测试近邻效力反转；
6. strict inactive test：实验确认非活性，不用简单属性诱饵替代；
7. external structural benchmark：评价早期富集和 pose/结构证据，不冒充功能效力验证。

### 5.2 泄漏红线

- 同一 canonical molecule 不得跨集合；
- 同一表示重复/measurement group 不得跨集合；
- scaffold/series/source 评估必须按相应组整体隔离；
- 预训练语料若包含测试分子的标签或文本，必须披露并做去重；
- 超参数选择不得查看最终测试集；
- 生成种子、QSAR 训练和候选评价的关系必须明确，禁止循环验证；
- 不能使用测试集定义药效团、hub、亚位点后再在同一测试集报告提升。

### 5.3 主要指标

- 回归：MAE、RMSE、R2、Spearman、pairwise accuracy；
- 分类：ROC-AUC、PR-AUC、MCC、balanced accuracy、sensitivity/specificity；
- 筛选：EF(0.5%)、EF(1%)、BEDROC、Top-k hit rate；
- 校准：Brier score、ECE、coverage/interval width；
- 所有主指标报告 bootstrap 95% CI，多随机种子报告均值与标准差。
- 跨来源/系列结果必须报告每系列指标、宏平均和最差系列；不同系列标签比例或分数尺度不同，禁止把 pooled AUC 作为主要域泛化结论。

## 6. 必须击败的 baseline

### Ligand-only

- Tanimoto kNN；
- ECFP + RF；
- ECFP + XGBoost/LightGBM；
- descriptor + SVM；
- MMP/局部 SAR baseline。

### 深度分子模型

- D-MPNN/Chemprop；
- GIN/GCN/GAT；
- 预训练分子编码器 + 小型预测头。

### 结构模型

- Vina 单结构；
- ensemble BEmin/BEavg；
- IFP similarity；
- pharmacophore；
- 几何接触；
- ligand descriptor + docking 的简单 late fusion。

### 消融

- 去掉 probe/assay 条件；
- 去掉结构耦合图；
- 去掉 MMP/activity-cliff loss；
- 去掉 mixture-of-experts；
- 去掉不确定性约束；
- 单结构 vs ensemble。

## 7. 对“真正有效”的表述边界

纯计算不能确保新分子真正具有 PAM 效力。没有湿实验时，最高等级证据是：严格时间盲测、整系列留出、activity-cliff、严格实验阴性、多 assay 一致性、外部结构基准、校准不确定性与可复现代码。

候选只能称为：

> 经严格回顾性/伪前瞻验证后得到的高置信度、可实验检验的 M4 PAM 假设。

不得称为“已发现/已验证的有效 PAM”。

## 8. 当前开发顺序

1. 冻结 Benchmark v1.0 和数据审计；
2. 重做强 baseline 与统计协议；
3. MMP/activity-cliff 模块；
4. probe/assay-aware 多任务模型；
5. 结构耦合图与亚位点专家；
6. 不确定性和适用域；
7. 双线条件生成；
8. 统一盲测、消融与候选输出。

生成模型不得早于 benchmark 和 baseline 完成而成为主开发任务。
