# PACER-M4 v1：功能效力预测、结构门控与分子生成闭环

> 冻结日期：2026-08-30
> 结论等级：完整、可复现的计算方法学工作；PACER-FS 是历史内部 benchmark 的冻结主方法，PACER-AssayConditional 是专利知情的下一代开发方法；**独立 Suven 外测未通过 SOTA 门槛，未宣称发现经实验验证的新 PAM**。

## 1. 最终回答

我们已经完成了从数据、强 baseline、跨系列与时间外推、activity cliff、三受体状态结构检验、双线生成、适用域审计，到 200 个候选结构门控和 24 个多样化候选假设的闭环。

最重要的结果不是一个漂亮但虚假的 AUC，而是明确定位了 M4 PAM 计算预测的瓶颈：

> M4 PAM 的功能效力不是静态结合能。它由 chemotype、正构探针、受体微状态、协同性、内在激动活性和实验 readout 共同决定。当前公开数据能够支持 ACh/calcium 条件下的回顾性系列排序，却不足以识别跨探针 cooperativity，也不足以让静态结构特征稳定提升新系列效力预测。

因此 PACER-M4 v1 不是“把 docking 分数、QSAR 概率和 QED 加权”的流水线，而是一个**可拒绝预测的多证据候选提出框架**：

1. 冻结数据与防泄漏 benchmark；
2. ligand-only 效力模型作为 null baseline；
3. 严格阴性风险只作安全标记；
4. ACh 状态结构只作独立口袋相容性门；
5. 生成器必须经过适用域审计；
6. 不用加权总分，使用 Pareto 与化学多样性选候选；
7. 动态对照 MD 和功能湿实验是确认 PAM 的必要末级验证。

## 2. 冻结 benchmark

- 连续效力：430 个唯一分子，12 个可评估来源/药化系列；
- 严格分类：463 PAM / 66 实验非活性；
- activity cliffs：694 个近邻分子对，其中 93 个 cliff；
- 时间字段：430/430 分子有年份；199 个用于 2018 年及以前训练，187 个用于 2019–2021 验证，44 个用于 2022 年及以后测试；
- 主指标：整来源留出后的系列内 Spearman、macro、worst-series 和 MAE；pooled 指标不能替代系列内泛化。

数据的关键限制是：精确效力主视图几乎都为 calcium/ACh，不能从该数据验证跨正构探针的 probe dependence。

## 3. 强 baseline 与被证伪的模块

### 3.1 效力系列外推

| 方法 | macro Spearman | worst-series | 正相关系列 | macro MAE |
|---|---:|---:|---:|---:|
| ECFP + RF | 0.170 | — | — | 0.453 |
| ChemBERTa + Ridge | 0.148 | -0.352 | 10/12 | — |
| Chemprop D-MPNN | 0.147 | -0.206 | 8/12 | 0.537 |
| LightGBM v1 | 0.199 | **-0.130** | 11/12 | **0.441** |
| LightGBM v2 | **0.207** | -0.138 | 11/12 | 0.458 |

LightGBM v1/v2 构成 ligand-only null baseline；v2 的排序略好，v1 的 MAE 与最差系列更好。没有深度模型稳定超过它们。

#### PACER-FS：系列偏移分解与三锚点适配

进一步诊断表明，绝对 pEC50 同时包含来源/assay/chemotype 的系列基线偏移与系列内 SAR。PACER-FS 在训练时对每个来源组件减去平均 pEC50，只学习可迁移的系列内相对效力；对新系列按零样本预测的 25%、50%、75% 分位选择 3 个功能实验锚点，再恢复绝对 pEC50 尺度。

在共同的 11 个 `n>=12` 完整系列、相同三锚点查询集上：

| 方法 | macro Spearman | worst-series | 正相关系列 | macro MAE |
|---|---:|---:|---:|---:|
| 同架构绝对 QSAR | 0.199 | -0.157 | 10/11 | 0.453 |
| **PACER-FS 系列中心化** | **0.263** | **-0.005** | 10/11 | **0.443** |
| 双通道 ΔSAR（探索性） | 0.268 | -0.041 | 10/11 | 0.441 |

中心化相对绝对 QSAR 的 Spearman 增益为 `+0.064`，按系列 bootstrap 95% CI `[+0.017,+0.113]`；12 组固定 LightGBM 参数中 11 组增益为正。仅未见 Murcko 骨架时仍为 `0.254→0.295`，同时要求 Tanimoto `<0.70` 时为 `0.218→0.249`。自由内层选择复杂 few-shot 适配器反而得到 `0.175 vs 0.184`，被淘汰。

在唯一 2022+ 系列中，中心化单独未改善；开放 3 个锚点后，均值守恒 ΔSAR 将排序从 0.148 提至 0.303/0.316，但该结果只有一个系列，只列为探索性证据。完整方法、拒绝规则和复现协议见 `docs/PACER_FS_METHOD.md`。

### 3.2 静态结构检验

- 7TRQ、7TRP、7TRS 三状态共 1290 个 pose，零失败；
- 7TRQ 共晶重对接质心约 1.61 Å、IFP Jaccard 0.923，说明 pose 流程具备基本质量；
- Vina affinity 与功能效力近随机，不能作为 PAM 效力分数；
- ACh 状态 7TRS 的低维 IFP-only 模型达到 macro Spearman 0.191、worst -0.085，是最可转移的独立结构信号；
- 但 7TRS 全特征融合降至约 0.110，三状态无条件拼接降至 0.112；
- StateResidualMoE：0.203 vs 0.199，delta 95% CI `[-0.040, +0.046]`，不晋级；
- 新增嵌套四专家稳健融合：0.151 vs 同轮 ligand-only 0.207，delta 95% CI `[-0.115, -0.003]`，被明确淘汰。

这说明结构质量门通过，并不等于结构特征可预测功能效力。静态 IFP 的作用被限定为候选结构相容性证据。

### 3.3 大型外部 M4R active/decoy baseline

在 Miao 2026 公开数据上重算 11.6–11.8 万分子的结构筛选：单 PDB Vina 的 ROC-AUC 为 72.53%，但 EF(0.5%) 仅 0.778；ensemble BEmin 的 AUC 为 73.38%，EF(0.5%) 反而只有 0.349；BEavg 的 AUC 为 69.58%，EF(0.5%) 为 0.959。也就是说，全局 AUC 看似高于随机，真正实验预算最关心的前 0.5% 却不优于随机。这一外部结果只评价 active/decoy enrichment，不评价 PAM 功能效力。

后续按 57,854 个 Murcko scaffold 严格五折训练早期富集模型时，9个理化描述符和 ECFP 已分别得到 EF(0.5%) 48.47 和 50.89（理论上限约51）；但在零外部骨架重叠的 463 PAM/66 实验非活性测试中，ECFP PU 与 DUD-E 监督模型的系列内 macro AUC 仅 0.530/0.527。该“假完美”确认人工 decoy 主要提供数据来源捷径，不能作为生物学 SOTA。完整审计见 `docs/PACER_M4_ENRICHMENT_ROUTE_AUDIT.md`。

### 3.4 holo-GaMD ensemble 的功能真值压力测试

复用 Miao 2026 开放的 ACh–M4R–G protein–MK-97 holo-GaMD 十簇与 PMF，按原文 Vina 协议重跑：Vina 1.2.7、30 Å、exhaustiveness 8、9 poses、energy range 3 kcal/mol，并用 molscrub 0.1.1 在 pH 7 枚举质子化/互变异构状态。主真值只保留 27 个 A-tier primary-source-confirmed inactive，并一对一匹配 27 个 PAM；全部 54 分子、810 个分子状态/构象组合零失败。

| 方法 | ROC-AUC | PR-AUC | EF10% |
|---|---:|---:|---:|
| GaMD cluster-00 | 0.499 | 0.511 | 1.00 |
| PMF-BEmin | 0.499 | 0.511 | 1.00 |
| PMF-BEavg | 0.457 | 0.460 | 0.33 |

BEmin 相对 cluster-00 的 ΔAUC 恰为 0；BEavg ΔAUC -0.044，paired bootstrap 95% CI `[-0.147,+0.045]`。M4R 自由能景观由 cluster-00 主导，PMF 惩罚使所有分子的 BEmin 都回到 cluster-00；对高能簇取算术平均则稀释信号。由此得到明确负结论：**受体 ensemble 能扩展姿势与微状态采样，但 Vina+PMF 仍只描述条件结合相容性，不能识别功能 PAM。**

同时，529 个 assay-specific 标签中 300 个与 Miao/ASD 库精确 SMILES 重叠；33 个重叠 functional negative 全部在 ASD 中至少被标为一次 active，6 个结构同时出现在官方 active/decoy 行。冲突阴性均为 database-text tier，因此它们需要原始表格复核；但标签语义审计已经证明 broad allosteric membership 不能替代指定 ACh/人M4功能 PAM endpoint。

### 3.5 适用域、时间与 cliff

- 固定相似度阈值的局部 kNN 在阈值 0.55 时覆盖 44.8%，macro Spearman 仅 0.143，最差 -0.252；
- 邻居效力方差与真实误差 Spearman 仅 0.024，不能充当校准不确定性；
- 原始 2022+ 零样本伪前瞻测试在主新系列上 Spearman -0.028，pooled -0.078，MAE 0.590；固定 PACER-FS 时间压力测试表明中心化单独仍失败，但三锚点双通道探索性恢复到 0.303/0.316；
- scaffold-held-out 的差分 SAR 能改善部分 cliff 指标，但 source-held-out DeltaRF 仅 delta Spearman 0.148，不能跨未见系列迁移；
- 静态状态切换对 cliff 的 source-fold AUC 约 0.44，不支持“activity cliff 可由单次 docking 状态切换解释”。

回顾性 LightGBM 绝对残差 90% 分位为 0.963 pEC50；留组覆盖率宏平均 0.898、最差 0.788。该区间只作为风险尺度，不是生成分子的前瞻覆盖保证。

## 4. 双线生成的真实结果

| 路线 | 有效输出 | 到已知集平均最大 Tanimoto | 适用域结论 |
|---|---:|---:|---|
| BRICS/受控片段生成 | 2954 | 0.405 | 210 个进入局部域，1574 个进入探索域 |
| SMILES LSTM | 306 | 0.248 | 全部在当前证据适用域外 |
| SMILES GPT | 16 | 0.237 | 全部在当前证据适用域外 |
| masked diffusion | 0 | — | 6400 次采样没有有效 SMILES，baseline 失败 |

这不是“diffusion 不先进”，而是当前字符级、少样本、CPU 版 masked diffusion 没有学会可逆化学语法。其失败已通过最大轮次和审计文件冻结，不能进入候选线。神经生成器的新颖性很高，但新颖性与可验证性发生冲突；最终 200 个入围分子全部来自受控片段路线。

## 5. 最终候选算法

### 5.1 预筛

3271 个去重生成分子经过：

- 排除已知分子的 exact canonical match；
- PAINS、分子量、logP、TPSA、转动键与 QED 审计；
- 与 430 个已知效力分子的 ECFP4 相似度和近邻数；
- 严格非活性分类器仅作风险标记；
- 局部效力均值仅作最近邻溯源，不作为已校准预测。

得到 160 个局部域和 40 个探索域分子进入 ACh 状态结构门控。

### 5.2 结构门控

200/200 分子在 7TRS ACh 状态完成 Meeko + Vina pose proposal，零失败。Vina 只产生姿势；Pareto 目标不包含 affinity。候选口袋残基覆盖率为 0.615–0.923，均值 0.835；已知效力集均值 0.807，说明候选总体可占据预定 ECV 口袋，但该差异不能解释为更高 PAM 效力。

### 5.3 Pareto 与最终组合

Pareto 目标为：ACh 状态口袋覆盖、QED、受限新颖性、低严格阴性风险；随后按 ECFP 多样性和角色配额选择。效力预测、Vina affinity 均未进入总分。

最终输出 24 个计算假设：

- 12 个 local exploitation；
- 6 个 local diversification；
- 6 个 exploratory hypothesis；
- 覆盖 20 个 Murcko 骨架；
- 到最近已知分子的 Tanimoto 0.356–0.821；
- 口袋覆盖率 0.769–0.923；
- 全部明确标注“需要三元复合物功能 assay”。

候选表位于 `results/pacer_candidates_v01/final/final_candidate_hypotheses.csv`；同时提供带全部属性的 24 分子 SDF 和较保守的 12 分子 local-exploitation SDF，可直接交给可视化、合成评估或送测环节。

冻结图件位于 `results/pacer_final_figures/`：系列外推 baseline、证据漏斗、最终候选组合和生成器审计均同时提供 PNG 与 SVG。

### 5.4 十构象动态证据层

24 个候选全部完成十簇 GaMD docking，共 300 个 pH7 状态/构象任务、零失败；15 个进入动态四目标 Pareto 集。目标仅包括 PMF-BEmin、平均原生口袋覆盖、效力参考下界和严格 inactive 风险，不拟合新的候选加权总分。

所有 24 个候选的 BEmin 最优构象均为 cluster-00；单 7TRS affinity 与 GaMD-BEmin 的 Spearman 为 0.689（p=0.00020），说明动态层与静态结构证据相关但不相同。跨状态 Vina range 中位数为 1.403 kcal/mol，平均口袋覆盖为 0.686–0.806。动态 Pareto 集给出 15 个结构稳健假设，其中效力参考下界较高且严格 inactive 风险低的首组六个为 PACER0076、PACER0057、PACER0024、PACER0026、PACER0058、PACER0145；该顺序是证据展示，不是经功能真值校准的 PAM 排名。

完整表位于 `results/m4_gamd_ensemble/candidate_ensemble_evidence.csv`。这些结果支持缩小后续显式膜 MD/湿实验组合，但不能把任何候选升级为确认 PAM。

## 6. 可支持的生物学结论

1. **M4 PAM affinity、efficacy、cooperativity 与 probe dependence 不能由一个静态或 PMF 重加权的结合能代替。** A-tier 真值上的十构象 GaMD-BEmin AUC 仍为 0.499，BEavg 反而降至 0.457；这支持“结合相容性与功能效力脱钩”，不是直接证明某个动态机制。
2. **正构探针状态影响结构证据的可转移性。** ACh 状态的低维口袋 IFP 比两个 iperoxo/PAM 共晶状态更接近功能效力 baseline，但这种信息不能通过高维拼接或简单门控稳定转成增益。
3. **M4 PAM SAR 强烈依赖 chemotype/系列上下文。** 新系列和后期系列的失败、activity-cliff 与局部模型结果相互一致。
4. **生成模型最大的风险是超出可验证化学域。** 当前神经生成器产生的分子更“新”，却没有足够近邻证据；受控片段生成更适合当前小数据阶段。
5. **ensemble docking 不是功能 MD 的替代品。** Published holo-GaMD 快照已用于受体柔性压力测试但没有带来功能富集；下一步必须比较阳性、A-tier 阴性、无 PAM 和不同 probe 下的三元复合物动态扰动，而非继续堆叠 docking 分数。
6. **跨文献的效力基线偏移会掩盖可迁移的系列内 SAR。** 去除训练系列固定效应后，新系列排序获得统计稳定提升，支持将相对 SAR 学习与绝对效力标定分开。

### 6.1 新增证伪：为什么不把“动态特征”直接做 PAM 分类器

我们从十个 GaMD cluster 构建了相对 cluster-00 的能量、口袋覆盖、位移和接触变化 state-response fingerprint。严格 matched-pair、scaffold-grouped 和 negative-component transfer 下，SRF ROC-AUC 分别为 0.500、0.411 和 0.173；加入 ECFP 后也低于 ECFP 单独结果。ECFP 在 scaffold split 的 0.723 又在 component transfer 降到 0.516，表明其主要利用系列捷径。该路线作为负结果冻结，不能称动态 PAM 识别器。

进一步测试了三锚点 Bayesian expert gating：每个外层系列用开发系列形成专家可靠性先验，再用 3 个功能锚点更新 LightGBM/RF/ExtraTrees/Ridge 权重。严格嵌套结果为 0.236，低于固定 PACER-FS 的 0.263，最差系列为 -0.033。结论不是“更复杂模型更好”，而是当前 11 个系列下，低方差的固定系列分解优于 few-shot 动态选模。

### 6.2 Anchor-guided trajectory metric：有效连续轴，但尚非稳定功能簇

我们按“不同 PAM/replica、相同耦合状态应接近”的设想完成了四级递进检验，而没有把单次漂亮聚类当作结论。首先，replica-invariant无标签triplet虽将held-out occupancy JSD从PCA的0.313降至0.285，但种子AMI仅0.234，未得到稳定簇。其次，用7TRQ/VU0467154和7TRP/LY2033298相对7TRS/ACh-only定义跨chemotype静态方向：21个C-alpha的210个内部距离中，141个在两个PAM中同向变化；第三PAM MK-97的六条GaMD轨迹有5/6沿该方向，唯一反向的replica 2在前100 ns为正、随后转入长期反向盆地。

冻结的无标签triplet cluster对这一事后外部注释的解释率eta-squared为0.321，低于PCA的0.410，证明普通时间/几何triplet没有自动学到PAM方向。随后预注册AnchorGuidedTriplet：训练五条replica中，以“跨replica、相同锚投影”为positive，以“几何相近但锚投影相反”为hard negative，第六条replica完全留出。其eta-squared达到0.695，对比同特征PCA 0.503，6/6折改善；但种子AMI仅0.449，时间持续率从0.635降至0.497，严格技术Go仍失败。

因此本工作建立了一个可检验的方法学增量：**静态跨PAM共识可作为弱监督重塑轨迹embedding，并跨replica恢复连续耦合轴；但离散cluster尚不稳定，不能命名为PAM功能微状态，更不能预测效力。** 下一步不是继续在单一MK-97体系调参，而是补齐VU0467154、LY2033298、ACh-only、弱PAM/无效物的匹配膜MD，以轨迹袋级cooperativity/efficacy端点训练StateMIL并做leave-chemotype-out验证。

后续独立压力测试进一步改变了最终方法。第三PAM结构7V68/LY2119620在冻结轴上为1.015，高于iperoxo-only与inactive结构；但allosteric agonist 7V6A/compound-110为1.177，证明单一结构轴会把内在激动误认成PAM。probe-matched `7TRQ/7TRP-7TRK`差分在独立`7V68-7V69`中复现约+0.54结构增量，仍无法单上下文排除compound-110。因此最终动态方法升级为PACER-DC：同一候选必须分别在`+orthosteric probe`与`without probe`上下文评价，联合报告coupled shift、orthosteric stabilization与intrinsic activation risk，不输出未经校准的单一PAM分数。

Wang 2022官方Source Data提供了首个外部动态端点：LY2119620在三条独立1 us轨迹中使iperoxo RMSD平均降低0.651 A，replica级精确单侧p=0.05，50 ns block-bootstrap 95% CI `[0.527,0.762]` A。MK-97六条公开GaMD中，完全不读取配体坐标的蛋白coupling coordinate与ACh稳定性在6/6副本同向，mean rho=0.196，循环移位p<0.0001。这支持“动态正构稳定化”而非结合能作为机制端点。hard cluster仍无增量；自相关校准的sticky soft prototype连续性达到0.928并保持6/6 ACh单调，但effective prototype count=2.979未过3.0冻结门槛，故只保留为attention表示头。

三锚点策略也做了独立 baseline。predicted-span 的系列偏移绝对误差为 0.219，随机选点为 0.226，差值 -0.007 的 95% CI 为 `[-0.108,+0.109]`；纯化学多样性数值最优 0.216。由于差异不显著，PACER-FS 的有效创新是“系列内相对 SAR 与绝对尺度分解”，不是“我们已找到最优主动选点算法”。

### 6.3 新增预注册模型：DeltaR与MetaMetric

PACER-DeltaR先用inner-LOSO得到跨系列absolute残差，再学习同系列分子对的残差变化，并强制反对称。其macro Spearman为0.216，相对Absolute 0.199的增益仅0.017，95% CI `[-0.015,+0.052]`，未晋级。

PACER-MetaMetric用episodic网络学习三个锚点之间的化学距离和残差传播，固定220轮、外层查询标签完全隐藏。它得到macro Spearman 0.254、9/11系列为正。查看该结果后只允许一次预注册的固定50:50互补性检验：MetaMetric与PACER-FS直接平均。PACER-Hybrid点估计0.266、MAE 0.414，略高于PACER-FS的0.263/0.428，但相对Absolute的95% CI `[-0.0067,+0.1623]`仍跨零，worst-series为-0.116而PACER-FS为+0.011。因此不把0.266包装成新SOTA，冻结PACER-FS为主方法、Hybrid为探索性候选。

### 6.4 系列控制的残基—效力生物学

对430个分子的7TRQ/7TRP/7TRS pose，在每个药化系列内中心化pEC50，并以接触数/重原子数消除直接大小偏差。104个残基×状态假设经过5000次系列内置换和BH-FDR后，7TRQ状态出现三个跨系列主信号：

- 更靠近Y89与更高效力相关：pooled rho -0.174，q=0.0208，11/11系列同向；
- 更靠近D432与更高效力相关：rho -0.182，q=0.0208，9/11系列同向；
- 更高Q184归一化接触与更高效力相关：rho +0.182，q=0.0208，9/11系列同向。

这些位点与预先定义的ECV anchor（Y89/Q184）及probe/species-sensitive区域（D432）一致。但进一步同时回归系列固定效应和MW、logP、HBD/HBA、TPSA、转动键、芳环、FractionCSP3、重原子数后，没有特征通过全局FDR；Q184仍最接近（adjusted rho +0.169，q=0.0624）。因此最强可支持结论是：**Y89–Q184–D432几何轴携带M4 PAM效力相关信息，但它与chemotype/理化性质纠缠，尚不是普适、因果的单残基药效团。**

24候选在7TRQ重对接中零失败，且没有一个同时高于已知集三项中位数。PACER0076满足Y89与D432几何检查但Q184不足；PACER0039满足Y89/Q184但D432不足。该结果支持将两者作为不同耦合几何的对照化实验假设，而不是继续寻找一个静态“全满分”分子。

### 6.2 冻结的前瞻闭环

24 个候选按最近已知系列归为 11/8/3/2 四组。由于 PACER-FS 的系列偏移不能跨 chemotype 共享，首轮只在最大 `SCOMP0003` 参考系列中选择 3 个确定性的 ECFP 多样性锚点，并同测 3 个封存查询；另加入 2 个跨系列探索分子及 4 个历史功能对照。多样性策略仅在回顾数据中数值略优，不能称显著最优。当前没有任何新湿实验结果，故 `candidate_pacer_fs_predictions.csv` 只含系列内相对 SAR，校准状态明确为等待三锚点。

批次、空白结果模板、其余 16 个候选和选择审计位于 `results/pacer_assay_closed_loop_v01/`。这完成的是“可执行且可回传的计算—实验闭环”，不是已经完成的湿实验验证。

### 6.5 独立 Suven 外测与 PACER-AssayConditional

为解决“只有同系列内部验证”的根本缺口，本项目从独立申请人 Suven Life Sciences 的 WO2025099660A1 恢复了 45 个 human M4 CRE-Luc/ACh-EC20 PAM 功能点、14 个 GloSensor 功能点和 8 个 M2 反筛点。45 个系统命名均经 OPSIN 转换、RDKit 规范化；与结构图/质谱交叉审计后 45/45 可用，与历史 430 分子零精确重叠，最近训练 ECFP4 Tanimoto 中位数 0.438。例 88 的报告质量与其名称/结构/同分异构体不一致，被保留为源文疑似笔误而非静默修正。

模型预测前冻结了主方法、基线和支持判据。零样本结果如下：

| 方法 | CRE-Luc Spearman | MAE | 50 nM ROC-AUC | Top-quartile recall |
|---|---:|---:|---:|---:|
| Absolute-QSAR | **0.257** | 1.193 | 0.569 | 0.583 |
| PACER-Centered | 0.071 | 0.837 | 0.528 | 0.250 |
| 1NN | -0.657 | 1.518 | 0.180 | 0.000 |
| Similarity-kNN | -0.274 | 1.079 | 0.334 | 0.250 |

Absolute-QSAR 的 95% bootstrap CI 为 `[-0.037,+0.521]`，因此未通过“CI 下界大于零”的外部支持条件；它相对 similarity-kNN 的 Δρ 为 `+0.530`，95% CI `[+0.239,+0.796]`。最关键的结论是：跨 chemotype 后，最近邻与相似度平滑不仅失效，而且方向反转。与此同时，14 个双测定化合物的实测 CRE-Luc–GloSensor Spearman 为 0.619，说明真实 M4 PAM 功能 SAR 可迁移但明显受细胞/读出平台影响。

在完整保留上述失败后，项目把专利数据降级为开发集，构建下一代 assay-conditional 模型：共享分支学习历史 ACh/calcium 功能先验，域特异分支学习新 assay/chemotype 的残差；历史样本权重与正则强度只在每个外层训练区内选择。严格分组结果为：

| 查询场景 | Local molecular Ridge | Pooled Ridge | PACER-AssayConditional |
|---|---:|---:|---:|
| 未见 core state | 0.733 | **0.749** | 0.738 |
| 未见 headgroup | **0.734** | 0.530 | 0.575 |
| 未见 core family | 0.683 | 0.705 | **0.747** |

这给出一个具体生物学/方法学分工：历史数据有助于迁移受体侧核心家族，但新 headgroup 的功能贡献主要依赖本系列数据；不能无条件池化。专利中 headgroup 的效力均值也呈稳定趋势：benzofuran/benzodioxin 类整体优于 methoxy/fluoropyridyl 类，但硬编码类别在未见因子上不稳，连续分子表示更可靠。

作为第二 campaign 回顾性复制，使用 VU6025733 第一轮 18 个功能点支持、预测后续 8 个 Generation-3 查询时，Pooled-Ridge 与 PACER-AssayConditional 的 ceiling-Spearman 均为 0.946（精确置换 p=0.00129），Pooled-Ridge MAE 0.360；旧 Augmented-Absolute 为 0.802/0.870。相对提升的 bootstrap CI 仍跨零，且这些查询在本项目早期已被查看，因此该结果只能称跨 campaign 复制，不能重新包装成独立 SOTA。

最后，Generation-4 四个非同位素立体异构体构成了机制压力测试：普通 2D QSAR 对四者输出完全相同；手性 QSAR Spearman -0.105；单构象、GaMD-BEmin、GaMD-BEavg 和平均 Vina 均为 0。该结果直接否定“加 MD ensemble 后 docking 就能预测 PAM 立体功能效力”的简化路线。

### 6.6 多维 operational-allostery：从失败中收敛的方法

项目从 ChEMBL 文献记录 `CHEMBL3769348` 恢复 Monash 2015 LY2033298 系列的 16 个结构：13 个具有完整功能 `pKB/logαβ/logτB`，3 个 N-acyl 类似物明确 Not Active；不把失活物伪造为连续数值。与历史 430 分子只有 1 个完全重叠，但整体仍属于 LY2033298 骨架邻域，因此这是 external endpoint-transfer stress test，不是 unseen-chemotype SOTA。

预注册 leave-family-out 结果中，主方法 PACER-Mechanism-PLS 对 `logαβ` 的 Spearman 为 0.555，Local Ridge 为 0.545；增益仅 +0.009，family-bootstrap 95% CI `[-0.375,+0.350]`，且最差系列为 -1。该轮明确失败。三个 N-acyl inactive 虽被 A-tier 分类器排在 active 之后，但其 PAM 概率仍约 0.84；只有三个同族阴性的描述性 AUC=1.0 不可晋升为有效分类性能。

该失败给出比单一分数更重要的生物学分解：13 个 active 中，`logαβ` 与 `logτB` 的 Spearman 为 0.918，而功能 `pKB` 与 `logαβ` 为 -0.330。O/N 链延长可保持甚至提高表观亲和力，却使协同性和内在激动性同步下降。因此静态 affinity、docking energy 或单终点 potency 都不能代表 PAM 功能质量。

随后从 ChEMBL 全 M4 靶点记录恢复两套多通路 operational-model 面板：Jörg 2023 的 9 个分子和 Liu 2024 的 18 个分子，覆盖 binding、GoB、cAMP、β-arrestin 的 `pKB/logα/logτB/logαβ`。两文献只有一个完全重叠结构。严格跨文献结果再次证伪零样本通用模型：PACER-Triad 在主 `consensus logαβ` 上为 0.124，1NN 为 0.120，增益 +0.005 的 95% CI `[-0.659,+0.667]`；cAMP 同名终点甚至出现负迁移。同一个重叠分子的 cAMP `logαβ` 在两文献为 2.34 与 0.53，β-arrestin `logτB` 为 0.99 与 0.02，表明 assay/protocol domain shift 与化学域偏移同样大。

最后执行三锚点 assay 校准。每个查询文献只揭示三个按 ECFP4 多样性、与标签无关选出的功能锚点；其余分子作为查询。对全部可用的 19 个 cAMP 查询，简单 endpoint offset 将跨文献 `logαβ` 从 -0.393 提升到 0.691，将 `logτB` 从 -0.123 提升到 0.411。1000 次随机三锚点中，`logαβ` 的中位 0.668、95% 范围 `[0.576,0.767]`；`logτB` 中位 0.514、范围 `[0.217,0.672]`。相反，`pKB` 仍为负，低秩/局部残差也没有稳定优于简单 offset。

由此最终方法不是“更深网络预测一切”，而是：**结构与动态只做结合/姿势门控；三个同 chemotype、同 probe、同 readout 的功能锚点校准 assay 域；分别预测 affinity、cooperativity 与 intrinsic agonism；用不确定性和 Pareto 规则选择实验分子。** 三锚点结果是在观察零样本失败后开发的回顾性证据，尚需新的未触碰 campaign 确认。

### 6.7 第四套 Acadia 外测：去重审计推翻了假阳性

项目在任何功能值进入模型前冻结了 WO2025122811A1 的来源、结构恢复和评估合同。49 个专利例中有 47 个结构通过 OPSIN/RDKit 与 LC-MS 质量差阈值；进一步发现 Example 12/13 与 32/33 分别是同一 active moiety 的游离形式/盐型或重复制备。若将其当独立样本，锚点可能看到查询分子的完全相同结构。主分析因此按 canonical active moiety 聚合为 45 个唯一结构，其中 43 个有精确 PAM pEC50，与历史训练集零精确重叠。

去重后的零样本最强方法是简单 similarity-kNN（ρ=0.398）。冻结的确定性三锚点 PACER-FS-DeltaSARHybrid 只有 ρ=0.055，较同样接收三锚点的 kNN 低 0.362，paired bootstrap 95% CI `[-0.721,+0.013]`；MAE 也未改善。该轮外部验证明确失败。1000 组标签盲随机锚点的 post-hoc 压力测试中，DeltaSARHybrid 中位 ρ=0.511，但仅 58.8% 锚点组在排序上优于 kNN、48.9% 在 MAE 上更好。结论是“局部差分 SAR 有潜力但锚点敏感”，不是 SOTA。

同一 campaign 还给出比单分数更强的功能约束。43 个唯一精确分子中，PAM pEC50 与候选单独给药的 intrinsic-agonist RE 呈强正相关（ρ=0.779，bootstrap 95% CI `[0.515,0.930]`）；PAM RE 与 intrinsic-agonist RE 的 ρ=0.828（`[0.688,0.905]`）。10 组结构匹配的 2-methyl→2,5-dimethyl 对中，10/10 PAM 效力提高，中位提高 5.1 倍（双侧 sign-test p=0.002）；但 8/10 的 intrinsic-agonist RE 也提高，中位 +19 个百分点（p=0.109），另有两个低激动性例外。因此这是该系列的 SAR trade-off，而非普适因果律。下一代 PACER 的优化对象必须从“最高 PAM pEC50”改为多目标功能空间：高 PAM potency/RE、低 intrinsic agonism、M2 反筛与不确定性共同 Pareto，结构/MD 仅作为机制与相容性证据。

### 9.4 PACER-DC可执行候选接口

为避免动态证据再次停留在叙述层，我们实现了replica级PACER-DC接口。输入合同显式要求
`candidate+probe`、`candidate-without-probe`、`probe-only`和`apo`四个上下文；输出分别为
coupled shift、正构探针稳定化、内在激动风险、结合相容性和不确定性，不训练未经功能
真值标定的总分。只有证据完整的候选才能进入Pareto比较。

真实公开证据冒烟测试中，LY2119620得到静态paired-delta `+0.539`和三副本MD正构稳定化
`+0.651 A`，但由于公开数据没有匹配的without-probe、apo和binding endpoint，程序返回
`evidence_complete=False`、`pareto_eligible=False`。这个“拒绝给答案”是设计要求，不是运行
失败：现有证据只支持LY2119620在该体系稳定正构配体，不能独立检验其是否存在内在激动。

首轮前瞻MD生产合同由一个已知PAM（LY2119620）、一个allosteric agonist反例
（compound-110）、两个计算候选（PACER0076/PACER0057）以及共享ACh-only/apo对照组成，
每个上下文三条配对种子，合计30条初始100 ns轨迹。该账本目前仅已冻结、尚未生产，故本
报告不声称候选通过了PACER-DC。

### 6.8 Function-Space v3 与第五套 US campaign：主方法再次失败，但标签问题被定位

Acadia 的 post-hoc 锚点诊断显示，确定性三锚点的实测 pEC50 仅跨 0.301 log，处于 1000 个随机三元组的 4.6 百分位；锚点响应跨度与查询排序的 Spearman 为 0.655。也就是说，化学上分散的锚点不必然覆盖功能响应空间。项目随后冻结 Function-Space v3：若初始锚点功能跨度不足则增加预测两端锚点，最多七个并允许拒绝预测；所有查询标签对选择与停止规则隐藏。

该策略首先在历史系列严格嵌套开发中失败：PACER-AdaptiveBracket macro Spearman 仅 0.031、平均使用 3.55 个锚点；固定七锚点为 0.156，也不足以晋级。原因是用尚不可靠的模型挑选“预测极端”会把模型误差反馈到实验设计。该自适应策略因此被淘汰，而不是带到外部数据继续追分。

第五套来源 US20260055116A1 披露人 M4 PAM pERK、大鼠 M4 PAM pERK 与人 M4 PAM GTPγS 的 A/B/C/D 效力档。111 个表格条目中只纳入 26 个可由完整系统命名经 OPSIN/RDKit 恢复的高置信结构；active-moiety 去重后为 25 个，与历史 430 分子零精确重叠，最近训练结构 Tanimoto 中位仅 0.216。由于来源资格审查阶段表格标签已经可见，该分析按冻结模型进行外部压力测试，但不是 pristine blind validation。

结果明确否定当前主方法：七锚点 DeltaSAR 在 18 个查询上的有序一致性为 0.472，低于常数 0.500、TanimotoGP 0.546 和 centered-LGBM+offset 0.560；相对所有 baseline 的 paired bootstrap 支持条件为 false。零样本最高为 centered-LGBM 0.584，与 similarity-kNN 0.581 接近，样本量不足以支持优越性。更重要的是，人 pERK 与大鼠 pERK 虽有 ρ=0.531，却只有 12% 完全同档；22 个非平局中人源效力档有 21 个更强（双侧 sign-test p=1.10×10⁻⁵）。人 pERK 与人 GTPγS 在 12 个双读出分子中为 ρ=0.787、同档率 50%，6 个非平局全部由 pERK 更强（p=0.031）。在冻结的 ECFP4≥0.55 且 |ΔMW|≤90 Da 定义下，62 个局部结构对中 23 对出现至少两档 endpoint cliff。这建立了一个更准确的生物学结论：**M4 PAM 不是单一分子标签，而是 ligand × species × probe × readout × receptor-state 的条件功能表型。**这些 ordinal 数据支持条件依赖，但不能单独证明 pathway bias 或具体分子机制。

因此候选闭环改为三联主测：M4+ACh EC20 的 PAM EC50/Emax、候选单独加药的 M4 agonism EC50/Emax、匹配 M2+ACh EC20 的 PAM EC50/Emax；ACh EC80 与 GTPγS 只作为 probe-window 和近端 G 蛋白确认。现有首轮 8 个生成候选与 Acadia 功能标签分子的最高 Tanimoto 仅 0.13–0.26，全部属于外推域；其数值先验只用于实验信息增益排序，不能判定 PAM、agonism 或选择性。当前确认 PAM 数仍为零。

作为该结论的首个算法实现，PACER-ContextKRR 将每条记录表示为“分子 × species × readout”，外层按完整分子留出，内层选择核与正则化，确保同一查询分子的其他 endpoint 不会泄漏。它的 macro 有序一致性为 0.714，略高于 endpoint-only KRR 的 0.706，却低于 pooled-no-context KRR 的 0.740，未通过冻结开发门槛。human-pERK 单终点由 0.678 提升至 0.696，但 human-GTPγS 未改善。由此保留“显式上下文建模”作为需要更多成对标签的数据结构，不宣称现有小样本模型成功；当前数据更支持一个共享活性轴叠加稀疏 context cliffs。

上述第五 campaign 的外测、物种转移、局部 endpoint cliff 与条件核 baseline 汇总于 `results/pacer_function_space_v3_figure/PACER_FUNCTION_SPACE_V3_EVIDENCE.png`，图中所有模型均来自冻结协议和完整样本，不筛选有利子集。

## 7. 尚不能声称的内容

- 不能声称 PACER-M4 已达到外部 target-specific SOTA；现已有独立 Suven 功能外测，但零样本主方法的 Spearman CI 跨零、50 nM AUC 仅 0.569。后续 PACER-AssayConditional 已使用该专利开发，必须由新的未触碰系列确认；
- 不能声称 24 个候选是真正 PAM；
- 不能声称对接、Glide、GNINA 或 MD 能单独给出 PAM 功能效力；
- 不能声称 probe-aware，因为现有精确标签几乎只有 ACh；
- 不能使用 pooled AUC 0.711、早期小样本 AUC 1.0 或旧几何 AUC 0.87 作为主成绩。
- 不能把 Monash 三个 N-acyl 阴性得到的描述性 AUC=1.0 当作分类 SOTA；模型对这些阴性的绝对 PAM 概率仍然过高。
- 不能把三锚点 cAMP 校准称为独立外部 SOTA；它是在跨文献失败后开发的回顾性方法，需要新的未触碰 assay/campaign 复现。
- 不能把 Acadia 的随机锚点中位结果称为独立支持；唯一 active-moiety 的预注册确定性外测已经失败，随机锚点分析只刻画敏感性。
- 不能把 Acadia Pareto 分子称为本项目发现；它们是已披露专利阳性对照，用于定义 potency–agonism–selectivity 功能空间。
- 不能把 US20260055116 命名子集称为盲外测；来源筛选时标签已可见，且只覆盖 25 个可可靠恢复结构的 active moieties。
- 不能把 Function-Space v3 的自适应锚点称为有效创新；其严格历史开发与第五 campaign 主方法均失败。
- 不能把候选的 Acadia 相似性加权 prior 称为 PAM/agonism 预测；8/8 候选均低于 0.50 的决策域阈值。
- 不能把 PACER-ContextKRR 称为新 SOTA；它未超过 frozen pooled baseline，且属于已查看来源上的方法开发。

## 8. 下一等级证据：实验设计

每个候选至少需要：

1. ACh 浓度—响应曲线，在多个 PAM 浓度下联合拟合 operational allosteric model，估计 affinity/cooperativity/efficacy；
2. 候选单独加药，排除 PAM-agonist 或直接激动；
3. 严格非活性近邻、VU0467154/LY2033298 阳性对照和 vehicle；
4. M1/M2/M3/M5 亚型、细胞毒性及早期安全性反筛；
5. 若条件允许，在 ACh 与 iperoxo 两种探针下重复，直接验证 probe dependence。

动态模拟按 `docs/PACER_M4_DYNAMIC_VALIDATION_PROTOCOL.md` 执行。当前已复用 published holo-GaMD ensemble 完成十构象压力测试，但没有自行生成生产级膜蛋白轨迹；候选特异的对照化 MD 仍需要 GPU/集群和一致的小分子参数化。

可直接交给合作实验室的两阶段筛选、全浓度矩阵、对照和预注册终点见 `docs/PACER_M4_FUNCTIONAL_VALIDATION_PLAN.md`。

## 9. 复现

### 9.1 Thompson–Miao 2026原文与PACER-XR

作者官方六路M4 docking score已经用作者分析代码精确复算；Figshare六条公共500 ns轨迹的3000个stride帧完成独立十簇重聚类。主簇占比82.77%（作者全流79.5%），十个代表口袋与作者发布状态的匹配RMSD均低于2 Å，均值1.12 Å。

在六路官方分数共同覆盖的114321个M4分子上，等权rank共识PACER-XR取得ROC-AUC 0.7710。进一步固定前1%为Glide-BEmin、剩余99%为XR的PACER-XR Cascade取得AUC 0.7775、PR-AUC 0.1067、logAUC 19.78%，同时保持EF0.5/1%为20.40/13.73；相对最佳单路AUC提升+0.0443，95% CI `[+0.0380,+0.0506]`。除EF2%外，它在作者公开M4 benchmark主指标上领先或并列领先。该结果仍是broad-AM retrieval，不是功能PAM预测。

进一步进行了exact-SMILES功能语义压力测试。299个可映射分子中有266个功能positive和33个functional negative；分子级平均rank得到Cascade表观AUC 0.783。但33个阴性全部来自B-tier database-text记录，A-tier primary-confirmed inactive为0，且33个又全部在ASD中被标作active。故这个数字只揭示broad-AM标签与指定ACh功能endpoint的语义冲突，不能晋升为功能PAM SOTA。更严格的27+27 A-tier平衡实验仍给出GaMD/Vina约随机的结果。候选域迁移同样不充分：PACER0076、PACER0057、PACER0026最近已知PAM在Cascade中仅位于top 10.03%、48.84%、38.22%。

为使公开基准创新能够被前瞻应用而不篡改协议，项目冻结了24候选的六通道score contract。任何通道缺失时程序直接拒绝排名；不会以Vina代替Glide，也不会把7TRS本地pose分数冒充原文PDB通道。当前模板保持空白，等待具有Schrödinger Glide许可的外部节点按原协议回填。

详细协议、论文—仓库版本差异和分子结论见 `docs/THOMPSON_MIAO_2026_EXACT_REPRODUCTION.md`。
功能标签语义审计见 `results/pacer_xr_functional_semantics_v01/VALIDATION_REPORT.md`。

```powershell
python project/scripts/run_pacer_complete.py
```

PACER-FS 方法证据单独复现：

```powershell
python project/scripts/run_pacer_fs_complete.py
```

GaMD ensemble 验证与候选动态证据复现：

```powershell
python project/scripts/run_m4_gamd_complete.py --skip-preparation
```

独立 Suven 外测、PACER-AssayConditional 嵌套开发、VU sequential 复制与立体化学压力测试：

```powershell
python project/scripts/run_pacer_external_frontier.py
```

多维 `pKB/logαβ/logτB` 数据构建、跨文献机制迁移与三锚点稳健性：

```powershell
python project/scripts/build_monash_ly2033298_allostery.py
python project/scripts/run_pacer_monash_allostery.py
python project/scripts/build_m4_mechanistic_multidomain.py
python project/scripts/run_pacer_mechanistic_multidomain.py
python project/scripts/run_pacer_mechanistic_fewshot.py
python project/scripts/run_pacer_assay_anchor_robustness.py
```

冻结首轮实验批次并生成回传前相对 SAR：

```powershell
python project/scripts/build_pacer_assay_closed_loop.py
python project/scripts/apply_pacer_fs_candidates.py
python project/scripts/apply_pacer_mechanism_calibration.py
```

昂贵对接已缓存；需要重跑时加 `--force-docking --workers 8`。所有失败模块与负结果均保留，不从结果目录删除。
