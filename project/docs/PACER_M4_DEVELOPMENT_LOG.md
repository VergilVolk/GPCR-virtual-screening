# PACER-M4 开发日志

## 2026-08-28：Benchmark v1.0 与第一组强基线

### 已完成

1. 冻结研究章程与摘要证据库；
2. 实现 `scripts/build_pacer_benchmark.py`；
3. 生成分子级回归、严格 PAM 分类、10 uM 阈值三个 benchmark 视图；
4. 实现零 scaffold/source-component 泄漏的预分配 folds；
5. 实现 `scripts/run_pacer_baselines.py`；
6. 完成 Tanimoto kNN、Random Forest、Extra Trees 的 OOF 预测与 1000 次 bootstrap CI。

### 数据审计发现

- 回归：430 分子、214 个 Murcko scaffolds、18 个主要来源；
- 严格分类：463 PAM / 66 实验非活性；
- 旧 QSAR 只得到 428 个回归分子，是因为在已冻结 modeling view 上再次要求 `identity_status == confirmed_structure`，误删了两个已经通过行级身份审计的分子；
- 来源年份在冻结 `sources.csv` 中基本缺失，时间外推测试暂时硬阻塞；
- 严格阴性只分布在四个可独立隔离的来源组件中，因此来源级分类最多做四折，而且有一折仅一个阴性，必须报告高方差和置信区间。

### Baseline OOF 结果

| 任务/划分 | 最佳模型 | 主要结果 |
|---|---|---|
| pEC50 / scaffold | Random Forest | Spearman 0.637；MAE 0.389；R2 0.412 |
| pEC50 / source | Random Forest | Spearman 0.472；MAE 0.452；R2 0.221 |
| PAM vs inactive / scaffold | Extra Trees | ROC-AUC 0.919；MCC 0.595；BalancedAcc 0.837 |
| PAM vs inactive / source | Random Forest | ROC-AUC 0.630；MCC 0.000；BalancedAcc 0.500 |
| 10 uM / scaffold | Extra Trees | ROC-AUC 0.857；MCC 0.682；BalancedAcc 0.810 |
| 10 uM / source | Extra Trees | ROC-AUC 0.600；MCC -0.024；BalancedAcc 0.498 |

### 科学解释

宽松的骨架级结果不能代表真实跨药化系列发现能力。整来源组件留出后，分类几乎退化到默认预测阳性，说明类别不平衡、系列偏置和来源特异 SAR 是当前核心瓶颈。PACER-M4 的创新模块必须优先改善 source/series-held-out、activity-cliff 和 calibration，而不是只提高 scaffold ROC-AUC。

### 下一开发门

1. 补齐并人工核验来源年份，冻结 temporal split；
2. 构建 matched molecular pair / activity-cliff benchmark；
3. 增加 LightGBM/XGBoost、Chemprop/D-MPNN 和预训练分子编码器 baseline；
4. 开发 probe/assay-aware 多任务模型；
5. 之后才接入结构耦合图和双线条件生成。

## 2026-08-28：Activity-cliff benchmark 与差分 SAR 基线

### Pair benchmark

冻结定义：ECFP4 Tanimoto >= 0.70 为近邻；`abs(delta pEC50) >= 1.0` 为 activity cliff；`>=1.5` 为 strong cliff；`<=0.3` 为 smooth pair。

- 430 分子形成 694 个近邻 pair；
- 93 个 activity cliffs，其中 24 个 strong cliffs；
- 46 个 cliff 具有相同 Murcko scaffold；
- 86/93 个 cliff 来自同一主要文献来源，提示 cliff/SAR 具有强烈系列结构；
- 严格 pair test 要求 pair 两侧分子同时位于同一 held-out fold；scaffold test 保留 376 pair/52 cliff，source test 保留 660 pair/93 cliff。

### 差分模型结果

| 划分/方法 | Delta Spearman | 方向准确率 | Cliff 方向准确率 | Cliff-vs-smooth AUC |
|---|---:|---:|---:|---:|
| scaffold / AbsoluteQSAR | 0.268 | 0.640 | 0.712 | 0.623 |
| scaffold / DeltaRidge | 0.282 | 0.671 | 0.731 | **0.705** |
| scaffold / DeltaRF | **0.309** | 0.626 | **0.808** | 0.649 |
| source / AbsoluteQSAR | 0.130 | 0.554 | **0.667** | 0.563 |
| source / DeltaRidge | -0.031 | 0.499 | 0.398 | 0.541 |
| source / DeltaRF | **0.148** | **0.574** | 0.559 | **0.589** |

### 结论

差分学习在同类化学上下文中能提高 cliff 检测或方向判断，但不能自动跨越整套未见药化系列。尤其线性结构变换在 source-held-out 上失效，说明“某个取代基变化提高效力”的规则高度依赖母核和系列上下文。

因此 PACER-M4 不采用独立 MMP 模型作为最终方案，而把 pairwise/activity-cliff loss 作为多任务正则，与预训练分子表示、assay 条件和结构耦合上下文联合训练。主选择指标仍以 source-held-out 为优先，scaffold 提升只能作为辅助证据。

### 下一步模型接口

每个训练 batch 同时支持：

1. 单分子样本：`molecule + assay/probe context -> pEC50 / class / uncertainty`；
2. 近邻 pair：`molecule_a, molecule_b, shared context -> delta pEC50 / direction / cliff`；
3. 后续结构样本：`pose/contact graph + orthosteric state -> coupling representation`。

下一开发模块为 probe/assay-aware 多任务数据层和强模型 baseline；pair loss 将在该框架内做开/关消融。

### 条件数据可识别性审计

扩展精确效力视图：472 个 measurement groups、460 个分子。

- assay family：calcium 437、cAMP 19、ERK 12、GIRK/thallium 2、unspecified 2；
- cell system 覆盖 69.1%：CHO 264、CHO-K1 43、HEK293 19；
- orthosteric agonist 覆盖 83.9%，但已知值全部为 acetylcholine；
- agonist concentration/ECxx 覆盖 24.8%，已知值全部为 EC20。

结论：当前数据支持 assay-aware 多任务和 source-domain generalization，不支持跨探针 probe-dependence 结论。下一模型优先开发 source-invariant / GroupDRO 类训练，因为首轮 baseline 的最大失败正是整来源外推；probe encoder 保留接口但不作为当前已验证创新点。

## 2026-08-28：来源域泛化与 PACER-M4 v0.1 多任务原型

### 更严格的阴性系列留出

66 个严格阴性主要集中在三个具有足够阴性的独立来源组件（16、16、33 个阴性）。新增 `negative-bearing-series holdout`：每次完整留出其中一个 PAM/Inactive 混合药化体系；其余 230 个 positive-only 或阴性不足分子只作为训练数据。最终测试覆盖 299 个分子。

### 传统与单任务 baseline（三系列留出）

| 方法 | ROC-AUC | PR-AUC | MCC | Brier |
|---|---:|---:|---:|---:|
| Tanimoto kNN | 0.575 | 0.810 | -0.005 | 0.196 |
| Random Forest | 0.622 | 0.870 | 0.000 | 0.181 |
| Extra Trees | 0.624 | **0.877** | 0.000 | 0.183 |
| Neural ClassifierOnly（3 seeds） | 0.594 | 0.819 | 0.110 | 0.203 |

### GroupDRO 结论

- 普通 source folds：3-seed ERM AUC 0.666，GroupDRO 0.678；配对 bootstrap delta AUC 95% CI `[-0.021, 0.046]`，不显著；
- 三套阴性系列留出（单 seed 探索）：ERM 0.660，GroupDRO 0.594；
- GroupDRO 改善部分 Brier/ECE，但损害 scaffold 或真实系列排序，因此不能作为核心创新单独保留。

### PACER-M4 v0.1

共享分子编码器联合：

1. PAM vs strict inactive 分类；
2. calcium/ACh pEC50 回归（masked loss）；
3. 近邻 pair 的 delta pEC50，cliff pair 加权。

三随机种子、三套阴性系列留出：

| 方法 | ROC-AUC | PR-AUC | Brier |
|---|---:|---:|---:|
| ClassifierOnly | 0.594 | 0.819 | 0.203 |
| MultiTask | 0.622 | 0.820 | 0.196 |
| **MultiTaskPair** | **0.651** | 0.833 | 0.190 |

MultiTaskPair 相对同架构 ClassifierOnly 的 delta AUC 为 **+0.056**，paired bootstrap 95% CI **[+0.016, +0.098]**，bootstrap `P(delta>0)=0.997`。这是第一个通过配对统计检验的模块增益，支持“功能效力 + 局部 SAR 联合学习有助于新系列 PAM 排序”。

但它尚未成为全面 SOTA：Extra Trees 的 ROC-AUC 较低（0.624），但 PR-AUC（0.877）和 Brier（0.183）仍优于 MultiTaskPair；固定 0.5 阈值下各模型仍明显偏向阳性。下一阶段必须改善早期阳性质量、校准和阈值决策，并通过 nested validation 确定多任务损失权重，不能使用最终三系列测试调参。

## 2026-08-28：PACER-M4 v0.2 嵌套验证与 Simpson 陷阱

v0.2 在每个外部阴性系列完全封存的条件下，用内部 scaffold CV 选择多任务权重、Platt 校准、MCC 阈值和神经/ExtraTrees 融合比例。聚合结果表面达到 fused ROC-AUC 0.711、PR-AUC 0.897、Brier 0.153。

但三个外部系列各自 ROC-AUC 仅为 **0.530、0.495、0.534**；macro AUC 约 **0.520**，最差系列 0.495。pooled AUC 的提升来自系列间阳性比例和预测尺度差异，而不是系列内 PAM/Inactive 排序，是 domain-confounding/Simpson 型假提升。

因此：

- pooled AUC 0.711 不作为算法成绩；
- 后续跨系列结果强制报告 per-series、macro 和 worst-group；
- 三个阴性承载系列不足以支撑高容量“通用 PAM 分类器”的方法开发；
- 主学习任务转为 11 个来源组件上的连续 pEC50 leave-one-series-out 排序，严格阴性分类作为独立安全门；
- 生物学目标从“脱离系列直接判定 PAM”调整为“在新 M4 PAM 化学系列内可靠排序功能效力，并识别已知严格非活性风险”。

### 12 系列零样本效力外推

对 n>=8 的 12 个来源组件逐一完整留出，以系列内 Spearman 为主指标：RF macro 0.170，NeuralERM 0.147，PACER-Rank 0.139；PACER-Rank 未超过 RF，淘汰固定 RankNet 作为独立核心模块。各系列表现高度异质（最佳约 0.58，部分为负相关）。

这支持一个药化层面的结论：M4 PAM 的局部 SAR 不是跨母核可直接转移的统一规则；activity-cliff/pair loss 必须由系列上下文条件化。下一步建立 0/1/3/5-shot 系列适配 benchmark，用少量已知功能锚点校准新系列候选排序。

## 2026-08-28：少样本适配、ChemBERTa 与结构证据线

### 简单少样本适配不能修复排序

在每个完整留出的药化系列中随机提供 0/1/3/5 个带 pEC50 标签的锚点，用 Tanimoto 加权的局部残差校准 RF：

| anchors | macro Spearman | macro MAE |
|---:|---:|---:|
| 0 | 0.170 | 0.453 |
| 1 | 0.167 | 0.436 |
| 3 | 0.184 | 0.423 |
| 5 | 0.176 | 0.414 |

少量锚点能改善绝对标尺（MAE），却几乎不能恢复系列内排序。这排除了“新系列只是整体 potency offset 不同”的简单解释：局部 SAR 次序本身会随母核和结合模式改变。

### 公共预训练表示不是充分解

冻结 44.1M 参数 `ChemBERTa-zinc-base-v1` 表示后做 LOSO：Ridge macro Spearman 0.148、最差 -0.352、10/12 系列为正；ExtraTrees macro 0.066。均未超过 ECFP+RF 的 0.170。现有证据说明瓶颈不是单纯模型容量或无监督预训练不足，而是功能标签稀疏、系列上下文和可能的变构结合模式差异。

### 7TRQ 结构线质量门

旧 OpenBabel 受体存在两个问题：零电荷，以及 Windows 环境下加氢会错误改写残基身份，故旧几何结果仅保留为探索性记录。新流程使用 Meeko 0.7.1 直接从 7TRQ chain R 制备 Gasteiger 受体，保留全部 13 个 VU0467154 共晶接触残基。共晶 PAM 在新受体上的最优 Vina affinity 为 -9.75 kcal/mol。

先导稳定性试验（12 分子 × 3 docking seeds，36/36 成功）显示平均跨 seed 标准差：Vina 0.23 kcal/mol、共晶 IFP Jaccard 0.029、配体质心 0.51 A、TYR439 最短距离 0.30 A。随后启动 430 分子全量 label-free pose 特征生成。

结构特征只作为待检验假设：原子对接触、逐残基最短距离、共晶 IFP 相似度、质心偏移和 TYR439 接触。TYR439 不预先定义为 PAM 标签；只有在整来源 LOSO 中改善 macro 与 worst-series Spearman，结构模块才可晋级。

### LightGBM 强 baseline

补充 ECFP4 + descriptors 的固定 LightGBM LOSO：macro Spearman **0.199**、median 0.139、worst -0.130、11/12 系列正相关、macro MAE 0.441。它超过 RF 的 0.170，成为当前效力系列外推的最佳传统 baseline。后续结构融合必须以 0.199 为最低比较线，不能只与较弱 RF 比较。

### 预注册 probe/state ensemble

为检验 7TRQ（VU0467154 + iperoxo）与主要 ACh 功能标签之间的 probe/state mismatch，预先准备三状态结构集：7TRQ、7TRP（LY2033298 + iperoxo）、7TRS（ACh，无 PAM）。Meeko 受体均制备成功。相对 7TRQ，7TRP/7TRS 的共同 CA 对齐 RMSD 分别为 0.368/0.601 A；7TRP 的原生 PAM 质心与映射的 7TRQ PAM 中心相距 3.855 A。后者提示两种共晶 PAM 在同一 extracellular vestibule 内采用不同微位点，可能构成 chemotype-dependent coupling 的结构基础。

分析顺序冻结为：先独立完成 7TRQ 单结构 LOSO；只有其结果封存后，再测试三状态 ensemble，防止根据测试标签事后挑结构。

### 7TRQ 单结构 LOSO（冻结结果）

430/430 分子完成 Meeko+Vina pose 和逐残基特征，无失败。主要消融：

| 方法 | macro Spearman | worst | 正相关系列 |
|---|---:|---:|---:|
| Vina only | 0.043 | -0.226 | 7/12 |
| IFP RF | 0.013 | -0.381 | 6/12 |
| TYR439 RF | 0.025 | -0.274 | 6/12 |
| literature coupling RF | 0.024 | -0.270 | 6/12 |
| full structure RF | 0.170 | -0.030 | 10/12 |
| 2D LightGBM | **0.199** | **-0.130** | **11/12** |
| 2D + structure LightGBM | 0.208 | -0.152 | 9/12 |

融合相对 2D LightGBM 的逐系列 delta 均值仅 +0.008，median -0.028；series-bootstrap 95% CI `[-0.106, +0.130]`，`P(delta>0)=0.547`。因此 7TRQ 单结构模块不晋级，不能声称有效提升。

关键解释：结构特征在部分系列明显有益、另一些系列明显有害，符合 M4 变构作用的 chemotype/probe/state dependence。当前静态结构是 VU0467154 + iperoxo + Gi1，而主要效力标签为 ACh 条件；下一项预注册检验是 7TRQ/7TRP/7TRS 状态 ensemble，而不是继续在 7TRQ 测试标签上挑残基或调融合权重。

事后机制审计（仅解释，不用于调模型）：系列平均/最大 ECFP Tanimoto 到 VU0467154 与结构融合 delta 的 Spearman 分别为 -0.140/-0.200（p=0.665/0.534）。因此“离共晶 PAM 越像，7TRQ 结构越有用”的简单相似度门控不成立，不能据此构造专家选择器。

### 7TRP 第二 PAM 微状态（冻结结果）

LY2033298 + iperoxo 状态全量 430/430 完成。`2D + 7TRP structure LightGBM` macro Spearman 0.148，低于 2D 的 0.199；但 worst-series 从 -0.130 改善到 -0.051。7TRP literature-coupling RF 本身 macro 0.112、9/12 系列正相关，但最差系列 -0.541，显示强烈的系列特异方向反转。

因此 7TRP 不能作为无条件融合模块晋级；它提供的是与 7TRQ 不同的稳健性/微状态信息。下一冻结对照为 ACh 状态 7TRS，用于区分 PAM chemotype mismatch 与 orthosteric probe-state mismatch。

### 7TRS ACh 状态与三状态 ensemble

7TRS 430/430 完成。最突出结果不是 docking affinity，而是低维 IFP-only RF：macro Spearman **0.191**、worst **-0.085**、10/12 系列正相关，接近 ligand-only LightGBM 0.199，并显著优于 7TRQ/7TRP 的 IFP 表现。Vina-only 为 -0.044；把全部 7TRS 结构特征直接融合则降到 0.110。

无条件拼接 7TRQ/7TRP/7TRS 全部结构特征同样失败：macro 0.112、worst -0.357。故“更多 receptor conformations + 高维特征拼接”不成立。当前可支持的生物学假设是：ACh 稳定的 ECV 口袋相容性比 iperoxo-PAM 共晶状态的静态能量更可转移，但其信息不能通过普通 feature concatenation 转成稳健效力增益。

### Chemprop D-MPNN baseline 与假完美拦截

Chemprop 2.3.1 严格 LOSO：macro Spearman 0.147、worst -0.206、8/12 系列正相关、macro MAE 0.537，未超过 LightGBM。首轮解析曾误读预测 CSV 中保留的真实 `pEC50` 列，产生 Spearman 1.0/MAE 0 的假完美值；在汇总前已拦截并修正为 `pred_0`，假结果未进入报告。

### PACER-StateMoE v0.1

以 LightGBM 为 null expert，使用三状态低维结构块预测其 OOF 残差；feature block、Ridge alpha 和残差强度 beta 均在每个 outer test 之外通过内层整来源验证选择，并采用来源平衡权重。

- Base macro Spearman 0.199，worst -0.130；
- StateResidualMoE macro 0.203，worst -0.165；
- paired delta 均值 +0.003，median +0.003；series bootstrap 95% CI `[-0.040, +0.046]`，`P(delta>0)=0.556`。

因此 v0.1 不晋级。低维 ACh-IFP 是独立结构证据，但其与 2D 模型误差不是简单线性残差关系。下一开发门转为：适用域/不确定性驱动的 abstention，以及 activity-cliff 中的状态切换检测；禁止继续在相同 12 个 outer series 上调融合权重。

## 2026-08-29：最终反证、生成审计与候选闭环

### 选择性预测与局部适用域

模型分歧、状态分歧和指纹距离与绝对误差的 Spearman 分别约 0.013、0.020、-0.031，不能支持可靠 abstention。固定 Tanimoto 门槛 0.35/0.45/0.55/0.65 的局部 kNN 覆盖率依次为 92.1%/70.4%/44.8%/16.1%，macro Spearman 依次为 0.121/0.079/0.143/0.103；邻居方差与误差相关仅 0.024。局部效力估计降级为近邻溯源，不作为校准预测。

### 时间外推与静态 cliff 状态切换

来源年份补齐后，430/430 分子具有时间字段。以 <=2018 训练、2019–2021 选模、>=2022 测试，后期主系列 Spearman -0.028、pooled -0.078、MAE 0.590。该测试是回顾性伪前瞻，且后期 chemotype 曾在其他 LOSO 中出现，仍显示未来系列泛化缺口。

静态 pose 的状态偏好在同数据上有弱信号，但整来源交叉验证中 SimilarityOnly/StateOnly/SimilarityPlusState 的 cliff AUC 约为 0.470/0.444/0.441，不支持用单次 docking 的“状态切换”解释 activity cliffs。

### 嵌套稳健结构融合

四专家（LGBM、来源平衡 LGBM、RF、ACh-IFP）在每个 outer source 外用 inner source folds 选择非负秩融合。结果 macro Spearman 0.151，低于同轮 ligand-only 0.207；delta -0.056，series bootstrap 95% CI `[-0.115,-0.003]`。静态结构融合正式淘汰，ACh-IFP 只保留为候选结构门。

### 生成器闭环

- 受控片段：2954；LSTM：306；GPT：16；masked diffusion：0；
- diffusion 在加入最大 100 轮终止条件后，6400 次采样得到 0 个有效 SMILES，确认采样器未收敛而非“运行太慢”；
- LSTM/GPT 到已知数据平均最大 Tanimoto 仅 0.248/0.237，全部落在当前证据域外；
- 3271 个去重生成分子中，210 个进入局部域、1574 个进入探索域；预筛 200 个全部来自片段路线；
- 200/200 在 7TRS ACh 状态完成 pose proposal；最终以结构覆盖、QED、受限新颖性、阴性风险标记做 Pareto，再按多样性得到 24 个计算假设、20 个 Murcko 骨架。

候选不使用 Vina affinity 或预测 pEC50 形成总分。生产级 MD 与功能 assay 仍是确认 PAM 的必要步骤。

## 2026-08-29：PACER-FS 新系列三锚点适配

旧 similarity-weighted residual 的 0/1/3/5-shot 结果主要改善 MAE，不能稳定改善系列内排序。自由组合 KernelResidual、TanimotoGP、DeltaSAR 与四种锚点策略后，开发集最好的 predicted-span 3-shot DeltaSARHybrid 达到 0.263，但严格嵌套选模降为 0.175、低于同查询基线 0.184；cliff 加权也不稳定。因此复杂适配器选择被判定为小样本过拟合。

可靠性门控 v0.2 将 MAE 从 0.457 降至 0.427，95% CI 完全低于零，但排序仅 0.199→0.205、CI 跨零。该实验定位了关键结构：锚点能校准绝对尺度，却不足以直接学习未见 chemotype 的排序。

据此提出系列固定效应分解：训练标签在每个来源组件内中心化，模型只学习系列内相对 SAR；3 个 predicted-span 功能锚点恢复新系列偏移。固定 LightGBM 下，11 个 `n>=12` 系列 macro Spearman 由 0.199 提至 0.263，增益 +0.064，按系列 bootstrap 95% CI `[+0.017,+0.113]`；worst 由 -0.157 提至 -0.005。12 组参数中 11 组正增益，严格新骨架/低相似度子域仍保持平均增益。

均值守恒双通道 ΔSAR 在总体上为 0.268，未显著超过中心化核心且 worst 略差，因此列为扩展。在唯一 2022+ 系列上，三锚点双通道把排序由 0.148 提至 0.303/0.316，但仅一个系列，不能作为外部 SOTA 证据。

冻结主方法为 PACER-FS CenteredCalibrated；完整协议、拒绝规则和复现见 `PACER_FS_METHOD.md`。当前结论是内部 benchmark-leading，不是外部 SOTA，也不提供未知分子的 PAM 身份保证。
