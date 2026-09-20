# PACER-M4：LY2033298 多维变构药理压力测试预注册

冻结日期：2026-08-30  
数据源：Szabo 等，MedChemComm 2015，DOI `10.1039/C5MD00334B`；ChEMBL 文献记录 `CHEMBL3769348`。

## 科学问题

本测试不再把 PAM 功能压缩成一个 docking energy 或 EC50。文章以同一批 LY2033298 类似物报告：

1. `functional_pKB`：功能实验估计的变构配体亲和力；
2. `log_alpha_beta`：ACh 条件下总体功能协同性，是本轮主连续终点；
3. `log_tauB`：变构配体自身的内在激动性；
4. `log_alpha`：结合实验中对 ACh 亲和力的协同性（仅部分分子）；
5. `active/inactive`：功能实验是否可检出活性。

核心检验是：分子模型能否区分“能结合”“能增强 ACh”“自身直接激动”三个不同机制维度，并识别小结构变化导致的功能悬崖。

## 数据与冻结边界

- 16 个结构：13 个 active，3 个明确 Not Active；不把 Not Active 伪造为连续数值。
- 与历史 430 分子做 canonical-SMILES 与 ECFP4 相似度审计；完全重叠分子不得计入严格 zero-shot 结构外推。
- 该文献团队和功能终点独立于当前模型开发，但 LY2033298 骨架并非完全新颖。因此本轮称为 **external endpoint-transfer / mechanistic stress test**，不称 unseen-chemotype SOTA。
- 在生成预测前冻结本文件、数据哈希、分组和指标。后续看到标签后的模型改动一律标记 post-hoc。

## 预先指定的评价任务

### T1：零样本机制迁移

查询集：排除历史完全重叠后的有数值分子。冻结历史 M4 PAM pEC50 模型不接触本论文标签，直接评估其排序与：

- 主终点 `log_alpha_beta`；
- 次终点 `functional_pKB`；
- 次终点 `log_tauB`。

基线：历史 Absolute-QSAR、1NN、similarity-kNN、常数中位数。主指标 Spearman；同时报告 MAE 仅在训练折内做线性尺度校准后使用。

### T2：系列内少样本机制学习

按 `modification_family` 留一族外推，所有超参数只在外层训练部分选择。比较：

1. Local Ridge；
2. Pooled Ridge（历史 pEC50 与本地机制标签直接共享）；
3. PACER-AssayConditional（共享分子支路 + 本终点域残差支路）；
4. PACER-Mechanism（联合预测 pKB、logαβ、logτB 的低秩多任务模型）；
5. 1NN 与 family-train mean。

主终点仍为 `log_alpha_beta`。主分组为 O-alkyl、5-halogen、N-alkyl、core-variant；单个 reference 不独立成测试折。N-acyl 三个 inactive 只进入分类任务。

### T3：活性悬崖分类

目标是识别 N-acyl 近邻从 active 到 inactive 的跳变。报告 leave-family-out 概率、Brier score 和逐分子结果；由于测试阴性仅 3 个，不用单个 ROC-AUC 宣称 SOTA。

## 统计规则

- 连续任务报告 pooled Spearman、每族 Spearman（族内至少 3 个时）、worst-family、MAE、预测跨度。
- 用按 family 重采样的 paired bootstrap 比较主方法与最强基线；小样本置信区间必须展示。
- 排名并列、缺失和固定值透明保留；不以删除困难分子改善指标。
- 若 95% CI 跨零，结论为“未建立增益”；若模型失败，失败结果仍进入总报告。

## 成功门槛与表述红线

只有同时满足以下条件，才能称“该数据集支持机制模型优于基线”：

1. `log_alpha_beta` pooled Spearman 高于所有预先指定基线；
2. paired family bootstrap 的增益 95% CI 下界大于 0；
3. worst-family 不低于最强基线；
4. pKB 与 logτB 的分解没有退化为同一预测排序。

无论计算结果多好，都不能把候选称为“已确认 PAM”。真正确认仍需要 ACh 浓度–反应矩阵、operational model 拟合及反筛实验。
