# PACER-M4 独立 Suven 专利外部验证预注册

冻结日期：2026-08-30  
数据源：WO2025099660A1，Suven Life Sciences  
状态：**模型预测前冻结**

## 科学问题

历史训练集主要是 Vanderbilt 系列的 hM4/Gqi5-CHO、ACh EC20 钙流功能数据。该外部集来自独立申请人、独立 chemotype，并使用两套不同的 cAMP 功能平台。因此它检验的不是同系列插值，而是：

> 仅凭二维结构学到的 M4 PAM 功能 SAR，能否跨研发组织、跨骨架和跨细胞读出排序新的 M4 PAM？

## 冻结数据

- 主终点：45 个化合物的 human M4 CRE-Luc PAM EC50（ACh EC20）。
- 次终点：其中 14 个的 human M4 GloSensor cAMP PAM EC50（ACh EC20）。
- 选择性终点：其中 8 个的 human M2 CRE-Luc EC50；`>10000 nM` 保留为右删失下界。
- 结构恢复：专利系统命名经 OPSIN 转换，RDKit 标准化；结构图和报告的 `(M+H)+` 用于独立 QC。
- 在任何预测前冻结结构、标签、去重结果及本文件。

## 模型与基线（禁止外测调参）

所有模型配置沿用 `run_pacer_external_2026.py` 的冻结配置：

1. `Absolute-QSAR`：Morgan/描述符 + 固定 LightGBM 回归。
2. `PACER-Centered`：训练系列内中心化目标，恢复到训练总体均值。
3. `1NN`：训练集中 ECFP4 最近邻活性。
4. `Similarity-kNN`：ECFP4 相似度三次方加权。
5. `Constant-median`：训练集活性中位数，仅检验 MAE，不参与排序比较。

主方法预先指定为 `Absolute-QSAR`。本轮不根据专利结果选择模型、不做 few-shot、不把专利分子加入训练。

## 主指标和支持判据

主指标：45 个 CRE-Luc 化合物上的 Spearman ρ。  
次指标：MAE、top-quartile recall、≤50 nM ROC-AUC、pairwise concordance。  
不确定性：分子级 bootstrap 50,000 次；同时报告相对 1NN 和 similarity-kNN 的配对差值置信区间。

预注册支持条件：

- `Absolute-QSAR` Spearman 的 95% bootstrap CI 下界 > 0；并且
- 相对 `Similarity-kNN` 的 Spearman 差值 95% CI 下界 > 0。

任何一项不满足，均不得宣称独立外部 SOTA。可以报告失败及其可诊断边界。

## 生物学一致性审计

- 在 14 个双测定化合物上报告 CRE-Luc 与 GloSensor 实测 pEC50 的 Spearman，作为跨平台可迁移性的经验上限线索，而不是模型成绩。
- 报告模型在两个平台上的排序方向是否一致。
- M2 数据只用于选择性风险审计；右删失值不作为精确回归标签。

## 禁止事项

- 不用专利标签调参、挑 seed、挑特征或决定汇报哪个模型。
- 不把同一专利内部随机拆分结果称为外部验证。
- 不把 docking/MD 结合能解释为 PAM 功能效力。
- 没有湿实验前，不把任何预测分子称为“已确认 PAM”。
