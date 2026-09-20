# PACER 2026 独立时间外验证预注册

冻结时间：2026-08-30，首次运行模型预测之前。

## 数据边界

- 来源：Gregro 等，ACS Chemical Neuroscience 2026，DOI 10.1021/acschemneuro.5c00963。
- 主数据：Table 1 中可由 Scheme 1 唯一恢复的 18 个 5/6 系列分子。
- 终点：hM4/Gqi5-CHO 钙流，ACh EC20 条件下 PAM EC50。
- 仅使用精确 EC50；模糊或无法唯一重建的结构不进入本轮。
- 与现有 430 分子精确结构重叠者在预测前剔除。

## 冻结模型与评价

训练集固定为 `m4_pam_v1/potency_molecules.csv`，外部标签不参与训练、特征选择或超参数选择。

零样本主要比较：

1. Absolute-QSAR：原始 pEC50 回归；
2. PACER-Centered：每个训练来源内中心化后的相对 SAR 回归；
3. 1-NN 与相似度加权 kNN：非参数基线。

少样本次要比较：先由零样本 Absolute-QSAR 的预测跨度选择 3 个锚点，再比较 offset-calibrated Absolute、Centered 与固定 PACER-FS DeltaSAR-Centered。锚点身份不由真实标签选择。

主要指标为全系列 Spearman；补充指标为 MAE、9 个 5→6 匹配分子对的方向准确率、强效分子富集。Bootstrap 仅描述不确定性，不改变模型。

## 成功判据与声明边界

- 外部支持：PACER 方法的 Spearman 高于 Absolute-QSAR，且配对 bootstrap 改善区间不跨 0。
- 若区间跨 0，只报告“方向性/探索性支持”或“未验证”，不得称 SOTA。
- 该验证证明的是 ACh 探针下同系列 PAM 功能效力排序，不证明新生成分子一定是 PAM，也不代替湿实验。
- 任何看过外部标签后的规则改动必须建立新版本，不能回写本轮结果。
