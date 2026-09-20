# PACER 顺序药化闭环验证预注册

冻结时间：2026-08-30，Generation-3 模型预测首次运行之前。

## 问题定义

模拟一轮真实 DMTA：历史 430 个 M4 PAM 数据构成先验，2026 论文 Table 1 的 18 个 5/6 系列实验值作为第一轮回传，随后预测同一项目下一代 29a-j。

主查询排除 29d/29e，因为其去同位素后的二维结构与 Table 1 的 5q/6q 相同。剩余 8 个分子为顺序查询；失活项仅作为 `pEC50 < 4.5` 的删失上界。

## 冻结方法

1. Historical-Absolute：仅历史 430 分子的绝对 QSAR；
2. Augmented-Absolute：历史 430 + 第一轮 18 分子的绝对 QSAR；
3. Table1-kNN：仅第一轮结构相似度加权局部基线；
4. Historical-DeltaSAR：历史数据训练 pairwise ΔpEC50，以第一轮 18 个实验锚点推断查询；
5. Augmented-DeltaSAR：在 pairwise 模型中加入 Table 1 作为一个新系列后再推断。

所有特征、LightGBM 参数、Morgan 半径和相似度权重沿用既有 PACER-FS，不用 Generation-3 标签选模型或参数。

## 评价

- 主指标：8 个查询的删失感知排序一致性；实现时报告保守 ceiling-Spearman 和 pairwise concordance；
- 补充：100 nM 强效分类 ROC-AUC、top-3 recall；
- 29d/29e 仅报告同位素归一化阳性对照误差，不并入主指标；
- 样本量仅 8，任何单点优势都不能称 SOTA。

声明边界：这是同论文、同项目的顺序药化模拟，不是独立实验室或前瞻湿实验验证。分析者已经读过论文标签，因此它是可审计的回顾性挑战，而非盲测。
