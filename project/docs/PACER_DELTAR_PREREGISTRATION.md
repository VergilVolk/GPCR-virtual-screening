# PACER-DeltaR v0.1 预注册

冻结日期：2026-08-30

## 问题

M4 PAM 功能效力数据包含明显的药化系列/文献基线漂移。绝对QSAR把系列偏移与系列内SAR混合；普通DeltaSAR又直接学习绝对效力差，容易重复基线模型已经学到的平滑信号。

## 冻结假设

PACER-DeltaR学习的是跨系列基线误差的局部变化：

`Delta residual(a,q) = [y(q)-base(q)] - [y(a)-base(a)]`

对每个外层留出系列：

1. 只用其他系列训练absolute LightGBM；
2. 对外层训练集再做inner leave-one-series-out，得到无泄漏base residual；
3. 仅用同系列分子对训练Delta-residual模型；
4. 新系列按absolute base预测分位选择3个predicted-span锚点；
5. 查询分子由`base(q) + anchor residual + predicted Delta residual(a,q)`得到；
6. 强制`Delta(a,q)=-Delta(q,a)`，并按ECFP4 Tanimoto三次方加权，固定shrinkage 1.5。

## 固定主终点

- 数据：`potency_molecules.csv`；
- 外层：`source_component` leave-one-series-out，系列大小至少12；
- 锚点：3个label-free predicted-span anchors；
- 主指标：equal-series macro Spearman；
- 次指标：worst-series Spearman、positive-series count、macro MAE；
- baseline：Absolute、Offset、PACER-FS DeltaSARHybrid；
- 不使用结构对接、候选标签或外层系列查询效力调参。

## 晋级规则

只有同时满足以下条件才进入主方法：

1. 相对同轮Absolute的series-bootstrap 95% CI下界大于0；
2. 相对冻结PACER-FS `0.263`不是明显退化；
3. worst-series不因平均值提升而显著恶化；
4. 时间外推、activity-cliff结果单独报告，不与主终点混合。

若失败，结果作为结构化负证据保留，不重新定义主终点。
