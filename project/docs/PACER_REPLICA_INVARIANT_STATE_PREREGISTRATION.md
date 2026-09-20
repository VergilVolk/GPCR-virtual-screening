# PACER Replica-Invariant StateMetric v0.2 预注册

日期：2026-09-01

## 科学问题

检验同一MK-97–M4R体系的六条GaMD replica中，是否存在可由无标签表示学习复现的共同动态状态。该实验不含不同PAM、inactive、ACh-only或功能效力标签，因此只能验证state representation，不识别PAM功能状态。

## 外层评价

完整leave-one-replica-out：每折只用5条replica拟合标准化、PCA、encoder与10个prototypes；第6条replica只用于最终投影和评价。

## v0.2表示学习

- 输入：训练折拟合的PCA-32口袋坐标；
- encoder：`32 -> 64 -> 16`，单位球归一化；decoder：`16 -> 32`；
- temporal positive：同replica相距1–5 ns；
- cross-replica positive：冻结PCA空间中其他训练replica的最近帧；
- negative：同replica相距至少50 ns且PCA距离最远帧；
- cosine triplet margin 0.2，temporal与cross-replica两项等权；
- reconstruction权重0.1，replica均值方差约束0.1，variance-floor 0.1；
- AdamW 1e-3、weight decay 1e-4、120 epochs、种子17/29/43；
- KMeans k=10只在训练replica embedding拟合，测试帧按最近prototype分配。

## Baseline

1. train-only PCA-10 + KMeans-10；
2. train-only tICA-5（lag 5 ns）+ KMeans-10；
3. ReplicaInvariantTriplet-16 + KMeans-10。

## 指标与Go条件

- heldout occupancy JSD：留出replica与训练池状态占比差异，越低越好；
- heldout temporal persistence：留出轨迹相邻1 ns帧同状态比例；
- common-state coverage：留出replica及每条训练replica均>=1%的状态数；
- heldout effective state count与max-state fraction防塌缩；
- 三种子heldout assignment AMI检验训练稳定性。

技术Go要求同时满足：mean seed AMI>=0.60；effective states>=3；max-state<=0.90；common coverage>=3；heldout JSD不高于PCA；temporal persistence不低于PCA。

即便通过，functional StateMIL仍为No-Go；只有补齐不同PAM chemotype、ACh-only、inactive/weak-PAM及matched replicas后，才可用轨迹级功能标签微调prototype attention。

