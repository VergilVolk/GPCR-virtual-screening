# PACER Unsupervised State Pilot v0.1 预注册

日期：2026-08-30

## 目的与边界

本 pilot 只检验：同一 MK-97–M4R 体系的六条 500 ns GaMD replica 是否支持可重复、非塌缩且具有时间连续性的 snapshot 状态表示。它不使用 PAM 效力标签，不判断 PAM/inactive，不预测候选效力，也不构成 StateMIL 或 SOTA 证据。

## 冻结数据

- Figshare 33283491 六条公开 1 ns-stride GaMD trajectory；
- 每条 500 帧，共 3000 帧；
- 使用作者口袋定义的 239 个重原子，经质量加权刚体对齐后的 717 维坐标；
- replica 身份只用于评估与构造时间邻近正样本，不作为模型输入。

## 冻结方法

全部方法固定 10 个 k-means 状态，比较：

1. `PCA-10`：坐标标准化后 PCA 10维；
2. `tICA-5`：先 PCA 32维，再按每条 replica 内 5 ns lag 求 time-lagged components；
3. `TemporalTriplet-16`：PCA-32 输入，64-32-16 MLP；同 replica 的 1–5 ns 邻帧为 positive，其他 replica 随机帧为 negative，cosine margin 0.2；不使用功能标签。

TemporalTriplet 固定 120 epochs、AdamW 1e-3、weight decay 1e-4、三个种子 17/29/43。结果出来后不在同一数据上调整 margin、维度或 epoch 追分。

## 冻结评估

- `replica occupancy JSD`：六条 replica 两两状态占比 Jensen–Shannon divergence，越低越可重复；
- `all-replica state coverage`：在六条 replica 均至少出现 1% 的状态数，越高越好；
- `temporal persistence`：相邻 1 ns 帧保持同状态的比例，仅作时间连续性指标；
- `effective state count`：`exp(entropy(global occupancy))`，用于识别塌缩；
- `max state fraction`：最大状态占比；
- `seed AMI`：TemporalTriplet 三个种子聚类的一致性。

## Go/No-Go

技术 Go 必须同时满足：

1. TemporalTriplet 三种子 mean AMI >= 0.60；
2. effective state count >= 3；
3. max state fraction <= 0.90；
4. 至少 3 个状态在全部六条 replica 中占比 >=1%；
5. mean replica JSD 不高于 PCA-10；
6. temporal persistence 不低于 PCA-10。

即使全部通过，也只晋级为“无标签动态表示可行”；功能 StateMIL 继续保持 No-Go，直到获得匹配的 ACh-only、PAM-only、inactive/weak-PAM、多 chemotype 对照轨迹。

