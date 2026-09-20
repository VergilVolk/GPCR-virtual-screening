# PACER Anchor-Guided Trajectory Triplet v0.3 预注册

日期：2026-09-01

## 假设

两个不同 PAM 静态复合物相对 ACh-only 结构定义的共同内部几何方向，可以作为
弱监督构造跨 replica triplet，使 snapshot embedding 按“共同耦合进度”组织，
而非按 replica 身份或单纯时间邻近组织。

## 严格边界

该标签来自静态结构签名，不是 PAM 效力标签。该实验检验的是第三 PAM 的结构方向
能否跨 replica 迁移，不检验 affinity、cooperativity 或 efficacy。

## 外层与输入

- 六折 leave-one-replica-out；所有缩放、PCA、triplet、encoder、KMeans 只拟合
  五条训练 replica；
- 输入为预注册 0.10 A 阈值下的 141 个共识内部距离，train-only 标准化后 PCA-32；
- held-out replica 的 signature 只用于最终评价。

## Triplet

- positive：其他训练 replica 中 signature 最接近的 20 帧里，PCA 距离最近者；
- hard negative：signature 差至少 0.50、但 PCA 距离最近的训练帧；
- cosine margin 0.20；另加时间邻近 triplet、重建与 variance floor；
- encoder `32 -> 64 -> 16`，三个冻结种子 17/29/43，120 epochs；
- 不增加直接 signature regression head，避免把评价退化成回归复制。

## Baseline 与评价

- 同一 141 特征的 train-only PCA-10 + KMeans-6；
- encoder embedding + KMeans-6；
- held-out eta-squared、相邻帧状态持续率、occupancy JSD、有效状态数、最大状态占比；
- 三种子 held-out assignment AMI。

## Go 条件

同时满足：Triplet eta-squared 比 PCA 平均提高 >=0.05 且至少 4/6 折提高；平均种子
AMI >=0.60；时间持续率不低于 PCA 0.05；有效状态数 >=3；最大状态占比 <=0.90。
不满足则冻结 v0.3，不在同一六条轨迹上调 margin/阈值追分。
