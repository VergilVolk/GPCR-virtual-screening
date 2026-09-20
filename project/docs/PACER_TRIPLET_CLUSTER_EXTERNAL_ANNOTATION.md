# PACER Triplet-Cluster External Annotation Test v0.1

日期：2026-09-01

## 目的与时间顺序

ReplicaInvariantTriplet v0.2 在本检验前已经训练并冻结，训练未读取任何静态
PAM 锚点、PAM signature 或功能标签。本检验在看到 Consensus-PAM signature
结果后定义，因此属于**回顾性外部注释测试**，不是预注册独立确认。

## 问题

检验无标签 triplet 轨迹 cluster 是否比 PCA 或 tICA cluster 更集中地划分由
7TRQ/7TRP 相对 7TRS 定义的跨 chemotype 结构方向。

## 冻结评价

- 对每个 held-out replica 单独计算 cluster 对连续 signature projection 的
  eta-squared（cluster 解释的方差比例）；
- PCA-10 与 tICA-5 严格复算 v0.2 的 train-only、leave-one-replica-out 分配；
- Triplet 使用已冻结的三个种子 held-out assignments，不重训、不选种子；
- Triplet 的折值为三个种子 eta-squared 的平均；
- 1000 次每 replica 独立循环移位 signature，保留其时间自相关，获得总体平均
  eta-squared 的经验 null 与 p 值；
- 报告六折 Triplet-PCA 差值，但 N=6 不宣称 SOTA。

## 解释边界

Triplet 若显著，只说明它捕获了与静态 PAM 共识方向相关的构象分区；不说明状态
导致 PAM 效力，也不证明可跨体系预测 PAM。若不优于 PCA，则 triplet 路线在现有
单体系数据上没有可证实增量，应冻结而不是调参。
