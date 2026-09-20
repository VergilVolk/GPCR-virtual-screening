# PACER Consensus-PAM Trajectory Signature v0.1 预注册

日期：2026-09-01

## 问题

检验两个不同化学型 PAM 复合物（7TRQ/VU0467154、7TRP/LY2033298）相对
无 PAM 的 ACh 复合物（7TRS）是否定义同向的 M4R 口袋内部几何变化，以及第三个
PAM（MK-97）六条公开 GaMD 轨迹是否复现该方向。

本实验只检验跨 PAM 的结构共识，不把结构接近解释为 PAM 效力、协同性或功能验证。

## 冻结特征与映射

- 使用 21 个预先指定的 M4R 口袋/耦合残基的 C-alpha 原子；
- public GaMD 拓扑在受体 ICL3 缺口前后使用两段编号。映射在读取结果前按完整
  受体序列同一性冻结：public 687--788 段为 `static = public - 598`，public
  852--875 段为 `static = public - 436`；程序同时检查每一对映射的残基名称；
- 特征为 21 个 C-alpha 的全部 210 个成对内部距离，不依赖整体平移旋转；
- 只保留在 `7TRQ-7TRS` 与 `7TRP-7TRS` 中方向相同、且两者变化绝对值均
  大于等于 0.10 A 的距离；0/0.05/0.20 A 仅作预注册敏感性分析，不替换主结果。

## 冻结 PAM 共识分数

令 `S` 为 7TRS 距离向量，`Q/P` 为两个 PAM 锚点，`d=(Q-S+P-S)/2`。
每一维按两个 PAM 变化幅度的一致性加权：
`w=min(|Q-S|,|P-S|)/max(|Q-S|,|P-S|)`。

轨迹帧的 signature projection 为其 `X-S` 在加权共识方向 `d` 上的投影，
并归一化使 ACh-only 锚点为 0、两个 PAM 锚点均值为 1。正值仅表示沿共同 PAM
结构方向移动。

## 聚类与外层检验

- 聚类不读取 signature projection；
- 全部 3000 帧进行一次不读取 signature 分数的无监督标准化、PCA（解释至少
  90% 方差且最多 10 维）与 KMeans，K 固定为 6，用于定义可比较的全局 cluster；
- 另做六折 leave-one-replica-out 聚类，held-out 帧只投影/分配；用 AMI 比较其
  与全局 cluster 的 held-out 划分，作为聚类稳定性而非生物学性能指标；
- 报告每条 replica 的分数分布，以及跨全部六条 replica、每条占比至少 1% 的
  全局 cluster 是否富集正 signature。

## Go 条件

同时满足才支持“存在第三 PAM 可采样的跨 chemotype 共同结构方向”：

1. 主阈值下至少 10 个同向距离特征；
2. 六条 replica 的 median projection 均大于 0；
3. 至少一个 cluster 在六条 replica 中均占比 >=1%，且 median projection >0；
4. 结论在 0.05 和 0.20 A 两个敏感性阈值方向一致。

即使 Go，也不允许称该 cluster 为已验证 PAM 功能状态；Triplet/StateMIL 的功能
训练仍要求不同 PAM、ACh-only、弱 PAM/无效物的匹配动态轨迹。
