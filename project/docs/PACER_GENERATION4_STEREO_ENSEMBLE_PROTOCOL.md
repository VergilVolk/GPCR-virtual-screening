# Generation-4 立体化学与 GaMD ensemble 压力测试

冻结时间：2026-08-30，33h-r 首次 docking 之前。

## 科学问题

四个非氘代 3-氟哌啶异构体 33h、33j、33l、33n 具有相同二维连接图，却有不同 hM4 PAM EC50。普通非手性 ECFP 必然给出相同表示。本测试判断 3D 构象与 M4 GaMD 受体 ensemble 是否提供可测的额外立体信息。

## 固定结构协议

- 受体：已经由公开 GaMD 轨迹复现的 10 个 M4 聚类代表构象；
- 口袋：各聚类共同模拟 MK-97 的质心，30 Å 立方盒；
- 配体：RDKit ETKDGv3 + MMFF，保留论文指定 (3S,4R)/(3R,4S)/(3S,4S)/(3R,4R)；
- 微状态：Molscrub pH 7.0；
- 对接：AutoDock Vina 1.2.7，exhaustiveness 8，每个 cluster 取 rank-1；
- ensemble 分数：文献复现定义 `BE_i = Vina_i + PMF_i`，主指标为 `-BEavg` 与 pEC50 的 Spearman；`-BEmin` 和单 cluster 00 为补充。

## 评价边界

- 主要集合仅四个非氘代立体异构体；N=4，只能作为机制压力测试；
- 氘代对和 33p-r 仅作一致性审计，不扩充有效样本量；
- 分数若不能恢复排序，结论是当前 docking/ensemble 不含足够功能效力信息，而不是“MD 无效”；
- 即使相关也不能证明 PAM 功能、协同性或探针依赖性。
