# PACER External Static PAM-Signature Validation v0.1

日期：2026-09-01

## 时间顺序

本协议在下载并查看外部结构坐标前，根据RCSB条目配体注释冻结。外部结构均未参与
7TRQ/7TRP/7TRS的Consensus-PAM signature构造或v0.3 triplet训练。

## 独立结构

- 7V68：iperoxo + LY2119620 PAM，独立第三PAM正对照；
- 7TRK：iperoxo-only，与7V68/7TRQ/7TRP共享正构探针的关键无PAM负对照；
- 5DSG：tiotropium-bound inactive M4，远端构象负对照。

## 冻结评分

沿用Consensus-PAM v0.1的21个C-alpha、210个成对内部距离、同向筛选、幅度一致性
权重与归一化公式；主阈值0.10 A，0.05/0.20 A为敏感性分析。不得重新选择残基、
旋转方向或根据结果改阈值。

## 支持条件

主阈值下必须同时满足：

1. `score(7V68) > 0`；
2. `score(7V68) > score(7TRK)`，排除仅由iperoxo/active-state造成的方向；
3. `score(7V68) > score(5DSG)`；
4. 三项顺序在0.05与0.20 A阈值保持。

通过仅支持“冻结结构签名对独立PAM结构有方向性外部效度”；不支持PAM效力预测，
也不将单个静态结构当作动态功能状态。
