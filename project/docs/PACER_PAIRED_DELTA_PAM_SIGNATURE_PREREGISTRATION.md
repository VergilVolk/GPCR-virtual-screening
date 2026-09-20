# PACER Paired-Delta PAM Signature v0.4 预注册

日期：2026-09-01

## 纠正的科学问题

旧Consensus-PAM signature以ACh-only 7TRS为基线，而7TRQ/7TRP使用iperoxo，因而
混入正构探针差异。v0.4只学习相同iperoxo背景下添加PAM后的结构增量。

## 开发与外测

- 开发差分1：`7TRQ(VU0467154+iperoxo) - 7TRK(iperoxo-only)`；
- 开发差分2：`7TRP(LY2033298+iperoxo) - 7TRK(iperoxo-only)`；
- 独立外测：`7V68(LY2119620+iperoxo) - 7V69(iperoxo-only)`；
- 特异性压力测试：7V6A compound-110 allosteric agonist。

使用相同21个C-alpha/210内部距离。保留两个开发差分同向、绝对值均>=0.10 A的
距离，按两差分幅度一致性加权；令7TRK为0，两个开发PAM均值为1。0.05/0.20 A
为冻结敏感性阈值。

## 条件

1. 独立配对增量 `score(7V68)-score(7V69) > 0`，三个阈值一致；
2. `score(7V68) > score(7V6A)`，三个阈值一致，才支持PAM-specific；
3. 若仅条件1通过，称probe-matched PAM structural increment，不称PAM-specific。

该分数不含效力、协同性标签。
