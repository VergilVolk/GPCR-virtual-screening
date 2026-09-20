# PACER Matched Static PAM-Specificity Extension v0.2

日期：2026-09-01

## 独立同研究结构

- 7V68：iperoxo + LY2119620，PAM；
- 7V69：iperoxo-only，同研究、同构建体、同正构探针的无PAM对照；
- 7V6A：compound-110 allosteric agonist，不是PAM，用于区分“PAM轴”与一般
  变构占位/激活轴。

三个PDB身份根据论文Data availability冻结；7V69/7V6A坐标尚未读取。沿用v0.1
全部141距离与评分，不允许重拟合。

## 判别

- `7V68 > 7V69` 且三个阈值一致：支持排除同研究active/probe构建差异；
- `7V68 > 7V6A` 且三个阈值一致：支持PAM结构特异性；
- 若7V6A同样高，则只能称allosteric/active structural axis，不得称PAM-specific。

静态特异性即使通过，也不等于效力或协同性预测。
