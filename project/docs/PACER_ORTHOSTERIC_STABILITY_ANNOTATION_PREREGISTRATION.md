# PACER Orthosteric-Stability Trajectory Annotation v0.1

日期：2026-09-01

## 问题

检验由静态PAM锚点学习、且训练时不含配体坐标的蛋白状态表示，是否对应第三PAM
MK-97轨迹中ACh结合姿势的动态稳定性。Wang 2022公开数据已独立证明PAM可降低
正构配体RMSD，因此ACh pose RMSD作为机制相关外部注释，而不是效力标签。

## ACh RMSD

- 每帧使用冻结21个M4R C-alpha对齐到run1第一帧；
- 计算ACh（resid 907）全部重原子相对run1第一帧的受体对齐pose RMSD；
- 同时报告相对各replica第一帧RMSD作敏感性分析；
- triplet训练从未读取ACh或MK-97坐标。

## 冻结评价

1. 各replica内Spearman：Consensus-PAM signature与`-ACh RMSD`；
2. 10000次各replica独立循环移位，保留时间自相关，检验六replica平均rho；
3. 对冻结leave-one-replica-out assignments计算cluster对ACh RMSD的eta-squared；
4. 比较AnchorGuidedTriplet与同输入PCA；三个triplet种子取平均，不选种子。

## 支持条件

- mean within-replica rho >0，至少4/6 replica同向，循环移位p<=0.05；
- triplet eta-squared比PCA平均提高>=0.05且至少4/6折提高。

通过只支持蛋白状态轴与正构配体稳定性耦合；不等于PAM效力、cooperativity数值或候选确认。
