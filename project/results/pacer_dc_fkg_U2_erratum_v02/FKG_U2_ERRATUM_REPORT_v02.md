# PACER-FKG 原 kernel_u2 标量筛选勘误 v02

**结论性质：新增估计目标的描述性复核，旧门禁已撤销；不授权新门禁、药理学判断或训练。**

## 确认的问题

- 历史 `kernel_u2` 是 U-statistic 的平方核范数估计，不是平方根；平方本身合法。
- 历史区域与 stable-core 使用独立的 per-call median 带宽，原样相减并据此设置门禁，缺乏共同核尺度。
- 旧图传播扩散的是非负 `sqrt(max(kernel_u2,0))`，并非有符号 RKHS 差分向量。
- 相邻 MD 帧并非独立；旧无偏表述和窗口 Spearman 不应作为统计显著性。

## v02 修复的估计目标

沿用冻结 G2-C 的 R2 通道预处理、共享节点 RBF 带宽、256D RFF 与种子。
先分别构造每个残基的有符号四上下文核均值差，再对区域内残基取平均；
该向量的范数在同一个 RFF 空间中与稳定核心范数比较。
**这是共享节点核下的 bag-of-residue 估计目标，不等于旧版 joint-region U-statistic。**
不根据结果设置或调节新门禁。窗口/块区间仅作敏感性分析。

## 旧筛选结果与修正描述量对照

| 区域 | 轴 | 历史筛选（撤销） | R2 中位范数差 | R3 中位范数差 | 新 pooled 余弦 |
|---|---|---|---:|---:|---:|
| compound110_extension | synergy_interaction | 曾通过 | -0.01092 | -0.03701 | -0.109 |
| cooperativity_mutagenesis | synergy_interaction | 未通过 | -0.03015 | -0.01155 | -0.143 |
| distal_control | synergy_interaction | 未通过 | 0.00581 | -0.00986 | -0.690 |
| intracellular_microswitches | synergy_interaction | 未通过 | -0.04085 | -0.03613 | -0.124 |
| orthosteric_activation_core | synergy_interaction | 未通过 | -0.06026 | -0.05669 | -0.013 |
| orthosteric_contact_union | synergy_interaction | 未通过 | -0.06444 | -0.05840 | 0.063 |
| pam_contact_consensus | synergy_interaction | 未通过 | -0.01063 | -0.00797 | -0.133 |
| pam_contact_union | synergy_interaction | 曾通过 | -0.01636 | -0.02611 | -0.134 |
| stable_core_control | synergy_interaction | 未通过 | 0.00000 | 0.00000 | 0.359 |
| compound110_extension | intrinsic_agonism | 未通过 | -0.00780 | 0.00501 | 0.321 |
| cooperativity_mutagenesis | intrinsic_agonism | 未通过 | -0.02266 | -0.01939 | 0.277 |
| distal_control | intrinsic_agonism | 未通过 | 0.01180 | -0.00722 | -0.774 |
| intracellular_microswitches | intrinsic_agonism | 未通过 | 0.02207 | -0.01600 | 0.053 |
| orthosteric_activation_core | intrinsic_agonism | 未通过 | -0.02291 | -0.03401 | 0.133 |
| orthosteric_contact_union | intrinsic_agonism | 未通过 | -0.02459 | -0.03764 | 0.100 |
| pam_contact_consensus | intrinsic_agonism | 未通过 | -0.01140 | -0.00793 | 0.118 |
| pam_contact_union | intrinsic_agonism | 曾通过 | -0.02682 | -0.01200 | 0.083 |
| stable_core_control | intrinsic_agonism | 未通过 | 0.00000 | 0.00000 | 0.134 |
| compound110_extension | conditional_pam_effect | 曾通过 | 0.01962 | 0.04554 | 0.176 |
| cooperativity_mutagenesis | conditional_pam_effect | 曾通过 | 0.04077 | 0.04442 | -0.022 |
| distal_control | conditional_pam_effect | 未通过 | 0.01346 | -0.00663 | 0.034 |
| intracellular_microswitches | conditional_pam_effect | 曾通过 | 0.04827 | 0.00222 | -0.205 |
| orthosteric_activation_core | conditional_pam_effect | 未通过 | 0.00387 | 0.00217 | 0.298 |
| orthosteric_contact_union | conditional_pam_effect | 未通过 | -0.00125 | -0.00060 | 0.300 |
| pam_contact_consensus | conditional_pam_effect | 未通过 | 0.05625 | 0.05650 | 0.262 |
| pam_contact_union | conditional_pam_effect | 曾通过 | 0.05483 | 0.03952 | 0.315 |
| stable_core_control | conditional_pam_effect | 未通过 | 0.00000 | 0.00000 | 0.393 |

## 校验与使用限制

- 输入 SHA256：40/40。
- G2-C 无图 27 项交叉回放：PASS。
- G2-B 的区域专用核与本次共享节点核不同，不能直接比较两者余弦大小。
- 原 FKG / G1 保留为 SUPERSEDED 历史数据；G2-B / G2-C 不因本次勘误而被改写。
- 长程 MD 尚未纳入；不能据此声称模型已验证或存在 PAM 功效。
