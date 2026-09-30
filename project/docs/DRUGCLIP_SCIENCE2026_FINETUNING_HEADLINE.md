# 头条主张：微调权重 vs 官方 Science-2026 checkpoint（13 靶严格 LOSO）

日期：2026-09-30。数据：`results/drugclip_science2026/migration_13target_v01/`
（6 个 LOSO 运行 + ECFP 基线 + 3 种子成簇 bootstrap，提交 0fbb9754/cd6248c4）。
可部署权重：`science2026_13target_finetuned.seed{20260925,26,27}.projection.pt`
（双侧全投影，13 靶 37,524 对全数据训练）。

## 主张（可写进报告）

在 13 靶整靶留出 + scaffold purge 协议下，我们的 10,240 参数双侧 rank-8
LoRA 微调（BCE + target retrieval，no-preservation）相对**官方 Science-2026
checkpoint 冻结基线**：

| 指标（宏） | 官方 raw | 我们微调 | Δ 95% CI | 判定 |
|---|---:|---:|---|---|
| ROC-AUC | 0.5442 | **0.6097** | **[+0.0522, +0.0668]** | ✅ 显著优于官方 |
| PR-AUC | — | — | **[+0.0247, +0.0381]** | ✅ |
| BEDROC20 | — | — | **[+0.0200, +0.0451]** | ✅ |
| EF1% | — | — | [−1.04, −0.20] | ❌ 官方略优 |
| 9 靶试点 (LOTO) | 0.5193 | **0.6080** | — | ✅ +0.089 |

随机标签控制 [+0.098, +0.115]（ROC）确认非退化。逐靶：ROC 9/13、
PR 10/13、BEDROC 10/13 改善。

## 必须保留的边界（不写即失实）

1. **vs ECFP4**：只有 ROC CI 为正 [+0.0187, +0.0433]；PR 跨零，
   BEDROC20/EF1%/EF5% 为负——ECFP 在早期富集上仍是更强的基线。
2. **EF1% vs 官方为负**：主张限定为"整体排序与检索质量显著改善"，
   不主张 Top-1% 富集改善。
3. **M4/GPCR 路由**：留出 M4 在新权重上退化（EF1% 1.56→0.64），
   M4/GPCR 结合线守旧 2023 权重（M4-LOSO EF1% 2.08 已冻结）；
   新微调权重用于 LIT-PCBA 型/非 GPCR 面板。
4. 微调回答的是结合筛选优先级，不是 PAM 功能（两阶段边界不变）。

## 支撑文件

- 全量结果与协议：`DRUGCLIP_SCIENCE2026_13TARGET_MIGRATION_RESULT.md`
- 15 靶全量官方 raw 基线：跑至 11/15（作官方参照，不再接 CGDA 链）
- CGDA 支线：9/30 团队决策关闭，反证见
  `DRUGCLIP_SCIENCE2026_CGDA_SUBSET_GATE_RESULT.md`
