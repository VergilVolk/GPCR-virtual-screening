# 13 靶 LOSO 迁移到 Science-2026 checkpoint：首轮结果

日期：2026-09-29 深夜。执行：`migrate_13target_to_science2026.ps1`（提交 11c3e5c7）。
协议与旧 checkpoint 完全一致，仅替换权重（`litpcba_identity_90.pt`，SHA 已核验）
与配套 embedding：GPCR 28,519 分子 + 40 GaMD 口袋用新 ckpt 重编码
（`science2026_ensemble_embeddings.npz`），外部 9 靶复用 pilot 的
`science2026_90.npz`。六跑 = 2 变体（preserve 0.2 / 0.0）× 3 种子，
ECFP4 在完全相同切分上重算。成簇 bootstrap CI 由
`summary_preserve*_3seed.json` 补充（生成中，本文数字先行为三种子原始均值）。

## 宏平均（preserve0，3 种子均值）

| 指标 | raw Science-2026 | LoRA 微调 | ECFP4 | 裁决 |
|---|---:|---:|---:|---|
| ROC-AUC | 0.5442 | **0.6097** | 0.5681 | ✅ 超 raw +0.066、超 ECFP +0.042，6/6 跑方向一致 |
| PR-AUC | 0.1780 | **0.2158** | 0.2014 | ✅ 超 raw +0.038、超 ECFP +0.014 |
| BEDROC20 | 0.1875 | 0.2305 | **0.2387** | ⚠️ 超 raw，略低于 ECFP（−0.008） |
| EF1% | 1.7878 | 1.3517 | **2.3768** | ❌ 微调后前端反而变差 |
| EF5% | 1.1269 | 1.3015 | **1.7080** | ⚠️ 超 raw，低于 ECFP |

随机标签对照 6/6 落在 0.49–0.51，增益来自监督信号。preserve0 稳定优于
preserve0.2（与旧 ckpt 消融方向一致）。

## 逐靶点（preserve0，3 种子均值，ROC-AUC）

改善 9/13：KAT2A +0.17、B2AR +0.16、VDR +0.16、FEN1 +0.15、M2R +0.13、
MTORC1 +0.11、GBA +0.11、MAPK1 +0.04、ESR1m +0.04。
恶化 4/13：**M4R 0.540→0.473（跌破 0.5，EF1% 1.56→0.64）**、PKM2 −0.11、
CCR2 −0.02、ALDH1 −0.00。

## 成簇 bootstrap 95% CI（preserve0，3 种子，已落地）

tuned − official（raw Science-2026）：

- ROC-AUC `[+0.0522, +0.0668, +0.0817]` ✅
- PR-AUC `[+0.0247, +0.0381, +0.0526]` ✅
- BEDROC20 `[+0.0200, +0.0451, +0.0728]` ✅
- EF1% `[−1.044, −0.203, +0.390]` 跨零；EF5% `[−0.122, +0.240, +0.586]` 跨零

tuned − ECFP4：

- ROC-AUC `[+0.0187, +0.0433, +0.0682]` ✅ **CI 为正**
- PR-AUC `[−0.007, +0.015, +0.036]` 跨零
- BEDROC20 `[−0.043, −0.006, +0.030]` 跨零
- EF1% `[−2.128, −1.144, −0.383]` ❌ **显著劣于 ECFP**
- EF5% `[−0.809, −0.376, +0.096]` 跨零

逐靶点改善计数（13 靶中）：ROC 9、PR 10、BEDROC 10、EF1% 4、EF5% 6。

**CI 层面的最终表述**：迁移后，相对 raw 的 ROC/PR/BEDROC 增益全部保持 CI 为正；
相对 ECFP 的优势从旧权重的四项（ROC/PR/BEDROC/EF5%）收缩为一项（ROC），
且 EF1% 显著更差。旧权重"全面优于 ECFP"的结论**不**随权重迁移。

## 与旧 checkpoint（2023 PDB-general）同协议对照

| | 旧 ckpt（2023） | 新 ckpt（Science-2026） |
|---|---|---|
| raw 宏 ROC | 0.562 | 0.5442 |
| 微调宏 ROC | 0.623 / 0.637(无保持) | 0.606 / 0.610(无保持) |
| 微调 vs ECFP（ROC） | +0.055/+0.069，CI 为正 | +0.042（CI 待出） |
| 微调 vs ECFP（EF1%/EF5%） | EF5% CI 为正、EF1% 跨零 | **两项均低于 ECFP** |
| 留出 M4 EF1% | 0.74 → **2.08**（改善） | 1.56 → **0.64**（恶化） |

## 结论（当前证据下的三条）

1. **宏增益可迁移**：参数高效双端适配在 Science-2026 权重上仍稳定提升宏
   ROC/PR（+0.066/+0.038，三种子方向一致），证明该方法不依赖旧权重的特定空间。
2. **M4/GPCR 优势不可迁移**：新权重的增益集中在 LIT-PCBA 型靶点；对留出 M4，
   13 靶适配在旧权重上是改善（EF1% 2.08）、在新权重上是恶化（0.64）。
   前端富集（EF1%）在新权重上全面弱于 ECFP。
3. **部署路由（证据支持）**：M4/GPCR 结合优先级线继续使用旧 2023 checkpoint
   + 其 13 靶适配权重；Science-2026 权重及其适配用于非 GPCR 靶点或作为
   跨家族对照。两条线的结果不可互换引用。

## 边界

- 全部为回顾性 LOSO，无 prospective 主张；
- 不改变"PAM 功能归四上下文"的模块边界；
- CCR2/PKM2 的失败模式与旧权重一致，非新权重特有；
- 本机为笔记本环境，长任务曾三次因系统休眠中断（看门狗已重启护航）。

## 复现

- 运行器：`scripts/migrate_13target_to_science2026.ps1`
- 输出：`results/drugclip_science2026/migration_13target_v01/`
- 逐种子 JSON 含完整 per-target 与 training audit（scaffold purge 记录）。
