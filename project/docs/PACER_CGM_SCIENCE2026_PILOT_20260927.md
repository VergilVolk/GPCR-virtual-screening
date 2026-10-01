# PACER-CGM × Science-2026：开发集审计

日期：2026-09-27。

## 方法

以经官方 SHA256 核验的 Science-2026 `litpcba_identity_90.pt` 为冻结 backbone。PACER-CGM 在原始 pocket--molecule cosine score 上增加低秩双线性关系残差；外层按整靶点留出，并从训练集清除留出靶点的全部 Murcko scaffold。三 seed 集成，并设置逐靶点 shuffled-label 对照。

## 九靶点开发集结果

| 方法 | ROC-AUC | PR-AUC | BEDROC20 | EF1% | EF5% |
|---|---:|---:|---:|---:|---:|
| Science-2026 raw | 0.5193 | 0.1978 | 0.2188 | 1.708 | 1.224 |
| PACER-CGM，三专家 | 0.5501 | 0.2248 | 0.2521 | 2.074 | 1.485 |
| PACER-CGM，单专家 | 0.5555 | 0.2230 | 0.2524 | 2.125 | 1.489 |
| Shuffled-label control | 0.5082 | 0.1922 | 0.1984 | 1.665 | 1.055 |
| ECFP4 ligand-only | **0.6075** | **0.2479** | **0.3019** | 1.589 | **1.695** |
| 固定 50:50 ECFP/PACER rank fusion | 0.6069 | 0.2489 | 0.2966 | **2.226** | 1.722 |

PACER-CGM 单专家相对 raw 的靶点级 bootstrap 95% CI：ROC-AUC `[+0.0159,+0.0568]`，PR-AUC `[+0.0098,+0.0458]`，BEDROC20 `[+0.0144,+0.0544]`，EF5% `[+0.0890,+0.4854]`；EF1% 区间跨零。

## 结论与边界

1. PACER 的结构关系适配在严格 target-LOTO/scaffold-purge 下具有可重复信号，且随机标签对照不能复现。
2. 三专家没有优于单专家，故不保留无证据的复杂 gating；当前支持的是单个共享低秩双线性关系。
3. ECFP ligand-only 明显更强，说明该九靶点开发集存在很强的配体/数据来源信号。PACER 尚未超过强基线，不能称 SOTA。
4. 固定融合改善前端排序，但相对 ECFP 的 bootstrap 区间跨零，不能称稳定提升。
5. 该数据集为约千分子/靶点的高活性比例开发子集，不代表正式虚拟筛选。下一闸门必须是官方完整 LIT-PCBA/DUD-E 与真实低阳性率评价。

本结果只支持：**Science-2026 backbone 上存在可迁移的结构关系适配信号，但它尚不足以取代强 ligand-only baseline。**

## 可部署 triplet adapter（2026-09-28）

已将现有排序目标明确冻结为 pocket-anchored hard triplet：口袋 embedding 为
anchor，同靶点实验 active 为 positive，原始 DrugCLIP 得分最高的 inactive 为 hard
negative；同时保留 target-balanced BCE 与有界残差约束。该目标仍属于结合排序，不是
PAM 功能效力学习。

在完全相同的 9 靶点 target-LOTO、held-target scaffold purge 协议下复跑，结果与原单专家
实验一致：AUROC `0.5193 -> 0.5555`，PR-AUC `0.1978 -> 0.2230`，
BEDROC20 `0.2188 -> 0.2524`，EF1% `1.708 -> 2.125`，EF5% `1.224 -> 1.489`。
AUROC、PR-AUC、BEDROC20 和 EF5% 的 target-bootstrap 95% CI 下界均大于零；EF1%
仍跨零。shuffled-label control AUROC 为 `0.4989`。

已用全部开发靶点训练并导出三个 seed checkpoint：

- `project/results/drugclip_science2026/litpcba_external_v01/deployable_triplet/pacer_cgm_triplet_seed20260925.pt`
- `project/results/drugclip_science2026/litpcba_external_v01/deployable_triplet/pacer_cgm_triplet_seed20260926.pt`
- `project/results/drugclip_science2026/litpcba_external_v01/deployable_triplet/pacer_cgm_triplet_seed20260927.pt`

三个 checkpoint 均已通过 strict state-dict reload 与 finite-value 测试。它们可以作为
PACER binding-ranking 模块的阶段性交付，但尚不能称为 SOTA，也不能把输出解释为 M4 PAM。

## 强版本：双侧投影 hard-triplet LoRA（2026-09-28）

将此前13靶点最强路线完整迁移到 Science-2026 权重：同时适配分子与口袋 projection，
使用 target-balanced BCE、口袋锚定 hard-triplet ranking 和跨靶点 retrieval，不使用在旧消融中
被证明限制迁移的 representation-preservation 项。仍采用严格 target-LOTO、held-target
scaffold purge、三随机种子和 shuffled-label control。

| 方法 | ROC-AUC | PR-AUC | BEDROC20 | EF1% | EF5% |
|---|---:|---:|---:|---:|---:|
| Science-2026 raw | 0.5193 | 0.1978 | 0.2188 | 1.708 | 1.224 |
| PACER-CGM relation adapter | 0.5555 | 0.2230 | 0.2524 | 2.125 | 1.489 |
| ECFP4 ligand-only | 0.6075 | 0.2479 | **0.3019** | 1.589 | 1.695 |
| **PACER dual-projection triplet LoRA** | **0.6080** | **0.2526** | 0.2959 | **2.184** | **1.740** |

相对 raw DrugCLIP，AUROC、PR-AUC、BEDROC20、EF5%的 target-bootstrap 95% CI 下界均
大于零；EF1%仍跨零。相对 ECFP，PACER 在四项指标上数值更高，但 BEDROC20 低 0.006；
配对 target-bootstrap 显示五项指标的95%区间均跨零。因此结论是“在最新权重上基本追平
强 ligand-only baseline，并在四项指标数值领先”，不是全面显著超越或 SOTA。完整比较见
`pacer_projection_vs_ecfp_bootstrap.json`。

结果文件：

- `project/results/drugclip_science2026/litpcba_external_v01/pacer_projection_triplet_no_preserve.json`
- `project/results/drugclip_science2026/litpcba_external_v01/deployable_projection_triplet/`

三个约2.63 MB的物化双侧投影 checkpoint 已通过 finite-value 与结构审计，可直接替换原始
DrugCLIP 的 molecule/pocket projection 用于后续 M4 binding-ranking。
