# DrugCLIP 大样本参数高效微调：结果与证据边界

更新日期：2026-09-26

## 1. 本阶段实际做了什么

本阶段直接微调 DrugCLIP 内部的分子与口袋投影层，不使用 ECFP 外接分类器冒充 DrugCLIP 微调。两个 Uni-Mol encoder 保持冻结，在 `mol_project` 和 `pocket_project` 的末层加入双侧 rank-8 LoRA，共训练 10,240 个参数。

训练数据包含 28,618 个 GPCR active/decoy 配对、28,483 个可编码分子和 B2AR、CCR2、M2R、M4R 四个别构口袋。目标函数包含：

1. 靶点与标签平衡的 active/decoy BCE；
2. 活性分子的靶点检索交叉熵；
3. 同分子毒蕈碱亚型 triplet；
4. 相对官方 DrugCLIP 的表示保持与分数蒸馏。

所有主结果使用全局 Murcko scaffold 五折 OOF；每折还会删除与测试骨架重叠的辅助 triplet。

## 2. 域内五折结果

| 方法 | 宏 ROC-AUC | 宏 PR-AUC | BEDROC20 | EF1% | EF5% | 靶点检索宏 Recall@1 |
|---|---:|---:|---:|---:|---:|---:|
| 官方 DrugCLIP | 0.629 | 0.195 | 0.245 | 4.31 | 2.39 | 0.413 |
| BCE-only LoRA | 0.956 | 0.847 | 0.895 | 10.89 | 9.99 | 0.799 |
| BCE + retrieval LoRA | 0.951 | 0.864 | 0.911 | 10.83 | 10.32 | **0.882** |
| 完整多任务 LoRA | 0.951 | **0.864** | **0.912** | 10.83 | **10.32** | 0.881 |
| 全随机标签对照 | 0.493 | 0.116 | 0.128 | 1.98 | 1.22 | 0.356 |

真标签模型与随机标签模型差异明显，说明主提升依赖监督信号，而不是训练代码自动产生的高分。

但 ECFP4-logistic 在同一域内五折上的宏 ROC-AUC/PR-AUC 为 0.986/0.979，仍高于 DrugCLIP LoRA。这一数据集存在强化学系列信号，不能据此声称通用 SOTA。

## 3. triplet 分支的消融结论

随机翻转 triplet 正负关系后，宏 PR-AUC 为 0.864、BEDROC20 为 0.912，和真实 triplet 几乎相同；检索宏 Recall@1 甚至为 0.887，略高于真实 triplet 的 0.881。

因此，当前 muscarinic triplet 辅助任务没有提供可重复的额外信息。现阶段有证据支持的是 **BCE + target retrieval + preservation**，而不是 triplet。triplet 仍可作为后续机制型任务，但不能作为当前性能提升的创新主张。

## 4. 严格留一整靶点外推

每轮完全删除一个靶点的全部分子、标签及重叠 scaffold；被留出口袋不参与任何训练损失。三种子平均结果如下：

| 方法 | 宏 ROC-AUC | 宏 PR-AUC | BEDROC20 | EF1% | EF5% |
|---|---:|---:|---:|---:|---:|
| 官方 DrugCLIP | 0.629 | 0.195 | 0.245 | 4.31 | 2.39 |
| pooled ECFP4 外推 | **0.663** | 0.193 | 0.229 | 1.69 | 2.32 |
| DrugCLIP BCE + retrieval 外推 | 0.655 | **0.218** | **0.313** | **4.48** | **3.57** |

DrugCLIP 微调相对官方模型的 scaffold-bootstrap 95% CI：

- ROC-AUC：`+0.008` 至 `+0.047`；
- PR-AUC：`+0.002` 至 `+0.049`；
- BEDROC20：`+0.038` 至 `+0.103`；
- EF5%：`+0.445` 至 `+1.683`；
- EF1%：区间跨零，尚不稳定。

这是当前最有价值的算法结果：普通 ROC-AUC 上 ECFP 略高，但 DrugCLIP 在未见靶点的 PR-AUC、BEDROC 和早期富集上明显更强，更符合虚拟筛选需求。

单独留出 M4 时，官方 DrugCLIP 的 ROC-AUC/EF1% 为 0.667/0.74，三靶点迁移模型为 0.692/2.08。它支持“GPCR 域微调可改善未见 M4 的 active/decoy 早期富集”，不等同于 PAM 功能预测。

## 5. 未通过的外部功能评价

冻结权重在 M4 外部药理面板上没有稳定改善：

- Acadia 功能集合：完整多任务与官方 DrugCLIP AUC 均为 0.837；
- Monash PAM/inactive：官方 AUC 0.576，完整多任务降至 0.364；
- Suven、VU6025733 和专利有序效力数据：Spearman 均未形成可靠正相关。

因此，新模型仍是结合/筛选模型，不是 PAM 身份、协同性或效力模型。PAM 功能必须由独立药理标签或四上下文动态模块处理。

## 6. 九靶点 LIT-PCBA 冻结外部测试

为检查 GPCR 微调是否破坏 DrugCLIP 的通用检索能力，额外构建了 9 个从未用于训练或选模的 LIT-PCBA 靶点，共 8,906 个配对和 8,764 个唯一分子。

| 方法 | 宏 ROC-AUC | 宏 PR-AUC | BEDROC20 | EF1% | EF5% |
|---|---:|---:|---:|---:|---:|
| 官方 DrugCLIP | 0.532 | 0.203 | 0.239 | 1.44 | 1.51 |
| BCE + retrieval | 0.535 | 0.202 | 0.238 | 1.79 | 1.54 |
| 完整多任务 | **0.541** | **0.204** | **0.239** | 1.73 | 1.47 |

完整多任务相对官方模型的 ROC-AUC 差值 95% CI 为 `+0.0007` 至 `+0.0168`；其余指标区间跨零。结果说明微调没有造成明显的跨家族灾难性遗忘，并有很小的 ROC-AUC 改善，但没有证明通用早期富集提升。

## 7. 扩展到 13 靶点的统一留靶点验证

进一步将四个 GPCR 与九个 LIT-PCBA 靶点统一建模，共 37,524 个配对。每轮完整留出一个靶点，并从其余 12 个靶点删除所有重叠 scaffold；被留出口袋不进入损失。结果为三个独立 seed 的平均。

| 方法 | 宏 ROC-AUC | 宏 PR-AUC | BEDROC20 | EF1% | EF5% |
|---|---:|---:|---:|---:|---:|
| 官方 DrugCLIP | 0.562 | 0.201 | 0.241 | 2.33 | 1.78 |
| 随机标签 LoRA | 0.557 | 0.198 | 0.228 | 1.58 | 1.78 |
| ECFP4 logistic | 0.568 | 0.201 | 0.239 | 2.38 | 1.71 |
| DrugCLIP BCE + retrieval | **0.623** | **0.247** | **0.300** | 2.39 | **2.28** |

相对官方 DrugCLIP 的 target-stratified scaffold-bootstrap 95% CI：ROC-AUC `[+0.046,+0.079]`、PR-AUC `[+0.029,+0.059]`、BEDROC20 `[+0.030,+0.089]`、EF5% `[+0.071,+0.870]`。相对 ECFP4 的对应区间为 `[+0.030,+0.081]`、`[+0.019,+0.065]`、`[+0.022,+0.093]`、`[+0.102,+0.902]`。EF1% 的两个区间均跨零。

三种子宏 ROC-AUC 为 `0.625/0.617/0.624`，方向稳定；13 个靶点中，ROC-AUC、PR-AUC、BEDROC20 分别有 10、11、10 个靶点优于官方模型。CCR2 仍明显失败，说明方法并非对所有口袋都有效。

这组结果比四 GPCR 域内五折更接近可发表证据：它同时包含完整未见靶点、全局 scaffold purge、随机标签、强二维基线、三种子和成簇置信区间。但它仍是回顾性 13 靶点结果，不应写成通用 SOTA。

## 8. docking 与动态构象消融

十个 GaMD 代表构象的 population mean/LSE/adaptive aggregation没有稳定超过单一 cluster0。严格嵌套选择 docking 方法与融合权重后，融合结果也低于单独 DrugCLIP target-LOSO。固定 50:50 融合的表面高分属于事后探索，不能作为冻结结论。

## 9. 当前可成立的创新点

当前可以严谨表述为：

> 通过对 DrugCLIP 双侧检索投影进行参数高效、能力保持的多靶点适配，并联合 active/decoy 判别与跨口袋检索监督，可在完整未见靶点上显著改善 ROC-AUC、PR-AUC、BEDROC 和 EF5%；该提升不能由随机标签、严格同切分的 ECFP4 或简单动态构象聚合解释。

它是一个真实的方法学正结果，但还不是“广泛 SOTA”。进一步发表需要增加独立靶点数，并在不使用外部标签选模的冻结测试上复现。

## 10. 可部署权重

- `results/gpcr_drugclip_screening_v01/large_multitask_controlled_result.bce_retrieval.projection.pt`
- `results/gpcr_drugclip_screening_v01/large_multitask_controlled_result.full_multitask.projection.pt`

推荐当前使用 `bce_retrieval` 权重；完整多任务的 triplet 增量尚未通过消融。

13 靶点正式部署使用三个全数据权重平均：

- `results/litpcba_drugclip_external_v01/drugclip_13target_bce_retrieval.seed20260925.projection.pt`
- `results/litpcba_drugclip_external_v01/drugclip_13target_bce_retrieval.seed20260926.projection.pt`
- `results/litpcba_drugclip_external_v01/drugclip_13target_bce_retrieval.seed20260927.projection.pt`

## 11. 关键复现文件

- `scripts/finetune_drugclip_large_multitask.py`
- `scripts/evaluate_drugclip_large_target_loso.py`
- `scripts/evaluate_drugclip_large_external_panel.py`
- `scripts/evaluate_drugclip_docking_nested_fusion.py`
- `results/gpcr_drugclip_screening_v01/large_multitask_controlled_result.json`
- `results/gpcr_drugclip_screening_v01/large_target_loso_3seed_result.json`
- `results/drugclip_m4_external_panel_v01/large_multitask_external_result.json`
- `scripts/build_litpcba_drugclip_external.py`
- `scripts/evaluate_litpcba_drugclip_external.py`
- `results/litpcba_drugclip_external_v01/frozen_external_result.json`
- `scripts/evaluate_drugclip_13target_loso.py`
- `scripts/evaluate_drugclip_13target_ecfp_loso.py`
- `scripts/summarize_drugclip_13target_multiseed.py`
- `results/litpcba_drugclip_external_v01/combined_13target_loso_3seed_summary.json`

## 12. 适配器与损失扩展消融

进一步的受控消融发现，当前性能主要来自口袋端重排与双侧联合适配，而不是分子端单独记忆。Target-DRO 明显有害，已否决。三种子严格 LOSO 中，去掉表示保持约束后，宏 ROC-AUC/PR-AUC/BEDROC20 从 `0.623/0.247/0.300` 提高到 `0.637/0.258/0.311`；三项 paired scaffold-cluster bootstrap 95% CI 均严格为正，EF1% 和 EF5% 差异仍跨零。

侧别锚定显示，保留分子空间、放开口袋端显著优于保留口袋空间、放开分子端。这支持“通用分子语义可保留，但口袋表示需要域适配”的解释。完整数据、失败方向和部署权重见 `docs/DRUGCLIP_TARGET_TRANSFER_ABLATION.md`。

同一协议下，更新完整末层投影的基线含 131,328 个训练参数，宏 ROC-AUC/PR-AUC/BEDROC20/EF5% 为 `0.591/0.221/0.266/2.00`，低于仅含 10,240 个参数的 rank-8 LoRA `0.626/0.251/0.307/2.34`。这支持低秩更新的参数效率，但 full-last 目前仅为预设学习率基线。
