# CGDA：上下文门控 DrugCLIP 未见靶点适配

更新日期：2026-09-26

## 问题

共享 DrugCLIP LoRA 在 13 靶点小面板中有效，但扩展到官方 Full LIT-PCBA 后会产生负迁移。不同靶点需要不同适配方向，单一共享投影无法同时处理异质口袋和化学空间。

## 方法

CGDA 使用目标口袋集合与共晶参考配体集合构造上下文：

`c = [mean(pocket), mean(reference), abs(p-r), p*r]`

门控网络根据 `c` 产生专家权重。每个专家是一个低秩残差映射；正式版本采用两个 rank-2 专家。每个参考配体独立适配，候选得分取与全部适配参考配体的最大余弦相似度。

训练目标包括：

1. 靶点内 active/decoy BCE；
2. 活性分子的跨靶点检索损失；
3. 高分已知阴性构成的困难负样本 pairwise ranking；
4. 相对原始参考配体表示的轻量保持约束。

DrugCLIP 编码器完全冻结。模型使用训练靶点标签学习共享专家；留出靶点只提供口袋、共晶参考配体及无标签候选，活性标签只用于最终评价。

## 数据与协议

- 官方 LIT-PCBA 15 靶点；
- 2,807,612 个候选；
- 129 个真实口袋构象；
- 每折完整留出一个靶点；
- 三个独立种子；
- 三种子候选分数平均；
- 15 个外层靶点配对 bootstrap 20,000 次。

这是 reference-ligand-assisted screening，不是 ligand-free zero-shot。

## 正式结果

| 方法 | ROC-AUC | PR-AUC | BEDROC80.5 | EF0.5% | EF1% | EF2% | EF5% |
|---|---:|---:|---:|---:|---:|---:|---:|
| Pocket DrugCLIP | 0.5672 | 0.02379 | 0.06109 | 7.98 | 5.36 | 3.45 | 2.08 |
| Reference retrieval | 0.5693 | 0.02594 | 0.07103 | 8.71 | 6.06 | 3.61 | **2.65** |
| CGDA-hardrank ensemble | **0.5788** | **0.02627** | **0.07378** | **8.79** | **6.31** | **4.06** | 2.38 |

相对 reference retrieval 的配对 target-bootstrap：

- ROC-AUC：`+0.00954`，95% CI `[+0.00134,+0.01840]`；
- BEDROC80.5：`+0.00274`，95% CI `[+0.00045,+0.00504]`；
- PR-AUC：CI 跨零；
- EF0.5%：CI 跨零；
- EF1%：`+0.247`，CI 下界为 `0`，但只有少数靶点出现离散改善；
- EF5%下降，CI包含零附近，不主张改善。

BEDROC 在 15 个靶点中的 12 个得到改善。ROC-AUC 与 BEDROC 是当前通过不确定性检验的正式正结果。

## 外部 GPCR 压力测试（DUD-E，冻结后一次性推理）

为检查结果是否仅限于 LIT-PCBA，使用未出现在源训练中的 AA2AR、ADRB1、CXCR4、DRD3 和 M3R（`mcr`）构成外部 GPCR 面板。源模型为三个只用 LIT-PCBA 训练的冻结 CGDA 权重；外部标签未用于拟合、调参或选模。共 91,311 个候选和 1,342 个 DUD-E 标注活性物。候选读取器会改变 LMDB 顺序，故在评价前以 SMILES 多重集合核验并逐 SMILES 显式回填标签；五靶点均通过。

相对 reference retrieval，CGDA 三种子均值在 5/5 靶点改善：宏 ROC-AUC `0.7233→0.7437`（配对靶点 bootstrap 95% CI `[+0.00587,+0.03480]`），PR-AUC `0.2658→0.2762`（`[+0.00329,+0.02175]`），BEDROC80.5 `0.3637→0.3738`（`[+0.00269,+0.02320]`）。EF0.5%、1%、2%、5% 的均值差也均为正，但只有 `3/5`、`3/5`、`3/5`、`4/5` 靶点改善，不能表述为全靶点早期富集胜出。

这不是对原始 pocket DrugCLIP 的全面胜利：其宏 ROC-AUC/BEDROC80.5 为 `0.8105/0.3936`，高于 CGDA 的 `0.7437/0.3738`，尤其 DRD3 与 M3R 明显退化。因此外部结果只支持“CGDA 相比参考配体检索具有跨 GPCR 的稳定增益”，不支持“CGDA 优于所有 DrugCLIP 打分方式”或通用 SOTA。

DUD-E 使用为活性物生成的 property-matched decoy，故该面板是外部迁移压力测试，不是实验 inactive、真实 HTS 命中率或 PAM 功能验证。原始逐靶点分数、审计和报告：`results/dude_gpcr_external_v01/cgda_external_dude_gpcr.json`、`data/external/drugclip_official/selected/dude_gpcr/assembled/dude_gpcr_assembly.audit.json`。

## 已否决：跨靶点监督 score router

为避免简单线性融合的局限，测试了一个低容量的固定 L2 logistic router。每一外层靶点只使用其余14个靶点标签训练；输入为 pocket 与CGDA的靶点内百分位、差值、乘积和绝对差值，不含靶点ID。该 router 的宏 ROC-AUC/PR-AUC/BEDROC80.5 仅为 `0.5605/0.00926/0.02109`，显著低于 pocket 与CGDA；相对CGDA的BEDROC target-bootstrap CI 为 `[-0.1102,-0.0146]`。这证明“分数与活性的跨靶点监督映射”本身不可迁移，后续不能沿此方向继续调参。结果：`results/litpcba_drugclip_external_v01/full15_cgda_score_router_loso.json`。

## 消融结论

- 四专家 rank-4 会过拟合；二专家 rank-2 更稳；
- 去除 retrieval 后 BEDROC下降；
- 无困难负样本时整体 ROC 提升，但 BEDROC 的 CI 跨零；
- 加入困难负样本排序后，BEDROC 的 CI 严格为正；
- 共享监督 LoRA、完整投影微调、口袋到共晶配体直接对齐、标签堆叠、多口袋聚合均未稳定超过强基线。

因此，CGDA 的有效部分不是简单增加参数，而是“目标上下文门控 + 小容量专家 + 困难负样本排序”。

## 创新与边界

当前可主张：

> 针对共享 DrugCLIP 适配在异质靶点上的负迁移，CGDA 使用口袋—参考配体上下文选择低秩专家，并在严格未见靶点的 Full LIT-PCBA 评价中显著改善 ROC-AUC 与 BEDROC。

当前不可主张：

- 通用虚拟筛选 SOTA；
- ligand-free zero-shot；
- 所有 EF 截断点都改善；
- M4 PAM 身份、协同性或效力；
- prospective 湿实验有效性。

CGDA 是 PACER 的 Binding 模块；PAM 功能仍由 PACER-DC 四上下文动态与独立药理证据处理。

## 复现入口

- `scripts/evaluate_drugclip_cgda_loso.py`
- `scripts/summarize_drugclip_cgda_multiseed.py`
- `scripts/evaluate_drugclip_cgda_nested_blend.py`
- `scripts/summarize_drugclip_cgda_nested.py`
- `results/litpcba_drugclip_external_v01/full15_cgda_small_hardrank_seed2026092*.json`
- `results/litpcba_drugclip_external_v01/full15_cgda_small_hardrank_3seed_summary.json`
