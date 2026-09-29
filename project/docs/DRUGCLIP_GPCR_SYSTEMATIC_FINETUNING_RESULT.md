# DrugCLIP 的 GPCR 参数高效微调：阶段结果与方法边界

更新日期：2026-09-25

## 1. 研究问题

本阶段不直接预测 M4 PAM 功能，而先回答一个更基础、可严格评价的问题：

> 在保留 DrugCLIP 预训练蛋白口袋—分子表征的前提下，能否用少量参数把模型特化到 GPCR 别构配体检索，并在严格留出数据上稳定优于官方 DrugCLIP？

DrugCLIP 的基础架构最初发表于 NeurIPS 2023；扩展后的基因组级虚拟筛选工作于 2026 年发表于 Science（DOI: `10.1126/science.ads9530`）。原模型将虚拟筛选表述为蛋白口袋与小分子的密集检索。本阶段沿用这一任务定义，不把检索分数解释为 PAM 效力。

## 2. 数据与切分

使用公开 GPCR 别构调节剂 benchmark 中四个靶点的阳性分子和相应 GaMD 代表构象：B2AR、CCR2、M2R、M4R。

- 原始阳性靶点—分子对：2,624；
- 删除 60 个跨靶点重复分子后：2,504 个单标签分子；
- DrugCLIP 三维构象成功：2,500 个；
- 独立 Murcko 骨架：968 个；
- 口袋直接取原 benchmark 的别构位点残基和 GaMD 代表构象，不使用正构口袋；
- 主评价采用五折 Murcko-scaffold OOF，每个分子只在未参与训练的折中评价。

样本极不平衡：B2AR 64、CCR2 16、M2R 170、M4R 2,250。因此总体 Recall 或总体 AUC 可能被 M4R 主导，主指标必须使用宏平均 Recall@1、逐靶点结果和强二维基线。

## 3. 方法

冻结官方 DrugCLIP 的两个 Uni-Mol encoder 及原始投影，在分子、口袋两个末端投影层分别加入 rank-4 LoRA：

\[
z_m=\operatorname{norm}(W_mh_m+B_mA_mh_m),\qquad
z_p=\operatorname{norm}(W_ph_p+B_pA_ph_p).
\]

双侧共训练 5,120 个参数。训练目标为靶点平衡的四分类交叉熵，并用保持项限制新旧 embedding 的余弦偏移。比较四种方法：

1. 官方冻结 DrugCLIP；
2. 随机标签 LoRA；
3. 靶点平衡 CE LoRA；
4. CE 加 hardest-wrong-pocket triplet。

外部强基线为 ECFP4 + class-balanced logistic regression，使用完全相同的五折骨架切分。

## 4. 已冻结的 OOF 结果

| 方法 | 宏 Recall@1 | 总 Recall@1 | MRR | pair ROC-AUC | pair PR-AUC |
|---|---:|---:|---:|---:|---:|
| 官方 DrugCLIP | 0.413 | 0.358 | 0.653 | 0.743 | 0.384 |
| 随机标签 LoRA | 0.235 | 0.840 | 0.912 | 0.915 | 0.717 |
| 平衡 CE LoRA | **0.925** | 0.990 | 0.994 | **0.9965** | **0.9860** |
| CE + hard triplet | 0.923 | 0.989 | 0.994 | 0.9965 | 0.9858 |
| ECFP4 logistic | **0.963** | **0.995** | **0.997** | **0.9996** | **0.9991** |

平衡 CE 相对官方模型的宏 Recall@1 提升约 51.2 个百分点。CE + hard triplet 相对 CE 的骨架 bootstrap 95% 区间为 -0.53 至 0.00 个百分点，因此 triplet 在这个四分类任务上没有增益，不进入最终模型。

随机标签模型的总 Recall@1 达到 0.840，但宏 Recall@1 只有 0.235，且 B2AR、CCR2 几乎完全失败。这说明不平衡数据上的总体 Recall 和 flattened AUC 可以产生严重误导。

## 5. 化学距离审计

每个测试分子只与其训练折计算最大 ECFP4 Tanimoto。结果表明 2,448/2,500 个分子的最大相似度不低于 0.50；真正远离训练化学空间（Tanimoto < 0.30）的分子只有 18 个。

| 与训练集的最大 Tanimoto | n | 官方宏 R@1 | CE 宏 R@1 | ECFP 宏 R@1 |
|---|---:|---:|---:|---:|
| < 0.30 | 18 | 0.425 | 0.417 | 0.467 |
| 0.30–0.50 | 34 | 0.220 | 0.813 | 0.958 |
| 0.50–0.70 | 822 | 0.403 | 0.926 | 0.985 |
| ≥ 0.70 | 1,626 | 0.433 | 0.999 | 1.000 |

因此，域内大幅提升主要来自 GPCR 化学系列特化；目前没有证据表明 CE adapter 在低相似新骨架上优于 ECFP。

## 6. 留一整靶点外推

每轮完全移除一个 GPCR 的全部分子，训练损失也不使用被留出口袋；四轮预测合并后评价：

| 方法 | 宏 Recall@1 | 总 Recall@1 | MRR | pair ROC-AUC |
|---|---:|---:|---:|---:|
| 官方 DrugCLIP | **0.413** | **0.358** | **0.653** | **0.743** |
| 三靶点 CE 迁移 | 0.370 | 0.292 | 0.610 | 0.720 |
| 随机标签迁移 | 0.077 | 0.173 | 0.561 | 0.704 |

CE 迁移相对官方模型的宏 Recall@1 差值 scaffold-bootstrap 95% CI 为 -9.29 至 +3.21 个百分点，中位数为负。专项 CE 微调不能取代官方模型的零样本能力。

## 7. 当前结论

### 可以成立

1. 双侧投影 LoRA 能以 5,120 个参数显著提高 DrugCLIP 在已知 GPCR 域内的靶点检索表现。
2. 类别平衡 CE 是当前最优微调目标；hard triplet 在该任务中无额外价值。
3. 宏平均、逐靶点和随机标签对照是必要评价；只报告总体 AUC/Recall 会得到错误结论。
4. 模型应采用任务路由：已知目标域使用 CE adapter；未见 GPCR 保留官方 DrugCLIP或已经验证的跨亚型关系微调。

### 不能成立

1. 不能称为通用 GPCR SOTA：ECFP4 在该回顾性数据上更强。
2. 不能声称改善未见靶点或真正新骨架泛化。
3. 不能把 target retrieval 分数解释为结合亲和力、PAM 身份、EC50、协同性或内在激动。

## 8. 下一评价闸门

下一阶段不继续堆叠损失函数，而完成与 DrugCLIP 原任务一致的四靶点 active/decoy 虚拟筛选：按靶点报告 ROC-AUC、PR-AUC、BEDROC、EF1%/EF5%，同时加入 ECFP、官方 DrugCLIP、静态 docking 和 GaMD ensemble docking。只有在外部靶点或低相似骨架上超过强基线，才能主张方法学提升。

## 9. 复现文件

- `project/scripts/build_gpcr_drugclip_retrieval_benchmark.py`
- `project/scripts/build_drugclip_inputs.py`
- `project/scripts/extract_drugclip_embeddings_cpu.py`
- `project/scripts/finetune_drugclip_gpcr_retrieval.py`
- `project/scripts/analyze_gpcr_drugclip_similarity_strata.py`
- `project/scripts/evaluate_gpcr_drugclip_target_loso.py`
- `project/scripts/fit_gpcr_drugclip_ce_checkpoint.py`
- `project/data/benchmarks/gpcr_drugclip_retrieval_v01/`
- `project/results/gpcr_drugclip_retrieval_v01/`
