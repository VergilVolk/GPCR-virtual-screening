# DrugCLIP 毒蕈碱受体 Triplet 微调：当前证据与边界

## 1. 要解决的问题

官方 DrugCLIP 主要学习一般性的蛋白口袋—小分子结合匹配，并不专门理解毒蕈碱受体亚型差异，更不能直接预测 M4 PAM 功能。本工作尝试用同分子跨亚型 triplet 微调其内部口袋投影：

\[
L=\max(0,m-s(M,P_{active})+s(M,P_{inactive})).
\]

同一个分子同时出现在正、负两侧，分子量、脂溶性和骨架等分子侧捷径在单个 triplet 内被抵消。训练冻结两个 encoder 和分子投影，只更新 DrugCLIP 的口袋投影；随机翻转正负口袋作为等规模对照。

## 2. 已确认的结果

### 2.1 已见亚型、未见分子的严格测试

训练集包含 716 个分子、1,383 个 M1/M2/M3/M5 triplet；测试集包含 33 个完全未见分子、54 个配对。

| 方法 | 配对准确率 |
|---|---:|
| 官方 DrugCLIP | 59.3% |
| 随机方向 triplet | 61.1% |
| 定向 GPCR triplet | **79.6%** |
| ECFP4-Ridge | **98.1%** |
| ECFP4-Random Forest | 96.3% |

定向 triplet 相对官方 DrugCLIP 的分子成簇 bootstrap 95% CI 为 **+5.9 至 +34.0 个百分点**，相对随机方向为 **+5.8 至 +32.1 个百分点**。但是二维模型明显更强；该测试主要是已知化学系列内插，不能据此声称 DrugCLIP 方法达到 SOTA。

### 2.2 未见 M4 口袋与未见分子的零样本测试

训练不含 M4 口袋，并删除全部训练分子后，冻结模型测试 153 个 M4—其他亚型配对、102 个分子。

| 方法 | 总准确率 | 方向平衡准确率 |
|---|---:|---:|
| 官方 DrugCLIP | 57.5% | 61.7% |
| 定向 GPCR triplet | **66.0%** | **69.0%** |

- 相对官方模型提升的分子成簇 95% CI：**+3.1 至 +14.3 个百分点**；
- M4 应更强时准确率 61.0%，M4 应更弱时 77.1%；
- 打乱分子与口袋配对后，零分布中位准确率 56.2%，实际 66.0%，单侧置换 `p=0.018`。

这是当前最有价值的阳性结果：微调学到了一部分分子依赖的跨亚型口袋匹配，而不只是把 M4 口袋整体抬高。该测试是开发过程中追加的冻结测试，仍属于回顾性证据。

### 2.3 五亚型 leave-one-target-out

进一步建立 1,976 对、841 个分子的五亚型 LOSO：每轮完全拿掉一个受体口袋，并删除该轮全部测试分子。

| 方法 | 五靶点宏平均准确率 | 宏平均方向平衡准确率 |
|---|---:|---:|
| 官方 DrugCLIP | 53.6% | 55.2% |
| 全口袋投影 triplet | 55.4% | 57.9% |
| rank-4 pocket LoRA（2,560 参数） | **55.6%** | **58.0%** |

结果具有明显靶点异质性：M2 改善最稳定，M3轻度改善，M1下降，M4/M5整体变化有限。因此 LoRA 只证明了参数效率，尚未证明普适的未见 GPCR 迁移能力。

## 3. M4 PAM 外部迁移

五套外部来源共 152 条记录、147 个唯一分子；每个药理终点独立评价，训练精确重叠被排除。

### PAM / inactive 判别

| 外部集 | 阴性数 | 官方 DrugCLIP AUC | GPCR-triplet AUC | M4 再微调 AUC |
|---|---:|---:|---:|---:|
| Acadia | 2 | 0.837 | **0.930** | 0.860 |
| Monash | 3 | 0.576 | **0.788** | 0.606 |

两个来源方向一致，但阴性总数只有 5 个，配对差值置信区间仍跨零。结果只支持“值得继续验证”，不能称为确定的外部 SOTA。

### 功能效力排序

| 外部集 | n | GPCR-triplet Spearman | M4 再微调 Spearman |
|---|---:|---:|---:|
| Suven pEC50 | 45 | -0.106 | -0.040 |
| VU6025733 pEC50 | 18 | -0.637 | -0.581 |
| US20260055116 pERK 分箱 | 25 | 0.107 | 0.103 |

当前 DrugCLIP 分数不能预测 PAM 效力。VU 系列出现显著反向关系，说明结合匹配、PAM 身份与同系列 EC50 排序不能混为一个任务。

## 4. 被否决的路线

1. M4R 人工 active/decoy AUC 接近 0.99，但 ECFP4 同样达到 0.9998，属于 decoy 来源捷径，不能作为模型证据。
2. 全分子投影 M4 功能微调在内部系列留出从 0.610 提到 0.655，但跨来源外部集退化，属于小样本过拟合。
3. rank-4 M4 分子 LoRA 的系列留出 AUC 为 0.630，随机配对对照为 0.641，未通过方向性门槛。
4. 54 对已见亚型测试上二维模型达到 98.1%，因此该集合不能作为结构模型的 SOTA 主证据。

## 5. 当前可成立的创新表述

> 通过同分子跨亚型药理 triplet 对 DrugCLIP 口袋空间进行参数高效特化，可显著改善官方模型对未见 M4 口袋和未见分子的零样本亚型排序；但这一结构迁移信号不能替代 PAM 功能与效力模型。

当前算法应拆成两个职责明确的专家：

1. **结构迁移专家**：DrugCLIP GPCR-triplet，回答候选与目标亚型口袋是否相容，重点用于未见靶点和低标签场景；
2. **功能药理专家**：二维 QSAR、实验上下文和后续四上下文动态特征，回答是否为 PAM、效力与内在激动风险。

只有在新的、位点匹配的外部数据上冻结验证这两个专家及其门控规则后，才可提出 SOTA 主张。任何计算候选仍只能称为待实验验证的 PAM 假设。

## 6. 对应证据文件

- `project/results/muscarinic_triplet_expansion_v01/external_result.json`
- `project/results/muscarinic_triplet_expansion_v01/strict_2d_baselines.json`
- `project/results/muscarinic_m4_zero_shot_v01/zero_shot_result.json`
- `project/results/muscarinic_target_loso_v01/target_loso_result.json`
- `project/results/muscarinic_target_loso_v01/lora_loso_result.json`
- `project/results/drugclip_m4_external_panel_v01/frozen_external_evaluation.json`
- `project/results/drugclip_m4_true_finetune_v01/matched_functional_lora_series.json`
