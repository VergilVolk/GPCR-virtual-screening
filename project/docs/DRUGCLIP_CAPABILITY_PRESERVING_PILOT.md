# DrugCLIP 能力保持式 GPCR 微调：单种子 pilot

更新日期：2026-09-25

## 研究目标

同时检验三件事：修复官方 DrugCLIP 在 GPCR 检索上的弱项；利用同分子跨毒蕈碱亚型 triplet 学习口袋差异；通过 embedding 与官方分数蒸馏限制灾难性遗忘。该任务不预测 PAM 功能或效力。

## 方法

- 冻结官方两个 encoder，仅在分子与口袋末端投影加入 rank-4 LoRA；
- 主损失：四 GPCR 靶点平衡交叉熵；
- 辅助损失：同一分子的 M1/M2/M3/M5 正负口袋 triplet；
- 保持项：新旧 embedding 余弦保持与官方分数 MSE 蒸馏；
- 五折全局 Murcko-scaffold OOF；每折删除与测试折同骨架的辅助 triplet；
- 对照：CE-only 与随机方向 triplet。

## 四 GPCR 检索结果

| 方法 | 宏 Recall@1 | 总 Recall@1 | MRR | Pair AUC | Pair PR-AUC |
|---|---:|---:|---:|---:|---:|
| 官方 DrugCLIP | 0.4130 | 0.3580 | 0.6525 | 0.7434 | 0.3837 |
| CE-only | 0.9289 | 0.9808 | 0.9894 | **0.9919** | **0.9783** |
| CE + 真 triplet | **0.9290** | **0.9812** | **0.9897** | 0.9916 | 0.9770 |
| CE + 随机方向 triplet | 0.9289 | 0.9808 | 0.9894 | 0.9919 | 0.9782 |

真实 triplet 只带来极小 Recall 增益，且随机方向对照几乎相同；目前不能声称辅助 triplet 改善四 GPCR 域内检索。

## 严格未见 M4 结果

为避免泄漏，另训练一个不含 M4R 靶点、且排除测试分子骨架的联合模型。训练只使用 B2AR/CCR2/M2R 与 M1/M2/M3/M5。

| 方法 | 配对准确率 | 方向平衡准确率 |
|---|---:|---:|
| 官方 DrugCLIP | 0.5752 | 0.6170 |
| 联合适配器 | 0.6209 | 0.6446 |
| 既有 muscarinic-triplet 专家 | **0.660** | **0.690** |

联合适配器相对官方模型的分子成簇 95% CI 为 `+0.0077` 至 `+0.0824`，分子置换检验 `p=0.039`。但联合训练弱于专门的 muscarinic-triplet 专家，说明存在任务负迁移。

## 决策

不推广单一联合 checkpoint。采用能力条件化路由：

1. 一般或未验证靶点：保留官方 DrugCLIP；
2. 有 GPCR 域标签的已知靶点：使用 target-balanced CE adapter；
3. 未见毒蕈碱亚型：使用 same-molecule muscarinic-triplet adapter；
4. M4 PAM 身份与效力：交给独立功能药理/动态模块，不由 DrugCLIP 分数代替。

下一阶段应在未参与开发的官方通用虚拟筛选集和新的 GPCR 外部集上冻结验证路由规则，而不是继续读取同一个 M4 测试集调权重。

## 复现文件

- `project/scripts/finetune_drugclip_capability_preserving.py`
- `project/results/gpcr_drugclip_retrieval_v01/capability_preserving_pilot.json`
- `project/results/gpcr_drugclip_retrieval_v01/capability_preserving_m4_zero_shot.json`

