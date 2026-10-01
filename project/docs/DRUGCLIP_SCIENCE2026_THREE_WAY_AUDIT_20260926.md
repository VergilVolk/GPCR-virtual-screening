# DrugCLIP 三者对照：首轮版本审计

审计日期：2026-09-26。此文件的目的不是宣称 SOTA，而是回答一个必要问题：项目此前在旧 DrugCLIP checkpoint 上观察到的适配结果，在已验证的 Science-2026 checkpoint 面前是否仍值得继续。

## 比较对象与共同条件

共同输入为冻结的九靶点 LIT-PCBA 派生外部子集（8,906 target--molecule pairs，8,763 unique molecules；固定口袋、固定 RDKit 构象、固定标签、相同 CPU encoder extraction）。该子集不是官方完整 LIT-PCBA，故不能与论文表格直接比较。

| 编号 | 方法 | 训练/权重身份 |
|---|---|---|
| A | Old DrugCLIP | 本地 2023 PDB-general task checkpoint |
| B | Previous PACER full-multitask | A 的冻结 embedding 上的既有适配结果 |
| C | Science-2026 DrugCLIP single | 官方 `litpcba_identity_90.pt`；SHA256 已核验 |

## 首轮结果（macro mean over 9 targets）

| 方法 | ROC-AUC | PR-AUC | BEDROC20 | EF1% | EF5% |
|---|---:|---:|---:|---:|---:|
| A: Old DrugCLIP | 0.5322 | 0.2029 | 0.2391 | 1.443 | 1.510 |
| B: Previous PACER | **0.5410** | **0.2042** | **0.2393** | **1.732** | 1.474 |
| C: Science-2026 single | 0.5193 | 0.1978 | 0.2188 | 1.708 | 1.224 |

## 严格解释

1. C 已确认为不同的、官方 hash 匹配的 2025-trained `pcba_90` checkpoint；它不是 A 的文件名变化。
2. 在**本项目这套九靶点、小规模、外部派生输入**上，C 未优于 A；这绝不推出 C 在官方全量 LIT-PCBA 或论文 protocol 上较差。训练分布、口袋定义、构象生成与论文官方输入都可能影响模型排序。
3. B 相对 A 的 ROC-AUC 增量为 +0.0088，但 BEDROC20 几乎不变、EF5% 下降；此前 bootstrap 也没有给出稳定的早期富集改善。因此 B 不是足以发表的强结论，只说明“旧权重上的低容量适配尚有微弱信号”。
4. 目前尚未在 C 上训练 PACER；将 B 的参数直接迁移到 C 是不合法的，因为 embedding 坐标系已变。故现阶段不能回答“PACER 是否超过 Science-2026”。

## 必须执行的下一步

1. 将 PACER training/evaluation 改为显式接收 checkpoint-specific embeddings，在 C 上从零训练；
2. 先完成本九靶点的 checkpoint-matched pilot（只判断是否有可重复增益）；
3. 下载并运行官方完整 LIT-PCBA/DUD-E 输入以及 6-fold 权重后，做预注册式的 full benchmark；
4. 主张边界：在第 3 步前仅可称为 *version-controlled pilot*，不可称 Science-2026 SOTA 超越或 PAM efficacy prediction。
