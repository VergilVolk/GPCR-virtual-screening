# DrugCLIP Science 2026 LIT-PCBA 协议审计（阶段性）

## 当前状态

- 完成：8/15 个 LIT-PCBA 靶点，506,984 个分子记录。
- 运行中：OPRK1；其余结果由断点续跑保留。
- 当前数值只能作阶段性诊断，禁止与论文 15 靶点宏平均直接比较。

## 已验证正确的部分

1. 权重为官方发布的 `litpcba_identity_90.pt`，checkpoint 架构为 `drugclip`，epoch 68。
2. checkpoint 严格载入，无 missing/unexpected keys；分子和口袋词表维度与 checkpoint 一致。
3. 数据编码、L2 归一化、口袋—分子余弦相似度及“多个口袋取最大值”与官方仓库评测代码一致。
4. 8 个已完成靶点均无重复 molecule ID、无 NaN/Inf。
5. 保存分数与 embedding 重算的最大绝对误差小于 `4e-7`。
6. AUC、BEDROC(α=80.5) 和 EF0.5/1/2/5% 与官方 RDKit 实现一致到 float32 误差范围；审计通过。
7. 论文方法说明 LIT-PCBA benchmark 使用单个、按 90% 序列同源性过滤训练数据的非集成模型；六折模型主要用于后续湿实验筛选，不是本基准的必要条件。

审计产物：

- `project/results/drugclip_science2026/full_litpcba/protocol_audit_interim.json`
- `project/scripts/audit_drugclip_science2026_protocol.py`

## 当前不能声称的内容

- 尚未完成 15 靶点，不能声称复现论文宏平均。
- 当前任务是官方 DrugCLIP baseline，不是 PACER/Triplet/GPCR 微调结果。
- LIT-PCBA 测量通用活性分子排序，不测 CHRM4 PAM 的协同性、探针依赖性或内在激动。
- 即使 LIT-PCBA 提升，也只能支持“通用虚拟筛选排序改善”，不能证明候选是功能性 PAM。
- 未与 ECFP、DrugCLIP、固定融合及我们的模型在完全相同 15 靶点协议上完成统计比较，不能声称 SOTA。

## 8 靶点阶段性数值（仅诊断）

| 指标 | 宏平均 |
|---|---:|
| AUROC | 0.5988 |
| PR-AUC | 0.0410 |
| BEDROC80.5 | 0.0761 |
| EF0.5% | 6.2397 |
| EF1% | 4.8633 |
| EF5% | 3.0452 |

这些数字受 PPARG 等少数靶点明显影响，不代表最终 15 靶点结果。

## 下一道发表级验收门

1. 完成官方 15 靶点 baseline，并与论文 DrugCLIP 指标核对。
2. 在完全相同输入、口袋聚合与指标协议下运行 PACER 改进模型。
3. 加入强 ECFP ligand-only baseline，报告每靶点结果而非只报宏平均。
4. 对 15 个靶点做配对 bootstrap/置换检验，重点检验 BEDROC 与 EF1%，同时报告 worst-target 与负迁移靶点。
5. 把通用 benchmark 证据与 CHRM4 PAM 功能证据分开陈述；后者仍需独立药理标签或实验验证。
