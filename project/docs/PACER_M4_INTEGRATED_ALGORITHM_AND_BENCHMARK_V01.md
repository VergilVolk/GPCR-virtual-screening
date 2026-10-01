# PACER-M4 整体算法与基准评估 v01

日期：2026-10-01

## 主线

PACER-M4 是一个有拒绝机制的证据级联，而不是加权总分：

1. CGDA / DrugCLIP 进行大规模结合检索；
2. docking、IFP 与 PACER-FS 检查姿势和口袋相容性；
3. A、P、C、CP 四上下文 MD 构造 `Delta_PAM`、`Delta_AGO` 和 `Delta_INT`；
4. 冻结 C1-BS256 + PACER-FKG 检查预定义区域的跨 replica 方向一致性；
5. 输出候选证据卡或拒绝功能结论。

## 整体结果

- Full LIT-PCBA 15-target LOTO：CGDA ROC-AUC 0.5788、BEDROC80.5 0.07378；相对 reference retrieval 的 paired target-bootstrap 增益分别为 +0.00954 `[+0.00134,+0.01840]` 和 +0.00274 `[+0.00045,+0.00504]`。
- 外部 DUD-E 5-GPCR：CGDA 相对 reference retrieval 的 ROC、PR、BEDROC 均改善，但未超过原始 pocket DrugCLIP。
- M4 外部功能面板：GPCR-triplet 在 Acadia 和 Monash 的 ROC-AUC 分别为 0.930 和 0.788；两组的负样本仅 2 和 3，paired improvement CI 均跨零。
- M4 外部效力：Suven 与 VU6025733 的 Spearman 为 -0.106 和 -0.637；当前 DrugCLIP 专化不能预测 PAM 效力。
- PACER-FKG：主区域 `STATE_MOTION / Delta_INT / compound110_extension` 的 R1/R3 cosine 为 compound110 +0.584、LY2119620 +0.634、CM00734 -0.327。该结果是冻结的 hard-negative partial specificity support，不是分类器性能。
- PACER-200 交接：同一 200 个候选经旧 DrugCLIP rank<=50 且 family-aug rank<=25 的固定规则保留 5 个候选；它们仍是待验证假设。

## SOTA 边界

当前不能宣称端到端 PAM SOTA。五类协议的样本、终点和统计单位不同，不能拼接成一个总 AUC。可主张的是：CGDA 在严格未见靶点检索上有显著但有限的增益；PACER-FKG 在三化合物 matched design 中获得部分 hard-negative 特异性支持；整个项目形成了从结合到动态功能复核的可审计算法链。

生成入口：

```powershell
$env:PACER_ARTIFACT_ROOT = 'D:\CLC'
python project\scripts\build_pacer_m4_algorithm_report_v01.py
```

输出：`output/pdf/PACER-M4_整体算法与基准评估_v01.pdf`。

