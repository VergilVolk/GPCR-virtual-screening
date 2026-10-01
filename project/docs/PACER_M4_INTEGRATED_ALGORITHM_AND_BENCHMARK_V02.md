# PACER-M4 整体算法与基准评估 v02

日期：2026-10-01

## 项目主线

PACER-M4 是有拒绝机制的证据级联，而不是加权总分：

1. CGDA / DrugCLIP 进行大规模结合检索；
2. docking、IFP 与 PACER-FS 检查姿势和口袋相容性；
3. A、P、C、CP 四上下文 MD 构造 `Delta_PAM`、`Delta_AGO` 和 `Delta_INT`；
4. 冻结 C1-BS256 + PACER-FKG 检查预定义区域的跨 replica 方向一致性；
5. 输出候选证据卡，证据不足时拒绝作出 PAM 结论。

## 当前主结果：预冻结 hard-negative 检验

- CM00717（PAM）与 CM00734（实验 inactive）来自同一系列，ECFP4 Tanimoto 约 0.955，是预先指定的困难负样本对。
- Stage B 沿用同一冻结 encoder、normalization、RFF、graph、region 和 contrast；未根据 CM00734 结果重拟合或选阈值。
- 主区域 `STATE_MOTION / Delta_INT / compound110_extension` 的 R1/R3 direction cosine 为 compound110 `+0.584`、LY2119620 `+0.634`、CM00734 `-0.327`；CM00734 在完整 9 个报告区域均为负。
- CM00734 的 Delta_INT 幅度并不为零：相对 LY2119620 的 magnitude 比例仍为 R1 `0.644`、R2 `0.798`、R3 `0.525`。有效信息来自跨 replica 的方向一致性，而不是变化幅度或结合强弱。
- 这是预冻结、apply-only、零结果驱动调参的回顾性计算验证；不是标签未知的前瞻盲测、湿实验或已验证 PAM 分类器。

静态证据必须准确表述：冻结 source-holdout 功能探针把 CM00734 排在 CM00717 之上（0.779 vs 0.642），ECFP-PU 对两者近乎饱和且 inactive 略高（0.99699 vs 0.99684）；raw DrugCLIP state-max 有轻微正确方向（0.0448 vs -0.0008），因此不能笼统说“所有 DrugCLIP 都分不开”。CM00734 可获得 ensemble docking 相容 pose（Vina `-9.02 kcal/mol`、口袋覆盖 `0.714`），但 CM00717 未进入同一冻结 docking 面板，不能宣称已完成严格成对 docking 失败证明。

## 其他冻结结果

- Full LIT-PCBA 15-target LOTO：CGDA ROC-AUC 0.5788、BEDROC80.5 0.07378；相对 reference retrieval 的 paired target-bootstrap 增益分别为 +0.00954 `[+0.00134,+0.01840]` 和 +0.00274 `[+0.00045,+0.00504]`。
- 外部 DUD-E 5-GPCR：CGDA 相对 reference retrieval 的 ROC、PR、BEDROC 均改善，但未超过原始 pocket DrugCLIP。
- M4 外部功能面板：GPCR-triplet 在 Acadia 和 Monash 的 ROC-AUC 分别为 0.930 和 0.788；负样本仅 2 和 3，差值区间均跨零。
- M4 外部效力：Suven 与 VU6025733 的 Spearman 为 -0.106 和 -0.637；当前 DrugCLIP 专化不能预测 PAM 效力。
- PACER-200 经固定双模型规则保留 5 个优先验证候选；尚未完成候选四上下文闭环，不能称为已确认 PAM。

## 可主张边界

当前不能宣称端到端 PAM SOTA，也不能把不同数据协议拼成一个总 AUC。最强且最贴合项目科学问题的结论是：在预冻结、零结果驱动调参的回顾性 hard-negative 测试中，PACER-FKG 区分了已知 PAM 的可重复协同方向与实验 inactive 高相似近邻的不一致方向，提供了初步功能特异性证据。

生成入口：

```powershell
$env:PACER_ARTIFACT_ROOT = 'D:\CLC'
python project\scripts\build_pacer_m4_algorithm_report_v01.py
```

输出：`output/pdf/PACER-M4_整体算法与基准评估_v02.pdf`。
