# Stage4 Temporal Metadata Note v01

归档日期：2026-10-03。状态：**DOCUMENTED — SOURCE RAW-FRAME METADATA DESCRIPTION ERROR**。

本说明登记源字段与实际数据的差异，引用已有核验记录，不生成第二套同义JSON，不改变任何原始文件、冻结结果或方法。

## 1. 源manifest的原始声明

`C:/projects/PACER_STAGE4_MD_backup/master_manifest.json`包含`trajectory_ps: 50`。若把该字段当作实际原始轨迹帧间隔，会与已下载DCD不符。保留其原始字节，不静默修订该文件。

## 2. DCD实际检查结果

最终科学摘要的`provenance.independent_dcd_header_checks`记录独立读取36/36条DCD头部：每条1000帧，间隔为`10.000000029814057 ps`（浮点表示下的10 ps）。`input_audit/INPUT_AUTHENTICATION_v01.json`记录相同1000帧/10 ps间隔、state.csv最终步数5,000,000及既有坐标有限性认证。2 fs生产步长对应10 ns/trajectory；36条共360 ns。

这些检查已在科学复核阶段完成。本次最终收尾只读取已有记录与文件哈希，不重新读取坐标进行科学分析，也不运行MD或PACER-FKG阶段。

## 3. PACER-FKG实际使用的stride

冻结Stage4 temporal contract为：

```text
1000 raw frames / trajectory, raw spacing 10 ps
→ frames[::5], retained indices 0, 5, ..., 995
→ 200 analysis frames, analysis spacing 50 ps
→ 20 frames / block
→ 10 correlated temporal blocks / trajectory
```

## 4. 分析间隔为何正确

实际分析帧间隔为`10 ps × 5 = 50 ps`。20帧组成名义1 ns block，19个内部lag覆盖950 ps。block内部计算时间差分，不跨block做差，不插值。不能把10个相关block当作独立统计样本。

## 5. 为什么不影响已完成的科学结果

适配器强制验证DCD的1000帧/10 ps实际契约；Phase2a使用固定`frames[::5]`构建200帧分析输入。它不根据源`master_manifest.json`中的`trajectory_ps`字段推导本次采样间隔。实际输入、冻结Stage4时间契约、回执和最终结果一致。

因此，问题属于源manifest的raw-frame metadata描述错误，而不是已经发现的时间采样错误或冻结状态污染。正式归档确认最终50 ps分析间隔正确，既有分析结果保持有效；不重新运行Phase1/2a/2b，不改DCD，不调整任何科学参数。历史数值freeze SHA仍为`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`。

## 6. 为什么保留源manifest

原manifest是来源与过程证据；静默修订会掩盖原始描述错误、破坏来源字节可追溯性。故在独立provenance文件中登记事实，原文件保留。

本说明不推测服务器启动参数、字段生成原因或覆盖记录。早期科学报告/摘要中的metadata caveat保留其复核时点的判断；本次项目收尾依据已核实实际时间契约，将影响定性登记为不影响既有分析的源字段描述错误，不重写早期记录。

## 7. 证据与claim boundary

- [最终科学摘要](../STAGE4_PROSPECTIVE_SCIENTIFIC_SUMMARY_v01.json)：`issues_requiring_confirmation`及`provenance.independent_dcd_header_checks`。
- [最终科学报告](../STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md)：Executive Summary、Dataset/Temporal Design及限制。
- [输入认证回执](../input_audit/INPUT_AUTHENTICATION_v01.json)：36条原始输入及其哈希。
- [Phase2a冻结回执](../phase2a_frozen_apply/PHASE2A_FREEZE_RECEIPT_v01.json)：实际temporal contract。

> Prospective frozen PACER-FKG dynamic differential evaluation only.
> No PAM/ago-PAM labels, probability, efficacy or independent-block inference.

metadata描述错误的登记不增加药理证据，不支持确定cooperative PAM mechanism，不产生新的功能标签、评分或预测概率。
