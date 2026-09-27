# PACER-FKG G1 · 冻结结果只读专项审计

## 范围与数据完整性

R2/R3 × W0–W4 × 四上下文；输入为已归档的 10 份窗口 JSON 及综合报告。
综合报告 SHA256：`b7dbcb11f07b543694d85a680125c25b1fac2edf16be92ad0f9ff84f174bd90d`；图 SHA256：`25593c11dac0f70ad6811c8c6aa7a4f06c1755f14397f1162eaa1157026be122`。

没有重新提取 Embedding、训练编码器、重新计算核统计或更改原有门禁。

## 原门禁及图传播描述性结果

| 轴 | 原门禁通过区域 | 图传播后窗口 rho 上升/下降/持平（七功能区域） |
|---|---|---|
| synergy_interaction | compound110_extension, pam_contact_union | 6/1/0 |
| intrinsic_agonism | pam_contact_union | 3/3/1 |
| conditional_pam_effect | compound110_extension, cooperativity_mutagenesis, intracellular_microswitches, pam_contact_union | 5/2/0 |

## 分区核带宽与跨 replica 描述性一致性

详细数据见 `G1_axis_region_comparison.csv`、`G1_replica_region_summary.csv` 和 `G1_window_diagnostics.csv`。

特别注意：当前区域核 U² 和先计算残基级核距离再图传播的区域均值不是同一估计量。

## 远端对照

- synergy_interaction：R2 中位差 0.244849，R3 0.135984；原 rho=-0.700，传播后 rho=0.000。
- intrinsic_agonism：R2 中位差 0.081301，R3 0.066316；原 rho=-1.000，传播后 rho=-0.700。
- conditional_pam_effect：R2 中位差 0.131091，R3 0.080531；原 rho=-0.400，传播后 rho=-0.300。

## 解释限制

- Five contiguous windows are not five independent replicates.
- Per-region/per-window and per-axis median-heuristic bandwidths limit cross-region raw kernel magnitude comparisons.
- Squared RKHS contrast norm is not the direction or pharmacological sign of a response.
- Graph-diffused residue kernel-root means and joint-region kernel U2 are distinct estimands; their numerical magnitudes are not directly comparable.
- No inference on graph causality from pre/post graph Spearman alone; no shuffled-graph control here.
- The frozen gate excludes distal control by definition, so passing it does not certify distal specificity.
- The prior-registration status of the frozen numerical gate thresholds has not been independently verified.

## 后续工作（不属于本次 G1 结果）

先核验 atom14/PBC 的独立 QC，再在新冻结方案下开展共同核空间方向比较、时间块重采样，以及无图/置乱图对照。
