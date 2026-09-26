# PACER-DC Progress Report
Date: 2026-09-23
Branch: audit/oneprot-pacer-dc-g0

## 1. Project objective

Evaluate reproducibility of frozen OneProt-MD trajectory
embeddings across four matched CHRM4 simulation contexts.

d_PAM = z_CA - z_A
d_AGO = z_C - z_0

These are descriptive differential features, not validated
biological activity predictions.

## 2. Completed four-context technical tests

LY2119620:
- Known PAM positive control.
- Replica 01, window 000 completed.
- Initial embedding-level differential calculated.

compound110:
- Allosteric agonist specificity control.
- Replica 01, window 000 completed.
- Four-context re-audit with extractor 498a2d2: PASS.
- Independent topology provenance verified.
- Checkpoint loading audit: PASS.
- All six old/new vectors exactly identical.

Initial descriptive differential results:

                  LY2119620    compound110
||d_PAM||2          0.630557       0.544822
||d_AGO||2          0.448729       0.744751
cosine             -0.564953       0.385785

These observations do not establish functional classification.

## 3. Apo replica 2

- Production: 5 ns complete.
- Recovery: binary checkpoint.
- 500 frames, 10 ps spacing.
- 227685 atoms.
- All trajectory coordinates finite.
- Mean receptor CA RMSD: 1.223 Angstrom.
- Mean temperature: 300.31 K.
- Preliminary structural and numerical QC: PASS.

## 4. Apo replica 3

- Production: 5 ns complete.
- Recovery from verified 50000-step binary checkpoint.
- Original DCD had 14 frames through 140 ps.
- Original XML/checkpoint state was at step 50000.
- Original files preserved and backed up.
- Recovery performed in an independent output directory.
- Recovery DCD: 490 frames.
- Final merged DCD: 500 frames, 10-5000 ps.
- Final merged CSV: 500 records, steps 5000-2500000.
- All merged coordinates match their source frames exactly.
- Maximum source-versus-merged coordinate difference: 0.
- Recovery-boundary adjacent CA RMSD: 0.566 Angstrom.
- Mean receptor CA RMSD: 1.181 Angstrom.
- Mean recovery-stage temperature: 300.32 K.
- Preliminary structural and numerical QC: PASS.

Canonical merged trajectory:
project/results/pacer_dc_recovery_merged_v01/apo/replica_03/trajectory_corrected.dcd

The recovery-stage thermodynamic statistics exclude
the original first 10 records.

## 5. Current limitations

- Five-nanosecond pilots do not prove conformational convergence.
- Replica 2 and 3 receptor RMSDs were computed relative
  to their respective first frames.
- The apo replica 3 recovery boundary passed CA RMSD
  continuity checks, not exhaustive PBC or per-atom checks.
- Remaining compound110 contexts require additional replicas.
- Four-molecule functional comparison is incomplete.
- LY2119620 re-audit compatibility remains to be confirmed.
- Training gate remains CLOSED.

## 6. Next execution steps

1. Freeze and commit small apo replica 2/3 QC manifests.
2. Confirm Windows compound110 re-audit assets.
3. Inventory remaining compound110 production start states.
4. Complete matched four-context, three-replica 5 ns pilots.
5. Audit cross-replica and cross-window differential stability.
6. Extend to additional predefined functional controls.

Large DCD files, checkpoints, membrane systems,
and temporary recovery files must not be committed.

Windows is the sole Git authority.

## 2026-09-23 最新里程碑：compound110 replicas 2 and 3

### Git milestone
Latest pushed commit: 8d51112. Branch: audit/oneprot-pacer-dc-g0. Two candidate_probe QC manifests committed and pushed.

### Production and QC
compound110 + ACh replicas 2 and 3: both completed 5 ns CUDA production, 500 frames each, 227742 atoms, 10 ps interval. Trajectory integrity, preliminary thermodynamic QC, receptor CA RMSD and centroid PBC audits PASS. Mean receptor CA RMSD: replica 2 = 1.102387 A; replica 3 = 1.200738 A. Mean temperature: 300.243 K and 300.277 K.

### Cross-replica observations
Mean compound110 RMSD: 4.379 A and 3.836 A. Final initial-pocket contacts: 26 and 28. Final core-centroid distances: 4.656 A and 4.485 A. Both replicas retain core contacts and lose contacts with GLY65, SER232, ILE234 and PRO235. THR148, VAL149 and PRO150 contact gains occur in replica 2 but not replica 3. A common final binding pose or biological mechanism is NOT established.

### Four-context progress
compound110 + ACh: replica 1 = 1 ns; replicas 2 and 3 = 5 ns + preliminary QC.
compound110 only: replica 1 = 1 ns; replicas 2 and 3 not started.
ACh only: replica 1 = 1 ns; replicas 2 and 3 not started.
apo: replica 1 = 1 ns; replicas 2 and 3 = 5 ns + preliminary QC.

### Next actions
1. Start probe_only replica 2 and 3 production after checking existing files.
2. Complete compound110 candidate_no_probe replica 2 and 3.
3. Audit replica 1 production status and complete matched 5 ns windows.
4. Perform matched four-context atom14 and frozen OneProt-MD analysis.
5. Expand the functional control matrix after QC.

### Limitations
Independent replicas share a common starting structure. Five nanoseconds does not establish conformational convergence. Simulation E-chain residue indices are not original CHRM4 sequence indices. Large DCD files, checkpoints and membrane systems remain local, not on GitHub. TRAINING_GATE = CLOSED.

## 2026-09-25 Progress Update

### Completed milestones
- compound110 R2 and R3: all four contexts completed 5 ns production, preliminary QC and provenance audit.
- R3 candidate-only QC archived in commit a189942.
- Extracted and verified 40 raw trajectory windows: 2 replicas x 5 windows x 4 contexts.
- All 90 transferred raw, topology and manifest files passed WSL/Windows SHA256 checks.
- R3 apo extraction uses the corrected recovery trajectory with verified provenance.
- Completed 40 frozen OneProt-MD embeddings (1024 dimensions) and 10 matched dPAM vectors.
- Added reusable Windows embedding batch runner and an apo trajectory override to the extractor.

### dPAM stability results
- Ten complete differential units and 25 differential pair comparisons.
- Matched R2/R3 window cosine similarities: 0.502745, 0.712310, -0.353644, 0.262047, 0.085430.
- R2 mean dPAM L2: 0.223579; R3 mean dPAM L2: 0.259751.
- Mean dPAM cross-replica cosine: 0.339248; L2 distance: 0.279383.
- Cancellation ratios: R2 0.353899; R3 0.399608.
- Window-level and cross-replica directional stability have not been established.

### Current decision
- TRAINING_GATE = CLOSED. No classifier training or fine-tuning.
- Descriptive results do not establish PAM classification or biological mechanism.
- Pause further computation pending the afternoon project meeting.
- Next scientific steps will be decided after discussion.


## 2026-09-25 统一QC、PCA与 dPAM 方向敏感性更新

### 一、compound110 三 replica × 四 context 统一 QC（共同 1 ns 窗口）

已对 compound110 的 12 组轨迹（4 context × 3 replica）完成统一 QC 与共同降维预处理。

共同分析范围与结论如下：

- 统一时间范围：每条轨迹取前 100 帧（约 1 ns，10 ps/frame）；
- 总帧数：1200；
- 统一受体选择：`chainID E and protein and name CA`；
- 共同受体映射检查通过：270 个残基、270 个 Cα、原子顺序一致；
- 全部 1200 帧坐标有限性检查通过；
- 12 条轨迹逐帧连续性检查通过；
- R3 apo 在共同 1 ns 分析中使用修复轨迹 `trajectory_corrected.dcd`。

当前结论：compound110 的三 replica 均可在统一 1 ns 范围内进行共同 PCA 与轨迹可视化分析。

### 二、共同 PCA 结果（compound110，12 组轨迹）

对 1200 帧共同对齐后的受体 Cα 坐标进行了确定性 PCA（full SVD）。

解释方差：

- PC1 = 0.286183
- PC2 = 0.256741
- PC3 = 0.184046
- PC1 + PC2 = 0.542924
- PC1 + PC2 + PC3 = 0.726970

阶段性观察：

- 二维 PCA（PC1/PC2）中，四个 context 与三个 replica 均存在分布差异；
- 12 面板时间演化图显示，部分轨迹存在明显时间漂移；
- R3 apo 在 PC3 方向上漂移最大，但前 1 ns 内未观察到异常逐帧跳变；
- R3 的存在会显著影响二维 PCA 子空间，但前三维子空间整体仍较稳定。

留一 replica 分析（相对于完整 1200 帧 PCA）：

- 排除 R1：PC1–PC2 overlap = 0.980830；PC1–PC3 overlap = 0.986651；
- 排除 R2：PC1–PC2 overlap = 0.973557；PC1–PC3 overlap = 0.989586；
- 排除 R3：PC1–PC2 overlap = 0.842806；PC1–PC3 overlap = 0.950084。

进一步的 R3 留一轨迹分析显示：

- 排除 R3 candidate + ACh：top2 overlap = 0.996180；top3 overlap = 0.998083；
- 排除 R3 candidate only：top2 overlap = 0.935915；top3 overlap = 0.994576；
- 排除 R3 ACh-only：top2 overlap = 0.567767；top3 overlap = 0.986661；
- 排除 R3 apo：top2 overlap = 0.962165；top3 overlap = 0.988111。

组件旋转诊断表明：排除 R3 ACh-only 后，二维 PCA 的变化主要来源于前三个主成分内部的旋转，尤其表现为原 PC2/PC3 的显著重排，而不是前三维子空间整体失稳。

当前结论：二维 PCA 图对个别轨迹（尤其 R3 ACh-only）较敏感；前三维子空间整体仍保持较高一致性。二维图适合做可视化提示，但不宜单独作为构象稳定性的定论依据。

### 三、R2/R3 五窗口 dPAM 方向敏感性审计

已基于现有冻结 OneProt-MD embedding，对 compound110 的 R2/R3 五个匹配窗口（W0–W4）完成 dPAM 方向敏感性审计。

定义：

- dPAM = candidate_probe - probe_only
- 跨 replica 变化定义为：Δ = R3 - R2

对每个窗口均验证：

- dPAM 数学一致性残差 < 1e-5；
- 差分分解平方残差 < 1e-5。

五窗口主要结果如下：

| Window | CA cosine | A cosine | dPAM cosine | dPAM angle (deg) | cos(ΔzCA, ΔzA) | cross term |
| --- | --- | --- | --- | --- | --- | --- |
| W0 | 0.895527 | 0.847766 | 0.502744 | 59.818 | -0.310129 | 0.156445 |
| W1 | 0.915691 | 0.918110 | 0.712309 | 44.577 | 0.169672 | -0.056393 |
| W2 | 0.832024 | 0.868390 | -0.353644 | 110.710 | -0.203894 | 0.121264 |
| W3 | 0.886168 | 0.614202 | 0.262047 | 74.808 | 0.476550 | -0.399468 |
| W4 | 0.882035 | 0.785358 | 0.085430 | 85.099 | -0.008958 | 0.005702 |

阶段性结论：

1. 原始 embedding 的跨 replica 相似度较高，并不保证 dPAM 方向一致；
2. 五个窗口中，两个原始变化向量既可能相互增强（W0/W2）、也可能相互抵消（W1/W3），W4 近乎正交；
3. W2 的 dPAM cosine 为负（-0.353644），显示最明显的方向反转；
4. W3 虽然存在较强抵消（cross term = -0.399468），但 dPAM 方向仍不稳定；
5. 当前结果支持“差分几何关系是 dPAM 方向波动的重要因素”，但尚不能区分其根源来自 MD 采样、编码器表示还是差分操作对两者相对方向的敏感性。

### 四、当前阶段状态

- 统一 QC：已完成（compound110，三 replica，四 context，共同 1 ns）；
- 轨迹降维可视化：已完成第一轮 PCA、12 面板时间图、PC1/PC3 诊断与稳定性分析；
- dPAM 方向敏感性：已完成 R2/R3 五窗口初步审计；
- TRAINING_GATE：继续保持 `CLOSED`。

### 五、已新增的本地结果与脚本

新增脚本：

- `project/pacer_dc_training/plot_dpam_direction_summary_R2_R3.py`
- `project/pacer_dc_training/audit_dpam_direction_R2_R3.py`

新增结果：

- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_decomposition_R2_R3_v01.json`
- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_decomposition_R2_R3_v01.csv`
- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_summary_R2_R3_v01.png`
- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_summary_R2_R3_v01.svg`
- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_summary_R2_R3_v01.json`

## 2026-09-26：G0 四上下文交互差分审计

### 一、研究目的

在现有 compound110 的 R2/R3 五窗口冻结 OneProt-MD embedding 基础上，检验从 dPAM 中扣除候选物单独作用后，四上下文交互差分 dINT 的跨 replica 稳定性是否改善。

本阶段不新增 MD、不重新运行 encoder，不修改已有 dPAM 历史结果。

### 二、差分定义

四个匹配 context：

- z_CA：candidate + ACh
- z_C：candidate only
- z_A：ACh only
- z_0：apo

定义：

- dPAM = z_CA - z_A
- dAGO = z_C - z_0
- dINT = dPAM - dAGO = z_CA - z_A - z_C + z_0

dINT 是 embedding 空间中的 candidate×ACh 交互差分，不直接等同于药理学协同效应。

### 三、输入和质量检查

分析对象为 compound110 的 R2/R3，每个 replica 包含五个匹配窗口，每个窗口包含四个 context。

共检查 40 个冻结 OneProt-MD embedding：

- 文件完整性：40/40，通过；
- 向量形状：全部为 (1024,)；
- 数值有限性：全部通过；
- L2 范数范围：0.9999996423721313 至 1.0000004768371582；
- 归一化检查：全部通过；
- dPAM 历史一致性：五窗口全部通过；
- dINT 差分数学一致性：全部通过。

### 四、五窗口跨 replica 比较

下表列出 R2/R3 匹配窗口差分向量的 cosine similarity。

| Window | dPAM cosine | dAGO cosine | dINT cosine |
|---|---:|---:|---:|
| W0 | 0.502744 | 0.531015 | 0.675959 |
| W1 | 0.712309 | -0.038464 | 0.362566 |
| W2 | -0.353644 | 0.222379 | 0.386707 |
| W3 | 0.262047 | 0.258136 | -0.426019 |
| W4 | 0.085430 | -0.456788 | -0.390655 |

相对于 dPAM，dINT 在 W0、W2 的方向一致性改善，在 W1、W3、W4 恶化。

尤其是 W2，dPAM cosine 为 -0.353644，而 dINT cosine 为 0.386707；相反，W3 和 W4 的 dINT 出现负 cosine。

因此，扣除 dAGO 没有产生跨窗口一致的稳定性改善。

### 五、五窗口平均向量比较

先在各 replica 内分别计算五个窗口的平均差分向量，再比较 R2/R3 平均向量。

| 指标 | dPAM | dAGO | dINT |
|---|---:|---:|---:|
| Cosine | 0.339247 | 0.330114 | 0.249590 |
| Angle (deg) | 70.169 | 70.724 | 75.547 |
| L2 distance | 0.279383 | 0.398976 | 0.430014 |
| R2 norm | 0.223579 | 0.241554 | 0.295126 |
| R3 norm | 0.259751 | 0.407143 | 0.394969 |

平均向量比较中，dINT 的 cosine 低于 dPAM，L2 distance 高于 dPAM。

上述指标不是五个窗口 cosine 的算术平均，也不构成统计显著性检验。

### 六、阶段性科学结论

1. 四上下文交互差分 dINT 已成功实现，且计算与历史 dPAM 结果一致。
2. dINT 未表现出一致优于 dPAM 的跨 replica 方向稳定性。
3. 二阶差分可能进一步累积四条独立轨迹的采样和表示波动，但当前分析不足以确定因果来源。
4. 当前结果仅涉及 compound110、两个 replica 和五个匹配窗口；窗口不是独立生物学重复。
5. 尚不能区分 MD sampling、隐藏表示、全局池化、MLP 投影和 L2 归一化的具体贡献。

### 七、脚本与结果

新增脚本：

- `project/pacer_dc_training/audit_dint_R2_R3.py`
- `project/pacer_dc_training/plot_dint_comparison_R2_R3.py`

结果目录：

`project/results/pacer_dc_four_context_v01/compound110/interaction_audit_R2_R3_v01/`

新增文件：

- `interaction_comparison_R2_R3_v01.csv`
- `interaction_comparison_R2_R3_v01.json`
- `dint_comparison_R2_R3_v01.png`
- `dint_comparison_R2_R3_v01.svg`

CSV 包含 15 条窗口比较记录。JSON 包含上述记录、3 条平均向量比较记录及全部 40 个输入 embedding 的 SHA256。

正式英文汇总图已生成，并通过人工视觉 QC。

### 八、下一阶段与训练门

G0 的计算、结果导出及绘图已经完成，GitHub 归档待完成。

下一阶段为 G1：审计 OneProt-MD 的实际 forward 路径，并在不改变现有正式输出的前提下，研究四层表示的提取：

1. MDGen 残基级时间隐藏表示；
2. 21 维全局池化表示；
3. MLP 投影后、L2 归一化前的 1024 维表示；
4. 当前归一化的 1024 维表示。

具体 tensor 语义与形状必须以实际代码和 checkpoint 审计为准。

TRAINING_GATE = CLOSED。不得根据本阶段结果宣称已验证 PAM 分类能力，也不启动分类器训练。
