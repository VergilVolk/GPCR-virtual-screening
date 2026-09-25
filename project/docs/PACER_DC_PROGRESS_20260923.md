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

