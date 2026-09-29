# PACER-DC / Geom2Vec Encoder 与 PACER-FKG v02 优化审计报告

**阶段性方法、结果、冻结状态与复现边界（截至 2026-09-29）**

生成目的：同学独立方法/结果审计。

## 1. 执行摘要

本项目完成了一条从表示层到 PACER-FKG 长程验证的受控优化路线。优化不是围绕最终 ΔPAM/ΔAGO/ΔINT 结果进行反向调参，而是先在短程历史数据上通过结构可访问性、数值稳定性、维度控制和 R2/R3 描述性复现逐级定位信息损失点，再冻结 encoder candidate，最后在 600 ns 长程 MD 上建立独立的 PACER-FKG v02 数值状态并进行第一次正式 outcome evaluation。

目前得到的核心结论是：

1. 原 production Geom2Vec 路线的主要信息损失点至少包括：ViSNet final readout、all-heavy-atom residue mean、以及简单 temporal averaging。
2. 最有证据支持的优化表示为：ViSNet 中间状态 C1 → backbone/side-chain 分离 → BS256。
3. 时间维度上，STATE_MOTION = [μZ, RMSΔ] 比单独 SIGNED_DRIFT = μΔ 更能在长程 MD 中形成跨 frozen application replicas（R1/R3）保持的 PACER-FKG directional structure。
4. PACER-FKG v02 与 v01 在输入维度、normalization、RFF/kernel numerical state、branch schema 上不兼容，因此 v02 已完全独立实现，没有覆盖或回写任何 v01 历史结果。
5. 当前结果支持“找到了一条 representation-level 优化路径”，但尚不能证明 v02 在药理学 efficacy 或整体方法学上全面优于 v01；正式 v01-v02 apples-to-apples comparative audit 尚未执行。

## 2. 审计边界与不变原则

本阶段所有实验均遵守以下边界：

- 长程 MD 在 encoder candidate 冻结前不用于 encoder 调优。
- 不以 ΔPAM、ΔAGO、ΔINT、历史成功区域、direction cosine 或 kernel_u2 magnitude gate 作为 encoder/readout 选择目标。
- 历史 kernel_u2 magnitude gate 已撤回，在新方法中未恢复。
- R3 在短程开发阶段仅作为已检查过的描述性工程 holdout，不作为新鲜生物学确认。
- PACER-FKG v01 的历史代码、calibration state 和结果目录均按 SHA256 保护；v02 使用独立目录、独立 numerical state 和独立 branch schema。
- v02 中 R2 仅用于 calibration；R1 与 R3 是 frozen application replicas。
- v02 Phase 3 之后不允许根据 outcome 回头修改 C1、B/S pooling、block size、branch 定义、normalization、bandwidth、RFF seed、graph 或 contrast formula。

## 3. 原始 production 路线与已识别的信息瓶颈

历史 production Geom2Vec 路线可概括为：

receptor coordinates
→ atom14
→ heavy-atom point cloud
→ frozen ViSNet
→ concat(scalar x, vector-channel norms ||v||)
→ all-heavy-atom mean per residue
→ residue feature 128D
→ downstream PACER-FKG v01

初始审计确认三个主要的不可逆压缩点：

- vector orientation 被 ||v|| 变换消除；
- atom-level heterogeneity 被 residue mean 消除；
- temporal averaging 会进一步消除时间顺序与短时动态。

此外，历史 G2-B normalization 与 PACER-FKG downstream 对 128D 输入存在显式/隐式绑定，因此优化后表示不能合法直接塞回 v01。

## 4. Encoder / readout 优化实验 A-H

### Experiment A — local/global decomposition

目标：检查现有 residue embedding 是否能在不改变信息的情况下分解为 local residual 与 global component。

定义：
H(t,i,:) = residue embedding
G_t = residue mean
R_ti = H_ti - G_t

结果：
- 40/40 历史 R2/R3 embedding 文件通过 hash、shape、dtype、frame ID、sequence 和 finite 检查。
- R + G → H 最大绝对重建误差约 1.11e-16。
- regional mean 重建最大误差约 2.22e-16。
- 109 个 region memberships，69 个 unique residues。

结论：local/global decomposition 可近乎无损保留已有 residue embedding 信息，但这本身不证明 G 或 R 的生物学优越性。

### Experiment B / B2 — region pooling 与 masked reconstruction

B 中，attention 在原始 masked-R 任务上相对 mean 有轻微改善，但简单的 R2 residue-specific temporal mean predictor 明显更强。进一步审计发现约 85.24% 的 R2 target energy 是 residue-stationary。

B2 去除 stationary component，并修正 full-frame centering 对 masked target 的代数泄漏。结果：
- 35 tests passed；
- 9/9 exact replay/reload；
- R3 上 0/9 learned seed/comparator 结果优于 zero-D baseline；
- attention 优势消失。

结论：停止 region-level attention/pooling complexity。该负结果只说明这一固定 masked reconstruction objective 未证明额外 dynamic information，不代表 ViSNet 本身无信息。

### Experiment C — ViSNet intermediate-state audit

固定 atom readout，仅比较三个 frozen ViSNet state：
- C0：production final normalized state；
- C1：前 3 个 message-passing layers 后的 accumulated state；
- C2：final accumulated pre-normalization state。

Phase 0-2 先验证：
- checkpoint SHA256：b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417；
- 152/152 parameters、437,670/437,670 elements strict load；
- installed inference source 与 Geom2Vec commit 371d642ec1061664f16e49fcac702d07fc8d0b51 一致；
- 40 个历史 atom14 输入全部匹配 archived hashes；
- smoke 中 C0/C1/C2 均为 [4,270,128] 且有限。

固定 torsion probe 结果（R3 combined backbone angular MAE）：
- C0：16.25°
- C1：15.39°
- C2：15.44°

C1 相对 C0 改善约 5.28%；C2 改善约 4.97%。

结论：ViSNet intermediate / pre-normalization states 保留更多线性可访问的局部结构信息。C1 作为后续 candidate state。

### Experiment D — atom-to-residue B/S separation

固定 C1 atomic descriptor [x, ||v||]，比较：
- M：all-heavy-atom mean，128D；
- B：backbone mean，128D；
- S：side-chain mean，128D；
- BS：concat(B,S)，256D；
- MC：[M,0]，256D dimension control。

R3：
- C1-M backbone combined MAE：15.39°
- C1-BS：12.80°
- C1-M chi1 MAE：43.44°
- C1-BS：35.87°

相对 M/MC：
- backbone angular MAE 改善 16.87%；
- chi1 MAE 改善 17.43%。

改善在 residue positions 上广泛存在：R3 backbone 约 80.4%，chi1 约 76.2% 的位置改善。

结论：all-heavy mean 会混合并损失 backbone 与 side-chain 的互补信息。BS256 是目前最强的 information-rich residue representation。

### Experiment E — compact 128D fusion

尝试 AVG、STATIC gate、DYNAMIC gate 将 B/S 压回 128D。

R3 中：
- M：backbone 15.392°，chi1 43.437°
- STATIC median：14.907°，41.400°
- DYNAMIC median：14.832°，41.660°
- BS256 reference：12.795°，35.868°

STATIC/DYNAMIC 虽优于 M，但仅恢复约 10–36% 的 BS 优势，低于预注册 50% criterion。DYNAMIC reconstruction 更好，但 chi1 并未跨 seeds 稳定优于 STATIC。

结论：停止 learned compact fusion。BS256 的信息优势不应强行压回单个 128D 表示。

### Experiment F — lag-1 temporal information

固定 M/MC/BS，使用相邻 stored frames 的 ΔZ 预测 torsion Δsin/Δcos，并与 orderless pair mean 和 zero-dynamic baseline 比较。

R3 上 BS 相对 zero baseline：
- phi：9.74% improvement；
- psi：9.71%；
- chi1：3.64%。

BS 在预定义 residue 上广泛优于 MC。20-frame averaging 仅保留约 13.34% 的 BS temporal variance，即约 86.66% 被 block averaging 丢弃。

结论：BS256 不仅保留更多静态结构信息，也保留额外的 signed short-lag dynamic information；简单 temporal mean 是新的信息瓶颈。

### Experiment G — branch-aware window representation

固定 20-frame block，定义：
- μZ：block state；
- μΔ：signed endpoint drift；
- RMSΔ：short-lag motion amplitude。

结果显示：
- μΔ 最适合 signed-dynamic target，加入 static/RMS 通常不改善；
- μZ 与 RMSΔ 对 static / motion-amplitude target 具有互补性；
- BS 在 static、signed、RMS 三类 matched task 上均优于 dimension-matched MC。

因此冻结 branch-aware temporal object：
- STATE_MOTION512 = [μZ, RMSΔ]
- SIGNED_DRIFT256 = μΔ

### Experiment H — PACER-FKG compatibility audit

H 的正式分类为 H-C：现有 v01 的代码与 frozen numerical state 都是 representation-specific，必须建立新版本并做实现适配。

历史 v01 输入是单一 (100,270,128) residue_features；frozen normalization 是 128D；G2-B RFF 输入宽度为 n_residues×128；G2-C node mapping 也绑定 128D；且 v01 没有 named branch interface。

可继承的科学定义：
- four-context formulation；
- graph topology；
- diffusion formula；
- RBF/kernel family；
- R2-only calibration principle；
- fixed 20-frame rule。

必须新建/重校准：
- branch-specific normalization；
- bandwidth / RFF numerical state；
- schema / routing；
- serialization / reporting。

## 5. PACER-FKG v02：长程 MD 独立实现

### Phase 0 — 数据与历史结果保护

长程 MD：
4 systems × 3 replicas × 50 ns = 600 ns
共 12 trajectories，每条 1000 stored frames，50 ps spacing，约 50 ps 至 50,000 ps。

系统：
- apo
- probe_only
- compound110__candidate_no_probe
- compound110__candidate_probe

replica seeds：
- R1 = 27101
- R2 = 38201
- R3 = 49301

完整性：
- 41 个历史 PACER-FKG v01 result files 重新 hash，保持不变；
- frozen baseline commit 0a4b4f2 的 630 个 tracked files 建立保护清单；
- 60 个 server-recorded production files 全部本地 checksum 匹配；
- 12 条 trajectory 均完成 50 ns；
- 12,000-frame receptor-CA scan 未发现 PBC split、severe chain break 或 temporal wrapping event；
- receptor 为有序 270 residues，backbone/atom14 mapping 完整。

Phase 0 gate：PASS_WITH_DOCUMENTED_CAVEAT。
非阻塞 caveats：
- backup 本身无 topology，使用经过认证的 short-data minimized PDB；
- state.csv 与 DCD header 时间约差 4 ps；
- 某些非 receptor PDB records 缺 residue identifiers，但 receptor mapping 不受影响。

### Phase 1 — frozen C1-BS256 extraction

对 12 条长程 trajectory 执行：
DCD → authenticated minimized PDB → chainID E and protein → historical atom14 → frozen C1 → [x,||v||] → B128/S128 → BS256。

每条 trajectory 输出：
[1000,270,256]

full extraction PASS；verify PASS。

### Phase 2 — R2-only calibration 与 numerical freeze

每条 1000-frame trajectory 分为 50 个连续 20-frame blocks，共 600 blocks。

R2-only calibration population：
4 systems × 50 blocks = 200 calibration blocks。

冻结数值状态：
- STATE_MOTION bandwidth = 32.30188361260893
- STATE_MOTION seed = 272340
- SIGNED_DRIFT bandwidth = 29.29934899925964
- SIGNED_DRIFT seed = 272084

freeze manifest SHA256：
b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd

17 tests passed；R1/R3 明确未进入 fitting。冻结后 --run 会主动拒绝再次生成 numerical state。

### Phase 3 — formal frozen long-MD evaluation

R2 明确标记为 CALIBRATION_REPLICA；R1/R3 为 FROZEN_APPLICATION_REPLICAS。两个 branch 分开报告，不做 branch selection / combination / magnitude gate / p-value / new bootstrap。

STATE_MOTION 的 R1/R3 agreement：
- ΔPAM：6/9 regions
- ΔAGO：7/9
- ΔINT：7/9

最强的预定义区域信号集中在 compound110_extension：
- ΔAGO R1/R3 direction cosine = 0.780998
- ΔINT R1/R3 direction cosine = 0.619660

但 ΔPAM 在同一区域为 -0.520858，说明 PAM-related directional structure 并未普遍稳定。

SIGNED_DRIFT 的 R1/R3 agreement：
- ΔPAM：3/9
- ΔAGO：4/9
- ΔINT：2/9

局部仍有正向结果，例如：
- ΔAGO / compound110_extension = 0.555629
- ΔAGO / pam_contact_consensus = 0.579508
- ΔAGO / pam_contact_union = 0.419299
- ΔINT / compound110_extension = 0.424821

整体上 SIGNED_DRIFT 的长程 directional reproducibility 明显弱于 STATE_MOTION。

## 6. 当前最重要的科学结论

1. 已找到一条有连续证据支持的 representation-level 优化路径：

production C0-M128
→ C1 intermediate ViSNet
→ backbone/side-chain separated BS256
→ branch-aware temporal representation
→ PACER-FKG v02

2. 三个主要信息瓶颈被逐级定位：
- final ViSNet readout；
- all-heavy-atom residue mean；
- temporal averaging。

3. BS256 的优势不是由简单维度增加解释。MC=[M,0] dimension control 在结构 probe 上复制 M，而 BS 显著更优。

4. 将 BS256 强行压回 128D 会重新损失大量已恢复的信息；因此 compact learned fusion 被停止。

5. 短时 signed dynamics 确实存在并可被 μΔ 访问，但这不自动转化为稳定的长程 PACER directional signal。

6. 目前长程正式结果主要支持 STATE_MOTION = [μZ,RMSΔ]。它在多个预定义 ΔAGO/ΔINT region 上具有较好的 R1/R3 directional consistency，尤其 compound110_extension。

7. 当前尚不能宣称：
- v02 全面优于 v01；
- v02 已证明 compound110 的药理 efficacy；
- 50 个 block 是独立 biological replicates；
- 所有 Agreement=True 都代表强复制（部分 cosine 非常接近 0）。

8. 尚未执行 frozen v01-v02 end-to-end apples-to-apples comparative audit，因此方法学 superiority 仍是开放问题。

## 7. 复现与冻结锚点

关键 provenance：

- Encoder optimization branch:
  experiment/pacer-encoder-opt-v01

- Frozen encoder candidate commit:
  0a4b4f2
  message: experiment: freeze encoder optimization A-H candidate

- Frozen encoder tag:
  encoder-candidate-ah-frozen-20260929

- PACER-FKG v02 branch:
  experiment/pacer-fkg-v02-longmd-v01

- Geom2Vec pinned source commit:
  371d642ec1061664f16e49fcac702d07fc8d0b51

- ViSNet checkpoint SHA256:
  b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417

- v02 Phase 2 freeze manifest SHA256:
  b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd

- Runtime used during encoder verification:
  Python 3.11.16
  PyTorch 2.6.0+cu126
  PyG 2.7.0
  NumPy 2.4.6
  CUDA 12.6
  torch-cluster 1.6.3+pt26cu126
  RTX 4070 Laptop GPU

审计时应优先验证以上 commit/tag/hash，再读取结果报告。

## 8. 建议同学重点审计的问题

建议独立审计按以下顺序进行：

1. Provenance
   - commit/tag 是否存在；
   - checkpoint/source commit 是否匹配；
   - v02 freeze manifest hash 是否匹配。

2. 数据完整性
   - 12 条长程 trajectory checksum；
   - 270-residue receptor mapping；
   - Phase 1 BS256 cache manifests。

3. 方法独立性
   - encoder candidate 是否确实在长程 PACER outcomes 之前冻结；
   - R1/R3 是否未参与 Phase 2 calibration；
   - Phase 3 是否未重新拟合 normalization/bandwidth/RFF。

4. Experiment D dimension control
   - MC=[M,0] 是否严格复制 M；
   - BS 优势是否仍存在于相同 256D dimensionality。

5. Experiment F/G temporal controls
   - ΔZ 是否优于 zero 与 orderless controls；
   - μΔ endpoint identity；
   - block boundary 是否没有跨 trajectory/window。

6. Phase 3
   - direction cosine 是否从完整 signed vectors 计算；
   - Agreement=True 是否仅表示方向符号/规则一致，而非强效应；
   - R2 是否始终标记为 calibration replica；
   - STATE_MOTION 与 SIGNED_DRIFT 是否从未组合成总分。

7. 历史保护
   - 41 个 v01 files 在 v02 全流程中 hash 是否不变；
   - kernel_u2 是否确实未恢复。

8. 开放问题
   - frozen v01-v02 end-to-end comparative audit 尚未执行；
   - historical v01 RFF numerical state 的 exact replay 能力仍需单独确认。

## 9. 主要结果文件索引

仓库内建议审计的核心文件/目录：

Encoder optimization：
- project/encoder_readout_v01/
- project/encoder_intermediate_v01/
- project/encoder_atom_readout_v01/
- project/encoder_bs_fusion_v01/
- project/encoder_temporal_v01/
- project/encoder_window_v01/
- project/encoder_candidate_v01/

PACER-FKG v02：
- project/pacer_fkg_v02/
- project/results/pacer_fkg_v02_longmd_v01/

重点结果：
- Phase 0 provenance / gate reports
- Phase 1 BS256 extraction reports/manifests
- calibration/V02_FREEZE_MANIFEST.json
- phase3/REPORT_PHASE3.md
- phase3 STATE_MOTION / SIGNED_DRIFT result JSONs
- phase3 block-level manifests
- integrity reports

注意：大体积 cache/results 按设计保持 Git-ignored；Git 仓库主要保存实现、方法定义、freeze receipts 与 provenance anchors。审计者若只拿到 GitHub 仓库而没有本地 result/cache，应同时取得对应结果目录或其压缩归档。

## 10. 当前停止点与下一步

当前建议停止继续调参。

已完成：
- representation optimization A-H；
- PACER-FKG v02 长程 Phase 0-3；
- frozen outcome evaluation。

尚未完成：
- historical PACER-FKG v01 exact replay qualification；
- frozen v01-v02 600 ns end-to-end comparative audit。

下一步若继续，应先确认 v01 numerical state（尤其历史 RFF）能否 exact replay；若不能，不应重新拟合后仍称为“历史 v01”，而应另立 matched-calibration comparator。
