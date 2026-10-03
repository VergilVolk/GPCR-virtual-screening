# Stage3 PACER-XR protocol deviation 事故记录

记录日期：2026-10-03；事故分支：`incident/pacer-xr-protocol-deviation-20261003`。
审计基线：`b78148e9c283aa834118f99ac5a4feb613885e28`。
范围：PACER-PROSPECTIVE-20261002 的 Stage1–Stage4 最终候选链；只读取已有代码、配置、manifest、receipt、provenance 和结果，新增本记录及归档索引。不重新执行科学计算。

## 1. 已确认的 protocol deviation

> 本次 Stage1–Stage4 最终候选链在 Stage3 没有按照既定的六路 PACER-XR 方案实行，而是另外采用了 Vina × 10 个 GaMD 受体构象的 ensemble docking 方案。

这是 **protocol deviation（协议偏离）**。十构象 GaMD ensemble 不能冒充六路 PACER-XR；引擎、受体来源、score definition、归一化和决策规则均须区分。Vina 任务完成不等于六路 PACER-XR 完成，也不等于 Vina_PDB / Vina_BEmin / Vina_BEavg 三路完成。

既定方案可核对 [PACER_XR_CASCADE_METHOD.md](PACER_XR_CASCADE_METHOD.md) 和 [apply_pacer_xr_candidates.py](../scripts/apply_pacer_xr_candidates.py)。后者明确要求六路完整、缺失即失败、禁止插补，以冻结官方 M4 分布计算各路 empirical mid-rank percentile，然后作六路等权融合及 Glide_BEmin 官方 top-1% cascade。

| 既定通道 | 本次最终候选链中的执行证据 |
|---|---|
| Glide_PDB | 未见本批候选的协议匹配分数或完成回执 |
| Glide_BEmin | 未见本批候选的协议匹配分数或完成回执 |
| Glide_BEavg | 未见本批候选的协议匹配分数或完成回执 |
| Vina_PDB | 实际十个 GaMD cluster 不能充作这一通道 |
| Vina_BEmin | 十构象最低 Vina 描述统计不能充作协议匹配 BEmin |
| Vina_BEavg | 十构象平均 Vina 描述统计不能充作协议匹配 BEavg |

这里的“既定六路”是本事故要求核对的协议基准，并非声称已发现本次生成批次的六路完成记录。历史公开 benchmark 存在六路数据和 PACER-XR 结果，不能移作这批新候选的 prospective evidence。

## 2. 实际最终候选链及结果

```text
M4-safe DrugCLIP
→ Vina × 10 GaMD receptor conformations
→ structural gate / diversity selection
→ PACER0010 / PACER0073 / PACER0027
→ four-context MD
→ PACER-FKG
```

以上是最终选择所用证据的逻辑链，箭头不额外声称两类 Stage3 作业的严格启动时间先后。以最终冻结 manifest 的输入哈希、身份合同和 selection rule 确认实际依赖；旧 handoff 中不同的摘要排列不覆盖这些证据。

| 阶段 | 本次实际结果 | 核心证据 |
|---|---|---|
| Stage1 | 12 个种子；BRICS fragment generation，17,146 次 join attempts，5,000 个唯一原始生成分子，2,605 通过，2,395 未通过生成规则；本批不是 LSTM 生成 | 外部执行目录 U 的 STAGE1_GENERATION_RECEIPT_v01.json、STAGE1_QC_v01.json、生成 CSV |
| Stage2 | 输入 2,605；留下 200（local 160 / exploratory 40），排除 2,405；排除的独占原因 13 exact-known、931 domain、15 risk、1,446 quota | U 的 STAGE2_FILTERING_RECEIPT_v01.json、STAGE2_QC_v01.json；predock_portfolio.csv |
| Stage3 | 同一 200 分子完成 M4-safe scoring；Vina 200 × 1 本批实际状态 × 10 cluster = 2,000 完成任务；200 个全部结构 gate PASS；140 个 Murcko scaffold 代表，最终选 3 个 | D 的 COMPLETION_RECEIPT.json；K 的 ledger/metadata/scores；S 的 selection manifest / provenance / validation |
| Stage4 | 3 候选 × 4 上下文 × 3 replicas = 36 条 10 ns 轨迹，合计 360 ns；36/36 输入被接受；PACER-FKG 最终冻结评估完成 | 当前 tree 的 MD master/job matrix；FKG 输入认证和三个阶段冻结 receipts、结果及科学报告 |

Stage1 的 2,395 是生成规则未通过数，不是已交付 SMILES 解析失败数。Stage2 重叠 risk 计数 55 不应与独占原因相加。旧训练、旧 benchmark、早期双模型 shortlist 均不替代最终 200 分子的 corrected scoring / Stage3D freeze。

### 2.1 DrugCLIP 与 docking 的实际处理

Corrected DrugCLIP 2023 backbone 重用本批已保存的 representation，并应用 M4-held-out GPCR-LOTO 投影；2026 使用独立 backbone / 三个晶体口袋，作为 second opinion。M4 的最终 routed 排名采用 2023 M4-safe 路由，不使用历史 DrugCLIP 候选排名。completion receipt 的 rows / graph_overlap 均为 200，missing / extra / duplicates 均为 0。该 receipt 自身的 stage3d_executed=false 仅说明 scoring 完成时的快照，后续 Stage3D manifest 已冻结选择。

Docking metadata 记录 Vina 1.2.7、cluster 0–9、exhaustiveness=8、num_modes=9，box 为各 cluster 的 MK-97 centroid 周围 30 Å 立方体。十构象分数保留在 framewise/statewise 表；best / median / mean / std 是 ensemble 描述统计，没有六路 empirical-rank fusion 或 PACER-XR cascade 输入。

结构门槛：cluster coverage ≥0.8、median native pocket coverage ≥0.5、median contacted residues ≥8、maximum centroid outlier distance ≤20 Å。所有 200 个通过；随后按冻结 M4-safe final_rank，每个 Murcko scaffold 保留一个代表，选出前三个 diverse candidates。不是以最低 Vina 能量筛出三候选。代表 pose 选取规则是接近 ensemble-median centroid，再按 native pocket coverage、Vina 和 cluster 决胜。

### 2.2 三个最终候选的冻结 lineage

| 候选 | M4-safe final_rank | diverse_rank | MD source cluster | 冻结 pose |
|---|---:|---:|---:|---|
| PACER0010 | 1 | 1 | 9 | PACER0010_c09_s00.pdbqt |
| PACER0073 | 2 | 2 | 0 | PACER0073_c00_s00.pdbqt |
| PACER0027 | 4 | 3 | 4 | PACER0027_c04_s00.pdbqt |

Stage2 portfolio SHA-256：`0dd7689860879eaae8abe806699a05a1c261a439e7355b99628e982abe46f52b`。
M4-safe ranking SHA-256：`bda22c9354a9ebf4a263ce92ccf1257c8b27e8e3e747404b0f4fa97c3fb34d89`。
Stage3D selected CSV SHA-256：`5672be147ee2d6b1cc4ece4d0b95b0e90512fe6b54dfacb88de976353a507dd3`。
MD pose manifest SHA-256：`8ef3c4be8ff3800543d09e3b37c11331671e961c7bd3f29cc70850b787e311a6`。

| 候选 | 冻结 pose SHA-256 |
|---|---|
| PACER0010 | 836a68d0cced23dc6e345a55d43e66119c27f26b9d49803635f97aed56e1e9de |
| PACER0073 | 3c503fc44a107cca2612dc71d9b0d529577bb5447eb2ea329dbe3adc32bcc37a |
| PACER0027 | 88cd9f7b228a25f65e82c9d82a49f65653cd67eb2abe53c3280146cc001d737e |

PACER0027 从 neutral parent 映射到冻结 +2 microstate；Stage4 保留该实际 docking state。没有更换候选、重新选 pose 或补算得分。

## 3. 未执行 PACER-XR 的证据与原因边界

1. XR 代码存在，接口可建立空六路模板，但接口存在不能证明执行。代码默认 roster 是旧 final_candidate_hypotheses.csv；该历史 roster 不构成本次三候选的完成回执。
2. 历史方法文档写明本机缺少 Schrödinger 许可，旧 24 候选尚无 Glide 分数；release manifest 也列出商业依赖。这说明已记录的前置限制，**不能仅据此断言 2026-10-02 每个执行环境的实时许可状态**。
3. 本批最终 Stage3D manifest 的 docking_inputs 指向十构象 Vina metadata / ledger / framewise / statewise scores，selection_rule 指向 structural gate + M4-safe final_rank + Murcko diversity；未指向六路 score table、XR application audit 或 cascade rank。因此实际 pipeline 使用了另一方案，未以既定六路作 prospective selection。
4. 已读记录没有足够证据判定是谁批准了替代、何时批准、是否有正式豁免，也不能虚构唯一运行时故障作为原因。已确认的是六路前置条件未在该候选链中完成、替代方案参与选择，以及相应 protocol deviation。

## 4. Stage4 结果与 claim boundary 保持不变

Stage4 使用各候选自身 cluster-matched A/P/C/CP（apo、ACh alone、candidate alone、candidate+ACh），保留全部 R1/R2/R3。最终回执状态是 `STAGE4_FULL_FROZEN_EVALUATION_COMPLETE`，不是药理确认。

[最终科学报告](../results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md) 与 [最终交接](PACER_STAGE4_FINAL_HANDOFF_v01.md) 仍是科学结论来源。本记录只引用其边界，不修改结论：PACER0010 有限的 probe/intrinsic dynamic pattern separation；PACER0027 的 STATE_MOTION candidate-alone perturbation 较一致；PACER0073 在 PAM-related regions 中 candidate-alone directional response 较稳定。没有候选在多数关键区域跨两 branch 展现稳健一致的 Delta_INT，不能确认 cooperative PAM mechanism。

> Prospective frozen PACER-FKG dynamic differential evaluation only.
> No PAM/ago-PAM labels, probability, efficacy or independent-block inference.

Delta_AGO / Delta_PAM / Delta_INT 是冻结方法术语；不能把计算动态特征写成实验确认的 PAM / ago-PAM 标签。Wet-lab functional validation 未实施。历史 calibration / training / benchmark 保持历史身份，不作为 Stage1–Stage4 新结果。

冲突处理：早期 builder 的 NOT_STARTED / NOT_RUN 是当时快照，最终 completion receipt 及输入认证决定最终状态。源 manifest 的 trajectory_ps=50 与实际 raw 10 ps 不同，按既有 temporal metadata note 说明，不修改源字段。PACER0010/CP/R3 的旧 running/9.9 ns progress 按既有认证特例接受，原 progress 不改。

## 5. 历史记录和后补分析的规则

- 不能修改历史记录把六路 PACER-XR 写成“已经执行”。
- 不能把十构象 GaMD ensemble 重新命名成六路 PACER-XR，也不能把十个 Vina 值压缩后冒称 Vina_PDB / Vina_BEmin / Vina_BEavg 已完成。
- 如果之后给现有三个候选补跑六路，只能叫 **post hoc supplementary PACER-XR analysis**。
- 不能把后补结果倒写成原 prospective selection evidence；原三候选选择时间、原证据、原分数和原冻结 manifest 必须保留。
- 后补工作若获授权，须单独建立有时间戳的输入身份、六路完整性、score protocol、receipt 和 provenance，并明确候选已预先选定。本次不启动任何补跑。
- 本次不 rerun docking / DrugCLIP / MD / PACER-FKG，不训练、不重新计算、不改 frozen results、不改 Stage4 scientific conclusions、不删除历史证据。

## 6. 证据位置与保存缺口

当前 incident tree 已包含 Stage4 轻量冻结结果、实现和正式报告。下列 U/D/K/S 属于**外部证据位置**，不等同于已经归档进本分支；本次仅读取、列索引，不复制或合并。

| 缩写 | 已读取位置 |
|---|---|
| U | C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01 |
| Stage1 / Stage2 输出 | C:/projects/GPCR-virtual-screening-final/project/results/generated；C:/projects/GPCR-virtual-screening-final/project/results/pacer_candidates_v01 |
| D | C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01 |
| K | C:/projects/GPCR-virtual-screening/project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01 |
| S | C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01 |
| Raw MD | C:/projects/PACER_STAGE4_MD_backup |

[完整归档索引](RECENT_EXECUTION_ARCHIVE_20261002_20261003.md) 记录两个指定 git log 的完整输出、全部 22 个提交和逐文件归属、339 个外部文本证据位置及哈希。首次基线实际 16/22 提交是 incident HEAD 祖先，另外 6 个位于其他 refs（原版本 15/7 为汇总错误，见新增 completion 节）；部分本批执行产物只在外部工作目录。故“最近两天所有代码/日志均已完整保存在该分支”首次核验结果为 **NO**；本次补档后原 22/22 已可达，轻量证据已复制，大型资产显式索引。不能改写历史或偷偷合并科学结果。

本事故收尾仅增加两份文档，是否修改任何科学结果：**NO**。

## Evidence preservation completion

补档日期：2026-10-03；工作分支：`incident/pacer-xr-protocol-deviation-20261003`。
补档前 HEAD：`d083b6eab27554ec96ec38c544b51f292e154df2`；ancestry 接入完成 HEAD：`9ca5356806b4683353aacae49e86fe636642ed39`。

### 数量更正与最终状态

原报告“15/22 已可达、7 个未可达”是汇总计数错误。重新逐条核验及原 22 行表本身均显示：**补档前 16/22 已可达，6 个未可达**。不能虚构第七个 SHA。原报告已保存在 d083b6e 的历史中，本节明确更正；不改写旧提交。

**原始 22/22 commits 现在全部 reachable from HEAD。** 时间窗重审还会显示第一次文档提交和本次四个 merge；这些收尾提交不计入“原 22 个科学/交付提交”的分母。全部原 SHA、补档前包含它们的 refs 和补档前/后可达性见 [COMMIT_REACHABILITY_AUDIT.csv](../archive/incident_20261003/COMMIT_REACHABILITY_AUDIT.csv)。

### 接入历史的方式

使用四个最少的 tip，依次执行 `git merge -s ours --no-ff --no-edit <SHA>`，没有 cherry-pick、没有重写源提交、没有用其内容覆盖当前 tree。

| 原先未可达的 SHA | 补档前来源 ref / branch | 接入 tip |
|---|---|---|
| 9dfaca3553ae5c175aad4b778ee4f58cea8d9669 | refs/remotes/origin/codex/drugclip-blind-benchmark-v01 | 9dfaca3553ae5c175aad4b778ee4f58cea8d9669 |
| 1efe33287b4ad0d0e1e3d69f38fb9835a1baf159 | refs/remotes/origin/codex/drugclip-blind-benchmark-v01 | 9dfaca3553ae5c175aad4b778ee4f58cea8d9669 |
| 492ce332f4727f441253b752a5711852657aa803 | refs/heads/codex/drugclip-blind-benchmark-v01 | 492ce332f4727f441253b752a5711852657aa803 |
| 270dd0eb5e88c0c369d493a58bd739cf118ebcdb | refs/remotes/origin/codex/drugclip-pacer-handoff | 270dd0eb5e88c0c369d493a58bd739cf118ebcdb |
| d948297b1011707cfe54ad9dbe47df446506aac1 | refs/remotes/origin/codex/drugclip-pacer-handoff | 270dd0eb5e88c0c369d493a58bd739cf118ebcdb |
| ddc11f970f508d89a75dac8ca30c41cfcd361731 | refs/heads/integration/final-freeze-v01；refs/tags/pacer-code-freeze-v01 | ddc11f970f508d89a75dac8ca30c41cfcd361731 |

来源 ref 指本次接入前包含该 commit 的 refs，不声称 Git commit 内记录了创建分支。9dfaca tip 覆盖 1efe332；270dd0 tip 覆盖 d948297；492ce332 与上述 remote blind-benchmark tip 已分叉，必须独立 merge；ddc11f97 独立接入。因此四个 tip 是覆盖实际六个缺失提交的最少集合。

四个 merge commit、parent/覆盖关系见 [ANCESTRY_PRESERVATION_RECEIPT.json](../archive/incident_20261003/ANCESTRY_PRESERVATION_RECEIPT.json)。接入前后 `git ls-files -s` 完全一致、status 都 clean、diff-tree 为空，tree hash 都是 `d6350f63cd92b7c462fdfc473d1ccb61208fe9a8`。索引原输出分别保存在 TRACKED_INDEX_BEFORE_MERGES.txt 和 TRACKED_INDEX_AFTER_MERGES.txt。这是保全 ancestry，不是把其他分支的旧代码应用到当前运行环境。

### 外部执行证据与大型资产

枚举 3420 个实际文件，按来源完整目录层级保存原文件名：
`project/archive/incident_20261003/external_evidence/C/projects/<source-directory>/...`。

- **2724 个轻量文件**已按原字节复制，合计 100,341,247 bytes；原 339 个已读文本证据全部包含。除日志、配置、manifest、receipt、provenance、CSV/MD/QC 外，也保留相关执行脚本、命令记录和轻量 docking pose。原已有 Git 文本的外部来源副本也可保留，避免混淆来源位置。
- [EXTERNAL_EVIDENCE_MANIFEST_SHA256.csv](../archive/incident_20261003/EXTERNAL_EVIDENCE_MANIFEST_SHA256.csv) 对每份复制件记录 source_absolute_path、archived_relative_path、size_bytes、sha256、category。
- **692 个大型/原始/缓存及辅助资产**仍不进入 Git，合计 77,087,365,186 bytes，逐项实际读取字节计算 SHA-256，见 [EXTERNAL_LARGE_ASSET_INDEX.csv](../archive/incident_20261003/EXTERNAL_LARGE_ASSET_INDEX.csv)。并非全部都是超过大小限制的文件；索引也显式保留小型 binary/cache/lock/raw-state 的未提交原因，不静默省略。
- 另有 **4 个 binary assets 的精确字节已经位于 reachable Git history**，在 ALREADY_PRESERVED_BINARY_ASSETS.csv 独立记录 blob ID，不重复列为未入 Git。
- 11 个历史依赖的记录相对路径在原工作目录不存在，但已在声明的 historical frozen root 找到大小与 SHA 完全匹配的文件，分别复制或列入大型索引；RECORDED_NONLOCAL_ASSET_REFERENCES.csv 保留原引用和匹配实际路径，未解决引用数 **0**。不修改原 manifest 的路径。
- 局部 .gitattributes 仅对 external_evidence/** 设置 `-text`，防止自动换行转换；无 LFS filter / policy 变更。局部 .gitignore 只允许已选择的 evidence，避免根目录日志忽略规则漏收。Windows 长路径仅在本次 Git 命令使用 `-c core.longpaths=true`，不改变全局配置。

大型/原始资产索引分类：

| 类型 | 未提交文件数 |
|---|---:|
| DCD | 120 |
| checkpoint (.chk) | 48 |
| NPY | 225 |
| NPZ | 22 |
| model/projection (.pt) | 5 |
| LMDB cache | 5 |
| LMDB lock | 5 |
| raw PDB / parse copies | 93 |
| raw MD System/State XML | 168 |
| runtime temporary auxiliary | 1 |
| 合计 | 692 |

索引包含本次运行及实际引用的历史冻结依赖，历史 calibration/model/benchmark 不转化为本次新科学结果。大型文件仍在源绝对位置保存；GitHub 上提交的是它们的身份索引，不声称 Git 内保存了全部原始轨迹/缓存。

### 验证与科学边界

复制件逐一与源文件、归档工作区文件、staged Git blob 核对大小和 SHA-256；receipt 见 PRESERVATION_COPY_RECEIPT.json 及 PRESERVATION_VERIFICATION_RECEIPT.json。已有 1,270 个 tracked 文件（排除本次允许更新的两份事故文档）验证原字节不变。最终新增归档内容与文档更新会改变最终 tree；“tree 保持一致”严格指 **四个 ancestry merges 前后**，不是声称添加归档后 tree 也不变。

scientific result files modified: **NO**。没有重跑、训练、docking、MD 或 PACER-FKG，没有修改 frozen assets 或 Stage4 scientific conclusions，没有删除 provenance。Stage3 protocol deviation 结论保持不变；六路 PACER-XR 仍未作为本批原 prospective selection 实施，后补只可称 post hoc supplementary PACER-XR analysis。
