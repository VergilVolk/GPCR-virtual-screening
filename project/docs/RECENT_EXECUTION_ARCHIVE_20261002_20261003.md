# Recent execution archive：2026-10-02 至 2026-10-03

记录日期：2026-10-03；分支：`incident/pacer-xr-protocol-deviation-20261003`。
审计基线 HEAD：`b78148e9c283aa834118f99ac5a4feb613885e28`。
工作目录：`C:/projects/GPCR-virtual-screening-incident`。

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

## 原始审计索引快照（2026-10-03 首次归档；下列状态对应 b78148e）

## 1. 范围、读取和完整性结论

本文件是现有代码、提交、manifest、receipt、provenance、日志及轻量结果的索引。只新增两份事故文档，不复制大型 DCD / cache / checkpoint，不重跑任何阶段，不修改科学结果。文件内容时间、run_id 和 commit 时间分别保留；例如 Stage1 回执 UTC 2026-10-01T18:15 对应北京时间 2026-10-02，并非不同批次。

对指定 --all 时间窗的 22 个提交逐个核对 HEAD ancestry 与当前 tree 路径；全部读取其 333 个文本文件版本（合计 40,687,065 字符，读取失败 0）。另完整读取下列 339 个外部文本证据，记录原字节 SHA-256。大型二进制原始轨迹、模型、缓存只核对现有 manifest/receipt 中的路径与身份，不加载推理、不解析重新分析轨迹、不重新计算科学指标。Git 对象及工作区字节哈希是保全核对，不是科学重算。

**首次审计时完整保存在 incident 分支：NO；本次补档后的 commits / 轻量证据状态见 completion 节。** 22 个时间窗提交中，16 个为基线 HEAD 祖先，6 个不在该分支 ancestry（原版本误记为 15/7，已在 completion 节更正）；--all 可见并不代表本分支已包含。另有本批 Stage1/Stage2 回执及 Stage3 docking/selection 文本在其他工作目录，不应把外部路径误写为 Git 内归档。当前分支保存了最终 Stage4 轻量结果和报告，但不是完整运行数据环境。本次不合并其他分支、不 force-add 科学产物、不删除外部证据。

检查基线已有 1,270 个 tracked 文件的工作区 SHA-256；两份文档是唯一授权新增文件。下表和附件的可达性、current tree 状态均对应文档加入前的上述基线，防止收尾提交日期落入 --all 时间窗后形成自引用。

## 2. 指定 Git 核验命令

在 incident worktree 完整执行，均 exit 0：

```powershell
git log --since="2026-10-02 00:00" --until="2026-10-04 00:00" --all --stat
git log --since="2026-10-02 00:00" --until="2026-10-04 00:00" --all --name-status
```

逐提交使用 git merge-base --is-ancestor COMMIT HEAD 核对可达性；逐路径使用 git cat-file -e HEAD:PATH 核对当前 tree。存在同名路径不代表外部分支的提交版本已合并，表中两个指标分别报告。时间窗按命令原样执行；commit 表保留各提交原始 +0800 日期。

## 3. 全部 commits

| Commit | 原日期 | Subject | HEAD 祖先 | 变更路径数 | 当前 tree 不含路径数 |
|---|---|---|---|---:|---:|
| b78148e9c283aa834118f99ac5a4feb613885e28 | 2026-10-03 21:17:36 +0800 | Finalize Stage4 prospective PACER-FKG evaluation | YES | 130 | 0 |
| 29fe47b5bdc18516e97989f5dabcff3cfebe4e1b | 2026-10-02 21:56:39 +0800 | Add prospective Stage4 membrane builder and bounded CUDA smoke | YES | 10 | 0 |
| 78ce1f30bc2e4855a4aa71dc6a1074a1e97a90e9 | 2026-10-03 13:16:30 +0800 | submission: standardized results.csv generator + Model Card (attachment-5 compliance) | YES | 4 | 0 |
| fb550656a16db20d0025ccbdbed23b59558ef527 | 2026-10-02 22:42:20 +0800 | fix: fresh-clone crash in --create-template + doc count drift | YES | 3 | 0 |
| 9dfaca3553ae5c175aad4b778ee4f58cea8d9669 | 2026-10-02 20:28:19 +0800 | artifacts: force-add frozen npz representations + README (gitignored by pattern) | NO | 3 | 3 |
| 1efe33287b4ad0d0e1e3d69f38fb9835a1baf159 | 2026-10-02 20:28:07 +0800 | artifacts: close Stage3 upstream gap - frozen 2026-side inputs | NO | 2 | 2 |
| 627825172102002363e54017258f820b1c8dd8e9 | 2026-10-02 17:55:09 +0800 | chore: record pre-Stage4 service and runtime audit baseline | YES | 3 | 0 |
| 05228ab1bf137629579d32e564531a6c43095991 | 2026-10-02 17:21:21 +0800 | benchmark: two-layer external M4 PAM functional benchmark v02 | YES | 18 | 0 |
| eb45810ad23f8691613ba64b20bee02c62283f42 | 2026-10-02 15:33:13 +0800 | build: align submission container with frozen runtime | YES | 10 | 0 |
| 492ce332f4727f441253b752a5711852657aa803 | 2026-10-02 15:25:24 +0800 | freeze PACER-M4 runtime dependencies | NO | 19 | 17 |
| c81cf4d26ecc4d7c7744887c7f381320a97c907f | 2026-10-02 15:11:29 +0800 | build: freeze Python 3.9 core requirements | YES | 8 | 0 |
| 85c08f74d2b903869da91a357d7c491a18cd5b75 | 2026-10-02 14:54:51 +0800 | feat: expose PACER-M4 through guarded local API | YES | 14 | 0 |
| bfbf94dadd34da0b76c540d8a3a696dd0e373150 | 2026-10-02 14:48:07 +0800 | feat: add audited PACER-M4 stage runner | YES | 17 | 0 |
| a8c4203c8b9cfb09f35bd8b88670d602aa266fa6 | 2026-10-02 14:39:19 +0800 | docs: establish PACER-M4 submission entry and evidence index | YES | 14 | 0 |
| d1f145832907c42850baf2e189c7142454f66519 | 2026-10-02 11:28:56 +0800 | ship M4 LOTO adapters and offline inference bundle | YES | 16 | 0 |
| f35acdd9d330e88ed8c690f81986e0305f5a515f | 2026-10-02 10:27:55 +0800 | fix router Python import example | YES | 1 | 0 |
| cd53f989ec22ef6df406b4924ceadd51567b3ea0 | 2026-10-02 10:27:27 +0800 | add deployable PACER DrugCLIP routing interface | YES | 3 | 0 |
| 512f8fa866a8c187a9a4923e7de4550c1c264c56 | 2026-10-02 10:19:51 +0800 | benchmark DrugCLIP generations with GPCR-safe routing | YES | 15 | 0 |
| 270dd0eb5e88c0c369d493a58bd739cf118ebcdb | 2026-10-02 09:38:34 +0800 | evidence: same-test-set M4R head-to-head 2023-triplet vs 2026-famaug | NO | 4 | 4 |
| d948297b1011707cfe54ad9dbe47df446506aac1 | 2026-10-02 09:24:27 +0800 | docs: repository organization -- mainline START_HERE entry + authoritative index | NO | 2 | 2 |
| 77f44dc6416f4d332bd2f0ff0bafbf4a9da9231c | 2026-10-02 09:22:56 +0800 | figs: per-module horizontal scoreboards + cascade funnel | YES | 3 | 0 |
| ddc11f970f508d89a75dac8ca30c41cfcd361731 | 2026-10-02 02:05:11 +0800 | freeze: finalize manifest-driven PACER production pipeline v01 | NO | 52 | 52 |

### 3.1 补档前未被 incident ancestry 包含的提交（本次均已接入）

以下六个提交在首次审计时未包含于本分支；现已按 completion 节以 ours merge 全部接入。下列清单和附件 A 保留首次审计的路径状态；没有 cherry-pick 或重写源提交：

- `9dfaca3553ae5c175aad4b778ee4f58cea8d9669`：artifacts: force-add frozen npz representations + README (gitignored by pattern)；3/3 个变更路径不在基线当前 tree。
- `1efe33287b4ad0d0e1e3d69f38fb9835a1baf159`：artifacts: close Stage3 upstream gap - frozen 2026-side inputs；2/2 个变更路径不在基线当前 tree。
- `492ce332f4727f441253b752a5711852657aa803`：freeze PACER-M4 runtime dependencies；17/19 个变更路径不在基线当前 tree。
- `270dd0eb5e88c0c369d493a58bd739cf118ebcdb`：evidence: same-test-set M4R head-to-head 2023-triplet vs 2026-famaug；4/4 个变更路径不在基线当前 tree。
- `d948297b1011707cfe54ad9dbe47df446506aac1`：docs: repository organization -- mainline START_HERE entry + authoritative index；2/2 个变更路径不在基线当前 tree。
- `ddc11f970f508d89a75dac8ca30c41cfcd361731`：freeze: finalize manifest-driven PACER production pipeline v01；52/52 个变更路径不在基线当前 tree。

ddc11f97 是 Stage1 run receipt 所指的 production_baseline，也是另一分支上的 manifest-driven implementation freeze；首次审计的 incident ancestry 未包含它，本次已接入历史但不覆盖当前 tree。不能因为当前 tree 中已有其他阶段代码就说该冻结实现已完整归档。492ce332 中有同名路径存在，也不证明该提交版本已合并。

## 4. Stage1 → Stage4 代码及执行证据索引

| 项目 | 核心文件 / 目录 | 保存状态和说明 |
|---|---|---|
| Stage1 实际实现 | project/scripts/fragment_generation.py | 基线 tree 中的 BRICS generator；实际运行身份以 U 的 Stage1 receipt 中 script hash 为准，本批非 LSTM |
| Stage1 运行配置/manifest/receipt/QC/logs | U/STAGE1_RUN_SPEC_v01.json；STAGE1_SEED_DERIVATION_RECEIPT_v01.json；STAGE1_GENERATION_RECEIPT_v01.json；STAGE1_QC_v01.json；logs_stage1_*；tools/derive_seeds.py；tools/stage1_qc.py | 外部工作目录；保留 authoritative pinned run 和历史 unpinned/rerun logs，不覆盖、不混算 |
| Stage1 实际输出 | G/generated_pam_analogs.csv | 外部：12 seeds→5,000 raw→2,605 pass；SHA a832c9cb2cfbed9b2e3db7902b48db5c384980cc70c551961c87f57db9e3be5d |
| Stage2 实际实现 | project/scripts/build_candidate_portfolio.py | 基线 tree 可读；receipt 是本批实施身份来源 |
| Stage2 manifests/receipts/logs | U/STAGE2_FILTERING_RECEIPT_v01.json；STAGE2_QC_v01.json；logs_stage2_filtering.txt；logs_stage2_qc.txt；tools/stage2_qc.py | 外部；2,605 input→200 pass /2,405 reject |
| Stage2 选择输入 | P/predock_portfolio.csv；all_generated_audit.csv；audit.json | 外部：160 local +40 exploratory；predock SHA 0dd7689860879eaae8abe806699a05a1c261a439e7355b99628e982abe46f52b |
| Stage3 初始 orchestrator / blocked receipt / runtime repair | U/tools/stage3_orchestrate.py；STAGE3_BLOCKED_RECEIPT_v01.json；STAGE3_DRUGCLIP_RECEIPT_v01.json；stage3/；STAGE3_OUTPUTS/ | 外部；历史阻断和初始 shortlist 是快照，不覆盖 corrected final ranking |
| Stage3 corrected scoring / provenance / logs | D/complete_apply_only.py；COMPLETION_RECEIPT.json；FINAL_VALIDATION.json；2023_provenance_audit.json；2026_provenance_audit.json；router_provenance_audit.json；2023_scoring_stdout.log；2026_extraction_stdout.log；router_stdout.log | 外部；200 分子同身份，无新增训练；2023 representations reuse，2026 second opinion，M4-safe route |
| Stage3 corrected scoring 输出 | D/prospective200_2023_gpcr_loto_scores.csv；prospective200_2026_family_aug_scores.csv；prospective200_two_model_scores.csv；prospective200_m4_safe_routed_ranking.csv | 外部；最终路线 SHA bda22c9354a9ebf4a263ce92ccf1257c8b27e8e3e747404b0f4fa97c3fb34d89 |
| Stage3 Vina ensemble | K/metadata.json；ledger.jsonl；framewise_scores.csv；statewise_scores.csv；poses/；ligands/ | 外部；200×10=2,000完成；不能称六路 PACER-XR；大体积原始 poses/ligands 保留在来源目录，不复制 |
| Stage3 selection 代码 | S/identity_contract.py；run_corrected_stage3d.py；stage3d_c_gate.py；stage3d_def_merge_diversity_pose.py | 外部实际 selection 实现 |
| Stage3 selection artifacts | S/stage3d_selection_manifest.json；stage3d_provenance_audit.json；stage3d_validation.json；stage3d_stdout.log；pacer200_structural_gate.csv；pacer200_stage3d_merged.csv；pacer_stage3d_diverse_shortlist.csv；stage3d_ranked_candidates.csv；stage3d_selected_candidates.csv；md_shortlist_pose_manifest.csv；md_shortlist_poses/ | 外部；200 gate PASS、140 scaffold reps→0010/0073/0027；selected SHA 5672be147ee2d6b1cc4ece4d0b95b0e90512fe6b54dfacb88de976353a507dd3 |
| XR 定义、接口和依赖 | project/docs/PACER_XR_CASCADE_METHOD.md；project/scripts/apply_pacer_xr_candidates.py；project/config/pacer_m4_release_manifest_v01.json | 基线 tree；六路必要条件与商业依赖，不是本批完成证据 |
| Stage4 MD 实现/配置 | project/scripts/pacer_stage4_prospective.py；build_pacer_stage4_prospective_membrane.py；run_pacer_stage4_membrane_smoke.py；M/campaign.py；master_manifest.json；job_matrix.csv；builder_input_validation.json；local_implementation_validation.json；BUILDER_HANDOFF.md | commit 29fe47b5，当前 tree；早期 NOT_STARTED 是实现交接快照 |
| Stage4 raw MD artifacts | B/systems/；B/production/；B/master_manifest.json；B/job_matrix.csv | 外部 raw MD backup；36×10ns=360ns；原 DCD、state/checkpoint、progress 保留原位置；未复制 |
| Stage4 PACER-FKG 实现 | project/pacer_fkg_v02/stage4_prospective_common_v01.py；run_stage4_prospective_phase1_bs256_v01.py；run_stage4_prospective_phase2a_frozen_apply_v01.py；run_stage4_prospective_phase2b_graph_region_v01.py；project/tests/test_stage4_prospective_fkg_v02_v01.py | commit b78148e9，当前 tree；本次读取既有测试记录，不重新执行测试或阶段 |
| Stage4 输入认证 | F/input_audit/INPUT_AUTHENTICATION_v01.json | 当前 tree；36_OF_36_ACCEPTED；含 raw-input hash / 原位路径；特例 0010/CP/R3 progress 保持原字节 |
| PACER-FKG Phase1 | F/phase1_bs256/PHASE1_FREEZE_RECEIPT_v01.json；SMOKE_RECEIPT_v01.json；manifests/full/；manifests/smoke/ | 当前 tree：36 full +1 smoke manifest，36 cache 记录；大型 cache 外置，不因回执存在误称已入 Git |
| PACER-FKG Phase2a | F/phase2a_frozen_apply/PHASE2A_FREEZE_RECEIPT_v01.json；manifests/full/ | 当前 tree：36 jobs、72 manifests；receipt SHA 04e95f7fa9a2e85f95ce25bcc01daa59ad84cdfaf016cecfd02fe3b0dd940f60 |
| PACER-FKG Phase2b/results | F/phase2b_graph_region/PHASE2B_FREEZE_RECEIPT_v01.json；STAGE4_PROSPECTIVE_RESULTS_v01.json | 当前 tree；STAGE4_FULL_FROZEN_EVALUATION_COMPLETE；结果 19,365,166 bytes，SHA ab224fa812588ccd89e5f988fd21d47a24f133a8c6cd25945564256f21be9a70；18 NPZ 外置 |
| Stage4 provenance | F/provenance/IMPLEMENTATION_VALIDATION_v01.json；TOPOLOGY_SELECTION_REPAIR_VALIDATION_v01.json；STAGE4_TEMPORAL_METADATA_NOTE_v01.md | 当前 tree；既有修复/验证和 temporal note，原 manifest 的 raw-spacing 字段不静默修改 |
| Final scientific reports | F/STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md；STAGE4_PROSPECTIVE_SCIENTIFIC_SUMMARY_v01.json；project/docs/PACER_STAGE4_FINAL_HANDOFF_v01.md | 当前 tree；科学结论不修改；No PAM/ago-PAM labels, probability, efficacy or independent-block inference |
| Stage1–4 已有只读审计报告 | C:/projects/GPCR-virtual-screening/project/results/STAGE1_STAGE4_RESULTS_AUDIT_REPORT_v01.md | 原工作目录外部报告；本次不移入、修改或将其当作 frozen receipt |
| 历史 PACER-FKG freeze/model/calibration | C:/projects/GPCR-virtual-screening-fkg-v02/project/results/pacer_fkg_v02_longmd_v01/calibration/V02_FREEZE_MANIFEST.json；历史 encoder/model assets | 外部历史依赖；freeze SHA b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd；Stage4 全部 EVALUATION_ONLY，不新训练/校准 |

缩写是索引别名，下面的全路径为实际读取位置；不是把外部目录映射成 incident tree：

| 别名 | 根路径 |
|---|---|
| U | C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01 |
| G | C:/projects/GPCR-virtual-screening-final/project/results/generated |
| P | C:/projects/GPCR-virtual-screening-final/project/results/pacer_candidates_v01 |
| D | C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01 |
| K | C:/projects/GPCR-virtual-screening/project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01 |
| S | C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01 |
| M | C:/projects/GPCR-virtual-screening-incident/project/results/pacer_stage4_prospective_10ns_v01 |
| F | C:/projects/GPCR-virtual-screening-incident/project/results/pacer_stage4_prospective_fkg_v02_v01 |
| B | C:/projects/PACER_STAGE4_MD_backup |

## 5. 历史 benchmark 与本次候选结果的区分

时间窗内包含 router / LOTO adapters、DrugCLIP generation benchmark、外部 M4 functional benchmark、scoreboards、服务 API、container/dependency freeze 和 submission packaging。它们是方法、历史验证或交付实现，不是本次 Stage1–Stage4 新候选实验结果。05228ab1、512f8fa8、270dd0eb、77f44dc6 的 benchmark/figs 文件仍按原提交索引，不把 benchmark AUC/EF 写成最终三候选性能。

历史 PACER-XR 官方六路 broad-AM benchmark 已执行，不证明本批候选六路完成。详细 protocol deviation 见 [事故记录](INCIDENT_PACER_XR_PROTOCOL_DEVIATION_20261003.md)。后补六路只能是 post hoc supplementary PACER-XR analysis，不能倒写原 prospective selection evidence。

## 附件 A：逐提交全部变更路径

“当前 tree=YES”仅表示同名 Git 路径存在于基线，不替代 commit ancestry 或字节版本验证。

### b78148e9c283aa834118f99ac5a4feb613885e28

2026-10-03 21:17:36 +0800 · Finalize Stage4 prospective PACER-FKG evaluation · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `.gitattributes` | YES |
| M | `README.md` | YES |
| M | `project/PROJECT_STATUS.md` | YES |
| A | `project/docs/PACER_STAGE4_FINAL_HANDOFF_v01.md` | YES |
| M | `project/docs/START_HERE_FOR_COLLABORATORS_V01.md` | YES |
| A | `project/pacer_fkg_v02/run_stage4_prospective_phase1_bs256_v01.py` | YES |
| A | `project/pacer_fkg_v02/run_stage4_prospective_phase2a_frozen_apply_v01.py` | YES |
| A | `project/pacer_fkg_v02/run_stage4_prospective_phase2b_graph_region_v01.py` | YES |
| A | `project/pacer_fkg_v02/stage4_prospective_common_v01.py` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_SUMMARY_v01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/input_audit/INPUT_AUTHENTICATION_v01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/PHASE1_FREEZE_RECEIPT_v01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/SMOKE_RECEIPT_v01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_no_probe__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_no_probe__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_no_probe__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_probe__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_probe__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_probe__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_no_probe__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_no_probe__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_no_probe__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_probe__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_probe__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_probe__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_no_probe__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_no_probe__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_no_probe__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_probe__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_probe__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_probe__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__apo__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__apo__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__apo__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__probe_only__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__probe_only__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__probe_only__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__apo__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__apo__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__apo__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__probe_only__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__probe_only__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__probe_only__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__apo__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__apo__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__apo__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__probe_only__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__probe_only__replica_02.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__probe_only__replica_03.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/smoke/PACER0073__candidate_probe__replica_01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/PHASE2A_FREEZE_RECEIPT_v01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_01__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_01__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_02__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_02__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_03__SIGNED_DRIFT.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_03__STATE_MOTION.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2b_graph_region/PHASE2B_FREEZE_RECEIPT_v01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/phase2b_graph_region/STAGE4_PROSPECTIVE_RESULTS_v01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/provenance/IMPLEMENTATION_VALIDATION_v01.json` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/provenance/STAGE4_TEMPORAL_METADATA_NOTE_v01.md` | YES |
| A | `project/results/pacer_stage4_prospective_fkg_v02_v01/provenance/TOPOLOGY_SELECTION_REPAIR_VALIDATION_v01.json` | YES |
| A | `project/tests/test_stage4_prospective_fkg_v02_v01.py` | YES |

### 29fe47b5bdc18516e97989f5dabcff3cfebe4e1b

2026-10-02 21:56:39 +0800 · Add prospective Stage4 membrane builder and bounded CUDA smoke · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/results/pacer_stage4_prospective_10ns_v01/BUILDER_HANDOFF.md` | YES |
| A | `project/results/pacer_stage4_prospective_10ns_v01/builder_input_validation.json` | YES |
| A | `project/results/pacer_stage4_prospective_10ns_v01/campaign.py` | YES |
| A | `project/results/pacer_stage4_prospective_10ns_v01/job_matrix.csv` | YES |
| A | `project/results/pacer_stage4_prospective_10ns_v01/local_implementation_validation.json` | YES |
| A | `project/results/pacer_stage4_prospective_10ns_v01/master_manifest.json` | YES |
| A | `project/scripts/build_pacer_stage4_prospective_membrane.py` | YES |
| A | `project/scripts/pacer_stage4_prospective.py` | YES |
| A | `project/scripts/run_pacer_stage4_membrane_smoke.py` | YES |
| A | `project/tests/test_pacer_stage4_prospective.py` | YES |

### 78ce1f30bc2e4855a4aa71dc6a1074a1e97a90e9

2026-10-03 13:16:30 +0800 · submission: standardized results.csv generator + Model Card (attachment-5 compliance) · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/docs/MODEL_CARD_V01.md` | YES |
| A | `project/results/pacer_submission_results_v01/results.csv` | YES |
| A | `project/results/pacer_submission_results_v01/results.meta.json` | YES |
| A | `project/scripts/generate_submission_results_v01.py` | YES |

### fb550656a16db20d0025ccbdbed23b59558ef527

2026-10-02 22:42:20 +0800 · fix: fresh-clone crash in --create-template + doc count drift · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| M | `project/docs/REPRODUCTION_QUICKSTART_V01.md` | YES |
| A | `project/results/pacer_candidates_v01/final/final_candidate_hypotheses.csv` | YES |
| M | `project/scripts/apply_pacer_xr_candidates.py` | YES |

### 9dfaca3553ae5c175aad4b778ee4f58cea8d9669

2026-10-02 20:28:19 +0800 · artifacts: force-add frozen npz representations + README (gitignored by pattern) · HEAD 祖先=NO

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/README.md` | NO |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/cands_science2026.npz` | NO |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/m4r_pocket_science2026_subsets.npz` | NO |

### 1efe33287b4ad0d0e1e3d69f38fb9835a1baf159

2026-10-02 20:28:07 +0800 · artifacts: close Stage3 upstream gap - frozen 2026-side inputs · HEAD 祖先=NO

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/cands200_famaug_scores.csv` | NO |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/predock_portfolio.csv` | NO |

### 627825172102002363e54017258f820b1c8dd8e9

2026-10-02 17:55:09 +0800 · chore: record pre-Stage4 service and runtime audit baseline · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/core_environment_audit.json` | YES |
| A | `reproduced_scores.audit.json` | YES |
| A | `reproduced_scores.csv` | YES |

### 05228ab1bf137629579d32e564531a6c43095991

2026-10-02 17:21:21 +0800 · benchmark: two-layer external M4 PAM functional benchmark v02 · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/data/benchmarks/pacer_m4_external_functional_v02/benchmark.audit.json` | YES |
| A | `project/data/benchmarks/pacer_m4_external_functional_v02/endpoint_records.csv` | YES |
| A | `project/data/benchmarks/pacer_m4_external_functional_v02/endpoint_summary.csv` | YES |
| A | `project/data/benchmarks/pacer_m4_external_functional_v02/molecules.csv` | YES |
| A | `project/data/benchmarks/pubchem_aid624126_m4_pam_v01/aid624126_screen_subset.csv` | YES |
| A | `project/data/benchmarks/pubchem_aid624126_m4_pam_v01/source.audit.json` | YES |
| A | `project/docs/PACER_M4_EXTERNAL_FUNCTIONAL_BENCHMARK_V02.md` | YES |
| A | `project/results/pacer_m4_external_functional_v02/benchmark_result.json` | YES |
| A | `project/results/pacer_m4_external_functional_v02/confirmatory_endpoint_skills.csv` | YES |
| A | `project/results/pacer_m4_external_functional_v02/figures/fig_extfunc_v02_two_layers.png` | YES |
| A | `project/results/pacer_m4_external_functional_v02/frozen_predictions_by_endpoint.csv` | YES |
| A | `project/results/pubchem_aid624126_m4_pam_v01/baseline_predictions.csv` | YES |
| A | `project/results/pubchem_aid624126_m4_pam_v01/baseline_result.json` | YES |
| A | `project/scripts/acquire_pubchem_m4_pam_screen.py` | YES |
| A | `project/scripts/build_pacer_m4_external_functional_benchmark_v02.py` | YES |
| A | `project/scripts/evaluate_pacer_m4_external_functional_benchmark_v02.py` | YES |
| A | `project/scripts/evaluate_pubchem_m4_pam_screen_baselines.py` | YES |
| A | `project/scripts/make_external_functional_benchmark_fig_v02.py` | YES |

### eb45810ad23f8691613ba64b20bee02c62283f42

2026-10-02 15:33:13 +0800 · build: align submission container with frozen runtime · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| M | `Dockerfile` | YES |
| M | `README.md` | YES |
| M | `docker-compose.yml` | YES |
| M | `project/config/pacer_m4_release_manifest_v01.json` | YES |
| M | `project/docs/LOCAL_API_AND_DOCKER_V01.md` | YES |
| M | `project/docs/RELEASE_VALIDATION_V01.md` | YES |
| A | `project/docs/SUBMISSION_CHECKLIST_V01.md` | YES |
| M | `project/release_audit_v01.json` | YES |
| M | `pyproject.toml` | YES |
| A | `requirements-api.txt` | YES |

### 492ce332f4727f441253b752a5711852657aa803

2026-10-02 15:25:24 +0800 · freeze PACER-M4 runtime dependencies · HEAD 祖先=NO

| Change | Path | 当前 tree |
|---|---|---|
| M | `.gitignore` | YES |
| A | `environment-drugclip.yml` | NO |
| A | `environment-md.yml` | NO |
| A | `environment-pacer-fkg.yml` | NO |
| A | `environment.yml` | NO |
| A | `project/config/runtime_dependency_map_v01.json` | NO |
| A | `project/docs/DEPENDENCY_MANIFEST.md` | NO |
| A | `project/results/dependency_audit_v01.json` | NO |
| A | `project/results/dependency_runtime_evidence_v01.json` | NO |
| A | `project/results/dependency_scientific_assets_v01.json` | NO |
| A | `project/results/dependency_smoke_core_v01.json` | NO |
| A | `project/results/dependency_smoke_drugclip_v01.json` | NO |
| A | `project/results/dependency_smoke_fkg_v01.json` | NO |
| A | `project/results/dependency_smoke_md_v01.json` | NO |
| A | `project/results/dependency_smoke_overview_v01.json` | NO |
| A | `project/scripts/audit_runtime_dependencies.py` | NO |
| A | `project/scripts/check_runtime_dependencies.py` | NO |
| A | `requirements-pacer-fkg.txt` | NO |
| A | `requirements.txt` | YES |

### c81cf4d26ecc4d7c7744887c7f381320a97c907f

2026-10-02 15:11:29 +0800 · build: freeze Python 3.9 core requirements · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| M | `README.md` | YES |
| M | `project/config/pacer_m4_release_manifest_v01.json` | YES |
| M | `project/docs/RELEASE_VALIDATION_V01.md` | YES |
| M | `project/docs/REPRODUCTION_QUICKSTART_V01.md` | YES |
| M | `project/release_audit_v01.json` | YES |
| A | `project/scripts/audit_core_environment.py` | YES |
| M | `pyproject.toml` | YES |
| A | `requirements.txt` | YES |

### 85c08f74d2b903869da91a357d7c491a18cd5b75

2026-10-02 14:54:51 +0800 · feat: expose PACER-M4 through guarded local API · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `.dockerignore` | YES |
| A | `Dockerfile` | YES |
| M | `README.md` | YES |
| A | `docker-compose.yml` | YES |
| A | `pacer_m4/api.py` | YES |
| M | `pacer_m4/runner.py` | YES |
| M | `pacer_m4/stages.py` | YES |
| M | `project/config/pacer_m4_release_manifest_v01.json` | YES |
| M | `project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md` | YES |
| A | `project/docs/LOCAL_API_AND_DOCKER_V01.md` | YES |
| M | `project/docs/RELEASE_VALIDATION_V01.md` | YES |
| M | `project/release_audit_v01.json` | YES |
| A | `project/tests/test_pacer_m4_api.py` | YES |
| M | `pyproject.toml` | YES |

### bfbf94dadd34da0b76c540d8a3a696dd0e373150

2026-10-02 14:48:07 +0800 · feat: add audited PACER-M4 stage runner · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| M | `README.md` | YES |
| A | `pacer_m4/__init__.py` | YES |
| A | `pacer_m4/__main__.py` | YES |
| A | `pacer_m4/cli.py` | YES |
| A | `pacer_m4/runner.py` | YES |
| A | `pacer_m4/stages.py` | YES |
| M | `project/config/pacer_m4_release_manifest_v01.json` | YES |
| M | `project/docs/CODE_AND_ARTIFACT_MAP_V01.md` | YES |
| M | `project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md` | YES |
| A | `project/docs/RELEASE_VALIDATION_V01.md` | YES |
| M | `project/docs/REPRODUCTION_QUICKSTART_V01.md` | YES |
| M | `project/docs/START_HERE_FOR_COLLABORATORS_V01.md` | YES |
| M | `project/release_audit_v01.json` | YES |
| M | `project/scripts/select_pareto_candidates.py` | YES |
| A | `project/tests/test_pacer_m4_cli.py` | YES |
| A | `project/tests/test_select_pareto_candidates.py` | YES |
| A | `pyproject.toml` | YES |

### a8c4203c8b9cfb09f35bd8b88670d602aa266fa6

2026-10-02 14:39:19 +0800 · docs: establish PACER-M4 submission entry and evidence index · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| M | `.gitignore` | YES |
| A | `README.md` | YES |
| M | `START_HERE_ONEPROT_PACER_DC.md` | YES |
| M | `project/PACER_M4_FINAL_REPORT.md` | YES |
| M | `project/PIPELINE.md` | YES |
| A | `project/config/pacer_m4_release_manifest_v01.json` | YES |
| A | `project/docs/CODE_AND_ARTIFACT_MAP_V01.md` | YES |
| A | `project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md` | YES |
| A | `project/docs/REPRODUCTION_QUICKSTART_V01.md` | YES |
| A | `project/docs/START_HERE_FOR_COLLABORATORS_V01.md` | YES |
| A | `project/release_audit_v01.json` | YES |
| A | `project/scripts/validate_pacer_m4_release.py` | YES |
| M | `project/tests/test_science2026_triplet_checkpoint.py` | YES |
| A | `project/tests/test_validate_pacer_m4_release.py` | YES |

### d1f145832907c42850baf2e189c7142454f66519

2026-10-02 11:28:56 +0800 · ship M4 LOTO adapters and offline inference bundle · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/README.md` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/bundle_manifest.json` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/drugclip2023_m4_loto_seed20260925.projection.pt` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/drugclip2023_m4_loto_seed20260926.projection.pt` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/drugclip2023_m4_loto_seed20260927.projection.pt` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_2023_gpcr_loto_scores.audit.json` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_2023_gpcr_loto_scores.csv` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_m4_2023_frozen_representations.npz` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_m4_safe_routed_ranking.audit.json` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_m4_safe_routed_ranking.csv` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_two_model_scores.csv` | YES |
| A | `project/artifacts/drugclip2023_m4_loto_pacer200_v01/requirements.txt` | YES |
| A | `project/scripts/build_drugclip2023_m4_loto_bundle.py` | YES |
| A | `project/scripts/export_drugclip2023_m4_loto_adapters.py` | YES |
| M | `project/scripts/pacer_drugclip_router.py` | YES |
| A | `project/scripts/score_pacer200_drugclip2023_m4_loto.py` | YES |

### f35acdd9d330e88ed8c690f81986e0305f5a515f

2026-10-02 10:27:55 +0800 · fix router Python import example · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| M | `project/docs/PACER_DRUGCLIP_ROUTER_INTERFACE_V01.md` | YES |

### cd53f989ec22ef6df406b4924ceadd51567b3ea0

2026-10-02 10:27:27 +0800 · add deployable PACER DrugCLIP routing interface · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/docs/PACER_DRUGCLIP_ROUTER_INTERFACE_V01.md` | YES |
| A | `project/scripts/pacer_drugclip_router.py` | YES |
| A | `project/tests/test_pacer_drugclip_router.py` | YES |

### 512f8fa866a8c187a9a4923e7de4550c1c264c56

2026-10-02 10:19:51 +0800 · benchmark DrugCLIP generations with GPCR-safe routing · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/docs/DRUGCLIP_GENERATION_BAKEOFF_PROTOCOL_V01.md` | YES |
| A | `project/docs/DRUGCLIP_GENERATION_BAKEOFF_RESULT_V01.md` | YES |
| A | `project/docs/DRUGCLIP_LOCKBOX_ACQUISITION_V01.md` | YES |
| A | `project/docs/PACER_M4_PROJECT_MAINLINE_DEVELOPMENT_AND_INFERENCE_V01.md` | YES |
| A | `project/docs/assets/drugclip_generation_bakeoff_v01/audit.json` | YES |
| A | `project/docs/assets/drugclip_generation_bakeoff_v01/generation_dumbbell.png` | YES |
| A | `project/docs/assets/drugclip_generation_bakeoff_v01/macro_metric_bars.png` | YES |
| A | `project/docs/assets/drugclip_generation_bakeoff_v01/macro_metrics.csv` | YES |
| A | `project/docs/assets/drugclip_generation_bakeoff_v01/matched_docking_comparison.csv` | YES |
| A | `project/docs/assets/drugclip_generation_bakeoff_v01/paired_delta_forest.png` | YES |
| A | `project/docs/assets/drugclip_generation_bakeoff_v01/per_target_metrics.csv` | YES |
| A | `project/docs/assets/drugclip_generation_bakeoff_v01/per_target_roc_heatmap.png` | YES |
| A | `project/docs/assets/drugclip_generation_bakeoff_v01/scaffold_bootstrap_ci.csv` | YES |
| A | `project/scripts/run_drugclip_generation_bakeoff.py` | YES |
| A | `project/tests/test_drugclip_generation_bakeoff.py` | YES |

### 270dd0eb5e88c0c369d493a58bd739cf118ebcdb

2026-10-02 09:38:34 +0800 · evidence: same-test-set M4R head-to-head 2023-triplet vs 2026-famaug · HEAD 祖先=NO

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/results/drugclip_science2026/family_aug_v01/m4_library_oldtriplet_ranking.csv` | NO |
| A | `project/results/drugclip_science2026/family_aug_v01/m4r_head2head_2023vs2026.json` | NO |
| A | `project/scripts/m4r_head2head_2023_vs_2026.py` | NO |
| A | `project/scripts/screen_m4_library_oldtriplet.py` | NO |

### d948297b1011707cfe54ad9dbe47df446506aac1

2026-10-02 09:24:27 +0800 · docs: repository organization -- mainline START_HERE entry + authoritative index · HEAD 祖先=NO

| Change | Path | 当前 tree |
|---|---|---|
| A | `START_HERE_MAINLINE.md` | NO |
| A | `docs/REPOSITORY_MAINLINE_INDEX.md` | NO |

### 77f44dc6416f4d332bd2f0ff0bafbf4a9da9231c

2026-10-02 09:22:56 +0800 · figs: per-module horizontal scoreboards + cascade funnel · HEAD 祖先=YES

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/results/project_wide_integration_benchmark_v02/figures/fig_v02_07_module_scoreboards.png` | YES |
| A | `project/results/project_wide_integration_benchmark_v02/figures/fig_v02_08_cascade_funnel.png` | YES |
| A | `project/scripts/make_module_scoreboards_v01.py` | YES |

### ddc11f970f508d89a75dac8ca30c41cfcd361731

2026-10-02 02:05:11 +0800 · freeze: finalize manifest-driven PACER production pipeline v01 · HEAD 祖先=NO

| Change | Path | 当前 tree |
|---|---|---|
| A | `project/drugclip_freeze/RUNBOOK.md` | NO |
| A | `project/drugclip_freeze/__init__.py` | NO |
| A | `project/drugclip_freeze/apply.py` | NO |
| A | `project/drugclip_freeze/dual_model.py` | NO |
| A | `project/drugclip_freeze/model_a.py` | NO |
| A | `project/drugclip_freeze/model_b.py` | NO |
| A | `project/drugclip_freeze/projection.py` | NO |
| A | `project/integration_freeze/FINAL_CODE_FREEZE_CANDIDATE_REPORT.md` | NO |
| A | `project/integration_freeze/FROZEN_MODULE_LEDGER_v01.json` | NO |
| A | `project/integration_freeze/RECEIPT_INDEX_v01.json` | NO |
| A | `project/integration_freeze/WSL_HISTORICAL_WORKSPACE_RECEIPT_v01.md` | NO |
| A | `project/integration_freeze/receipts/CM00734_PAYLOAD_REGRESSION_v01.json` | NO |
| A | `project/integration_freeze/receipts/FROZEN_EOL_NORMALISATION_APPLIED_v01.json` | NO |
| A | `project/integration_freeze/receipts/FROZEN_EOL_NORMALISATION_DRYRUN_v01.json` | NO |
| A | `project/integration_freeze/receipts/FROZEN_EOL_PREFLIGHT_AT_FREEZE_v01.json` | NO |
| A | `project/integration_freeze/receipts/PHASE1_ENCODER_REPLAY_v01.json` | NO |
| A | `project/integration_freeze/receipts/WSL_DIRTY_CODE_DIFF_v01.patch` | NO |
| A | `project/integration_freeze/receipts/WSL_DIRTY_FILES_BLOB_HASHES_v01.tsv` | NO |
| A | `project/integration_freeze/receipts/WSL_GIT_BRANCHES_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_GIT_BRANCH_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_GIT_DIFF_NUMSTAT_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_GIT_DIFF_STAT_IGNORE_CR_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_GIT_HEAD_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_GIT_LOG1_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_GIT_REMOTES_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_GIT_STASH_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_GIT_STATUS_PORCELAIN_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_RESULTS_FILES_v01.tsv` | NO |
| A | `project/integration_freeze/receipts/WSL_RESULTS_SHA256_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_RESULTS_TOTAL_BYTES_v01.txt` | NO |
| A | `project/integration_freeze/receipts/WSL_TOPFILES_v01.tsv` | NO |
| A | `project/integration_freeze/receipts/WSL_WORK_EXTRAS_SHA256_v01.txt` | NO |
| A | `project/integration_freeze/run_specs/CM00734_REGRESSION_v01.json` | NO |
| A | `project/integration_freeze/schemas/final_pipeline_run_spec_v1.schema.json` | NO |
| A | `project/integration_freeze/tools/_shim_helper.py` | NO |
| A | `project/integration_freeze/tools/diagnose_regression.py` | NO |
| A | `project/integration_freeze/tools/normalize_frozen_eol.py` | NO |
| A | `project/integration_freeze/tools/run_regression.py` | NO |
| A | `project/integration_freeze/tools/verify_phase1_replay.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/__init__.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/blocks.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/cli.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/engines.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/fkg_phase1.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/fkg_phase2a.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/fkg_phase2b.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/io_utils.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/paths.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/regression.py` | NO |
| A | `project/pacer_fkg_v02/frozen_runner/spec.py` | NO |
| A | `project/tests/test_drugclip_freeze.py` | NO |
| A | `project/tests/test_frozen_runner.py` | NO |

## 附件 B：外部文本证据读取索引

以下 339 个文本文件已完整读取；bytes 和 SHA-256 以外部原始文件为准，未复制其内容。原目录中的重跑/修复日志是历史记录，读取它们不代表本次执行了这些操作。二进制 DCD / cache / checkpoint 不加入本索引的文本读取计数。

| 实际原位路径 | Bytes | SHA-256 |
|---|---:|---|
| `C:/projects/GPCR-virtual-screening-final/project/results/generated/generated_pam_analogs.csv` | 262621 | a832c9cb2cfbed9b2e3db7902b48db5c384980cc70c551961c87f57db9e3be5d |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_candidates_v01/all_generated_audit.csv` | 625626 | 085c59bc874a68036ae93312c688a6badbc57c7072e7e23b0ca6914fce88da53 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_candidates_v01/audit.json` | 453 | cc1645d10632630f76415572fa1ec190b4dad5f999181e53e7e5ce5bb657bf98 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_candidates_v01/predock_portfolio.csv` | 50012 | 0dd7689860879eaae8abe806699a05a1c261a439e7355b99628e982abe46f52b |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/_stage1_first_run_copy.csv` | 331000 | 077aad8d8da2cbd89c994b9ea38d12faa05dcd75528d7d02f09bdbbb71a2e738 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/_stage1_pinned_A.csv` | 262621 | a832c9cb2cfbed9b2e3db7902b48db5c384980cc70c551961c87f57db9e3be5d |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/_stage1_pinned_B.csv` | 262621 | a832c9cb2cfbed9b2e3db7902b48db5c384980cc70c551961c87f57db9e3be5d |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/FINAL_RUN_REPORT.md` | 9276 | 6bc3319d740d00268c88c8e625075014dc1723a3a41a7df76bd6b3ecae2d189d |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/logs_derive_seeds.txt` | 672 | 4275f3791bcec1db1ee6215ad58b2fc82b37d77ed3e09e72480850017b2ca253 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/logs_stage1_generation_rerun.txt` | 5730 | ccb443b83276cd3dd4a5249e92dda1b932e494dfdb9ced3e8e35aa1b2b5e509e |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/logs_stage1_generation.txt` | 5734 | a8d0f4e849e48420e21ddffa5e448fd1adf78c0b9c758d2768f955cd20834936 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/logs_stage1_pinned_A.txt` | 5734 | 852c1da333463a3fcacf8a83b65db18bfc12bbfc819cce085585c8529a1fda64 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/logs_stage1_pinned_B.txt` | 5734 | 852c1da333463a3fcacf8a83b65db18bfc12bbfc819cce085585c8529a1fda64 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/logs_stage2_filtering.txt` | 2422 | b01726e3680edc4e9ed8e1e47cc561fb7bed49d07ba8c54d6950413c1be9bedb |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/logs_stage2_qc.txt` | 6370 | 113b14fffad4da7126818c3f2672a8a1f4406435f318034eceddb342f8965d0f |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/logs_stage3_orchestrate.txt` | 5174 | 6a4b0de888118ac8ba19cc2989533a0c4c5ff7bd271c1a23fb07a764c541fa9d |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/README.md` | 2438 | 51d02d64854f4eccc6d920cd128fec4bed93ede3f368ebd42569cf7de69d0bf2 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE1_GENERATION_RECEIPT_v01.json` | 5971 | 0833a5a79221fdcff41d41558cfc39fca0fd5dc003b21a6ad5b955001c5eba25 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE1_QC_v01.json` | 941 | 0cc0db38220ff9d637654482c7bd039f7cc00806f1db5c6679b14a44cb459193 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE1_RUN_SPEC_v01.json` | 4583 | bb5b497f938865a8a90a1cccd2f17278527b3cacb4e75e59fcf7416621aa685a |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE1_SEED_DERIVATION_RECEIPT_v01.json` | 3674 | 09758020f6b25d2a3c3f34f18d1257bcd953333ebd2a039a54d8db3d8f3b6b4a |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE2_FILTERING_RECEIPT_v01.json` | 4915 | 56bda891052a4ef1377d8eb954eee4d5b40b4139c0a9e401a67a8ee92f7dfc91 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE2_QC_v01.json` | 2631 | 206acbc3c855c02bcbf102357c6997aa62d4c83cf4f703b4e87d5f6dc12074f9 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE3_BLOCKED_RECEIPT_v01.json` | 4504 | eeb2300956ed4d6818fc035080a1177bf1af275ad443805af35874dd2e13e1ae |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE3_DRUGCLIP_RECEIPT_v01.json` | 128010 | 502cd42d13820317c324d3bd602e645e3501a8ebdaed6e11563195a5d6c78474 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE3_OUTPUTS/model_a_scores.csv` | 34286 | a263a866db1e4ccdf6ac05d83f81cae3d29704d2ba00cb0f13ee67d538ee05ba |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE3_OUTPUTS/model_b_scores.csv` | 27815 | 4323130ee12c2774c8ba217bcd5235360e9e7f0ad288116b15db0350867c3364 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE3_OUTPUTS/model_b_scores.json` | 37990 | faac31c6a87b9979000590bc274bedbf64cfd648a0ba47eb755bf4975a4b1ca6 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE3_OUTPUTS/STAGE3_DUAL_MODEL_DECISION_v01.json` | 1017 | d8a68aeedd8a50038f1acacea821468d9019573f40a128a2ffaad6e5b1c631c2 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE3_OUTPUTS/STAGE3_DUAL_MODEL_SHORTLIST_v01.csv` | 858 | 10fc0c8ed474a309271b0817bff38e2080847d3ca0a7272f6f9313a8fe7265ac |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/STAGE3_SHA256_MANIFEST_v01.txt` | 3233 | 6dacc46a66e1809399609ea60cc7f65286545d94e8e4025622ec8e33067a9a0f |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/embeddings.audit.json` | 1226 | e42c095d38c722301ed1c54b7c8ecb113fcf2c14c38bd56406f092ab6dbd628e |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/embeddings.json` | 1227 | 1eb4fc087273e02a48c78b616ec7b73cdefecde78ea8b51e34ee49041eb829de |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/ENV_EXPLICIT_lock.txt` | 15285 | 8913585e08f61fb267eaaba498a3d138aab609eff5673cf7ee9135aec00ad691 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/ENV_EXPORT_final.yml` | 6227 | e19092389e545963621bebb76bba0949a4a16f36ab441851bc0d9e5132e6f309 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/full_extract_stdout.txt` | 1420 | a51bbd203d826c770adb8a59d2209137b743563f14859f75d77933238aa8e0b6 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/molecules.lmdb.audit.json` | 74 | 257610ffd0d1db97787bc57650620e92c77272e584b59afdc7471d6558ec3a2f |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/smoke_embeddings.audit.json` | 1222 | dfe8eb58f258ae83cdd0e5adcf40bb64bcb944276ea74a5a53dbd2cc36c4701c |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/smoke_extract_stdout.txt` | 1416 | 3c185fd856e4a79f021765c171a43fe37ee068f5fb972871e5859aa80e467c61 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/smoke_molecules.csv` | 101 | 8f1db10f03fee6cf37c501655208ea68f285745a6e6fff63c97b1ef0962562cc |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/smoke_molecules.lmdb.audit.json` | 70 | a7510bbd3d389c36317c6d1a0cc941f9ea211f79ab043952076b38bb500da836 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/stage3/STAGE3_RUNTIME_COMPATIBILITY_RECEIPT_v01.json` | 1341 | b39816f9156c82d3b9b234793b55276c4e287a9d46de102ea0d2dc1cc32e039c |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/tools/derive_seeds.py` | 3599 | 30e74a59b479e7c2e8a45261dedb93ae5a503e6e299c8cfc11884c5663258c06 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/tools/stage1_qc.py` | 4110 | 2eae0076bd8a2600753ae8196b3f98c454e138a8271dab38d3362f8048f15a2c |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/tools/stage2_qc.py` | 6173 | 82db1e6f12805519a3603e969eed0388ae1a154ebe05dd23ae1b7b892cc08141 |
| `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01/tools/stage3_orchestrate.py` | 7922 | 95b6f1f57085bd8f05a4fcc497b3cfd47eeed45be3e134c1aebbc28537d03318 |
| `C:/projects/GPCR-virtual-screening/project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01/framewise_scores.csv` | 397915 | 7eab001a2bb579e964cb6f93f21b89236f3092991bd449581bc3b74b88682d93 |
| `C:/projects/GPCR-virtual-screening/project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01/ledger.jsonl` | 857741 | 829f32c2f8c9f77cee977196b275d52f28129d395d7a28a746fc8e4268ba037a |
| `C:/projects/GPCR-virtual-screening/project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01/metadata.json` | 605 | 8d1703f4f3a40b7b73cc6fb74fb7b04f63ade20ba040608e6624607c5b0746d3 |
| `C:/projects/GPCR-virtual-screening/project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01/statewise_scores.csv` | 397915 | c9959c05a319f9e7813400381b4179f56cfe0ee380a6c4f3d2205592eb6e7f9e |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/2023_provenance_audit.json` | 2554 | 2e6f0a3e209a7b2c2d594ca8501fd3226ebc7e45744bf0855b53fd32dd9b7a00 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/2023_scoring_stdout.log` | 1047 | dd419df5065a5d89e6d43da8142713e4fcac09a0674bf86d0312323fe996aa09 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/2026_extraction_stdout.log` | 1578 | 07920d35c2e1cebcfe4aa529fd836e72d68945e15cb6e91891069721efe4d142 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/2026_provenance_audit.json` | 5633 | c964ee42d0619bafe4417a4a20bd62107e76f0701eee3353bb4dbd313344c113 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/bundle_manifest.json` | 1942 | 890f645d7c795d189834db7d20ff4f2df5334c0f616211e0d40ad6a835653198 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/complete_apply_only.py` | 13620 | e49d9de72d851d77a466d476ff1413effcf3a284b4275cd8b232072a50a35883 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/COMPLETION_RECEIPT.json` | 1557 | 186ccd108367f1c69175c2ae8c81971d2d6bc73bcdc32bf453c025d205d1172c |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/FINAL_VALIDATION.json` | 1370 | 07c6cfdd9ce3b0261acc699e17cb8f4acaee2370e619e81f292ee6fc14994d07 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_2023_gpcr_loto_scores.audit.json` | 1269 | 13fffa1523d7172a6b6e9463cd92aad473003d72e3197d02ff1095a3969e6de3 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_2023_gpcr_loto_scores.csv` | 25928 | 4763c4283b93c257ee348db0488d8e381977e6d6df60d3bf13df4074bbf28c15 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_2026_family_aug_scores.csv` | 24986 | 813ac675668b3b36c17ed0d3752a76b2fabee9ebd2d21bf32877dcbd0ab464ca |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.audit.json` | 971 | 63da29cceb0788405fbc51af2c4eab105b0c2eb734077d81b2664c5697e103fd |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_m4_safe_routed_ranking.csv` | 51455 | bda22c9354a9ebf4a263ce92ccf1257c8b27e8e3e747404b0f4fa97c3fb34d89 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_science2026_representations.audit.json` | 1217 | 1eb29d3fdce11f33cfa4210006f877418a25f6f4bb9fdf8d62676d8f6a6754d8 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/prospective200_two_model_scores.csv` | 22339 | 493a886353129502dcd3a72978634e7bb5bfe74fe4ad0cdb75162fe96c1c3034 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/router_provenance_audit.json` | 1853 | 8751cc2648877ee1bc38a497247a35ecc9a3e4d67d37c9d74270fd6a6e9c1f35 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01/router_stdout.log` | 972 | 74361edb8c9a9895e99076d71829468980bd9b5acb5408a725e615558d6cdc6e |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/identity_contract.py` | 3389 | ab121b2f346fe5a00b15c6efb82ab79ec423c713395b864f2a6966615eb2c78d |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/md_shortlist_pose_manifest.csv` | 3823 | 8ef3c4be8ff3800543d09e3b37c11331671e961c7bd3f29cc70850b787e311a6 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/pacer_stage3d_diverse_shortlist.csv` | 97794 | 2116b0bc3cdf26a2f88ebdfd19e67d7410e7e8d4de4c8fde3e5ba24f7f377cbb |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/pacer200_stage3d_merged.csv` | 113681 | d2a679786ef3e39277476c4278ea569ccdd3450d2f581cb0f4de93cdeb45c6a9 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/pacer200_structural_gate.csv` | 60772 | 34a9f08e62902f63606cbf3a69972f5d712d2d06135374cda3ba67c62cc67809 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/run_corrected_stage3d.py` | 9804 | 1c659a3d14613cea71f1dd6088a4c5c4bbb300bf6d929636e536bebd6f6aa5ed |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/stage3d_c_gate.py` | 6045 | c86a804703dd043e6847c677d98624e4f54fd9ce51a2d4027c561d574c1b838d |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/stage3d_def_merge_diversity_pose.py` | 7886 | f2a804da6d5187d6051e227c3cc3144cf91db513929072c5572da06a530e8987 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/stage3d_provenance_audit.json` | 7125 | 7ea6037409e95f6c604a35c3081351fc20be5de9c9aca6ce9925d763e4cce0e1 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/stage3d_ranked_candidates.csv` | 124799 | acc52dbd22fce024be6d42e6912ad37b68dec7c509cea94db87fc0b9229ba9c6 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/stage3d_selected_candidates.csv` | 3061 | 5672be147ee2d6b1cc4ece4d0b95b0e90512fe6b54dfacb88de976353a507dd3 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/stage3d_selection_manifest.json` | 6034 | 430986b24e6c63152bf6ea8f3846cd013969e6c321d3741ea4b22c493ed5c162 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/stage3d_stdout.log` | 2024 | d7d4ed075d4aaaac0f9c1d3320315ba79c2cf87a18e23a34747d744e9d396e23 |
| `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01/stage3d_validation.json` | 1128 | b2da207ab4baffe4c5a0d4b31ebe86a3143974a7907346ce9d455845803f01ca |
| `C:/projects/PACER_STAGE4_MD_backup/builder_input_validation.json` | 1105 | 9cc962540650dd4c83bdb1ea0705eb67944b04c073718969c562bf9722eb644f |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster0__apo/replica_01/audit.json` | 1418 | aab27bb54fdd1536337059c5d2a97e36ed8522b5bcc8b611672bb91110d6d9d1 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster0__apo/replica_01/npt_state.csv` | 1363 | 2be7ab524d62876a5d7babbeec68c4306318b7e54f2b2f969e7fc520f2878662 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster0__apo/replica_01/nvt_state.csv` | 1366 | 6a2f7fe1ba59406da58e0755ae6f2097eb3e43831687c9c90486673fbbaf886f |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster0__probe_only/replica_01/audit.json` | 1429 | 2e307150a1666fb5a2d1a95f5d4e9fd2e015072e92e3f62ffc26f57a3672de1a |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster0__probe_only/replica_01/npt_state.csv` | 1365 | 4eb8acdab499ee54a141a14f776abbcf6e52d64c393043f3cb27549271b4c4b2 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster0__probe_only/replica_01/nvt_state.csv` | 1365 | 16b1a254ba237e2270d69ad360c54c729a86b6ae6fd61ff56b58c4e278337aea |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster4__apo/replica_01/audit.json` | 1426 | 195d3bfefc304db47be7447e8bde8a8a7d53a467f8b15b6c15377268e42314c7 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster4__apo/replica_01/npt_state.csv` | 1380 | 6540e9afc98f85ad2da3592ec5fc1d8d9a43547eeff73fb1916f3678682b96f2 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster4__apo/replica_01/nvt_state.csv` | 1382 | cf687dd0f6913d7cf9052dfd38dc565147a16a6ddd6c19e2d9cbb2d3e7ff5215 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster4__probe_only/replica_01/audit.json` | 1432 | c698cdc0442a7b99c93834282c9fe8177b08fcdcebe995a7f47628e82e848e8f |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster4__probe_only/replica_01/npt_state.csv` | 1382 | 62782a9d878db087efc6b123f96672fdce53c2829b1e5ccf164161159f3c2e74 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster4__probe_only/replica_01/nvt_state.csv` | 1378 | 311bd001cb84773889f2bb9777ba26b697a4c192558d2ecc6d000fdc6471567f |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster9__apo/replica_01/audit.json` | 1425 | 5319384ff2c4397fe0f6d59c376a0be7c81faaec07abeb3d501f97e9e44a8610 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster9__apo/replica_01/npt_state.csv` | 1382 | a949146dd94452ef38177aebcf6da9d18685540404cba8afefb9ea27a671a21d |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster9__apo/replica_01/nvt_state.csv` | 1385 | cf99d2208c8d2e2d2311fc8b94626afd3a10da7c5e71354366fbfef55ea156e9 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster9__probe_only/replica_01/audit.json` | 1433 | e8418ce6d2f8bc30b4cfd543bc57982b5d9a4c6b61ef1ff119ca43414f028ca6 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster9__probe_only/replica_01/npt_state.csv` | 1378 | fa622589ce8b638971fc02f6ba051fb7f568ee251d9c6fc15dc2e010a1df6972 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/cluster9__probe_only/replica_01/nvt_state.csv` | 1395 | 3e9ee95d1b764c4556da04701a38fece391a29186ec2c1f2e75960cfb738e92f |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0010__candidate_no_probe/replica_01/audit.json` | 1442 | 6c7d2026b972c24a763baad0f0dfb17926b47e1283a5b876a069b5fd21d78c02 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0010__candidate_no_probe/replica_01/npt_state.csv` | 1382 | df552c23014c59b0ff708dc1d06d6822f7ceea92174bd2e8835b605fc2219a75 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0010__candidate_no_probe/replica_01/nvt_state.csv` | 1382 | 6f9d1ab45945d7ee0d1ffd04117ad953f200ff2c24e319022126a8891a65edb7 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0010__candidate_probe/replica_01/audit.json` | 1439 | 4b7691adbdcf37665de4ed01bc29f83a2dffc3bd86902e2709788ca5a320ee2b |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0010__candidate_probe/replica_01/npt_state.csv` | 1385 | e3950f1a33f84f35dfaac82290745350bf3b17ebcf93a933b79466937b151af6 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0010__candidate_probe/replica_01/nvt_state.csv` | 1382 | fe5bca579dd5e1cda5419f65a23066d49c41fb92a3c16780fa4e0b890d8a9214 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0027__candidate_no_probe/replica_01/audit.json` | 1446 | a049d66830d69949767c089619edb7b688ba3ed67d8c86fd80fcc2d9ecc8ffde |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0027__candidate_no_probe/replica_01/npt_state.csv` | 1389 | 09125fa4d14bc51d18f7b0f6a7fd0aa99d951cd75e319ec0496152d218c201d2 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0027__candidate_no_probe/replica_01/nvt_state.csv` | 1375 | a3f96ea5d0aa98f7f84216d9d7156d91448887bbcfbb5527b99550f0845c9d95 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0027__candidate_probe/replica_01/audit.json` | 1442 | 7ce0bb612c435264a30707787722f1858fcebe56a0c4fe86742316d35bf9244c |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0027__candidate_probe/replica_01/npt_state.csv` | 1389 | 5dd35f54ce71ab1d59044278edf93cdfd14e88940756b6f9810107990db7439d |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0027__candidate_probe/replica_01/nvt_state.csv` | 1385 | a019f69cd603b2a27b6d81924b60e9c44d35965cea8fc456eacb263b21266500 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0073__candidate_no_probe/replica_01/audit.json` | 1438 | 3df25bc2111518405a036905ab9825fedc20a341a41275d4f7310014ff4cfd71 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0073__candidate_no_probe/replica_01/npt_state.csv` | 1362 | 128c3d6834717c62f13465e4630a578a8e2b718fe2c3b1e22be35be1b8518ed6 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0073__candidate_no_probe/replica_01/nvt_state.csv` | 1367 | efa134a9f66ab8a61a7e7fff7fc9a1ff3a1c4938af9687a20fbbd5909b7e6af5 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0073__candidate_probe/replica_01/audit.json` | 1433 | 10e54f953431c8f2bc17ca3d4f1f02d2ef9d1ae351bb67a0d2d3793e09749af3 |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0073__candidate_probe/replica_01/npt_state.csv` | 1358 | ae5f5d0b49cc9699f1e65fb9b7d935dfdfd03ee1edb264d3c49d63878540ba3c |
| `C:/projects/PACER_STAGE4_MD_backup/equilibration/PACER0073__candidate_probe/replica_01/nvt_state.csv` | 1363 | 31d21f55021550aa638ae8b796ba0d2662a6c03f6c8532a86d335af41f5d59a3 |
| `C:/projects/PACER_STAGE4_MD_backup/job_matrix.csv` | 9498 | 9749068bd82ac80a7b88961ead6dbb7bef2b053f3feb5f7f237c7965cc2134f2 |
| `C:/projects/PACER_STAGE4_MD_backup/master_manifest.json` | 16382 | a9b5a208accd467326c9b31c0262eff3d96d8e4e6fc326ac8a64752e2cef2329 |
| `C:/projects/PACER_STAGE4_MD_backup/orientation/cluster0/alignment_audit.json` | 84098 | 82a774eef5a595f2f86c2b7826033b6c7326b827b2c8b31f3309abc76654ba3f |
| `C:/projects/PACER_STAGE4_MD_backup/orientation/cluster4/alignment_audit.json` | 84089 | 0d3e09702105b59b6ce8a0cc62c4310389198780defc96534209df5990f00f14 |
| `C:/projects/PACER_STAGE4_MD_backup/orientation/cluster9/alignment_audit.json` | 84096 | 1ffb24ed074dd6e6f97b51e65fad4929ddcff4f2ae0db1e209cc6ac0bdb6076a |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__apo/replica_01/progress.json` | 513 | 3bb35c4195db898dfeeb70efa60cdb6b0f354ff6335dab352397d666b067d151 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__apo/replica_01/state.csv` | 128407 | e44d86a3f94f05162a9f80fa99aad28778fd4cd952f792b553e13f488c4e3250 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__apo/replica_02/progress.json` | 514 | 323e13c8b434c0cd4fe8403f30ffeb8fdecd9425ce8e90ec6c4216f8f9a7d0d0 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__apo/replica_02/state.csv` | 128407 | bed07b53f0d96a03ecb041c501d34f76e97aa1e6944d5f79746d936474a7fb4f |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__apo/replica_03/progress.json` | 514 | e1839258f23dce64b803e3380714a584b0a9edfc2230adaae57469c6d364d6ec |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__apo/replica_03/state.csv` | 128367 | 2517402bd9c85eaae39235da341fb507da5325f233adb2b205074781768c9d81 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__probe_only/replica_01/progress.json` | 521 | 6a1c3de53521538d4ba500d83906e9d859574b9bdb367ad355122ac6a5e7ebe2 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__probe_only/replica_01/state.csv` | 128378 | 50dec77d963e3d4247e3ca21b009081cf7d71980a5aa275abe3216c4e6940215 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__probe_only/replica_02/progress.json` | 520 | 6748f31c3d3279688cf37d4e165f741e74d7d37c44e8efd881a48210905207e1 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__probe_only/replica_02/state.csv` | 128347 | 8ed77aad7fc5b6b2d71e7a6d73ed55b9155dc8adb2aaabc62859ca6b51e17f00 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__probe_only/replica_03/progress.json` | 521 | 69d12328b0692d7de5f5e5a342333ff362aea83714966369c4d3a1488e1ee329 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster0__probe_only/replica_03/state.csv` | 128353 | 9932cca60a3c11d7dc1e6a67a64191ad1c769c76993bead2b98593de7235179c |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__apo/replica_01/progress.json` | 513 | 1b65b6540a52a177ea14d29672300cc6e6808309f3f7aefe19759be8509c583a |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__apo/replica_01/state.csv` | 130033 | 3f3c93ded57df9df075ec58bb884c4c688e479bc6ce6e5fd8652353ac1a9eced |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__apo/replica_02/progress.json` | 513 | 3b00d406f89a7e7f8266bd1f0bcccb83a4fc520e6e51cc4a457cbf27966cdaeb |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__apo/replica_02/state.csv` | 130075 | 1d6ea0e054623df2c5daa770aa26b291d56b5eb577dd4631abb431ef3040c3d9 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__apo/replica_03/progress.json` | 514 | c3f4264712a8320bf6195af75be5e7d70f77e071479cadbc05e4ce7912ecfc62 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__apo/replica_03/state.csv` | 130018 | 38a1bce9b7d9c3ffa854b145856d78f68c0aedc001901e8855d6514f94729eba |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__probe_only/replica_01/progress.json` | 521 | b0e280bae575178d79e71cda18ae715839c60e63413a92b8bd8e186e14c57353 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__probe_only/replica_01/state.csv` | 130061 | d61ef3f213fbd688636890ab85edcc5074be16ae8a6feea7bb10b9e5f5e5eecc |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__probe_only/replica_02/progress.json` | 521 | d7b54927d31a0b95e7a9fc20f8d57a1b7f21a709fd236cb015f0e1b6bdcb7e04 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__probe_only/replica_02/state.csv` | 130012 | 332dbd4856699baeeb87268dad437d5c4e43fee9611246c0d5f297d7b07ed53e |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__probe_only/replica_03/progress.json` | 521 | 1d801fa6a9ffb6663b2b7d47a1868633bec6cd312681b7cd4d1575f300776874 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster4__probe_only/replica_03/state.csv` | 130012 | 665697a22d60d472ac82d7f3bc8707f5efed58e2fe66d912dc16216df44204e9 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__apo/replica_01/progress.json` | 514 | a371fbbd35bd98b616a2effc8182cdd58d00ad0aff6b26f2c97a160c3f86a98b |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__apo/replica_01/state.csv` | 129884 | 83b33cb4c25b78961483f8e2b2161c40e6f8ab7edf7cf104a162f89e109b4dde |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__apo/replica_02/progress.json` | 514 | 5ddfd3d9cb07c327904c907122ef7ac5eff670e064df87e0d4a90d7a29b18d42 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__apo/replica_02/state.csv` | 129930 | 8fcacf2cbc020e19b30e25b346a081b43fc7acd72211af7576a6526a1d221d1a |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__apo/replica_03/progress.json` | 514 | 63f053f522741933c9de6f36457007cf4f8659c3b197fc3f6333323c60e71397 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__apo/replica_03/state.csv` | 129953 | dc468b4527ff83d82e9c39fd10800a4ee24b58ccc7d96564a3ea716b86aee705 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__probe_only/replica_01/progress.json` | 521 | b389077c23c42ddb120e39137d34a04f94172dac6e78bd5746e59649f69c97f1 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__probe_only/replica_01/state.csv` | 130004 | 625b4fed1beb4e66f1f47d311720f6f7a21dbe3548c82a19b78bb4259b39950c |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__probe_only/replica_02/progress.json` | 521 | 699f5c2b54016e8410b269a6ae6093b57c71c3c082f49a1a363f671585194163 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__probe_only/replica_02/state.csv` | 129948 | ad572970010a2df1a4bafe811be2b3d6e810b236809000433d22c106ba2bc39e |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__probe_only/replica_03/progress.json` | 521 | a28358eb90bbb139786abe8eaf77487f72688e8b06ed0b062d52d86dd7c9de29 |
| `C:/projects/PACER_STAGE4_MD_backup/production/cluster9__probe_only/replica_03/state.csv` | 129982 | df3084c512b0a36bb28f7a15156e3c255a0881b9489fe7de696ab6b183feff41 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_no_probe/replica_01/progress.json` | 530 | 01ae7cb624c63f4dedf2533c9d89fa5e16f15fd4f9e39fb384f86379eac6b1cd |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_no_probe/replica_01/state.csv` | 129957 | 08fc7d5ec59e6bda82b2a934a9045f2c15996fc93e1646776b280b266f73704d |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_no_probe/replica_02/progress.json` | 530 | 05d3491c97304896ce6e08fa1fff3bbdbde4d410b9fe4c7b8fbcff688b83dc6b |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_no_probe/replica_02/state.csv` | 130041 | 834e795418aea4754f517d476afe1903ad7afd11eee74f36ed279b5349a977fd |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_no_probe/replica_03/progress.json` | 529 | 2e86abcc086c4914af1f54488c3d5f7ed312a029e1ca153d56548f08bbe72266 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_no_probe/replica_03/state.csv` | 129943 | f1072a6134c054dbcd48ca66cb3b9b02fc18eb31e1b7b095cca16cfd0dca3110 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_probe/replica_01/progress.json` | 526 | 699c8189bcd7f8399d09440d6851917920d2edfc3fb4ed0f04c8db65c5afab12 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_probe/replica_01/state.csv` | 129897 | 91d8b8c4f97f896a954bec163f03de144cf217c96438dbf1cc162b3c053c507a |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_probe/replica_02/progress.json` | 527 | 98cb6c30ad67b8f10fd790e9b3af01da93d4c2ba595c2e0b34b89891f48debf0 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_probe/replica_02/state.csv` | 129988 | 03b5eff07ec590a60472dc868e01b5019dba42263730a9ff77a0b25ac977de36 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_probe/replica_03/progress.json` | 431 | cf4518f467a8fe35f65cbea1474ba7e6ff6637408bc77d66a1a3a2fe6c50cff4 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0010__candidate_probe/replica_03/state.csv` | 129925 | 2685f1427b34fc1c62b9449ceebeea8349736249620a5f15074089c7f7af546c |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_no_probe/replica_01/progress.json` | 530 | 5274cb71040d060575b5ebaa3d8281abfe11eb728b6459bb202792f8a5b6366b |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_no_probe/replica_01/state.csv` | 129911 | b006b24cb8734e0a8b64774ffe8d0b03df4cdf33995a15b4313d99d9941b09af |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_no_probe/replica_02/progress.json` | 529 | af48cdd2a4e81a61634524b9a07b4d26d5f753ff6d2d2c66add5d412b59aee96 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_no_probe/replica_02/state.csv` | 130079 | 8e99faccec017b61857d6b5fc33d141b2978f96ef953c277a63da33ee5a8c2f7 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_no_probe/replica_03/progress.json` | 529 | 44370c587cabcffbbc253524f5f6f371612257271c8a43e2ec9d49d82fc7b891 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_no_probe/replica_03/state.csv` | 130033 | d2d1da75904c1336f67e7ea8b85a518627cc0715e9a793980b363d505fd6b845 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_probe/replica_01/progress.json` | 526 | 7b39a8b9a2bf4cf999a59a4f196799e06855d40415ec4d062a887844dcb3060c |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_probe/replica_01/state.csv` | 130011 | d2287f4d7940d3454dbf24f1e9e27395698b1f551deb9f63142b501641fc61d4 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_probe/replica_02/progress.json` | 526 | 39aa62c60437f0e91885be797991f4109ee605df3cce30d629ac567d0a9d5a46 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_probe/replica_02/state.csv` | 129956 | 4c4965df86bca9f9746223309bf84d8d706d4eef631d10378607f5e48e12bd84 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_probe/replica_03/progress.json` | 527 | c1b5933e94e20d3ad351a0c1571d9de58eedce938e37dda87def6f26b36b2ba6 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0027__candidate_probe/replica_03/state.csv` | 130035 | 5e4d0b554a1bc55ba4e2af18888400c80dc1c4cc53e0453f7a87cb54f2c9bd1c |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_no_probe/replica_01/progress.json` | 530 | 259aa06164dd2a034e1fc56e51a26a8cf59c9e712d318b96210545e7769194ad |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_no_probe/replica_01/state.csv` | 128410 | 1bc33411b5909672aade4ac28ac39cf58f343a0341ea5bb4fa914a6400147216 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_no_probe/replica_02/progress.json` | 528 | ed44e3f4dadb21de2de495a162d2e4f50b8cfff0b413290b2a274b4197a07969 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_no_probe/replica_02/state.csv` | 128399 | 74631f104f3c4aed8a14bffbc762747c0af9843a97cadf280b895e003aff598f |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_no_probe/replica_03/progress.json` | 530 | 800003d1c0781ad7eb47abb904e3bdf61d5977af6cf4c3d78ce4b11421fdee39 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_no_probe/replica_03/state.csv` | 128377 | 533c1ddfc2367e2b3daf388f3ac5782a0e4f887315aad620541a9f7037616176 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_probe/replica_01/progress.json` | 526 | 57b9393a01db6886e3413c881eda7d4f22eb169a05d0558d418fa6750f01b196 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_probe/replica_01/state.csv` | 128367 | 81f43da7d379c87c61c81c9c6603e75c60f1103b6a5a8698f484a399b0b847a9 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_probe/replica_02/progress.json` | 526 | b73063a911b7207957036085ddb16bf30e34c581a0f268378f22248545f00050 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_probe/replica_02/state.csv` | 128333 | 04f7c68b9e09df0c6bdacf0e5713c74e8eccc2603f7afb8e428a43bf635d52b7 |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_probe/replica_03/progress.json` | 527 | bdd8eb02109eda6e4c6e91fea06f5921bc3c0f2553402aaa79d845d80c8fb1ad |
| `C:/projects/PACER_STAGE4_MD_backup/production/PACER0073__candidate_probe/replica_03/state.csv` | 128387 | 1291ba27b0e1d4b4bd19ffc303e7326c4d4a7abe971fb9fd0761100fe63e2d16 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__apo/replica_01/audit.json` | 1547 | 5012a18e04e6db3e28814a37ad3aa13632b6db3fc2a65442e654654102f72a17 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__apo/replica_01/stage_01_k_500.csv` | 1354 | e455f09d714e9db6ea3d692c584cb43279e3d27335bc0473d9ed30007daa0d22 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__apo/replica_01/stage_02_k_100.csv` | 1352 | a41cdc26cbaaae2fcba34f1b84a1b5382141d8eb6c46d382e7db3aa0a99a2996 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__apo/replica_01/stage_03_k_10.csv` | 1355 | 2fc03b76a62c259457a955889f5e8e2b447807807f92ee072edd1702bafd96fc |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__apo/replica_01/stage_04_k_0.csv` | 1360 | 8e7718dd1ef50bc809aef39f71012f19df257c84abd2d97685f03d37bfd984ef |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__probe_only/replica_01/audit.json` | 1557 | dbd27eeca54a6f3aaa0b38a7af235d15ee9a866c1ff248b548ab3f92e4c45487 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__probe_only/replica_01/stage_01_k_500.csv` | 1360 | 645b8c2b3d46cb22ad449caf01ea0bc010563ffe454c5505f5f2e04803d562a8 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__probe_only/replica_01/stage_02_k_100.csv` | 1357 | 0678bd6b6685aa2918117b126ba020c4ad821e3e6bf85d19418c5cf0e95e6180 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__probe_only/replica_01/stage_03_k_10.csv` | 1357 | 2c6f0a2dbb7ec10020121e330629f0f22b375d0d916f37664b630827b6af0eb5 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster0__probe_only/replica_01/stage_04_k_0.csv` | 1351 | 256c3402960533c71c67a47a5715b5c562dd91c7464ce51f065c4e83aae7cc83 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__apo/replica_01/audit.json` | 1561 | f3491a35be2104729baad0736a17ed2ccc3e7e45aa6ae1975b9bdef2b7f921f4 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__apo/replica_01/stage_01_k_500.csv` | 1370 | 1de2268c089522009cbd5fb14c1403ff7b4799b4203b677248400b4f70011aa4 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__apo/replica_01/stage_02_k_100.csv` | 1372 | 68da973e44eca9a2b18c1463c0e1b4e93ee93e64aba99a222857b02140a8b71e |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__apo/replica_01/stage_03_k_10.csv` | 1364 | 8e92e6fa12654bc5b4ddf2a47fccc4a5aa1476625a3e08390a03b337673b8f2f |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__apo/replica_01/stage_04_k_0.csv` | 1372 | a7e6de38bb7888abeb46d24fe16180a747bceab1b90c7830bd793dd1cefa6f10 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__probe_only/replica_01/audit.json` | 1567 | 6b07ab3fdf15ed125db1f5479fba5b76855c802876f84b677eda2e00a33202a7 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__probe_only/replica_01/stage_01_k_500.csv` | 1368 | 965cc4c10d4b8bb9423265862e85513efe0804a3c4c866bba196d382aed3b43a |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__probe_only/replica_01/stage_02_k_100.csv` | 1362 | cce4e448191f5b17b4346fd4500d32da93f545c340a4532ca5c0de6b7d130b30 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__probe_only/replica_01/stage_03_k_10.csv` | 1374 | 7ae455ed36d8e1d34e4a4df48b40bb74c4e366995c8a7b5979743720d663a835 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster4__probe_only/replica_01/stage_04_k_0.csv` | 1371 | 0d5afb2f28e547e22523ba40316d5e592a1a46b71ffe4620c80da80e4aa498a8 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__apo/replica_01/audit.json` | 1560 | d6f121a56944cd325e3a950d4dcd1f5eb95877acc01e420d2dc54c4d06af0739 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__apo/replica_01/stage_01_k_500.csv` | 1369 | a1104d70ae7d1631e623c4f5916dd62dece856c903a9e99d95ff709bf1d22f45 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__apo/replica_01/stage_02_k_100.csv` | 1367 | ebed69b5360370eb60d545226714c927af593f2ed8de59993e3ae4d9b3d36a08 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__apo/replica_01/stage_03_k_10.csv` | 1374 | c09164ad36fa20b9ee09980e17e2fed36d288e50c32a5d7d417f2cd043b8a7b3 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__apo/replica_01/stage_04_k_0.csv` | 1372 | 90aed0bed7d2c705e4e032ac0471bf2a14c92ebeaf87ff972385b88bbcfb6885 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__probe_only/replica_01/audit.json` | 1568 | 1db9ec0d39c03268c6eca5f5c6c89e9a92c49d260c496ae44ea17b13a0905147 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__probe_only/replica_01/stage_01_k_500.csv` | 1374 | cdd939d5e664c488d81adc76cae5976e043ebf8ceeda38430247cecfa054dbc0 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__probe_only/replica_01/stage_02_k_100.csv` | 1366 | 6a57e1763b2fe066bee2a017672cc1d374cdd023911b3272e2d7c15c13dfb3fd |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__probe_only/replica_01/stage_03_k_10.csv` | 1363 | a13d6a1692a214c75e488cfa2ccf7608260a62bbc0ab8421c16b2fe0777738b2 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/cluster9__probe_only/replica_01/stage_04_k_0.csv` | 1364 | 39b0a3f3a06bd0fa83f5e8ee3e6209d7502fc051489de82a816d79369d11d715 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_no_probe/replica_01/audit.json` | 1577 | be4f6fe43c8b8de87f39cbce578d13433910eff604a4fc7e1001eb16a428e7a9 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_no_probe/replica_01/stage_01_k_500.csv` | 1362 | 5f089b2802277bbc7b7bd68c41ba7f27cf2667ea914aa269a6134a57c3a1ea7e |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_no_probe/replica_01/stage_02_k_100.csv` | 1369 | 118678a3fcec298bd9064fb5630a9e12546d07e359a62877daded32236033177 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_no_probe/replica_01/stage_03_k_10.csv` | 1372 | efe20fb2f8915e4191fd6139c80ef18fb93fdf6c5fe30aa5d7f35d657f3c69b9 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_no_probe/replica_01/stage_04_k_0.csv` | 1368 | bd10f1d7259f9a8a37409002399a0b6953a8e85359637629ab4cce455d5cf619 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_probe/replica_01/audit.json` | 1572 | 00103fbfb55aaac1fdab550a4d3b22ce079cb4fbb6ea36c64cd6ae3ab6e3288d |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_probe/replica_01/stage_01_k_500.csv` | 1368 | 2a1428a3b550c34c8705665644a6677c5035113f53cf7ff895cad06f0ff292b8 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_probe/replica_01/stage_02_k_100.csv` | 1373 | ccccb3043523d8a771b1c5cb6924659feef195484c88c69b0f0d303d3dc1c2dd |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_probe/replica_01/stage_03_k_10.csv` | 1369 | 9c0a7e57c72cf61b27127ad6eaa35ea14f1024c2d53adb751b445f55d34c9702 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0010__candidate_probe/replica_01/stage_04_k_0.csv` | 1379 | 118bbd9b70334bc5ea4a561c42256eb159ee1834ff98bc331e4e28b21107b99e |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_no_probe/replica_01/audit.json` | 1575 | 762b07d351120173863c9da901aa9a5a6f630e4bcbd8ee77177ef8bd0efd7028 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_no_probe/replica_01/stage_01_k_500.csv` | 1369 | 7124e9bddea6d26d3cbf0d89970ce26eaeb437d9ed8e215ac0b25fe5931e9d34 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_no_probe/replica_01/stage_02_k_100.csv` | 1375 | ef4236278b7e46e5e9fa8ceebb784cd369cba709dbc0019b787f499bdf160c2f |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_no_probe/replica_01/stage_03_k_10.csv` | 1368 | cd399e37aed5c463b7eee4aae5d4bf8da24920a5dd4f12c83551af9343c1ad7f |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_no_probe/replica_01/stage_04_k_0.csv` | 1366 | e88ee696e4bece875b211b9cf7a7eae4962fd6980c47552c626cdbccd9091667 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_probe/replica_01/audit.json` | 1571 | 901b12522b50524775c51023e90fbca9d853538178f8a8f2f041ccb600c6c971 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_probe/replica_01/stage_01_k_500.csv` | 1366 | 98861df97ee55eb0db8f42e8c3a98477507565272994a8d187fcf5d58ecb9111 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_probe/replica_01/stage_02_k_100.csv` | 1377 | b2136e63e0b3fbbbe27673ee2d1071f8b35baa12770d4e013f70b866ae7578e5 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_probe/replica_01/stage_03_k_10.csv` | 1359 | db2a2cc038f8492df0b756b9e52c391d56e2b59eba27658539eddd0446ac18a5 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0027__candidate_probe/replica_01/stage_04_k_0.csv` | 1371 | 4617a61af814beaa82721d8f5d7a713db9200ec3136714f03a24e1c0f62166e7 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_no_probe/replica_01/audit.json` | 1565 | ac74feb69e9d22f786952f4fbd3ef08f97976827b66dafd1835fcb44d50e0418 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_no_probe/replica_01/stage_01_k_500.csv` | 1359 | 97ac42abb856633fe1daf9319abbd365021c925ae02d4b82d30cbf5b3ab3480f |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_no_probe/replica_01/stage_02_k_100.csv` | 1350 | 525bf3ff3fa5c3866ce05b03914f92e9c004bc74ff60e8895a738fb1b296f058 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_no_probe/replica_01/stage_03_k_10.csv` | 1350 | 4d6276e938d998724b00ceb15b84669e313f80368acca68269fbf19ad1c6ce99 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_no_probe/replica_01/stage_04_k_0.csv` | 1358 | 031ac636988bbc9d1086c603da2e3401f9c10fa7438fb25f05dd13a2920b154c |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_probe/replica_01/audit.json` | 1563 | 7486c70aea4893925e009e4afec8919138c67f6d72bcd10d7019b4d8bb0a9e99 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_probe/replica_01/stage_01_k_500.csv` | 1353 | 41aa42810c4e020dceb95f6649f3e6414b6eca5745a3cb797af9cb0353683661 |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_probe/replica_01/stage_02_k_100.csv` | 1349 | 7dcfb3b10c7fba889faffddcd5daa18b3a2c603c2c4084b46e75d9d28c5ee0be |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_probe/replica_01/stage_03_k_10.csv` | 1347 | d8ffacb5e427ed15e84c964fa2e76e4f72a2681ed4c3f07159662e380a30aede |
| `C:/projects/PACER_STAGE4_MD_backup/restraint_release/PACER0073__candidate_probe/replica_01/stage_04_k_0.csv` | 1353 | 2fa61cd3c756947d2c1aee0dd9e85d60db938bce79c32e913a0ff03af481d325 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster0__apo/audit.json` | 1187 | dc2fe02848bc47f35a23b652f1c0c88acf9a92ecdc83aff1450a55bedd1a4a6a |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster0__apo/preflight.json` | 576 | cab3aad6b011b12786fb74b75e51e1b173229e07de5ee7c01ad5984f6b29197e |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster0__apo/refined_source_sha256.txt` | 193 | 4e06cd49e18fb0d680ee547e186682642835436e8e6b181a43ca12930237bfe7 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster0__apo/state.csv` | 1153 | a43b7694e866636811e2ab378fe5d5daba1a260c0a5221f50a9546fadc3681c4 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster0__probe_only/audit.json` | 1194 | ca2b3514a1b28a0163423f13e4c0352cba747a27e46f4a647b530eae5600369c |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster0__probe_only/preflight.json` | 585 | d2b5876aad3783536a150a9acca7a832d6cd7c4d6e37b893a422b1797a4b97a0 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster0__probe_only/refined_source_sha256.txt` | 193 | 34b13850fd7401386ae4fc8a346266861d3534168637c1740db86e7d7551f122 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster0__probe_only/state.csv` | 1152 | be29d3328410a0856188a1d9263a8db21972040181424d3bcc1bfff35319d51b |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster4__apo/audit.json` | 1187 | 3229664f8756311f1b18395365a3c6ad2fc45f48e6c7487855d6fb6dbb9235a6 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster4__apo/preflight.json` | 577 | 9b3e74a2ef28669d8e04e370eb56c68e42d1f4f0d163d6e269f2980f060831e5 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster4__apo/refined_source_sha256.txt` | 193 | ca93ec5eb26c8dbfe7e0baa4620e661226a1581c8734e03e56bbbc1fac1385c3 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster4__apo/state.csv` | 1162 | d3ac85905f11c5313910fb517f0e10000e01f67e51eb3c84c9fd4dfcfa7573b5 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster4__probe_only/audit.json` | 1195 | a72e8532274a0395a332daeafe666046da68823229fa2e9c7f6c28536dbdcb74 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster4__probe_only/preflight.json` | 584 | 19b4f29fbf79766f0f14030924a656444d45196885cf7d299762bc1b77098473 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster4__probe_only/refined_source_sha256.txt` | 193 | f09607a2a21a6655738d56f12024f99886943a388f16089e648ec0f15639f662 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster4__probe_only/state.csv` | 1164 | 7aa02b45f495fb2d5283a04c73b7faff58a995355f2802ec0e648ca344258ca4 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster9__apo/audit.json` | 1186 | b8e552c102aa905dbf8ac34abac12457baf8ba25062a19b0210b829079ee59d2 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster9__apo/preflight.json` | 576 | 58c3c4b60d889f09f4fae8616d438ad0f9c3f9be211833ade44cf704c50dbc18 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster9__apo/refined_source_sha256.txt` | 193 | 68c9997facf96b0c6d503c31d08b82ca0a7955863686494286fe1a2edf3b6354 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster9__apo/state.csv` | 1164 | 80e470b100e009cf95a0adff2f4c259bfb07a04ff3a9e71b68a0ebe3568ac144 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster9__probe_only/audit.json` | 1195 | 5fedb98bb97750e85f85cfcc523cbcf244b614e531a77dc2c75372f2594e0465 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster9__probe_only/preflight.json` | 584 | 934546d0569e3387498f36930a01b4556f175e2b2575ef288e733b0f79e4299c |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster9__probe_only/refined_source_sha256.txt` | 193 | 008178683bd68279752c1c525bbddd1924181e3cd1e8f22fae544cdb30e38cfb |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/cluster9__probe_only/state.csv` | 1164 | 231adb9cc73e65d9f8c15c4a4969e6962eeb9110548df6255fcf47a34107f856 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0010__candidate_no_probe/audit.json` | 1203 | 53a1bb2e10b241217765e982e559cd50429d3aac400083f56937b21d77037abc |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0010__candidate_no_probe/preflight.json` | 593 | 5787d2f10dc2f368ab1e0eabe0b2daa2b71fc7c2c2fdcc868fca0ac629fdd3a3 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0010__candidate_no_probe/refined_source_sha256.txt` | 193 | af28dc8bced44855bf4747111cb5272158970852f953c60a29acb81578144a91 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0010__candidate_no_probe/state.csv` | 1153 | fb06a72f4ea408bf1826a723413865dd68e89f667924523d53ef864e7e4c7b4d |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0010__candidate_probe/audit.json` | 1201 | 5fa8454a3fdbbf3a849246332235ea080b2e395a10486375b1cde92a92b3d2b4 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0010__candidate_probe/preflight.json` | 591 | 098f2909854ba428ad0bd3b459d388941fe03dff6554091ed2dabbba1a28461b |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0010__candidate_probe/refined_source_sha256.txt` | 193 | cb3377ad64f23b7fa7d085b8fb4c7b05134466a314864f7a834a6d896e023ebb |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0010__candidate_probe/state.csv` | 1150 | 31e9e53054dd1abf0ff744b0a9ab86644f371e938ded4a469b95f1897ebe2f6c |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0027__candidate_no_probe/audit.json` | 1204 | 0d334cb51030b85070fd1559cfe5831ece159e0f67e8d24155ba83bdb5eebc28 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0027__candidate_no_probe/preflight.json` | 593 | 25a9d374e6c946348c0ed615884d8913bc5d5fd64efc1e8c8ab51ff6785850d1 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0027__candidate_no_probe/refined_source_sha256.txt` | 193 | 2b3b3bd086ec87a35f46943d84008f72c8b039f291b475eb977e955b53ad4f8b |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0027__candidate_no_probe/state.csv` | 1165 | c898014e51608b68bbccd2e0332bc892948124b116c38dc2034889e30405ce15 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0027__candidate_probe/audit.json` | 1199 | 4f524428ac346742f7c2120f9ac327aee5fc4b5ad3a0af8d9c0fbcfddc4c24ab |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0027__candidate_probe/preflight.json` | 590 | 0e1a264dafa90f6bb187ee338b666163bdd58477cb5c2bda185260c8351d99f0 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0027__candidate_probe/refined_source_sha256.txt` | 193 | cedcb90e3845a344107d97458cd3cfcb3428a5a55d60080a376ae773f2272cac |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0027__candidate_probe/state.csv` | 1173 | be117e7f6c8632ef11aa45089c3f3b7544f294b107b8e1559702319870819f32 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0073__candidate_no_probe/audit.json` | 1200 | 58d69c5c5129d5448b7b06708984d4656f1530915d8055c878a408d165255614 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0073__candidate_no_probe/preflight.json` | 594 | 52825e2e9ad86c063e0071d85de2e0f037f8d95d72b3c8f452c97212e7ce0034 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0073__candidate_no_probe/refined_source_sha256.txt` | 193 | 319bc6a88166bbc373505bb7f10b1e74fd444c6b12cb6c5f41504fd51d7ab417 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0073__candidate_no_probe/state.csv` | 1158 | 4eb06e44c0eab65e5e1d0d25cb0716e743c0fe3a3275193c96a53dcf2e70ad33 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0073__candidate_probe/audit.json` | 1199 | f1564971efb1ddc9792183db9916556cb8cb2abded7a9fb97e976887ac031b5d |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0073__candidate_probe/preflight.json` | 591 | 6a91e136a5f2787594484f0129b17876ccc043c16f83accad0c5c8e9500d1ad8 |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0073__candidate_probe/refined_source_sha256.txt` | 193 | 23274178422627a57c1cbdfc3e26e03d9be4357d5317306d5c4917794d82eedb |
| `C:/projects/PACER_STAGE4_MD_backup/smoke/PACER0073__candidate_probe/state.csv` | 1148 | 775dc405acb9ab09d95399d327f5ac7eca55c2d608b6125f945f73d25130c006 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster0__apo/build_audit.json` | 92992 | 3bbeef4c47c8647b610e3727ca956895b76f3eedb38438f7b56faf4fbb8513aa |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster0__apo/build.log` | 131 | bbcd0505973d5d4adf2dee71af1094f656b947c0f3e6490946b7f1988f8249ab |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster0__apo/parameters.json` | 3 | ca3d163bab055381827226140568f3bef7eaac187cebd76878e0b63e9e442356 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster0__probe_only/build_audit.json` | 93014 | 90a52f3dc4fbd102b3af8c086fc090b3ad83a81fc5e91961a5d84ac77b17be8d |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster0__probe_only/build.log` | 132 | 4d40e955c0245f091054932d18f60bb4a3e0d0bd7e796d6c2d1f61696d3ce3c9 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster0__probe_only/parameters.json` | 1243 | ef830153b01c838d2d4d9f961e5a6cfccae09c190ef089d14571db9a93c413f9 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster4__apo/build_audit.json` | 92993 | dcbef8a437f615538499a453fa73be76e147cc756e9a2ecf9c6c01bbc92df42d |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster4__apo/build.log` | 133 | f340873f597ea2f6cae1dab003df5c1ff63b2e01d423d55737ff7077f7fb9c48 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster4__apo/parameters.json` | 3 | ca3d163bab055381827226140568f3bef7eaac187cebd76878e0b63e9e442356 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster4__probe_only/build_audit.json` | 93014 | 9b7926ea91484ba768fb25585b445b759b2d66213056de1e6306543145fbb478 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster4__probe_only/build.log` | 133 | a061ad9bbfe55a4555f4ff69b6eaff1ff9e7db037a695eb2489a53dfd2d84b81 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster4__probe_only/parameters.json` | 1243 | ef830153b01c838d2d4d9f961e5a6cfccae09c190ef089d14571db9a93c413f9 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster9__apo/build_audit.json` | 92995 | 4651b04d4202dde692c35e9637200f6e7c9091ab79ba55136c787ba579308613 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster9__apo/build.log` | 133 | 5e2e0e2e6d0a7aebc5e6fec7d719579de0e9c53f6dc7af30c7fc85a140661f42 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster9__apo/parameters.json` | 3 | ca3d163bab055381827226140568f3bef7eaac187cebd76878e0b63e9e442356 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster9__probe_only/build_audit.json` | 93014 | e22812aaf1559436cc3ffae933827d45a420e196c8e47b0a82618f843052db94 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster9__probe_only/build.log` | 131 | c6532c42df8d29c1f3fb3979bb0d38a6c80655b9bef693b75f21be98fc74ec4c |
| `C:/projects/PACER_STAGE4_MD_backup/systems/cluster9__probe_only/parameters.json` | 1243 | ef830153b01c838d2d4d9f961e5a6cfccae09c190ef089d14571db9a93c413f9 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0010__candidate_no_probe/build_audit.json` | 93030 | bc300a5620e62bf570157d15e80613d744cf5958b31204dc9181ac434ee02bb7 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0010__candidate_no_probe/build.log` | 133 | ee7e87be44dbd5bf7b367549d9c9ba63d4f9b9716862a2d8089b9282a996e886 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0010__candidate_no_probe/parameters.json` | 2369 | c10f2e550f171e6254c5bded5c6940736bcb50b3e59539b54969435abe895f23 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0010__candidate_probe/build_audit.json` | 93039 | 45ba6b2cd47984e4b07387f6cfde36c1cea621d1ecbd6a1da7073be3b8989185 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0010__candidate_probe/build.log` | 133 | bcbca26deafbd625ccb4aab750fdb995423441fa9f02a37abfac08f6c65cc6a7 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0010__candidate_probe/parameters.json` | 3609 | c8d89482ef61a5c23042f6d8e510328f9ec999909c4edf8d51be40f661a9759c |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0027__candidate_no_probe/build_audit.json` | 93027 | 4ac77bcbbaa51a0872de07c20d8d96895da7ffc706ab2fe8916a190f35365260 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0027__candidate_no_probe/build.log` | 132 | 88c98dbffa1bafa3e8b7b91112b451b8c2d818dee15d50b60920d547fa1b29d7 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0027__candidate_no_probe/parameters.json` | 2696 | 81b0b12eb66ce350d6d082f1da84b14dfa4f6de4b2975aa7db9daac2ec72a150 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0027__candidate_probe/build_audit.json` | 93036 | 5644b2fbf1036cb31603e0593352a6f49eb850f1f16e01120dc14bffa5df574a |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0027__candidate_probe/build.log` | 132 | d24e982399ad0bd21fdf509ce1175ebc94da69b831553cccf7a12191dbc99493 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0027__candidate_probe/parameters.json` | 3936 | 72b3c7d7ee2594badf40fdfaf21bf326912c122c1e4fe18f311a4847ad2ecdb0 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0073__candidate_no_probe/build_audit.json` | 93028 | 384b801a388cd6fd36ed0521b5d1d8cb71ac8548e1255e7aad3f72b7a1a4d030 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0073__candidate_no_probe/build.log` | 132 | e492498edfcec46030336fcd146dad62649817df449b978964f055bc58dfac77 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0073__candidate_no_probe/parameters.json` | 2105 | f2ec9985b37d59a3a2b83fa42f0bb6fb27cc91be747bc196c87a19a3d6f2c0d3 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0073__candidate_probe/build_audit.json` | 93037 | e0d9758368a4d51b0e5e50735573157bce527e69ac7e70fc2e5fd24b7b278617 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0073__candidate_probe/build.log` | 132 | 0bfbb90f3e01eac8423e255d754518af08751d944f31c3602ad38838b794f324 |
| `C:/projects/PACER_STAGE4_MD_backup/systems/PACER0073__candidate_probe/parameters.json` | 3345 | c5a1b8593e41cdda259d9ae8496e0fd80f2818b931ad8262eafc9a073e73895b |

## 附件 C：指定 git log --all --stat 完整原输出

为 Markdown/Git 空白检查仅移除行尾空格；提交、消息、stat 和 name-status 条目完整保留。

```text
commit b78148e9c283aa834118f99ac5a4feb613885e28
Author: VergilVolk <barneygong329@gmail.com>
Date:   Sat Oct 3 21:17:36 2026 +0800

    Finalize Stage4 prospective PACER-FKG evaluation

 .gitattributes                                     |      7 +
 README.md                                          |     24 +
 project/PROJECT_STATUS.md                          |     38 +-
 project/docs/PACER_STAGE4_FINAL_HANDOFF_v01.md     |    180 +
 project/docs/START_HERE_FOR_COLLABORATORS_V01.md   |     14 +-
 .../run_stage4_prospective_phase1_bs256_v01.py     |    240 +
 ..._stage4_prospective_phase2a_frozen_apply_v01.py |    149 +
 ..._stage4_prospective_phase2b_graph_region_v01.py |    186 +
 .../pacer_fkg_v02/stage4_prospective_common_v01.py |    505 +
 .../STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md    |    671 +
 .../STAGE4_PROSPECTIVE_SCIENTIFIC_SUMMARY_v01.json |  11700 +
 .../input_audit/INPUT_AUTHENTICATION_v01.json      |  80423 ++++
 .../phase1_bs256/PHASE1_FREEZE_RECEIPT_v01.json    |    420 +
 .../phase1_bs256/SMOKE_RECEIPT_v01.json            |    245 +
 .../PACER0010__candidate_no_probe__replica_01.json |   6197 +
 .../PACER0010__candidate_no_probe__replica_02.json |   6197 +
 .../PACER0010__candidate_no_probe__replica_03.json |   6197 +
 .../PACER0010__candidate_probe__replica_01.json    |   6197 +
 .../PACER0010__candidate_probe__replica_02.json    |   6197 +
 .../PACER0010__candidate_probe__replica_03.json    |   6197 +
 .../PACER0027__candidate_no_probe__replica_01.json |   6197 +
 .../PACER0027__candidate_no_probe__replica_02.json |   6197 +
 .../PACER0027__candidate_no_probe__replica_03.json |   6197 +
 .../PACER0027__candidate_probe__replica_01.json    |   6197 +
 .../PACER0027__candidate_probe__replica_02.json    |   6197 +
 .../PACER0027__candidate_probe__replica_03.json    |   6197 +
 .../PACER0073__candidate_no_probe__replica_01.json |   6197 +
 .../PACER0073__candidate_no_probe__replica_02.json |   6197 +
 .../PACER0073__candidate_no_probe__replica_03.json |   6197 +
 .../PACER0073__candidate_probe__replica_01.json    |   6197 +
 .../PACER0073__candidate_probe__replica_02.json    |   6197 +
 .../PACER0073__candidate_probe__replica_03.json    |   6197 +
 .../manifests/full/cluster0__apo__replica_01.json  |   6197 +
 .../manifests/full/cluster0__apo__replica_02.json  |   6197 +
 .../manifests/full/cluster0__apo__replica_03.json  |   6197 +
 .../full/cluster0__probe_only__replica_01.json     |   6197 +
 .../full/cluster0__probe_only__replica_02.json     |   6197 +
 .../full/cluster0__probe_only__replica_03.json     |   6197 +
 .../manifests/full/cluster4__apo__replica_01.json  |   6197 +
 .../manifests/full/cluster4__apo__replica_02.json  |   6197 +
 .../manifests/full/cluster4__apo__replica_03.json  |   6197 +
 .../full/cluster4__probe_only__replica_01.json     |   6197 +
 .../full/cluster4__probe_only__replica_02.json     |   6197 +
 .../full/cluster4__probe_only__replica_03.json     |   6197 +
 .../manifests/full/cluster9__apo__replica_01.json  |   6197 +
 .../manifests/full/cluster9__apo__replica_02.json  |   6197 +
 .../manifests/full/cluster9__apo__replica_03.json  |   6197 +
 .../full/cluster9__probe_only__replica_01.json     |   6197 +
 .../full/cluster9__probe_only__replica_02.json     |   6197 +
 .../full/cluster9__probe_only__replica_03.json     |   6197 +
 .../PACER0073__candidate_probe__replica_01.json    |   5201 +
 .../PHASE2A_FREEZE_RECEIPT_v01.json                |    733 +
 ...ndidate_no_probe__replica_01__SIGNED_DRIFT.json |    415 +
 ...ndidate_no_probe__replica_01__STATE_MOTION.json |    415 +
 ...ndidate_no_probe__replica_02__SIGNED_DRIFT.json |    415 +
 ...ndidate_no_probe__replica_02__STATE_MOTION.json |    415 +
 ...ndidate_no_probe__replica_03__SIGNED_DRIFT.json |    415 +
 ...ndidate_no_probe__replica_03__STATE_MOTION.json |    415 +
 ..._candidate_probe__replica_01__SIGNED_DRIFT.json |    415 +
 ..._candidate_probe__replica_01__STATE_MOTION.json |    415 +
 ..._candidate_probe__replica_02__SIGNED_DRIFT.json |    415 +
 ..._candidate_probe__replica_02__STATE_MOTION.json |    415 +
 ..._candidate_probe__replica_03__SIGNED_DRIFT.json |    415 +
 ..._candidate_probe__replica_03__STATE_MOTION.json |    415 +
 ...ndidate_no_probe__replica_01__SIGNED_DRIFT.json |    415 +
 ...ndidate_no_probe__replica_01__STATE_MOTION.json |    415 +
 ...ndidate_no_probe__replica_02__SIGNED_DRIFT.json |    415 +
 ...ndidate_no_probe__replica_02__STATE_MOTION.json |    415 +
 ...ndidate_no_probe__replica_03__SIGNED_DRIFT.json |    415 +
 ...ndidate_no_probe__replica_03__STATE_MOTION.json |    415 +
 ..._candidate_probe__replica_01__SIGNED_DRIFT.json |    415 +
 ..._candidate_probe__replica_01__STATE_MOTION.json |    415 +
 ..._candidate_probe__replica_02__SIGNED_DRIFT.json |    415 +
 ..._candidate_probe__replica_02__STATE_MOTION.json |    415 +
 ..._candidate_probe__replica_03__SIGNED_DRIFT.json |    415 +
 ..._candidate_probe__replica_03__STATE_MOTION.json |    415 +
 ...ndidate_no_probe__replica_01__SIGNED_DRIFT.json |    415 +
 ...ndidate_no_probe__replica_01__STATE_MOTION.json |    415 +
 ...ndidate_no_probe__replica_02__SIGNED_DRIFT.json |    415 +
 ...ndidate_no_probe__replica_02__STATE_MOTION.json |    415 +
 ...ndidate_no_probe__replica_03__SIGNED_DRIFT.json |    415 +
 ...ndidate_no_probe__replica_03__STATE_MOTION.json |    415 +
 ..._candidate_probe__replica_01__SIGNED_DRIFT.json |    415 +
 ..._candidate_probe__replica_01__STATE_MOTION.json |    415 +
 ..._candidate_probe__replica_02__SIGNED_DRIFT.json |    415 +
 ..._candidate_probe__replica_02__STATE_MOTION.json |    415 +
 ..._candidate_probe__replica_03__SIGNED_DRIFT.json |    415 +
 ..._candidate_probe__replica_03__STATE_MOTION.json |    415 +
 .../cluster0__apo__replica_01__SIGNED_DRIFT.json   |    415 +
 .../cluster0__apo__replica_01__STATE_MOTION.json   |    415 +
 .../cluster0__apo__replica_02__SIGNED_DRIFT.json   |    415 +
 .../cluster0__apo__replica_02__STATE_MOTION.json   |    415 +
 .../cluster0__apo__replica_03__SIGNED_DRIFT.json   |    415 +
 .../cluster0__apo__replica_03__STATE_MOTION.json   |    415 +
 ...ter0__probe_only__replica_01__SIGNED_DRIFT.json |    415 +
 ...ter0__probe_only__replica_01__STATE_MOTION.json |    415 +
 ...ter0__probe_only__replica_02__SIGNED_DRIFT.json |    415 +
 ...ter0__probe_only__replica_02__STATE_MOTION.json |    415 +
 ...ter0__probe_only__replica_03__SIGNED_DRIFT.json |    415 +
 ...ter0__probe_only__replica_03__STATE_MOTION.json |    415 +
 .../cluster4__apo__replica_01__SIGNED_DRIFT.json   |    415 +
 .../cluster4__apo__replica_01__STATE_MOTION.json   |    415 +
 .../cluster4__apo__replica_02__SIGNED_DRIFT.json   |    415 +
 .../cluster4__apo__replica_02__STATE_MOTION.json   |    415 +
 .../cluster4__apo__replica_03__SIGNED_DRIFT.json   |    415 +
 .../cluster4__apo__replica_03__STATE_MOTION.json   |    415 +
 ...ter4__probe_only__replica_01__SIGNED_DRIFT.json |    415 +
 ...ter4__probe_only__replica_01__STATE_MOTION.json |    415 +
 ...ter4__probe_only__replica_02__SIGNED_DRIFT.json |    415 +
 ...ter4__probe_only__replica_02__STATE_MOTION.json |    415 +
 ...ter4__probe_only__replica_03__SIGNED_DRIFT.json |    415 +
 ...ter4__probe_only__replica_03__STATE_MOTION.json |    415 +
 .../cluster9__apo__replica_01__SIGNED_DRIFT.json   |    415 +
 .../cluster9__apo__replica_01__STATE_MOTION.json   |    415 +
 .../cluster9__apo__replica_02__SIGNED_DRIFT.json   |    415 +
 .../cluster9__apo__replica_02__STATE_MOTION.json   |    415 +
 .../cluster9__apo__replica_03__SIGNED_DRIFT.json   |    415 +
 .../cluster9__apo__replica_03__STATE_MOTION.json   |    415 +
 ...ter9__probe_only__replica_01__SIGNED_DRIFT.json |    415 +
 ...ter9__probe_only__replica_01__STATE_MOTION.json |    415 +
 ...ter9__probe_only__replica_02__SIGNED_DRIFT.json |    415 +
 ...ter9__probe_only__replica_02__STATE_MOTION.json |    415 +
 ...ter9__probe_only__replica_03__SIGNED_DRIFT.json |    415 +
 ...ter9__probe_only__replica_03__STATE_MOTION.json |    415 +
 .../PHASE2B_FREEZE_RECEIPT_v01.json                |     49 +
 .../STAGE4_PROSPECTIVE_RESULTS_v01.json            | 427895 ++++++++++++++++++
 .../provenance/IMPLEMENTATION_VALIDATION_v01.json  |    817 +
 .../STAGE4_TEMPORAL_METADATA_NOTE_v01.md           |     55 +
 .../TOPOLOGY_SELECTION_REPAIR_VALIDATION_v01.json  |   2902 +
 .../tests/test_stage4_prospective_fkg_v02_v01.py   |    319 +
 130 files changed, 785733 insertions(+), 12 deletions(-)

commit 29fe47b5bdc18516e97989f5dabcff3cfebe4e1b
Author: VergilVolk <barneygong329@gmail.com>
Date:   Fri Oct 2 21:56:39 2026 +0800

    Add prospective Stage4 membrane builder and bounded CUDA smoke

 .../BUILDER_HANDOFF.md                             |  86 +++++
 .../builder_input_validation.json                  |  49 +++
 .../pacer_stage4_prospective_10ns_v01/campaign.py  | 199 +++++++++++
 .../job_matrix.csv                                 |  37 ++
 .../local_implementation_validation.json           |  25 ++
 .../master_manifest.json                           | 381 +++++++++++++++++++++
 .../build_pacer_stage4_prospective_membrane.py     | 205 +++++++++++
 project/scripts/pacer_stage4_prospective.py        | 205 +++++++++++
 project/scripts/run_pacer_stage4_membrane_smoke.py | 107 ++++++
 project/tests/test_pacer_stage4_prospective.py     |  77 +++++
 10 files changed, 1371 insertions(+)

commit 78ce1f30bc2e4855a4aa71dc6a1074a1e97a90e9
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Sat Oct 3 13:16:30 2026 +0800

    submission: standardized results.csv generator + Model Card (attachment-5 compliance)

    - generate_submission_results_v01.py: main result entry per section 4;
      emits UTF-8 results.csv with required fields (candidate_id, track,
      SMILES, key metrics, model+version, structure file, remarks);
      Stage3D-manifest mode outputs the frozen 3-candidate MD shortlist,
      default mode the 200-candidate pre-MD ranking (tested, 200 rows)
    - MODEL_CARD_V01.md: model inventory with open-source disclosure,
      finetune provenance, scope and known limitations per section 3
    - claim boundary ('confirmed PAM count = 0') embedded in every row

 project/docs/MODEL_CARD_V01.md                     |  48 +++++
 .../pacer_submission_results_v01/results.csv       | 201 +++++++++++++++++++++
 .../pacer_submission_results_v01/results.meta.json |  11 ++
 project/scripts/generate_submission_results_v01.py | 111 ++++++++++++
 4 files changed, 371 insertions(+)

commit fb550656a16db20d0025ccbdbed23b59558ef527
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 22:42:20 2026 +0800

    fix: fresh-clone crash in --create-template + doc count drift

    - final_candidate_hypotheses.csv (24-row frozen candidate roster, the
      only input of make_template) was never committed on any branch -
      force-added so a fresh clone can build the six-channel contract
    - make_template() now exits with a targeted message + recovery command
      instead of a raw FileNotFoundError traceback
    - quickstart test count corrected 25 -> 37 passed (2 skipped), with a
      note that pytest output is authoritative

    Verified both paths: missing file -> friendly SystemExit(1); present
    file -> contract written.

 project/docs/REPRODUCTION_QUICKSTART_V01.md        |  2 +-
 .../final/final_candidate_hypotheses.csv           | 25 ++++++++++++++++++++++
 project/scripts/apply_pacer_xr_candidates.py       |  9 ++++++++
 3 files changed, 35 insertions(+), 1 deletion(-)

commit 9dfaca3553ae5c175aad4b778ee4f58cea8d9669
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 20:28:19 2026 +0800

    artifacts: force-add frozen npz representations + README (gitignored by pattern)

 .../upstream_frozen_inputs/README.md               |  35 +++++++++++++++++++++
 .../upstream_frozen_inputs/cands_science2026.npz   | Bin 0 -> 478833 bytes
 .../m4r_pocket_science2026_subsets.npz             | Bin 0 -> 19491 bytes
 3 files changed, 35 insertions(+)

commit 1efe33287b4ad0d0e1e3d69f38fb9835a1baf159
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 20:28:07 2026 +0800

    artifacts: close Stage3 upstream gap - frozen 2026-side inputs

    Adds the local-only assets the router chain depended on: famaug
    scores CSV (sha-matched to handoff audit), frozen Science2026
    molecule representations for the 200 candidates, M4R-only pocket
    representation subset (10 pockets, extracted from the 69MB full npz),
    and the predock portfolio. With these, Stage3 + score recomputation
    close on any clone with CPU only - no base checkpoint, no lmdb,
    no forward inference required (README documents the minimal path).

 .../cands200_famaug_scores.csv                     | 201 +++++++++++++++++++++
 .../upstream_frozen_inputs/predock_portfolio.csv   | 201 +++++++++++++++++++++
 2 files changed, 402 insertions(+)

commit 627825172102002363e54017258f820b1c8dd8e9
Author: VergilVolk <barneygong329@gmail.com>
Date:   Fri Oct 2 17:55:09 2026 +0800

    chore: record pre-Stage4 service and runtime audit baseline

 project/core_environment_audit.json | 162 +++++++++++++++++++++++++++++
 reproduced_scores.audit.json        |  19 ++++
 reproduced_scores.csv               | 201 ++++++++++++++++++++++++++++++++++++
 3 files changed, 382 insertions(+)

commit 05228ab1bf137629579d32e564531a6c43095991
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 17:21:21 2026 +0800

    benchmark: two-layer external M4 PAM functional benchmark v02

    Layer 1 (confirmatory): 145 unique molecules / 384 endpoint units /
    13 endpoints (6 enter macro) from 5 independent sources; all static
    representations near-random macro directional skill (official +0.008,
    triplet -0.054, m4-matched -0.068, ECFP-RF -0.136, similarity +0.056)
    - honest negative evidence that binding retrieval cannot replace
    functional PAM prediction.

    Layer 2 (large screen): PubChem AID 624126 M4 PAM primary screen,
    frozen deterministic subset 21,445 molecules (1,445 active + 20,000
    hash-sampled inactive), 15,421 scaffolds, 0.76% active-scaffold overlap
    with training, median max-train-sim 0.234. 2D baselines: ROC~0.51-0.54
    but EF@0.5% up to 4.1 - a benchmark that similarity alone cannot
    saturate; DrugCLIP/docking scoring handed off to GPU host.

    Includes acquisition/build/eval scripts, two-layer figure, audit JSONs.

 .../benchmark.audit.json                           |    19 +
 .../endpoint_records.csv                           |   399 +
 .../endpoint_summary.csv                           |    14 +
 .../pacer_m4_external_functional_v02/molecules.csv |   148 +
 .../aid624126_screen_subset.csv                    | 21451 +++++++++++++++++++
 .../pubchem_aid624126_m4_pam_v01/source.audit.json |    14 +
 .../PACER_M4_EXTERNAL_FUNCTIONAL_BENCHMARK_V02.md  |   103 +
 .../benchmark_result.json                          |   844 +
 .../confirmatory_endpoint_skills.csv               |    31 +
 .../figures/fig_extfunc_v02_two_layers.png         |   Bin 0 -> 118943 bytes
 .../frozen_predictions_by_endpoint.csv             |   385 +
 .../baseline_predictions.csv                       | 21446 ++++++++++++++++++
 .../baseline_result.json                           |   106 +
 project/scripts/acquire_pubchem_m4_pam_screen.py   |   110 +
 ...d_pacer_m4_external_functional_benchmark_v02.py |   178 +
 ...e_pacer_m4_external_functional_benchmark_v02.py |   163 +
 .../evaluate_pubchem_m4_pam_screen_baselines.py    |   130 +
 .../make_external_functional_benchmark_fig_v02.py  |    68 +
 18 files changed, 45609 insertions(+)

commit eb45810ad23f8691613ba64b20bee02c62283f42
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 15:33:13 2026 +0800

    build: align submission container with frozen runtime

 Dockerfile                                        | 16 ++++++--
 README.md                                         |  8 +++-
 docker-compose.yml                                |  3 +-
 project/config/pacer_m4_release_manifest_v01.json |  2 +
 project/docs/LOCAL_API_AND_DOCKER_V01.md          | 11 +++--
 project/docs/RELEASE_VALIDATION_V01.md            | 21 +++++++++-
 project/docs/SUBMISSION_CHECKLIST_V01.md          | 49 +++++++++++++++++++++++
 project/release_audit_v01.json                    | 40 ++++++++++++------
 pyproject.toml                                    |  4 +-
 requirements-api.txt                              |  4 ++
 10 files changed, 134 insertions(+), 24 deletions(-)

commit 492ce332f4727f441253b752a5711852657aa803
Author: VergilVolk <barneygong329@gmail.com>
Date:   Fri Oct 2 15:25:24 2026 +0800

    freeze PACER-M4 runtime dependencies

 .gitignore                                         |   12 +
 environment-drugclip.yml                           |   33 +
 environment-md.yml                                 |   27 +
 environment-pacer-fkg.yml                          |   10 +
 environment.yml                                    |   20 +
 project/config/runtime_dependency_map_v01.json     | 1630 +++++
 project/docs/DEPENDENCY_MANIFEST.md                |  166 +
 project/results/dependency_audit_v01.json          | 6619 ++++++++++++++++++++
 .../results/dependency_runtime_evidence_v01.json   | 1293 ++++
 .../results/dependency_scientific_assets_v01.json  |  514 ++
 project/results/dependency_smoke_core_v01.json     |  108 +
 project/results/dependency_smoke_drugclip_v01.json |  156 +
 project/results/dependency_smoke_fkg_v01.json      |  100 +
 project/results/dependency_smoke_md_v01.json       |  128 +
 project/results/dependency_smoke_overview_v01.json |  308 +
 project/scripts/audit_runtime_dependencies.py      |  252 +
 project/scripts/check_runtime_dependencies.py      |  111 +
 requirements-pacer-fkg.txt                         |   17 +
 requirements.txt                                   |   18 +
 19 files changed, 11522 insertions(+)

commit c81cf4d26ecc4d7c7744887c7f381320a97c907f
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 15:11:29 2026 +0800

    build: freeze Python 3.9 core requirements

 README.md                                         |  10 ++
 project/config/pacer_m4_release_manifest_v01.json |   2 +
 project/docs/RELEASE_VALIDATION_V01.md            |   2 +
 project/docs/REPRODUCTION_QUICKSTART_V01.md       |  10 +-
 project/release_audit_v01.json                    |  32 ++++--
 project/scripts/audit_core_environment.py         | 120 ++++++++++++++++++++++
 pyproject.toml                                    |   2 +-
 requirements.txt                                  |  32 ++++++
 8 files changed, 197 insertions(+), 13 deletions(-)

commit 85c08f74d2b903869da91a357d7c491a18cd5b75
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 14:54:51 2026 +0800

    feat: expose PACER-M4 through guarded local API

 .dockerignore                                     |  25 +++++
 Dockerfile                                        |  22 ++++
 README.md                                         |  10 ++
 docker-compose.yml                                |  10 ++
 pacer_m4/api.py                                   | 129 ++++++++++++++++++++++
 pacer_m4/runner.py                                |  23 +++-
 pacer_m4/stages.py                                |  37 ++++++-
 project/config/pacer_m4_release_manifest_v01.json |   8 +-
 project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md   |   1 +
 project/docs/LOCAL_API_AND_DOCKER_V01.md          |  91 +++++++++++++++
 project/docs/RELEASE_VALIDATION_V01.md            |  15 ++-
 project/release_audit_v01.json                    |  72 ++++++++++--
 project/tests/test_pacer_m4_api.py                |  47 ++++++++
 pyproject.toml                                    |  18 +++
 14 files changed, 489 insertions(+), 19 deletions(-)

commit bfbf94dadd34da0b76c540d8a3a696dd0e373150
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 14:48:07 2026 +0800

    feat: add audited PACER-M4 stage runner

 README.md                                         |  11 +
 pacer_m4/__init__.py                              |   8 +
 pacer_m4/__main__.py                              |   6 +
 pacer_m4/cli.py                                   |  63 +++++
 pacer_m4/runner.py                                |  88 +++++++
 pacer_m4/stages.py                                | 165 ++++++++++++
 project/config/pacer_m4_release_manifest_v01.json |   8 +-
 project/docs/CODE_AND_ARTIFACT_MAP_V01.md         |  11 +-
 project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md   |   2 +-
 project/docs/RELEASE_VALIDATION_V01.md            |  78 ++++++
 project/docs/REPRODUCTION_QUICKSTART_V01.md       |  16 +-
 project/docs/START_HERE_FOR_COLLABORATORS_V01.md  |   2 +-
 project/release_audit_v01.json                    |  72 +++++-
 project/scripts/select_pareto_candidates.py       | 302 +++++++++++++++++-----
 project/tests/test_pacer_m4_cli.py                |  39 +++
 project/tests/test_select_pareto_candidates.py    |  42 +++
 pyproject.toml                                    |  19 ++
 17 files changed, 847 insertions(+), 85 deletions(-)

commit a8c4203c8b9cfb09f35bd8b88670d602aa266fa6
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 14:39:19 2026 +0800

    docs: establish PACER-M4 submission entry and evidence index

 .gitignore                                         |   3 +-
 README.md                                          |  86 ++++++++++
 START_HERE_ONEPROT_PACER_DC.md                     |   4 +
 project/PACER_M4_FINAL_REPORT.md                   |   5 +
 project/PIPELINE.md                                |   4 +
 project/config/pacer_m4_release_manifest_v01.json  |  35 ++++
 project/docs/CODE_AND_ARTIFACT_MAP_V01.md          | 127 +++++++++++++++
 project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md    | 155 ++++++++++++++++++
 project/docs/REPRODUCTION_QUICKSTART_V01.md        | 125 ++++++++++++++
 project/docs/START_HERE_FOR_COLLABORATORS_V01.md   | 114 +++++++++++++
 project/release_audit_v01.json                     | 179 +++++++++++++++++++++
 project/scripts/validate_pacer_m4_release.py       | 114 +++++++++++++
 .../tests/test_science2026_triplet_checkpoint.py   |  21 ++-
 project/tests/test_validate_pacer_m4_release.py    |  48 ++++++
 14 files changed, 1015 insertions(+), 5 deletions(-)

commit d1f145832907c42850baf2e189c7142454f66519
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 11:28:56 2026 +0800

    ship M4 LOTO adapters and offline inference bundle

 .../drugclip2023_m4_loto_pacer200_v01/README.md    |  59 ++++++
 .../bundle_manifest.json                           |  46 +++++
 ...drugclip2023_m4_loto_seed20260925.projection.pt | Bin 0 -> 2631059 bytes
 ...drugclip2023_m4_loto_seed20260926.projection.pt | Bin 0 -> 2631059 bytes
 ...drugclip2023_m4_loto_seed20260927.projection.pt | Bin 0 -> 2631059 bytes
 .../pacer200_2023_gpcr_loto_scores.audit.json      |  19 ++
 .../pacer200_2023_gpcr_loto_scores.csv             | 201 +++++++++++++++++++++
 .../pacer200_m4_2023_frozen_representations.npz    | Bin 0 -> 386211 bytes
 .../pacer200_m4_safe_routed_ranking.audit.json     |  22 +++
 .../pacer200_m4_safe_routed_ranking.csv            | 201 +++++++++++++++++++++
 .../pacer200_two_model_scores.csv                  | 201 +++++++++++++++++++++
 .../requirements.txt                               |   3 +
 .../scripts/build_drugclip2023_m4_loto_bundle.py   |  88 +++++++++
 .../export_drugclip2023_m4_loto_adapters.py        | 132 ++++++++++++++
 project/scripts/pacer_drugclip_router.py           |   2 +-
 .../scripts/score_pacer200_drugclip2023_m4_loto.py |  92 ++++++++++
 16 files changed, 1065 insertions(+), 1 deletion(-)

commit f35acdd9d330e88ed8c690f81986e0305f5a515f
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 10:27:55 2026 +0800

    fix router Python import example

 project/docs/PACER_DRUGCLIP_ROUTER_INTERFACE_V01.md | 6 +++++-
 1 file changed, 5 insertions(+), 1 deletion(-)

commit cd53f989ec22ef6df406b4924ceadd51567b3ea0
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 10:27:27 2026 +0800

    add deployable PACER DrugCLIP routing interface

 .../docs/PACER_DRUGCLIP_ROUTER_INTERFACE_V01.md    | 211 +++++++++++++++++++++
 project/scripts/pacer_drugclip_router.py           | 164 ++++++++++++++++
 project/tests/test_pacer_drugclip_router.py        |  58 ++++++
 3 files changed, 433 insertions(+)

commit 512f8fa866a8c187a9a4923e7de4550c1c264c56
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 10:19:51 2026 +0800

    benchmark DrugCLIP generations with GPCR-safe routing

 .../DRUGCLIP_GENERATION_BAKEOFF_PROTOCOL_V01.md    |  95 ++++++
 .../docs/DRUGCLIP_GENERATION_BAKEOFF_RESULT_V01.md |  80 +++++
 project/docs/DRUGCLIP_LOCKBOX_ACQUISITION_V01.md   |  42 +++
 ...OJECT_MAINLINE_DEVELOPMENT_AND_INFERENCE_V01.md | 182 +++++++++++
 .../drugclip_generation_bakeoff_v01/audit.json     |  44 +++
 .../generation_dumbbell.png                        | Bin 0 -> 49518 bytes
 .../macro_metric_bars.png                          | Bin 0 -> 157075 bytes
 .../macro_metrics.csv                              |  19 ++
 .../matched_docking_comparison.csv                 |  37 +++
 .../paired_delta_forest.png                        | Bin 0 -> 71967 bytes
 .../per_target_metrics.csv                         |  73 +++++
 .../per_target_roc_heatmap.png                     | Bin 0 -> 172714 bytes
 .../scaffold_bootstrap_ci.csv                      |  91 ++++++
 project/scripts/run_drugclip_generation_bakeoff.py | 351 +++++++++++++++++++++
 project/tests/test_drugclip_generation_bakeoff.py  |  21 ++
 15 files changed, 1035 insertions(+)

commit 270dd0eb5e88c0c369d493a58bd739cf118ebcdb
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 09:38:34 2026 +0800

    evidence: same-test-set M4R head-to-head 2023-triplet vs 2026-famaug

    Answering teammate challenge: the 0.69-vs-0.51 routing comparison was
    cross-protocol (illegitimate). New fair fight on the identical GaMD M4R
    screening pairs (25,249 rows, 2,301 actives), both sides never trained on
    these pairs (2023: ChEMBL triplet external_result; 2026: LOSO held-out
    predictions):

      config                  ROC    EF1%   BEDROC20  PR
      2023-triplet 10pock max 0.6673 0.39   0.1263   0.1295
      2023-triplet 10pock MEAN 0.7136 0.87  0.1802   0.1547  <- best
      2026-famaug LOSO 3-seed 0.5095 1.65   0.1118   0.0926
      2026 official raw       0.5401 1.56   0.1245   0.1070

    Verdict: 2023 base is genuinely the stronger M4 representation (ROC/PR/
    BEDROC all favor 2023-mean); 2026 edges only the noisy top-1% slice EF.
    Screening config frozen: 2023 external_result + MEAN pooling over 10 GaMD
    M4 pockets. Library ranking regenerated (top-1% overlap with famaug
    ranking only 8/285 -- genuinely independent second opinion).

 .../m4_library_oldtriplet_ranking.csv              | 28520 +++++++++++++++++++
 .../family_aug_v01/m4r_head2head_2023vs2026.json   |   101 +
 project/scripts/m4r_head2head_2023_vs_2026.py      |    54 +
 project/scripts/screen_m4_library_oldtriplet.py    |    41 +
 4 files changed, 28716 insertions(+)

commit d948297b1011707cfe54ad9dbe47df446506aac1
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 09:24:27 2026 +0800

    docs: repository organization -- mainline START_HERE entry + authoritative index

    Survey result: codex/project-integration-benchmark-v02 is the single
    superset (verified: every business branch has zero commits not in v02).
    Primary worktree fast-forwarded to v02 tip (77f44dc6).

    START_HERE_MAINLINE.md (root entry, tracked per repo convention; local
    README.md stays untracked per the upstream .gitignore) points to
    docs/REPOSITORY_MAINLINE_INDEX.md, which maps: the evidence-cascade
    mainline definition, per-line frozen conclusions (DrugCLIP fine-tuning /
    PACER-FKG Stage A+B / encoder C1-BS256 / OneProt retirement / CGDA
    closure / PACER-MCV / four-context MD), worktree statuses (geom2vec_pilot
    superseded), entry documents, next execution points, and no-retuning red
    lines. Marks the two older START_HERE_* files as archived sub-line entries.

 START_HERE_MAINLINE.md            | 24 +++++++++++++
 docs/REPOSITORY_MAINLINE_INDEX.md | 73 +++++++++++++++++++++++++++++++++++++++
 2 files changed, 97 insertions(+)

commit 77f44dc6416f4d332bd2f0ff0bafbf4a9da9231c
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 09:22:56 2026 +0800

    figs: per-module horizontal scoreboards + cascade funnel

    fig_v02_07: four same-protocol module leaderboards side by side
    (binding 13T, static gating all-near-zero, dynamic three-class,
    QSAR collapse) with explicit 'different rulers' caveat in title;
    fig_v02_08: the only legitimate overall comparison - the candidate
    funnel 200/24/5/(1-3 pending, log scale) from handoff audit v02

 .../figures/fig_v02_07_module_scoreboards.png      | Bin 0 -> 198286 bytes
 .../figures/fig_v02_08_cascade_funnel.png          | Bin 0 -> 66620 bytes
 project/scripts/make_module_scoreboards_v01.py     | 125 +++++++++++++++++++++
 3 files changed, 125 insertions(+)

commit ddc11f970f508d89a75dac8ca30c41cfcd361731
Author: VergilVolk <barneygong329@gmail.com>
Date:   Fri Oct 2 02:05:11 2026 +0800

    freeze: finalize manifest-driven PACER production pipeline v01

    Integration base: 18100f9b6b42f27e43760daa5b3dfdd028b330d3
    PACER-FKG freeze anchor: b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd
    Encoder checkpoint SHA256: b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417

    No frozen scientific numerical asset is modified by this commit. The four frozen
    v02 modules (run_phase1_bs256.py, run_phase2_calibration.py, run_phase3_fkg.py,
    path_guard.py) are unchanged; V02_FREEZE_MANIFEST.json is unchanged.

    Added:
    - project/pacer_fkg_v02/frozen_runner/  manifest-driven, apply-only Phase 1/2A/2B drivers
                                            with fail-closed spec validation and an
                                            import-only frozen re-export layer
    - project/drugclip_freeze/             frozen DrugCLIP Model A (official checkpoint)
                                            and Model B (projection-only) apply paths;
                                            dual-model decision old_top50 INTERSECT new_top25
    - project/integration_freeze/          run-spec JSON Schema, CM00734 regression spec,
                                            receipts, and the code-freeze tooling
    - project/tests/                       test_frozen_runner.py, test_drugclip_freeze.py

    Evidence:
    - frozen anchor verifies strictly: VERIFIED_FROZEN, 12 pairs, 600 blocks,
      200 calibration blocks, R2-only fitting, 41 protected v01 files
    - block constructor bit-identical to both frozen implementations (1000x50 and 400x20)
    - CM00734 payload regression: 1632 exact fields, 0 mismatched; 108 cosines with
      max |delta| 6.4e-8 and 0 sign changes; primary signal -0.326839 reproduced
    - frozen encoder replay on CUDA: relative RMSE 1.77e-7 vs limit 1e-4
    - tests: 16/16 frozen_runner, 9/9 drugclip_freeze

    Mandatory fresh-checkout preflight: run
    project/integration_freeze/tools/normalize_frozen_eol.py --apply
    before using the frozen chain from a new checkout. The v02 freeze manifest pins a
    MIXED line-ending state (run_phase2_calibration.py as LF, analyze_pacer_fkg_g2b/g2c.py
    as CRLF), so a CRLF Windows checkout fails verify_frozen_anchor() until normalised.

    Claim boundary: prospective representation-space mechanistic evidence only. No PAM
    classifier, no potency predictor, no PAM probability, no magnitude threshold.

 project/drugclip_freeze/RUNBOOK.md                 | 132 ++++++
 project/drugclip_freeze/__init__.py                |  21 +
 project/drugclip_freeze/apply.py                   |  81 ++++
 project/drugclip_freeze/dual_model.py              |  96 ++++
 project/drugclip_freeze/model_a.py                 | 117 +++++
 project/drugclip_freeze/model_b.py                 |  98 ++++
 project/drugclip_freeze/projection.py              |  92 ++++
 .../FINAL_CODE_FREEZE_CANDIDATE_REPORT.md          | 397 ++++++++++++++++
 .../FROZEN_MODULE_LEDGER_v01.json                  |  29 ++
 project/integration_freeze/RECEIPT_INDEX_v01.json  | 153 ++++++
 .../WSL_HISTORICAL_WORKSPACE_RECEIPT_v01.md        | 155 +++++++
 .../receipts/CM00734_PAYLOAD_REGRESSION_v01.json   |  88 ++++
 .../FROZEN_EOL_NORMALISATION_APPLIED_v01.json      |  85 ++++
 .../FROZEN_EOL_NORMALISATION_DRYRUN_v01.json       |  80 ++++
 .../FROZEN_EOL_PREFLIGHT_AT_FREEZE_v01.json        |  79 ++++
 .../receipts/PHASE1_ENCODER_REPLAY_v01.json        |  52 +++
 .../receipts/WSL_DIRTY_CODE_DIFF_v01.patch         | 273 +++++++++++
 .../receipts/WSL_DIRTY_FILES_BLOB_HASHES_v01.tsv   |   4 +
 .../receipts/WSL_GIT_BRANCHES_v01.txt              |   3 +
 .../receipts/WSL_GIT_BRANCH_v01.txt                |   1 +
 .../receipts/WSL_GIT_DIFF_NUMSTAT_v01.txt          | 327 +++++++++++++
 .../receipts/WSL_GIT_DIFF_STAT_IGNORE_CR_v01.txt   |   6 +
 .../receipts/WSL_GIT_HEAD_v01.txt                  |   1 +
 .../receipts/WSL_GIT_LOG1_v01.txt                  |   3 +
 .../receipts/WSL_GIT_REMOTES_v01.txt               |   2 +
 .../receipts/WSL_GIT_STASH_v01.txt                 |   0
 .../receipts/WSL_GIT_STATUS_PORCELAIN_v01.txt      | 333 +++++++++++++
 .../receipts/WSL_RESULTS_FILES_v01.tsv             | 514 +++++++++++++++++++++
 .../receipts/WSL_RESULTS_SHA256_v01.txt            | 514 +++++++++++++++++++++
 .../receipts/WSL_RESULTS_TOTAL_BYTES_v01.txt       |   2 +
 .../receipts/WSL_TOPFILES_v01.tsv                  |   6 +
 .../receipts/WSL_WORK_EXTRAS_SHA256_v01.txt        |   7 +
 .../run_specs/CM00734_REGRESSION_v01.json          | 116 +++++
 .../schemas/final_pipeline_run_spec_v1.schema.json | 309 +++++++++++++
 project/integration_freeze/tools/_shim_helper.py   |  10 +
 .../tools/diagnose_regression.py                   |  76 +++
 .../tools/normalize_frozen_eol.py                  | 109 +++++
 project/integration_freeze/tools/run_regression.py |  36 ++
 .../tools/verify_phase1_replay.py                  | 105 +++++
 project/pacer_fkg_v02/frozen_runner/__init__.py    |  21 +
 project/pacer_fkg_v02/frozen_runner/blocks.py      | 104 +++++
 project/pacer_fkg_v02/frozen_runner/cli.py         |  73 +++
 project/pacer_fkg_v02/frozen_runner/engines.py     | 181 ++++++++
 project/pacer_fkg_v02/frozen_runner/fkg_phase1.py  | 231 +++++++++
 project/pacer_fkg_v02/frozen_runner/fkg_phase2a.py | 139 ++++++
 project/pacer_fkg_v02/frozen_runner/fkg_phase2b.py | 228 +++++++++
 project/pacer_fkg_v02/frozen_runner/io_utils.py    |  36 ++
 project/pacer_fkg_v02/frozen_runner/paths.py       |  55 +++
 project/pacer_fkg_v02/frozen_runner/regression.py  | 235 ++++++++++
 project/pacer_fkg_v02/frozen_runner/spec.py        | 268 +++++++++++
 project/tests/test_drugclip_freeze.py              | 127 +++++
 project/tests/test_frozen_runner.py                | 180 ++++++++
 52 files changed, 6390 insertions(+)
```

## 附件 D：指定 git log --all --name-status 完整原输出

```text
commit b78148e9c283aa834118f99ac5a4feb613885e28
Author: VergilVolk <barneygong329@gmail.com>
Date:   Sat Oct 3 21:17:36 2026 +0800

    Finalize Stage4 prospective PACER-FKG evaluation

A	.gitattributes
M	README.md
M	project/PROJECT_STATUS.md
A	project/docs/PACER_STAGE4_FINAL_HANDOFF_v01.md
M	project/docs/START_HERE_FOR_COLLABORATORS_V01.md
A	project/pacer_fkg_v02/run_stage4_prospective_phase1_bs256_v01.py
A	project/pacer_fkg_v02/run_stage4_prospective_phase2a_frozen_apply_v01.py
A	project/pacer_fkg_v02/run_stage4_prospective_phase2b_graph_region_v01.py
A	project/pacer_fkg_v02/stage4_prospective_common_v01.py
A	project/results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md
A	project/results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_SUMMARY_v01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/input_audit/INPUT_AUTHENTICATION_v01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/PHASE1_FREEZE_RECEIPT_v01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/SMOKE_RECEIPT_v01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_no_probe__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_no_probe__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_no_probe__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_probe__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_probe__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0010__candidate_probe__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_no_probe__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_no_probe__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_no_probe__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_probe__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_probe__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0027__candidate_probe__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_no_probe__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_no_probe__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_no_probe__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_probe__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_probe__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/PACER0073__candidate_probe__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__apo__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__apo__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__apo__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__probe_only__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__probe_only__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster0__probe_only__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__apo__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__apo__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__apo__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__probe_only__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__probe_only__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster4__probe_only__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__apo__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__apo__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__apo__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__probe_only__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__probe_only__replica_02.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/full/cluster9__probe_only__replica_03.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/manifests/smoke/PACER0073__candidate_probe__replica_01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/PHASE2A_FREEZE_RECEIPT_v01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_no_probe__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0010__candidate_probe__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_no_probe__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0027__candidate_probe__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_no_probe__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/PACER0073__candidate_probe__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__apo__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster0__probe_only__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__apo__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster4__probe_only__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__apo__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_01__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_01__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_02__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_02__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_03__SIGNED_DRIFT.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/manifests/full/cluster9__probe_only__replica_03__STATE_MOTION.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2b_graph_region/PHASE2B_FREEZE_RECEIPT_v01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/phase2b_graph_region/STAGE4_PROSPECTIVE_RESULTS_v01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/provenance/IMPLEMENTATION_VALIDATION_v01.json
A	project/results/pacer_stage4_prospective_fkg_v02_v01/provenance/STAGE4_TEMPORAL_METADATA_NOTE_v01.md
A	project/results/pacer_stage4_prospective_fkg_v02_v01/provenance/TOPOLOGY_SELECTION_REPAIR_VALIDATION_v01.json
A	project/tests/test_stage4_prospective_fkg_v02_v01.py

commit 29fe47b5bdc18516e97989f5dabcff3cfebe4e1b
Author: VergilVolk <barneygong329@gmail.com>
Date:   Fri Oct 2 21:56:39 2026 +0800

    Add prospective Stage4 membrane builder and bounded CUDA smoke

A	project/results/pacer_stage4_prospective_10ns_v01/BUILDER_HANDOFF.md
A	project/results/pacer_stage4_prospective_10ns_v01/builder_input_validation.json
A	project/results/pacer_stage4_prospective_10ns_v01/campaign.py
A	project/results/pacer_stage4_prospective_10ns_v01/job_matrix.csv
A	project/results/pacer_stage4_prospective_10ns_v01/local_implementation_validation.json
A	project/results/pacer_stage4_prospective_10ns_v01/master_manifest.json
A	project/scripts/build_pacer_stage4_prospective_membrane.py
A	project/scripts/pacer_stage4_prospective.py
A	project/scripts/run_pacer_stage4_membrane_smoke.py
A	project/tests/test_pacer_stage4_prospective.py

commit 78ce1f30bc2e4855a4aa71dc6a1074a1e97a90e9
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Sat Oct 3 13:16:30 2026 +0800

    submission: standardized results.csv generator + Model Card (attachment-5 compliance)

    - generate_submission_results_v01.py: main result entry per section 4;
      emits UTF-8 results.csv with required fields (candidate_id, track,
      SMILES, key metrics, model+version, structure file, remarks);
      Stage3D-manifest mode outputs the frozen 3-candidate MD shortlist,
      default mode the 200-candidate pre-MD ranking (tested, 200 rows)
    - MODEL_CARD_V01.md: model inventory with open-source disclosure,
      finetune provenance, scope and known limitations per section 3
    - claim boundary ('confirmed PAM count = 0') embedded in every row

A	project/docs/MODEL_CARD_V01.md
A	project/results/pacer_submission_results_v01/results.csv
A	project/results/pacer_submission_results_v01/results.meta.json
A	project/scripts/generate_submission_results_v01.py

commit fb550656a16db20d0025ccbdbed23b59558ef527
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 22:42:20 2026 +0800

    fix: fresh-clone crash in --create-template + doc count drift

    - final_candidate_hypotheses.csv (24-row frozen candidate roster, the
      only input of make_template) was never committed on any branch -
      force-added so a fresh clone can build the six-channel contract
    - make_template() now exits with a targeted message + recovery command
      instead of a raw FileNotFoundError traceback
    - quickstart test count corrected 25 -> 37 passed (2 skipped), with a
      note that pytest output is authoritative

    Verified both paths: missing file -> friendly SystemExit(1); present
    file -> contract written.

M	project/docs/REPRODUCTION_QUICKSTART_V01.md
A	project/results/pacer_candidates_v01/final/final_candidate_hypotheses.csv
M	project/scripts/apply_pacer_xr_candidates.py

commit 9dfaca3553ae5c175aad4b778ee4f58cea8d9669
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 20:28:19 2026 +0800

    artifacts: force-add frozen npz representations + README (gitignored by pattern)

A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/README.md
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/cands_science2026.npz
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/m4r_pocket_science2026_subsets.npz

commit 1efe33287b4ad0d0e1e3d69f38fb9835a1baf159
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 20:28:07 2026 +0800

    artifacts: close Stage3 upstream gap - frozen 2026-side inputs

    Adds the local-only assets the router chain depended on: famaug
    scores CSV (sha-matched to handoff audit), frozen Science2026
    molecule representations for the 200 candidates, M4R-only pocket
    representation subset (10 pockets, extracted from the 69MB full npz),
    and the predock portfolio. With these, Stage3 + score recomputation
    close on any clone with CPU only - no base checkpoint, no lmdb,
    no forward inference required (README documents the minimal path).

A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/cands200_famaug_scores.csv
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/upstream_frozen_inputs/predock_portfolio.csv

commit 627825172102002363e54017258f820b1c8dd8e9
Author: VergilVolk <barneygong329@gmail.com>
Date:   Fri Oct 2 17:55:09 2026 +0800

    chore: record pre-Stage4 service and runtime audit baseline

A	project/core_environment_audit.json
A	reproduced_scores.audit.json
A	reproduced_scores.csv

commit 05228ab1bf137629579d32e564531a6c43095991
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 17:21:21 2026 +0800

    benchmark: two-layer external M4 PAM functional benchmark v02

    Layer 1 (confirmatory): 145 unique molecules / 384 endpoint units /
    13 endpoints (6 enter macro) from 5 independent sources; all static
    representations near-random macro directional skill (official +0.008,
    triplet -0.054, m4-matched -0.068, ECFP-RF -0.136, similarity +0.056)
    - honest negative evidence that binding retrieval cannot replace
    functional PAM prediction.

    Layer 2 (large screen): PubChem AID 624126 M4 PAM primary screen,
    frozen deterministic subset 21,445 molecules (1,445 active + 20,000
    hash-sampled inactive), 15,421 scaffolds, 0.76% active-scaffold overlap
    with training, median max-train-sim 0.234. 2D baselines: ROC~0.51-0.54
    but EF@0.5% up to 4.1 - a benchmark that similarity alone cannot
    saturate; DrugCLIP/docking scoring handed off to GPU host.

    Includes acquisition/build/eval scripts, two-layer figure, audit JSONs.

A	project/data/benchmarks/pacer_m4_external_functional_v02/benchmark.audit.json
A	project/data/benchmarks/pacer_m4_external_functional_v02/endpoint_records.csv
A	project/data/benchmarks/pacer_m4_external_functional_v02/endpoint_summary.csv
A	project/data/benchmarks/pacer_m4_external_functional_v02/molecules.csv
A	project/data/benchmarks/pubchem_aid624126_m4_pam_v01/aid624126_screen_subset.csv
A	project/data/benchmarks/pubchem_aid624126_m4_pam_v01/source.audit.json
A	project/docs/PACER_M4_EXTERNAL_FUNCTIONAL_BENCHMARK_V02.md
A	project/results/pacer_m4_external_functional_v02/benchmark_result.json
A	project/results/pacer_m4_external_functional_v02/confirmatory_endpoint_skills.csv
A	project/results/pacer_m4_external_functional_v02/figures/fig_extfunc_v02_two_layers.png
A	project/results/pacer_m4_external_functional_v02/frozen_predictions_by_endpoint.csv
A	project/results/pubchem_aid624126_m4_pam_v01/baseline_predictions.csv
A	project/results/pubchem_aid624126_m4_pam_v01/baseline_result.json
A	project/scripts/acquire_pubchem_m4_pam_screen.py
A	project/scripts/build_pacer_m4_external_functional_benchmark_v02.py
A	project/scripts/evaluate_pacer_m4_external_functional_benchmark_v02.py
A	project/scripts/evaluate_pubchem_m4_pam_screen_baselines.py
A	project/scripts/make_external_functional_benchmark_fig_v02.py

commit eb45810ad23f8691613ba64b20bee02c62283f42
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 15:33:13 2026 +0800

    build: align submission container with frozen runtime

M	Dockerfile
M	README.md
M	docker-compose.yml
M	project/config/pacer_m4_release_manifest_v01.json
M	project/docs/LOCAL_API_AND_DOCKER_V01.md
M	project/docs/RELEASE_VALIDATION_V01.md
A	project/docs/SUBMISSION_CHECKLIST_V01.md
M	project/release_audit_v01.json
M	pyproject.toml
A	requirements-api.txt

commit 492ce332f4727f441253b752a5711852657aa803
Author: VergilVolk <barneygong329@gmail.com>
Date:   Fri Oct 2 15:25:24 2026 +0800

    freeze PACER-M4 runtime dependencies

M	.gitignore
A	environment-drugclip.yml
A	environment-md.yml
A	environment-pacer-fkg.yml
A	environment.yml
A	project/config/runtime_dependency_map_v01.json
A	project/docs/DEPENDENCY_MANIFEST.md
A	project/results/dependency_audit_v01.json
A	project/results/dependency_runtime_evidence_v01.json
A	project/results/dependency_scientific_assets_v01.json
A	project/results/dependency_smoke_core_v01.json
A	project/results/dependency_smoke_drugclip_v01.json
A	project/results/dependency_smoke_fkg_v01.json
A	project/results/dependency_smoke_md_v01.json
A	project/results/dependency_smoke_overview_v01.json
A	project/scripts/audit_runtime_dependencies.py
A	project/scripts/check_runtime_dependencies.py
A	requirements-pacer-fkg.txt
A	requirements.txt

commit c81cf4d26ecc4d7c7744887c7f381320a97c907f
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 15:11:29 2026 +0800

    build: freeze Python 3.9 core requirements

M	README.md
M	project/config/pacer_m4_release_manifest_v01.json
M	project/docs/RELEASE_VALIDATION_V01.md
M	project/docs/REPRODUCTION_QUICKSTART_V01.md
M	project/release_audit_v01.json
A	project/scripts/audit_core_environment.py
M	pyproject.toml
A	requirements.txt

commit 85c08f74d2b903869da91a357d7c491a18cd5b75
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 14:54:51 2026 +0800

    feat: expose PACER-M4 through guarded local API

A	.dockerignore
A	Dockerfile
M	README.md
A	docker-compose.yml
A	pacer_m4/api.py
M	pacer_m4/runner.py
M	pacer_m4/stages.py
M	project/config/pacer_m4_release_manifest_v01.json
M	project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md
A	project/docs/LOCAL_API_AND_DOCKER_V01.md
M	project/docs/RELEASE_VALIDATION_V01.md
M	project/release_audit_v01.json
A	project/tests/test_pacer_m4_api.py
M	pyproject.toml

commit bfbf94dadd34da0b76c540d8a3a696dd0e373150
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 14:48:07 2026 +0800

    feat: add audited PACER-M4 stage runner

M	README.md
A	pacer_m4/__init__.py
A	pacer_m4/__main__.py
A	pacer_m4/cli.py
A	pacer_m4/runner.py
A	pacer_m4/stages.py
M	project/config/pacer_m4_release_manifest_v01.json
M	project/docs/CODE_AND_ARTIFACT_MAP_V01.md
M	project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md
A	project/docs/RELEASE_VALIDATION_V01.md
M	project/docs/REPRODUCTION_QUICKSTART_V01.md
M	project/docs/START_HERE_FOR_COLLABORATORS_V01.md
M	project/release_audit_v01.json
M	project/scripts/select_pareto_candidates.py
A	project/tests/test_pacer_m4_cli.py
A	project/tests/test_select_pareto_candidates.py
A	pyproject.toml

commit a8c4203c8b9cfb09f35bd8b88670d602aa266fa6
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 14:39:19 2026 +0800

    docs: establish PACER-M4 submission entry and evidence index

M	.gitignore
A	README.md
M	START_HERE_ONEPROT_PACER_DC.md
M	project/PACER_M4_FINAL_REPORT.md
M	project/PIPELINE.md
A	project/config/pacer_m4_release_manifest_v01.json
A	project/docs/CODE_AND_ARTIFACT_MAP_V01.md
A	project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md
A	project/docs/REPRODUCTION_QUICKSTART_V01.md
A	project/docs/START_HERE_FOR_COLLABORATORS_V01.md
A	project/release_audit_v01.json
A	project/scripts/validate_pacer_m4_release.py
M	project/tests/test_science2026_triplet_checkpoint.py
A	project/tests/test_validate_pacer_m4_release.py

commit d1f145832907c42850baf2e189c7142454f66519
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 11:28:56 2026 +0800

    ship M4 LOTO adapters and offline inference bundle

A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/README.md
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/bundle_manifest.json
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/drugclip2023_m4_loto_seed20260925.projection.pt
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/drugclip2023_m4_loto_seed20260926.projection.pt
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/drugclip2023_m4_loto_seed20260927.projection.pt
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_2023_gpcr_loto_scores.audit.json
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_2023_gpcr_loto_scores.csv
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_m4_2023_frozen_representations.npz
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_m4_safe_routed_ranking.audit.json
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_m4_safe_routed_ranking.csv
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_two_model_scores.csv
A	project/artifacts/drugclip2023_m4_loto_pacer200_v01/requirements.txt
A	project/scripts/build_drugclip2023_m4_loto_bundle.py
A	project/scripts/export_drugclip2023_m4_loto_adapters.py
M	project/scripts/pacer_drugclip_router.py
A	project/scripts/score_pacer200_drugclip2023_m4_loto.py

commit f35acdd9d330e88ed8c690f81986e0305f5a515f
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 10:27:55 2026 +0800

    fix router Python import example

M	project/docs/PACER_DRUGCLIP_ROUTER_INTERFACE_V01.md

commit cd53f989ec22ef6df406b4924ceadd51567b3ea0
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 10:27:27 2026 +0800

    add deployable PACER DrugCLIP routing interface

A	project/docs/PACER_DRUGCLIP_ROUTER_INTERFACE_V01.md
A	project/scripts/pacer_drugclip_router.py
A	project/tests/test_pacer_drugclip_router.py

commit 512f8fa866a8c187a9a4923e7de4550c1c264c56
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 10:19:51 2026 +0800

    benchmark DrugCLIP generations with GPCR-safe routing

A	project/docs/DRUGCLIP_GENERATION_BAKEOFF_PROTOCOL_V01.md
A	project/docs/DRUGCLIP_GENERATION_BAKEOFF_RESULT_V01.md
A	project/docs/DRUGCLIP_LOCKBOX_ACQUISITION_V01.md
A	project/docs/PACER_M4_PROJECT_MAINLINE_DEVELOPMENT_AND_INFERENCE_V01.md
A	project/docs/assets/drugclip_generation_bakeoff_v01/audit.json
A	project/docs/assets/drugclip_generation_bakeoff_v01/generation_dumbbell.png
A	project/docs/assets/drugclip_generation_bakeoff_v01/macro_metric_bars.png
A	project/docs/assets/drugclip_generation_bakeoff_v01/macro_metrics.csv
A	project/docs/assets/drugclip_generation_bakeoff_v01/matched_docking_comparison.csv
A	project/docs/assets/drugclip_generation_bakeoff_v01/paired_delta_forest.png
A	project/docs/assets/drugclip_generation_bakeoff_v01/per_target_metrics.csv
A	project/docs/assets/drugclip_generation_bakeoff_v01/per_target_roc_heatmap.png
A	project/docs/assets/drugclip_generation_bakeoff_v01/scaffold_bootstrap_ci.csv
A	project/scripts/run_drugclip_generation_bakeoff.py
A	project/tests/test_drugclip_generation_bakeoff.py

commit 270dd0eb5e88c0c369d493a58bd739cf118ebcdb
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 09:38:34 2026 +0800

    evidence: same-test-set M4R head-to-head 2023-triplet vs 2026-famaug

    Answering teammate challenge: the 0.69-vs-0.51 routing comparison was
    cross-protocol (illegitimate). New fair fight on the identical GaMD M4R
    screening pairs (25,249 rows, 2,301 actives), both sides never trained on
    these pairs (2023: ChEMBL triplet external_result; 2026: LOSO held-out
    predictions):

      config                  ROC    EF1%   BEDROC20  PR
      2023-triplet 10pock max 0.6673 0.39   0.1263   0.1295
      2023-triplet 10pock MEAN 0.7136 0.87  0.1802   0.1547  <- best
      2026-famaug LOSO 3-seed 0.5095 1.65   0.1118   0.0926
      2026 official raw       0.5401 1.56   0.1245   0.1070

    Verdict: 2023 base is genuinely the stronger M4 representation (ROC/PR/
    BEDROC all favor 2023-mean); 2026 edges only the noisy top-1% slice EF.
    Screening config frozen: 2023 external_result + MEAN pooling over 10 GaMD
    M4 pockets. Library ranking regenerated (top-1% overlap with famaug
    ranking only 8/285 -- genuinely independent second opinion).

A	project/results/drugclip_science2026/family_aug_v01/m4_library_oldtriplet_ranking.csv
A	project/results/drugclip_science2026/family_aug_v01/m4r_head2head_2023vs2026.json
A	project/scripts/m4r_head2head_2023_vs_2026.py
A	project/scripts/screen_m4_library_oldtriplet.py

commit d948297b1011707cfe54ad9dbe47df446506aac1
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 09:24:27 2026 +0800

    docs: repository organization -- mainline START_HERE entry + authoritative index

    Survey result: codex/project-integration-benchmark-v02 is the single
    superset (verified: every business branch has zero commits not in v02).
    Primary worktree fast-forwarded to v02 tip (77f44dc6).

    START_HERE_MAINLINE.md (root entry, tracked per repo convention; local
    README.md stays untracked per the upstream .gitignore) points to
    docs/REPOSITORY_MAINLINE_INDEX.md, which maps: the evidence-cascade
    mainline definition, per-line frozen conclusions (DrugCLIP fine-tuning /
    PACER-FKG Stage A+B / encoder C1-BS256 / OneProt retirement / CGDA
    closure / PACER-MCV / four-context MD), worktree statuses (geom2vec_pilot
    superseded), entry documents, next execution points, and no-retuning red
    lines. Marks the two older START_HERE_* files as archived sub-line entries.

A	START_HERE_MAINLINE.md
A	docs/REPOSITORY_MAINLINE_INDEX.md

commit 77f44dc6416f4d332bd2f0ff0bafbf4a9da9231c
Author: 2501_93004763 <2501_93004763@noreply.gitcode.com>
Date:   Fri Oct 2 09:22:56 2026 +0800

    figs: per-module horizontal scoreboards + cascade funnel

    fig_v02_07: four same-protocol module leaderboards side by side
    (binding 13T, static gating all-near-zero, dynamic three-class,
    QSAR collapse) with explicit 'different rulers' caveat in title;
    fig_v02_08: the only legitimate overall comparison - the candidate
    funnel 200/24/5/(1-3 pending, log scale) from handoff audit v02

A	project/results/project_wide_integration_benchmark_v02/figures/fig_v02_07_module_scoreboards.png
A	project/results/project_wide_integration_benchmark_v02/figures/fig_v02_08_cascade_funnel.png
A	project/scripts/make_module_scoreboards_v01.py

commit ddc11f970f508d89a75dac8ca30c41cfcd361731
Author: VergilVolk <barneygong329@gmail.com>
Date:   Fri Oct 2 02:05:11 2026 +0800

    freeze: finalize manifest-driven PACER production pipeline v01

    Integration base: 18100f9b6b42f27e43760daa5b3dfdd028b330d3
    PACER-FKG freeze anchor: b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd
    Encoder checkpoint SHA256: b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417

    No frozen scientific numerical asset is modified by this commit. The four frozen
    v02 modules (run_phase1_bs256.py, run_phase2_calibration.py, run_phase3_fkg.py,
    path_guard.py) are unchanged; V02_FREEZE_MANIFEST.json is unchanged.

    Added:
    - project/pacer_fkg_v02/frozen_runner/  manifest-driven, apply-only Phase 1/2A/2B drivers
                                            with fail-closed spec validation and an
                                            import-only frozen re-export layer
    - project/drugclip_freeze/             frozen DrugCLIP Model A (official checkpoint)
                                            and Model B (projection-only) apply paths;
                                            dual-model decision old_top50 INTERSECT new_top25
    - project/integration_freeze/          run-spec JSON Schema, CM00734 regression spec,
                                            receipts, and the code-freeze tooling
    - project/tests/                       test_frozen_runner.py, test_drugclip_freeze.py

    Evidence:
    - frozen anchor verifies strictly: VERIFIED_FROZEN, 12 pairs, 600 blocks,
      200 calibration blocks, R2-only fitting, 41 protected v01 files
    - block constructor bit-identical to both frozen implementations (1000x50 and 400x20)
    - CM00734 payload regression: 1632 exact fields, 0 mismatched; 108 cosines with
      max |delta| 6.4e-8 and 0 sign changes; primary signal -0.326839 reproduced
    - frozen encoder replay on CUDA: relative RMSE 1.77e-7 vs limit 1e-4
    - tests: 16/16 frozen_runner, 9/9 drugclip_freeze

    Mandatory fresh-checkout preflight: run
    project/integration_freeze/tools/normalize_frozen_eol.py --apply
    before using the frozen chain from a new checkout. The v02 freeze manifest pins a
    MIXED line-ending state (run_phase2_calibration.py as LF, analyze_pacer_fkg_g2b/g2c.py
    as CRLF), so a CRLF Windows checkout fails verify_frozen_anchor() until normalised.

    Claim boundary: prospective representation-space mechanistic evidence only. No PAM
    classifier, no potency predictor, no PAM probability, no magnitude threshold.

A	project/drugclip_freeze/RUNBOOK.md
A	project/drugclip_freeze/__init__.py
A	project/drugclip_freeze/apply.py
A	project/drugclip_freeze/dual_model.py
A	project/drugclip_freeze/model_a.py
A	project/drugclip_freeze/model_b.py
A	project/drugclip_freeze/projection.py
A	project/integration_freeze/FINAL_CODE_FREEZE_CANDIDATE_REPORT.md
A	project/integration_freeze/FROZEN_MODULE_LEDGER_v01.json
A	project/integration_freeze/RECEIPT_INDEX_v01.json
A	project/integration_freeze/WSL_HISTORICAL_WORKSPACE_RECEIPT_v01.md
A	project/integration_freeze/receipts/CM00734_PAYLOAD_REGRESSION_v01.json
A	project/integration_freeze/receipts/FROZEN_EOL_NORMALISATION_APPLIED_v01.json
A	project/integration_freeze/receipts/FROZEN_EOL_NORMALISATION_DRYRUN_v01.json
A	project/integration_freeze/receipts/FROZEN_EOL_PREFLIGHT_AT_FREEZE_v01.json
A	project/integration_freeze/receipts/PHASE1_ENCODER_REPLAY_v01.json
A	project/integration_freeze/receipts/WSL_DIRTY_CODE_DIFF_v01.patch
A	project/integration_freeze/receipts/WSL_DIRTY_FILES_BLOB_HASHES_v01.tsv
A	project/integration_freeze/receipts/WSL_GIT_BRANCHES_v01.txt
A	project/integration_freeze/receipts/WSL_GIT_BRANCH_v01.txt
A	project/integration_freeze/receipts/WSL_GIT_DIFF_NUMSTAT_v01.txt
A	project/integration_freeze/receipts/WSL_GIT_DIFF_STAT_IGNORE_CR_v01.txt
A	project/integration_freeze/receipts/WSL_GIT_HEAD_v01.txt
A	project/integration_freeze/receipts/WSL_GIT_LOG1_v01.txt
A	project/integration_freeze/receipts/WSL_GIT_REMOTES_v01.txt
A	project/integration_freeze/receipts/WSL_GIT_STASH_v01.txt
A	project/integration_freeze/receipts/WSL_GIT_STATUS_PORCELAIN_v01.txt
A	project/integration_freeze/receipts/WSL_RESULTS_FILES_v01.tsv
A	project/integration_freeze/receipts/WSL_RESULTS_SHA256_v01.txt
A	project/integration_freeze/receipts/WSL_RESULTS_TOTAL_BYTES_v01.txt
A	project/integration_freeze/receipts/WSL_TOPFILES_v01.tsv
A	project/integration_freeze/receipts/WSL_WORK_EXTRAS_SHA256_v01.txt
A	project/integration_freeze/run_specs/CM00734_REGRESSION_v01.json
A	project/integration_freeze/schemas/final_pipeline_run_spec_v1.schema.json
A	project/integration_freeze/tools/_shim_helper.py
A	project/integration_freeze/tools/diagnose_regression.py
A	project/integration_freeze/tools/normalize_frozen_eol.py
A	project/integration_freeze/tools/run_regression.py
A	project/integration_freeze/tools/verify_phase1_replay.py
A	project/pacer_fkg_v02/frozen_runner/__init__.py
A	project/pacer_fkg_v02/frozen_runner/blocks.py
A	project/pacer_fkg_v02/frozen_runner/cli.py
A	project/pacer_fkg_v02/frozen_runner/engines.py
A	project/pacer_fkg_v02/frozen_runner/fkg_phase1.py
A	project/pacer_fkg_v02/frozen_runner/fkg_phase2a.py
A	project/pacer_fkg_v02/frozen_runner/fkg_phase2b.py
A	project/pacer_fkg_v02/frozen_runner/io_utils.py
A	project/pacer_fkg_v02/frozen_runner/paths.py
A	project/pacer_fkg_v02/frozen_runner/regression.py
A	project/pacer_fkg_v02/frozen_runner/spec.py
A	project/tests/test_drugclip_freeze.py
A	project/tests/test_frozen_runner.py
```

本次修改任何科学结果：**NO**。唯一新增内容为本索引和事故记录；Git 元数据操作仅建立/提交/推送用户指定 incident 分支，不改变既有科学文件。
