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
