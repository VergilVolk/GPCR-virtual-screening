# PACER Stage4 PBCFIX01 → 冻结 PACER-FKG v02 交接

适用生产轮次：`pacer_stage4_rerun_v01_pbcfix01`。分析接口版本：
`PACER_STAGE4_PBCFIX01_FKG_V02_v01`。从仓库根目录执行本文所有命令。
这是独立分析接口，没有修改或调用 MD 生产入口，也没有提交或推送 GitHub。

## 下载完成后一条命令

在本机既有 `pacer_dc_geom2vec` 环境及冻结资产路径保持不变时：

```powershell
conda activate pacer_dc_geom2vec
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --all --source-root D:/data/pacer_stage4_rerun_v01_pbcfix01
```

`--all` 顺序执行全量输入认证、Phase1 四帧 Encoder smoke、Phase1 全量、
Phase2a、Phase2b 和最终校验。任何失败都会停止，退出码非零。
重复执行时只有身份和哈希完全一致的结果才能复用。不会合并或重写 DCD。
输入目录必须是下载后的 **campaign 根目录**，不是 production 子目录或外层部署目录。
请对完成下载的稳定副本运行；仍在变动的 progress 或分段会导致认证失败。

换机器时，先指定冻结资产位置和 Encoder 依赖位置，再执行相同命令：

```powershell
$env:PACER_STAGE4_FROZEN_ROOT = 'D:/assets/GPCR-virtual-screening-fkg-v02'
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --all --source-root D:/data/pacer_stage4_rerun_v01_pbcfix01 --geom2vec-source D:/assets/geom2vec-source --checkpoint D:/assets/visnet_l6_h64_rbf64_r75.pth --device cuda
```

默认冻结目录为 `C:/projects/GPCR-virtual-screening-fkg-v02`，geom2vec 源码为
`C:/projects/geom2vec-source`，checkpoint 为仓库内
`project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth`。
这些默认路径不是下载 MD 时需要重新生成的资产。
`PACER_STAGE4_PBCFIX01_SOURCE_ROOT` 可设置各阶段默认数据根目录，显式 `--source-root` 优先。

本地验证环境：Python 3.11、MDAnalysis 2.10.0、NumPy 2.4.6、PyTorch 2.6.0+cu126，
以及既有 geom2vec、torch_geometric、冻结 Encoder 依赖。
纯 `--contract` 不需要 MDAnalysis/torch。预检和分段测试需要 MDAnalysis；
Encoder 推理需要模型及其依赖，`--device cpu` 可用于小规模推理检查。
不需要 OpenMM 来恢复 checkpoint，也不会启动 MD。

## 必须保留的下载输入

建议保留完整 campaign 原目录结构。校验所需的文件如下：

```text
<source-root>/
  SHA256SUMS_PBCFIX01.txt
  manifests/STAGE4_RERUN_JOB_MATRIX.csv
  manifests/STAGE4_RERUN_SYSTEM_MANIFEST.csv
  config/pacer_stage4_rerun_v01_pbcfix01.json
  scripts/...                     # 完整性清单中的全部静态文件
  source.json
  DEPLOYMENT_GUIDE.md              # 若在完整性清单中
  authorization.json
  evidence/source_snapshot.json
  evidence/...                    # 完整性清单中的文件
  logs/server_preflight_pacer_stage4_rerun_v01_pbcfix01.json
  systems/<system_id>/
    build_audit.json
    minimized.pdb
    system.xml
    state.xml
    ...                           # build_audit.outputs 中全部文件
  smoke/<system_id>/audit.json
  equilibration/<system_id>/replica_01/
    ready_pacer_stage4_rerun_v01_pbcfix01.json
    nvt/receipt.json, state.xml, state.csv, trajectory.dcd
    npt/...
    release_01/... through release_04/...
  production/<system_id>/replica_01/
    progress.json
    segments/00050000/trajectory.dcd, state.csv, checkpoint.chk, receipt.json
    segments/00100000/...
    ...
    segments/05000000/...
```

每个体系均需 replica_01、replica_02、replica_03。每条生产轨迹必须包含
100 个规范分段；总计 36 条轨迹、3600 个生产 checkpoint、360 ns。
checkpoint 仅按字节校验 SHA256，不反序列化、不恢复、不修改。
历史失败目录 `*.interrupted.*` 不作为输入；额外的规范数值分段目录会报错。
无需下载旧轮 MD 轨迹或旧轮绝对路径 `source.json.source_root` 指向的目录。
需要的是已冻结 PACER-FKG 数值资产及其有序映射。

生产身份按照实际 runner 的 fingerprint 重建，逐项关联包清单、config、build、
topology/system/state、smoke、authorization、preflight、source_snapshot、六个准备阶段回执及最终起始状态。
不会用旧轮 progress 的字段或“已知陈旧 progress 例外”放行本轮。
job matrix 的静态 status 可以仍是 `AWAITING_SERVER_BUILD`；完成条件取自生产回执和 progress。

## 36 任务身份契约

| 候选 | cluster | 独立情境 | 副本/种子 |
|---|---|---|---|
| PACERGEN02432 | 01 | A、P、C、CP | R1/27101、R2/38201、R3/49301 |
| PACERGEN00123 | 01 | A、P、C、CP | R1/27101、R2/38201、R3/49301 |
| PACERGEN00462 | 01 | A、P、C、CP | R1/27101、R2/38201、R3/49301 |

system_id 格式为
`pacer_stage4_rerun_v01_pbcfix01__<candidate>__cluster_01__<context>`；
context 分别为 apo、probe_only、candidate_no_probe、candidate_probe。
job_id 为 `<system_id>__R1` 等。A/P 按候选独立，不共享轨迹；相同 seed 不改变独立任务身份。
所有 R1/R2/R3 都只作冻结态评估，本轮 R2 不参与校准。

## 分阶段命令

```powershell
$data = 'D:/data/pacer_stage4_rerun_v01_pbcfix01'

# 可在下载生产轨迹之前运行：36任务、12体系、配置和包完整性清单
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --contract --source-root $data

# 下载期间检查冻结资产及36条轨迹文件可用性；不代表数据通过认证
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --preflight --source-root $data

# 少量真实分段：只需该任务前2段及其准备/身份资产，允许progress为RUNNING
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --sample --source-root $data --job-id pacer_stage4_rerun_v01_pbcfix01__PACERGEN02432__cluster_01__apo__R1 --sample-segments 2

# 下载完成后：认证36条轨迹全部3600段，发布正式输入回执
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_inputs_v01 --run --source-root $data

# Phase1：正式smoke认证全部输入，只嵌入PACERGEN02432/CP/R1的4帧
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_phase1_bs256_v01 --smoke --source-root $data
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_phase1_bs256_v01 --run --source-root $data
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_phase1_bs256_v01 --verify --source-root $data

# Phase2a：原始1000帧 → stride5 → 200帧 → 10个20帧块 → 冻结归一化/RFF
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_phase2a_frozen_apply_v01 --run --source-root $data
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_phase2a_frozen_apply_v01 --verify --source-root $data

# Phase2b：逐候选、逐副本匹配四情境 → 冻结图扩散和区域读出
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_phase2b_graph_region_v01 --run --source-root $data
python -m project.pacer_fkg_v02.run_stage4_pbcfix01_phase2b_graph_region_v01 --verify --source-root $data
```

单独运行 Phase1 时可同样传入 `--checkpoint`、`--geom2vec-source` 和 `--device`。
`--sample` 始终标记 `SAMPLE_ONLY_NOT_FULL_ACCEPTANCE`，不发布全量输入或 Phase1 完成回执。
Phase2a/2b 的 `--smoke` 依赖已完成的全量上游缓存，不替代少量 MD 分段预检。

## 时间、映射和冻结数值约定

每段 50000 步、2 fs/步、100 ps。state.csv 每 5000 步报告一次：
首段 5000–50000 步/10–100 ps，第二段 55000–100000 步/110–200 ps，
末段结束于 5000000 步/10000 ps。每段必须正好 10 帧、DCD dt=10 ps。
新版 reporter 不请求 Density 列；接口只要求实际生产代码写出的能量、温度和体积列。

OpenMM 每段创建新的 DCDReporter，文件头时间可能每段回到 10–100 ps。
接口接受经过验证的局部时钟或正确绝对时钟；其他时钟报错。
绝对时间由回执 start/end step 和 CSV 建立。
ChainReader 使用 `continuous=False` 保留全部帧，再只设置读取时间为 `(global_frame+1)*10 ps`。
不使用会按重叠头时间丢帧的连续时间筛选，也不拼接生成大 DCD。
坐标不 unwrap、居中、对齐或施加 PBC 变换，沿用冻结 Encoder 预处理约定。

1000 帧覆盖 10–10000 ps 的采样点，对应 10 ns 生产区间；不存在虚构的 t=0 帧。
`frames[::5]` 保留索引 0、5、…、995，物理时间 10–9960 ps、间隔50 ps。
每20帧块有19个内部 lag，跨度950 ps，名义区间1 ns；差分不跨块。
STATE_MOTION 和 SIGNED_DRIFT 计算式与既有 Stage4 完全相同。

冻结清单 SHA256：
`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`。
冻结 Encoder 代码/权重、normalization center/scale、bandwidth、RFF weights/bias/seed、
graph、regions 都按原回执认证并复用；没有训练、重新拟合或阈值选择入口。

270 个节点顺序来自冻结 `G2_RESIDUE_MAPPING_v01.csv`，先核对冻结图节点身份。
实际膜体系用拓扑顺序中的唯一序列子序列投影，不按可变 chain/resid 猜配。
缺失或多义映射、缺少任一 ATOM14 原子或重复原子名必须报错。
每个体系保存具体 source_resindex/resid/chain/segid 和选择表达式。
仅在独立输出目录保存去除 CONECT 的小 PDB 解析副本，原始拓扑不改动。

## 输出与恢复

独立结果目录：`project/results/pacer_stage4_pbcfix01_fkg_v02_v01/`。
独立缓存目录：`project/cache/pacer_stage4_pbcfix01_fkg_v02_v01/`。
输出根与输入根重叠时直接拒绝。写入保护禁止落入上一轮目录，禁止 junction 重定向。

主要交付物：

- `input_audit/INPUT_AUTHENTICATION_v01.json`：全量身份、分段时间、SHA256、实际残基映射。
- `input_audit/samples/*.json`：真实少量分段预检，明确非全量认证。
- `phase1_bs256/PHASE1_FREEZE_RECEIPT_v01.json`：36份 `[1000,270,256]` float32缓存和确定性重放。
- `phase2a_frozen_apply/PHASE2A_FREEZE_RECEIPT_v01.json`：72份双分支描述符与冻结RFF。
- `phase2b_graph_region/STAGE4_PROSPECTIVE_RESULTS_v01.json`：3候选、2分支、3副本描述性结果。
- `phase2b_graph_region/PHASE2B_FREEZE_RECEIPT_v01.json` 和18份区域块向量 npz。

BS256 全量缓存约9.27 GiB，Phase2a约0.65 GiB，另需单任务临时缓存和小型清单；
建议分析输出盘至少预留12 GiB。输入数据及checkpoint空间另计。
独立校验会重新读取源SHA256，大量文件的磁盘I/O是预期行为。

已发布文件不可覆盖；身份和内容相同可以跳过。
输入移动、代码变更或数据变动后，旧结果不会静默复用。
发现部分缓存/清单或残留 `.partial.*` 时先检查失败原因，再由操作者归档本轮分析目录并重新运行。
不删除或更改生产数据来“修复”预检失败。结果只表示冻结动态差异，不能当作 PAM/ago-PAM 标签、
药效概率或独立时间块统计推断。

## 测试与仍缺失的真实验证

```powershell
python -m unittest project.tests.test_stage4_pbcfix01_fkg_v02_v01 -v

# 下载前2段后启用真实分段集成测试（默认PACERGEN02432/A/R1）
$env:PACER_STAGE4_PBCFIX01_TEST_ROOT = $data
python -m unittest project.tests.test_stage4_pbcfix01_fkg_v02_v01.OptionalRealSegments -v
```

可用 `PACER_STAGE4_PBCFIX01_TEST_JOB` 选择其他矩阵 job_id。
测试使用少量合成 DCD/CSV/checkpoint 字节验证边界读取、完整性失败、时间断裂、错误帧数/dt、
错误身份、ATOM14不兼容等；合成数据不作为真实 MD 验证证据。
100分段/1000帧索引测试重复读取两个小测试文件，不宣称已验证100个真实生产分段。
实际36任务矩阵、包完整性、历史冻结资产及真实 oriented_cluster_01 受体基线可在本机验证。

截至本次实现，本地没有新版 committed segment、最终膜体系 minimized.pdb、production progress
或服务器准备链。真实两分段测试必须跳过，不能报成 PASS。
因此尚缺：12个最终膜体系的270残基/ATOM14校验、真实分段头时钟及CSV一致性、
3600份checkpoint/receipt完整性、全量36条轨迹/360ns连续性、真实Encoder确定性重放和
Phase1→Phase2a→Phase2b全量端到端执行。代码已经提供这些检查，需下载后执行上面的命令完成。
本次本地验证记录位于独立结果目录的 `local_validation/`。
