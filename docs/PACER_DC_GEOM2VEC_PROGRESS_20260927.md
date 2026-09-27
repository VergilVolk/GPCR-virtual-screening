# PACER-DC Geom2Vec 阶段进度与科学结果（2026-09-27）

**状态：阶段性暂停；TRAINING_GATE = CLOSED。**

## 1. 冻结范围与可追溯性

- 项目：人源 M4 受体、compound110 四上下文 MD 差分框架。
- 当前工作分支（用户最后确认）：`codex/pacer-dc-geom2vec-pilot`；HEAD `ccd5b33`（尚未包含本进度归档）。
- OneProt-MD G0–G4 已冻结，不因 Geom2Vec 结果调整历史分析。
- 冻结 Geom2Vec 官方源码提交：`371d642ec1061664f16e49fcac702d07fc8d0b51`。
- 冻结 ViSNet checkpoint：`project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth`；SHA256 `b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`。
- 本地部署验证：Python 3.11.16、PyTorch 2.6.0+cu126、CUDA 可用、RTX 4070 Laptop GPU；PyG CUDA radius graph PASS；Geom2Vec 1.0.0 官方模块导入 PASS；`pip check` 无依赖冲突。
- 区域映射：旧 G2 归档中的 `G2_REGION_MAP_v02.json`（提交 `fde6e9d`），包括新增的 12 残基 `distal_control`；其余预设区域 ACh_pocket 12、compound110_pocket 12、ECV_anchor 4、W435_gate 1、species_probe 2、global 270。

## 2. 已完成数据处理

- 第一阶段正式数据限定 compound110 的 R2/R3：两个独立 replica × 5 个连续窗口（W0–W4）× 4 个上下文（apo、probe_only、candidate_no_probe、candidate_probe），共 40 份 atom14 输入。
- 每份 NPY 原始形状 `(100,270,14,3)`；此前 atom14 输入 QC 40/40 PASS，但**atom14 原子顺序与 PBC 处理的独立复核仍未完成**。
- 冻结 ViSNet 采用 `stride=1`、`batch_size=1`、CUDA 提取，40 份全部完成。`batch_audit.json` 第一层审计 40/40 PASS，批处理总耗时 288.65 秒。
- 独立读取 40 份 `.geom2vec.npz`：共 4,000 帧，每份残基表示 `(100,270,128)`，有限性、帧号 `[0..99]`、序列长度及全局表示与残基表示均值的一致性全部 PASS。
- 本地数据目录（**非 Git 提交证明**）：`project/results/pacer_dc_geom2vec_R2R3_full_v01/`。本进度文档根据终端提供的执行和审计输出编写，尚未在本环境重新读取该 Windows 目录。
- R1 只有 W0 四上下文输入（4/20），曾用于方法开发和 CPU 小试，**不纳入当前 R2/R3 独立验证**。云端新增 MD 暂未纳入本次提取，其运行总规模及完成状态未核实。

## 3. 当前审计方法及结果

脚本：`project/pacer_dc_training/audit_geom2vec_cross_replica.py`。对于每个 replica × 窗口 × 上下文，先取 100 帧的 Geom2Vec 残基表示均值，再对预定义区域内残基平均，得 128 维区域向量 `μ`。定义：

- `dAGO = μ(candidate_no_probe) − μ(apo)`；此处名称是代码中的差分标签，不能仅凭命名认定已证实激动剂药理效应。
- `dINT = μ(candidate_probe) − μ(probe_only) − μ(candidate_no_probe) + μ(apo)`。

计算 R2/R3 对应窗口（W0 对 W0，…，W4 对 W4）区域差分向量余弦，报告五窗口中位数、10000 次窗口 bootstrap 的 95% 区间和预设门禁（中位数 > 0.5 且 bootstrap 下界 > 0）。这些窗口是相邻轨迹切片，**不是五个独立重复**；窗口区间不应解读为跨 replica 置信区间。

| 预定义区域 | dAGO 余弦中位数 | dINT 余弦中位数 | dINT 预设门禁 |
|---|---:|---:|---|
| ACh_pocket | 0.2962 | 0.2793 | 未通过 |
| compound110_pocket | -0.2624 | 0.1188 | 未通过 |
| ECV_anchor | 0.2475 | -0.2116 | 未通过 |
| W435_gate | 0.3730 | -0.0802 | 未通过 |
| species_probe | 0.6313 | -0.3049 | 未通过 |
| global | 0.7673 | -0.1315 | 未通过 |
| distal_control | -0.7761 | -0.7099 | 未通过 |

特别观察：

- **全局 dAGO** 五窗口余弦 `[0.607,0.917,0.767,0.721,0.799]`（据终端四舍五入显示），中位数 `0.7673`，窗口 bootstrap 95% 区间 `[0.607,0.917]`，按脚本门禁通过。
- **全局 dINT** 中位数 `-0.1315`，窗口 bootstrap 区间 `[-0.571,0.347]`，未通过。
- **功能区域 dINT**：5 个预设功能区域均未通过。`pilot_pass=false`，`functional_regions_passing_dINT=[]`，`distal_dINT_gate=false`，`full_replacement_authorized=false`。
- 远端对照 `dINT` 中位数 `-0.7099`；W2–W4 余弦 `[-0.710,-0.916,-0.861]`。负余弦可能反映方向反转，也可能受弱向量模长影响；尚未做差分幅度诊断，不能据此认定可靠的阴性对照表现。
- 输出原始审计文件（Windows 本地，归档时需实际添加）：`project/results/pacer_dc_geom2vec_R2R3_cross_replica_v01.json`。

## 4. 当前科学结论与严格边界

**已支持：** 冻结 Geom2Vec ViSNet 能在真实 R2/R3 MD 数据上完成可复核的 GPU 特征提取；线性时间均值/区域均值方法下，全局 dAGO 具有正向跨 replica 方向一致性。

**未支持：** 目前 5 个预定义功能区域的四上下文 dINT 均未达到 R2/R3 预设稳定性门禁；不能宣布 Geom2Vec 已解决 OneProt-MD 暴露的局部交互不稳定问题，也不能把本次数据视为 PAM 效力或效价验证。

**尚未测试：** 本次脚本采用的是 **Geom2Vec + 线性均值四上下文双差分**，并非 PACER-FKG 核分布均值嵌入或冻结多结构残基图传播。故本次 `pilot_pass=false` 不能外推为 PACER-FKG 已失败，也不能保证 FKG 会成功。

**替换门禁：** 即使两 replica pilot 后续通过，也仍需第二药理学对照 LY2119620；当前 `full_replacement_authorized=false`。`TRAINING_GATE=CLOSED`。

## 5. 暂停节点及续接任务

1. 将本进度文档与原始批处理审计 JSON、跨 replica 审计 JSON 一并归档；提交前核实文件、哈希和 `git status`。不重写 G0–G4、不覆盖既有结果。
2. 不默认将 40 份大型 NPZ 直接提交 Git；先生成独立文件清单与 SHA256，并另行确定大文件保存策略（Git LFS/外部归档），确保审计数据可追溯。
3. 继续时先做冻结线性审计的差分向量模长/近零向量诊断，然后按事先定义的参数运行 PACER-FKG 分布差分和冻结图分析；不得依据本次未通过结果临时挑选有利区域、参数或阈值。
4. 云端长程 MD 的 replica 总数及实际进度需要回到服务器核实，当前不纳入结论。
5. 后续独立核验 atom14 原子顺序及周期性边界处理；审视只有两个独立 replica 带来的统计限制。

## 6. 建议的 Git 归档文件

- `project/docs/PACER_DC_GEOM2VEC_PROGRESS_20260927.md`（本文件）
- `project/results/pacer_dc_geom2vec_R2R3_full_v01/batch_audit.json`（用户 Windows 工作区现有文件，须现场核对后添加）
- `project/results/pacer_dc_geom2vec_R2R3_cross_replica_v01.json`（用户 Windows 工作区现有文件，须现场核对后添加）

**注意：本文件的生成本身不代表以上 Windows 文件已被提交或 GitHub 已完成 push。**
