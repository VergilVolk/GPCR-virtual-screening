
# PACER-DC Geom2Vec 阶段进度与科学结果（2026-09-27；G2 完成后更新）

**状态：PACER-FKG G2 已完成并由项目负责人确认统一归档；阶段性暂停；`TRAINING_GATE = CLOSED`。**

本文件在原 Geom2Vec 提取与线性 pilot 进度报告基础上追加 PACER-FKG 方法学验证及与原 OneProt-MD PACER-DC 的比较。历史 pilot 结果保留，不因后续分析而更改其门禁或重新解释为正式 FKG 结果。

## 1. 冻结范围与可追溯性

- 项目：人源 M4 受体、compound110 四上下文 MD 差分框架。上下文分别为 CA（compound110 + ACh）、A（ACh）、C（compound110）、0（apo）。
- 当前 Geom2Vec / PACER-FKG 工作分支：`codex/pacer-dc-geom2vec-pilot`。原进度文档创建时的历史基准为 `ccd5b33`，首次进度归档提交为 `0447da4`；这些均非本次更新后的 HEAD 声明。
- OneProt-MD PACER-DC G0–G4 已冻结，不改写历史分析、区域定义或门禁。旧分析的正式交接文件：`project/pacer_dc_training/G0_G4_FREEZE_AND_H_HANDOFF_v01.md`。
- 冻结 Geom2Vec 官方源码提交：`371d642ec1061664f16e49fcac702d07fc8d0b51`。
- 冻结 ViSNet checkpoint：`project/tools/geom2vec_checkpoints/visnet_l6_h64_rbf64_r75.pth`；SHA256 `b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`。
- 部署时验证：Python 3.11.16、PyTorch 2.6.0+cu126、CUDA / RTX 4070 Laptop GPU；PyG CUDA radius graph、Geom2Vec 1.0.0 导入及 `pip check` 均通过。
- 旧线性 pilot 使用冻结的 `G2_REGION_MAP_v02.json`；后续 PACER-FKG 使用另外冻结的 M4 多结构残基图及九区域定义。两个阶段区域和表示空间不完全相同，不能直接比较其余弦数值或通过区域的数量。

## 2. Geom2Vec 提取与原线性 pilot（历史冻结结果）

- R2/R3 两个 replica × W0–W4 五个连续 1 ns 窗口 × 四上下文，共 40 份 atom14 输入；每份原始形状 `(100,270,14,3)`。每个窗口 100 帧，所有 40 份均已提取。
- 冻结 ViSNet 使用 `stride=1`、`batch_size=1` 在 GPU 上提取。原 `batch_audit.json` 40/40 PASS，批处理总耗时 288.65 秒；独立读取 40 份 `.geom2vec.npz` 得到 4,000 帧，各文件残基表示 `(100,270,128)`，有限值、帧号及全局均值一致性检查通过。
- 原始完整输入目录：`project/results/pacer_dc_geom2vec_R2R3_full_v01/`；相关批处理及跨 replica 审计已单独归档。R1 仅有 W0 的四上下文输入，未纳入 R2/R3 分析。长程云端 MD 不属于本批数据。
- 原线性 pilot 脚本：`project/pacer_dc_training/audit_geom2vec_cross_replica.py`。对帧与区域直接求 128 维均值后，计算 `dAGO = μ_C − μ_0`、`dINT = μ_CA − μ_A − μ_C + μ_0` 的匹配窗口 R2/R3 余弦。原预设门禁为五窗口中位数 > 0.5 且窗口 bootstrap 下界 > 0；连续窗口不是独立重复，区间仅作描述性用途。

| 原预定义区域 | 线性 pilot dAGO 余弦中位数 | 线性 pilot dINT 余弦中位数 | dINT 原门禁 |
|---|---:|---:|---|
| ACh_pocket | 0.2962 | 0.2793 | 未通过 |
| compound110_pocket | -0.2624 | 0.1188 | 未通过 |
| ECV_anchor | 0.2475 | -0.2116 | 未通过 |
| W435_gate | 0.3730 | -0.0802 | 未通过 |
| species_probe | 0.6313 | -0.3049 | 未通过 |
| global | 0.7673 | -0.1315 | 未通过 |
| distal_control | -0.7761 | -0.7099 | 未通过 |

原 pilot 的全局 dAGO 门禁通过；五个原预定义功能区域的 dINT 均未通过，`pilot_pass=false`、`full_replacement_authorized=false`。这只是**线性 Geom2Vec 均值分析**的结果，不代表随后完成的 PACER-FKG 核分布分析。

原输出：`project/results/pacer_dc_geom2vec_R2R3_cross_replica_v01.json`。

## 3. PACER-FKG 主分析与 G1 诊断

PACER-FKG 在已冻结的 40 份 Geom2Vec 残基表示上研究四上下文分布差分，定义：

- `ΔPAM = μ_CA − μ_A`
- `ΔAGO = μ_C − μ_0`
- `ΔINT = ΔPAM − ΔAGO`

这里的 `μ` 在核分析中指核均值嵌入；`PAM`、`AGO` 与 `INT` 是操作性对比名称，不是已经确定的药理学功效或构象机制。

- 冻结 M4 多结构图：270 个节点、1,304 条边、九个冻结区域。主分析包括非负的局部 RBF 核 `U²` 响应幅度以及图传播描述量；其非负响应幅度不能被解释为差分方向。
- 主分析共完成 10 个 R2/R3 窗口 JSON 及汇总审计。原区域特异性门禁使用相对稳定核心的 `U²`、两 replica 中位响应方向、至少 4/5 正向窗口和非负匹配窗口 Spearman，同时检查远端与稳定核心对照。
- 在原先的**幅度与窗口排序门禁**下：`ΔINT` 的 `compound110_extension`、`pam_contact_union` 通过；`ΔAGO` 的 `pam_contact_union` 通过；`ΔPAM` 的 `compound110_extension`、`cooperativity_mutagenesis`、`intracellular_microswitches`、`pam_contact_union` 通过。
- G1 的 270 个窗口记录、54 个分 replica 区域汇总、27 个跨 replica 区域—差分轴比较揭示：远端对照存在较强幅度响应及不稳定的跨 replica 排序；不同区域的带宽和非负图扩散会影响结果解释。因此原门禁保持冻结，新增 G2 检验方向、数值稳健性和图结构贡献。

主分析：`project/results/pacer_dc_fkg_R2R3_v01/PACER_FKG_R2R3_AUDIT.json`。

G1：`project/results/pacer_dc_fkg_G1_diagnostics_v01/`。

## 4. G2-A：atom14 映射、几何与 PBC 独立审计

- 40/40 个 atom14 输入通过几何和原子映射筛查；使用上游残基常量核对了残基—原子定义。此结果纠正了原进度文档中“atom14 原子顺序尚未独立复核”的历史待办。
- PBC 检查覆盖八条原始 DCD 轨迹、4,000/4,000 帧有效周期盒，与 40/40 个冻结原始窗口匹配；定义内未检出 Cα 周期盒拆分或时间跳跃。
- 原 R3 apo 来源为恢复后的 `pacer_dc_recovery_merged_v01/apo/replica_03/trajectory_corrected.dcd`，不应误用普通 production DCD。
- 正式状态：`PASS_PBC_SCREEN_NOT_FULL_PROOF`。这些检查支持当前 atom14 输入的可用性，**不是对每个原子的 PBC 无误性或 MD 收敛性的完整证明**。最初 CSV 报告字段异常已基于完成的审计 JSON 修复，无须重跑八条轨迹。

正式审计：`project/results/pacer_dc_fkg_G2A_PBC_v01/G2A_PBC_AUDIT.json`；对应脚本位于 `project/pacer_dc_training/`。

## 5. G2-B：共同核空间差分方向与补充稳健性

G2-B 冻结 R2 拟合的通道 median/MAD 和区域核带宽，使用共同随机傅里叶特征（RFF，512 维）来比较 R2/R3 的核均值**差分方向**。40/40 输入哈希通过，九个区域 × 三种差分得到 27 行，135 个匹配窗口结果，400 次时间块重采样仅作描述性敏感性分析。

| 区域 | ΔINT 方向余弦 | ΔAGO 方向余弦 | ΔPAM 方向余弦 | ΔPAM 描述性区间 |
|---|---:|---:|---:|---|
| pam_contact_consensus | -0.135 | -0.168 | **0.514** | 0.397～0.547 |
| cooperativity_mutagenesis | 0.011 | 0.357 | **0.442** | 0.314～0.528 |
| orthosteric_contact_union | 0.423 | 0.184 | 0.433 | 0.112～0.543 |
| pam_contact_union | **-0.228** | -0.241 | **0.414** | 0.278～0.456 |
| orthosteric_activation_core | 0.467 | 0.179 | 0.334 | 0.051～0.468 |
| compound110_extension | **-0.171** | -0.021 | 0.122 | -0.101～0.288 |
| distal_control | **-0.604** | **-0.810** | 0.134 | -0.069～0.326 |

原 ΔINT 幅度门禁通过的 `compound110_extension` 与 `pam_contact_union`，在新的共同核**方向**分析中余弦分别为 -0.171、-0.228：幅度或窗口排序重复出现并不自动意味着分布变化方向一致。

相反，PAM 接触及协同突变区域的 ΔPAM 呈现正向方向一致性。不过各区域可能重叠，不能视为独立证据。

补充稳健性研究比较 RFF 维度 256/512/1024、三个随机种子及同一确定性平衡子样本的精确 RBF 核：

| 区域／轴 | 原完整采样 512D RFF 余弦 | 子样本精确核余弦 |
|---|---:|---:|
| pam_contact_consensus / ΔPAM | 0.514 | 0.368 |
| cooperativity_mutagenesis / ΔPAM | 0.442 | 0.357 |
| pam_contact_union / ΔPAM | 0.414 | 0.311 |
| compound110_extension / ΔINT | -0.171 | -0.192 |
| pam_contact_union / ΔINT | -0.228 | -0.102 |
| distal_control / ΔINT | -0.604 | -0.322 |

上述组合在两类计算中符号一致，但**全量数据与子样本结果不是同一估计量**，差值不可全归因于 RFF 近似。

所报告的 1024 维最大绝对误差按区域—轴为 0.014～0.080；不能据此宣布所有近零方向余弦均已收敛。没有事后挑选随机种子，也没有修改原门禁。

主报告：`project/results/pacer_dc_fkg_G2B_v01/G2B_DIRECTION_REPORT.md`。

补充报告：`project/results/pacer_dc_fkg_G2B_robustness_v01/G2B_ROBUSTNESS_REPORT.md`。

## 6. G2-C：带符号共同核图传播消融

G2-C 使用同一套 R2 拟合的通道预处理，另行校准全残基共享的核空间，使用 256 维 RFF 进行有符号差分节点传播。

它**不是** G2-B 的分区域核估计量，数值大小不得跨阶段直接比较。

实验比较无图、冻结真实图和 12 个接触边重连的置乱图（保持无权度数及骨架连接，不保持加权节点强度）。40/40 输入通过；27 个区域—轴结果、135 个匹配窗口比较和 324 条置乱图记录已经生成。

| 区域／轴 | 无图余弦 | 真实图余弦 | 置乱图中位数 | 真实图减无图 | 描述性变化区间 |
|---|---:|---:|---:|---:|---|
| pam_contact_consensus / ΔPAM | 0.262 | -0.056 | 0.223 | **-0.318** | -0.389～-0.191 |
| pam_contact_union / ΔPAM | 0.315 | -0.028 | 0.281 | **-0.343** | -0.417～-0.129 |
| compound110_extension / ΔINT | -0.109 | 0.035 | -0.182 | 0.144 | 0.003～0.209 |
| intracellular_microswitches / ΔAGO | 0.053 | 0.388 | 0.141 | 0.335 | 0.167～0.417 |
| distal_control / ΔINT | -0.690 | -0.691 | -0.669 | -0.001 | -0.031～0.046 |

真实图效应依赖区域及差分轴。在两个 PAM 接触区域的 ΔPAM 分析中，图传播明显降低正向方向一致性，而置乱图中位数保持正向；因而目前**没有证据支持将冻结图作为普遍有效的方向一致性增强组件**。

少数其他区域存在描述性改善，但不能据此断言结构图具有因果贡献。12 次置乱和时间块重采样均不构成显著性推断。

正式报告：`project/results/pacer_dc_fkg_G2C_v01/G2C_GRAPH_ABLATION_REPORT.md`。

## 7. PACER-FKG 相比原 PACER-DC 的新进展（独立科学比较）

本节比较**新增的科学证据与可识别的问题**，而非宣称跨模型数值指标可直接对齐。

原 PACER-DC 使用 OneProt-MD 的 L0/L1 等表示；Geom2Vec 线性 pilot 以及 PACER-FKG 的 RBF 均值嵌入又属于不同的统计对象。

共享的是 compound110 四上下文、短程 R2/R3 与操作性差分定义，不共享全部表示空间及区域定义。

### 7.1 原 PACER-DC 已经建立的基线

原 PACER-DC G4 的五窗口平均 ΔINT R2/R3 余弦：

- ACh pocket：L0 为 0.305，L1 为 **0.548**；
- compound110 pocket：L1 为 -0.070；
- 远端对照：L1 为 -0.042；
- 全局：L1 为 0.011。

ACh pocket 的 L1 ΔINT 在五个匹配窗口全部为正，说明原方法并非完全没有局部重复信号。

不过其表现受 FinalLayer 区域依赖性、replica 偏差抵消及有限短程采样影响；独立 Cα 几何特征未能稳健确认对应构象机制。

上述值是旧 L1 表示中的余弦，**不得与 FKG 核余弦直接比较大小**。

### 7.2 FKG 的实质性新增能力与证据

1. **超越线性均值：** 从区域均值向量的线性差分，扩展为 RBF 核分布响应与核均值嵌入，可表征均值变化之外的分布差异。原 Geom2Vec 线性 pilot 五个预定义功能区域 ΔINT 全未通过旧门禁，并不等于核分析没有信息；但 FKG 尚未通过独立长程或药理学对照证明整体优越性。

2. **把响应幅度和响应方向分开审计：** 原 FKG 的 ΔINT 幅度门禁在 compound110 extension 和 PAM contact union 通过；G2-B 却发现其跨 replica 方向余弦分别为 -0.171、-0.228。这一诊断能力防止将较强的非负核响应误报为同向协同交互机制。

3. **得到可独立验证的条件性差分线索：** G2-B 在 PAM contact consensus（0.514）、PAM contact union（0.414）和 cooperativity mutagenesis（0.442）观察到正向 ΔPAM 方向一致性；同一确定性子样本的精确核计算保留这三个方向的正号。这些是两 replica 的描述性结果，不构成化合物 PAM 效力或跨分子泛化证据。

4. **提供输入质量和数值稳健性的可追溯证据：** G2-A 完成 atom14/几何/PBC 独立筛查；G2-B 补充了固定样本的精确 RBF 核对照和不同 RFF 维度/种子复核。数据质量筛查通过不等于构象采样收敛，RFF 近似敏感性也应继续显式报告。

5. **首次检验冻结图是否真正提供附加信息：** G2-C 的真实图与无图及置乱图消融显示图效应不统一，甚至削弱核心 PAM 接触区域的 ΔPAM 方向一致性。这是否定性但重要的方法学结果：目前不应将真实图的传播作用包装为已经验证的模型增益。

### 7.3 尚未克服的限制

- 原 PACER-DC 的 ACh pocket ΔINT 与 PACER-FKG 的 PAM 区域 ΔPAM **不是一对一的正面性能比较**。尚未在相同区域、相同估计目标和预先冻结的评估准则下完成公平 benchmark。
- G2-B 与 G2-C 使用不同核空间及不同统计对象，数值不得直接比较；区域重叠、单一化合物及仅两条 replica 均限制可推广性。
- 五个连续短程窗口与每窗口时间块存在序列相关性，描述性区间不是总体置信区间，也没有有效 p 值、药理学 PAM 效应或因果图机制结论。
- 远端对照可出现较强非负核响应，且 G2-B 的远端 ΔINT、ΔAGO 余弦分别为 -0.604 和 -0.810；区域特异性及真实结构对应关系仍需新增独立采样核验。

**阶段性判断：** PACER-FKG 相比原 PACER-DC，已经在非线性分布分析、核方向复现诊断和图传播对照上扩展了方法学能力，并发现部分值得独立检验的 ΔPAM 区域信号。现有证据尚不足以证明 PACER-FKG 在稳定性、PAM 特异性、药理学效能或模型总体性能上优于 PACER-DC。

## 8. 归档状态、冻结边界与续接事项

- 截至本次文档更新，项目负责人已确认 G2-A、G2-B、G2-B 稳健性复核、G2-C 的本地统一归档工作已完成。**本文件的此次更新仅为待提交的文档改动，不代表它已自动写回 GitHub。** 远端 HEAD 和本次文档更新的提交号需在用户本地 Git 检查后记录。
- 原 Geom2Vec 数据与旧 PACER-DC G0–G4 保持冻结。原 FKG 标量门禁与 G2-B、G2-C 报告不因本阶段科学解释而重算或修改。
- `TRAINING_GATE = CLOSED`；没有授权分类器训练、药理学效力推断或根据目前结果重新挑选图结构与有利超参数。
- 长程独立 MD 归入后续独立验证阶段；开始新数据的机制性分析前，应冻结检验指标、区域、采样质量与自相关检查方案。新增 MD 的完成状态、收敛性及结果不能从此次短程审计推定。
- 核心证据路径：
  - `project/pacer_dc_training/G0_G4_FREEZE_AND_H_HANDOFF_v01.md`
  - `project/results/pacer_dc_geom2vec_R2R3_cross_replica_v01.json`
  - `project/results/pacer_dc_fkg_R2R3_v01/PACER_FKG_R2R3_AUDIT.json`
  - `project/results/pacer_dc_fkg_G1_diagnostics_v01/`
  - `project/results/pacer_dc_fkg_G2A_PBC_v01/G2A_PBC_AUDIT.json`
  - `project/results/pacer_dc_fkg_G2B_v01/G2B_DIRECTION_REPORT.md`
  - `project/results/pacer_dc_fkg_G2B_robustness_v01/G2B_ROBUSTNESS_REPORT.md`
  - `project/results/pacer_dc_fkg_G2C_v01/G2C_GRAPH_ABLATION_REPORT.md`

**记录边界：** 本文仅记录已完成和经核实的短程审计，不把尚未验证的长程采样、未来方法设想或后续方案调整写成当前正式研究方案。
