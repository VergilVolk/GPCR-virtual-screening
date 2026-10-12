# Stage4 PBCFIX01 PACER-FKG v02 独立科学分析 v01

## 核心结论与结论等级

**直接数值事实：** 三候选共同具有 STATE_MOTION 分支 PAM 接触区域的 ΔAGO、ΔPAM 重复方向一致性；该结论没有延伸到全受体、所有胞内微开关或两个分支。PAM 接触区域的 ΔINT 均出现至少一对负余弦。SIGNED_DRIFT 的 ΔPAM、ΔINT 没有呈现与上述 STATE_MOTION 接触信号相当的重复方向证据。幅度稳定和方向稳定必须分别判断。

**探索性机制解释：** 冻结特征中存在候选相关、在无 probe 和有 probe 背景下均出现的 PAM 区域动态差异；它可能反映局部构象/波动分布的改变。STATE_MOTION 合并平均状态与运动强度，现有区域 RFF 向量不能分离两者贡献。没有形成三个重复都支持的、从接触区到胞内输出的完整一致链条。

**待实验功能假说：** 三者可作为检验“候选自身效应、probe 背景调制、二者交互”的实验对象，但冻结 ΔAGO/ΔPAM 名称不等于 agonist/PAM 活性标签。本文不推断 PAM/ago-PAM 活性、效力、概率或因果机制，不给出综合药效排名。

## 1. 输入身份、冻结协议及审计

最终 JSON：`../pacer_stage4_pbcfix01_fkg_v02_v01/phase2b_graph_region/STAGE4_PROSPECTIVE_RESULTS_v01.json`。其 SHA256 为 `ee3a6621b6df039c18cd762a5d2448cd60b1c69642e7131203df281a11558236`，与 Phase2b 最终冻结回执一致。全部18份 NPZ 的字节数和 SHA256 与 JSON 清单一致；候选×分支×R1/R2/R3 为完整 3×2×3 网格，路径、seed、cluster=1、A/P/C/CP 身份及上游 manifest 哈希均核对。未剔除、填补、调换任何重复。

独立复算了 1458 个有限的 (10,512) 数组、810 个 JSON 重复/区域摘要（包括两个辅助 probe 差分）及主要差分的 486 个重复对余弦。JSON 摘要复核容差 rtol=1e-10、atol=1e-12。上下文线性重构与冻结 contrast 的最大绝对误差 9.55e-09；因分别扩散和 float32 存储产生舍入差异，按 rtol=2e-5、atol=2e-7 核对，不替换冻结数组。

全量输入认证 `INPUT_AUTHENTICATION_v01.json` 状态为 `36_OF_36_ACCEPTED`，36条轨迹各100段，总计3600段。本文读取并记录该回执，不重新验证原始 DCD/checkpoint 内容。早期 `LOCAL_VALIDATION_v01.json` 和 `DEPLOYMENT_ACCEPTANCE_v01.json` 以及交接文档中的“未完成”属于先前部署/采样快照，不能代替最终输入认证与 Phase2b 回执，也不能将早期快照称为最终验证。本文验证范围是冻结结果身份与数值复算，不声称重新执行了原始 MD 到 Encoder 的最终流水线验证。

冻结历史 manifest SHA256：`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`；图 SHA256：`25593c11dac0f70ad6811c8c6aa7a4f06c1755f14397f1162eaa1157026be122`。历史校准来自历史长 MD 的 R2；本次 Stage4 R2 与 R1、R3 一样仅用于评价。历史定义和数值资产只做读取/哈希核对。

360 ns = 3候选×4上下文×3重复×10 ns；不是单条360 ns。每条1000原始帧，以 frames[::5] 取200帧；采样点10–9960 ps，间隔50 ps。每20帧形成1个名义1 ns块，块内首末跨度950 ps，每条10块。没有丢弃前段、重新划块或额外平衡期选择。10块彼此相关，仅作时间分布描述，独立重复数为3；全部区域亦非独立观测。

## 2. 差分与统计量的解释

| 名称 | 冻结公式 | 可作的描述 |
| --- | --- | --- |
| ΔAGO | C − A | 无 probe 背景的候选相关动态差异 |
| ΔPAM | CP − P | 有 probe 背景的候选相关动态差异 |
| ΔINT | CP − C − P + A = ΔPAM − ΔAGO | 两个背景的候选差异之差 |

A=apo、P=probe_only、C=candidate_no_probe、CP=candidate_probe，全部取同一候选、同一重复匹配上下文。配对种子和时间块编号不代表两条随机轨迹的构象同步；ΔINT 的块幅度还包含上下文采样波动。ΔINT 非零不直接证明生物协同或非加和药效：RBF-RFF 本身非线性，其差分发生在冻结表示空间。

STATE_MOTION = concat(mean(Z), sqrt(mean(diff(Z)^2)))；SIGNED_DRIFT = mean(diff(Z))，后者是块内50 ps采样步的平均有符号增量，经独立冻结归一化与 RFF 后的对比，不等于受体整体运动速度。两分支使用各自冻结尺度和随机特征，原始幅度不能跨分支比较，RFF 方向亦没有通用“激活方向”。正余弦表示同一差分在不同重复的特征方向对齐；负余弦表示方向冲突，不等于抑制/反向药理活性。

令 v[r,b] 为一个重复的块向量，m[r]=mean_b v[r,b]。报告 ||m[r]||（时间平均向量幅度）、||v[r,b]|| 的 min/q10/median/q90/max，以及 R1–R2、R1–R3、R2–R3 的 cos(m[i],m[j])。幅度 CV 用三个 ||m[r]|| 的样本标准差(ddof=1)/均值。时间相干度 = ||mean(v)||/mean(||v||)，用于观察块间抵消；early/late 为前五块和后五块平均向量余弦，不能充当收敛检验。

不做 p 值、置信区间、独立块 bootstrap、显著性宣称、数据驱动阈值、模型训练或综合权重。统计三个余弦是否全正仅作为几何符号事实（零点不是药效阈值）；接近零的正值不据此宣称强一致性。图中全部区域使用固定顺序。幅度图每个分支独立设置统一面板尺度。

## 3. 区域依赖和特异性限制

保留9个冻结区域：PAM consensus/union、cooperativity_mutagenesis、compound110_extension；正交 activation_core/contact_union；intracellular_microswitches；stable_core_control/distal_control。区域大小、每个残基的结构编号与 embedding_index 见 `REGION_MEMBERSHIP_v01.csv`。区域名称是结构先验，不能因名称将其动态变化认作该功能已发生。

| 区域对 | 共同残基数 | Jaccard | 图扩散共享质量 |
| --- | --- | --- | --- |
| pam_contact_consensus / pam_contact_union | 12 | 0.750 | 0.871 |
| pam_contact_consensus / cooperativity_mutagenesis | 9 | 0.562 | 0.802 |
| pam_contact_consensus / orthosteric_activation_core | 1 | 0.042 | 0.240 |
| pam_contact_consensus / stable_core_control | 0 | 0.000 | 0.021 |
| pam_contact_consensus / distal_control | 0 | 0.000 | 0.005 |
| orthosteric_activation_core / orthosteric_contact_union | 13 | 0.929 | 0.953 |
| orthosteric_activation_core / intracellular_microswitches | 1 | 0.048 | 0.282 |

图采用冻结 alpha=0.65、20步重启扩散，每步 H[t+1]=0.35I+0.65 T H[t]，H[0]=I。仅对结构算子计算区域平均影响权重 w，未对 MD/Encoder/RFF 数值重新运行扩散；这是评估共享信息的解析性诊断。共享质量=sum_j min(w_region1[j],w_region2[j])，范围0–1；它不是实测信号的独立性校正，也没有用于扣减、重加权或重新排名。

PAM consensus 是 union 的子集；正交 activation core 也基本嵌在 contact union 内，cooperativity 与 PAM 区域重叠。因此多个此类区域一起重复，不是多份独立支持。即使直接无交集，图扩散仍共享信息；对照区域同样经历扩散。区域取均值降低残基数的简单求和效应，但区域大小、结构平滑程度与局部变异仍不同，区域范数比不能自动解释为选择性。

特异性同时审阅两个对照：每个重复的 region/stable-control 与 region/distal-control 范数比均输出到区域 CSV，不跨分支、不先合并重复。远端对照并非静止：某些 ΔINT 的远端幅度超过 PAM 接触区域，限制“局部特异交互”的解释。`CROSS_REGION_DIRECTION_v01.csv` 另外给出同重复全部区域对的方向余弦，防止将平滑后的共同方向当作信号传播路径。

## 4. 全部主要差分的候选横向比较

| 候选 | 分支 | 差分 | 9区域中三对余弦均正的数量 | 最小余弦 | 三对最小余弦最高的区域 |
| --- | --- | --- | --- | --- | --- |
| PACERGEN02432 | STATE_MOTION | Delta_AGO | 4 | -0.314 | pam_contact_consensus |
| PACERGEN02432 | STATE_MOTION | Delta_PAM | 6 | -0.359 | pam_contact_consensus |
| PACERGEN02432 | STATE_MOTION | Delta_INT | 0 | -0.635 | pam_contact_union |
| PACERGEN02432 | SIGNED_DRIFT | Delta_AGO | 1 | -0.235 | pam_contact_consensus |
| PACERGEN02432 | SIGNED_DRIFT | Delta_PAM | 0 | -0.438 | distal_control |
| PACERGEN02432 | SIGNED_DRIFT | Delta_INT | 0 | -0.329 | distal_control |
| PACERGEN00123 | STATE_MOTION | Delta_AGO | 3 | -0.519 | cooperativity_mutagenesis |
| PACERGEN00123 | STATE_MOTION | Delta_PAM | 5 | -0.290 | pam_contact_consensus |
| PACERGEN00123 | STATE_MOTION | Delta_INT | 1 | -0.689 | intracellular_microswitches |
| PACERGEN00123 | SIGNED_DRIFT | Delta_AGO | 1 | -0.309 | intracellular_microswitches |
| PACERGEN00123 | SIGNED_DRIFT | Delta_PAM | 1 | -0.181 | compound110_extension |
| PACERGEN00123 | SIGNED_DRIFT | Delta_INT | 0 | -0.359 | intracellular_microswitches |
| PACERGEN00462 | STATE_MOTION | Delta_AGO | 6 | -0.287 | cooperativity_mutagenesis |
| PACERGEN00462 | STATE_MOTION | Delta_PAM | 3 | -0.258 | pam_contact_consensus |
| PACERGEN00462 | STATE_MOTION | Delta_INT | 3 | -0.606 | stable_core_control |
| PACERGEN00462 | SIGNED_DRIFT | Delta_AGO | 5 | -0.401 | cooperativity_mutagenesis |
| PACERGEN00462 | SIGNED_DRIFT | Delta_PAM | 0 | -0.126 | orthosteric_activation_core |
| PACERGEN00462 | SIGNED_DRIFT | Delta_INT | 1 | -0.251 | orthosteric_activation_core |

该表展示空间分布，没有将重叠区域计数转化为证据强度分数。候选比较仅限相同分支、相同差分和区域；不按幅度、计数或任一余弦组合生成药效排序。各候选使用各自冻结 A/P 背景，不混用对照。

### PACERGEN02432

**支持（直接事实）：** STATE_MOTION 的 PAM consensus ΔAGO 余弦为0.562/0.424/0.744，ΔPAM 为0.742/0.729/0.583（顺序均R1–R2/R1–R3/R2–R3）。ΔPAM 在 consensus、union、cooperativity 的三对均正，orthosteric_activation_core 亦为0.225/0.349/0.274，但后者对齐较弱。ΔAGO 的 PAM consensus 幅度 CV=0.071左右，说明该局部差异在三个重复中幅度也较稳定；精确值见附表。

**冲突（直接事实）：** PAM consensus ΔINT 为0.287/−0.182/−0.003，extension ΔINT 的R1–R3为−0.635。胞内 ΔPAM 为0.322/−0.158/−0.072，ΔINT 为0.426/−0.445/−0.052；局部差异未伴随可重复胞内方向。SIGNED_DRIFT ΔPAM 的 PAM consensus 为−0.191/0.121/−0.329。SIGNED_DRIFT 胞内 ΔINT 的三个幅度几乎相等（0.06578/0.06606/0.06620），方向却为0.306/−0.106/−0.043，是“幅度一致不等于方向一致”的具体例子。

**探索性解释与边界：** 更符合可重复的 PAM 接触区状态/波动表示变化，并有较弱的正交区域背景依赖差异；不能证明协同交互、持续定向漂移或胞内激活。

### PACERGEN00123

**支持（直接事实）：** STATE_MOTION 的 PAM consensus ΔAGO=0.704/0.608/0.430，ΔPAM=0.624/0.315/0.692。cooperativity 的 ΔAGO=0.749/0.761/0.715，是本候选最一致的主要区域信号。PAM consensus ΔPAM 幅度 CV约0.094。胞内 ΔPAM 三对为0.122/0.189/0.413，虽全正但对齐较弱；胞内 ΔINT 的R1–R3仅0.003，不能用全正计数强化它。

**冲突（直接事实）：** 正交 activation core 的 ΔAGO 三对均负（−0.443/−0.032/−0.046），ΔPAM 亦有R2–R3=−0.144。PAM consensus ΔINT=0.377/−0.407/−0.292，union R1–R3=−0.559。extension ΔINT=0.668/−0.689/−0.525，若只看R1–R2会遗漏严重冲突。distal ΔINT R2幅度0.1802，高于R1=0.06665及R3=0.09934，CV约0.507，方向又为−0.332/0.249/−0.637。SIGNED_DRIFT 的 PAM consensus ΔPAM=0.083/−0.085/0.035，缺乏重复对齐。

**探索性解释与边界：** 接触/协同性结构先验区域的无 probe 差异较稳定；正交激活区和交互差分与之不一致。胞内弱正对齐不足以支持完整激活机制，更不能跳过远端R2差异。

### PACERGEN00462

**支持（直接事实）：** STATE_MOTION PAM consensus ΔAGO=0.820/0.640/0.814；cooperativity ΔAGO=0.868/0.749/0.819，其三对最小余弦高于其他两候选同区域同差分。ΔPAM PAM consensus=0.746/0.572/0.739，cooperativity=0.732/0.571/0.784。正交 activation core 的 ΔAGO=0.520/0.322/0.308，ΔINT=0.260/0.444/0.164，均全正但明显弱于局部接触区。SIGNED_DRIFT 正交 ΔAGO=0.134/0.177/0.237，提供弱方向一致事实，不能与 STATE_MOTION 幅度相加。

**冲突（直接事实）：** PAM consensus ΔINT=0.052/−0.606/0.077、cooperativity R1–R3=−0.430。胞内 ΔAGO=−0.287/−0.141/0.197、ΔPAM=0.206/−0.258/−0.197、ΔINT=−0.278/−0.153/0.070，均混合。ΔPAM stable-control 的R3幅度0.06658约为R1=0.03290的2倍，CV约0.423；该重复的背景变化削弱选择性。ΔINT stable-control 三对为0.281/0.183/0.208，接近甚至超过部分正交交互对齐，进一步限制独特激活解释。SIGNED_DRIFT PAM consensus ΔPAM=0.080/0.205/−0.100。

**探索性解释与边界：** 在同分支同区域比较中，它具有较整齐的接触/协同性区域 ΔAGO、ΔPAM 方向，并出现较弱的正交 ΔAGO/ΔINT 一致性；但 PAM 区域交互方向严重不稳定、胞内方向冲突、对照也变化，不能称其更可能有功能活性或已建立变构通路。

## 5. 幅度、时间分布和交互几何复核

附表列出全部三个重复的时间平均范数和三对余弦；CSV 进一步保留每个重复10块的分位数、极值、时间相干度和 early/late 余弦。范数是非负数，噪声或抵消也可产生非零范数。块范数 q10–q90 是观测分布，不是均值置信区间；三个重复的CV低也不能补偿向量方向相反。

ΔINT=ΔPAM−ΔAGO 已逐块核对。`CONTRAST_RELATIONS_v01.csv` 输出同重复 AGO/PAM 夹角以及 ||INT||/(||AGO||+||PAM||)，用于区分“两个背景同向且近似”与“两个背景方向不同”。由于共享C/A/P/CP项，三差分不是独立证据；ΔINT 的更大幅度可能来自两个方向相减。现有10 ns重复和有限区域均值不能排除起始姿势、短时松弛、背景取样或随机特征表示的影响。

### PAM consensus：幅度分布与对照审阅

以下每个范围是三个重复各自统计量的最小至最大值，不将30个相关块合并成独立样本。q10最小与q90最大给出三个重复描述区间的外包络，不是置信区间。

| 候选 | 分支 | 差分 | 块q10最小–q90最大 | 时间相干度范围 | 前/后半段余弦范围 | 相对stable范数比范围 | 相对distal范数比范围 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PACERGEN02432 | STATE_MOTION | Delta_AGO | 0.109–0.147 | 0.761–0.847 | 0.655–0.855 | 1.972–2.573 | 1.130–1.730 |
| PACERGEN02432 | STATE_MOTION | Delta_PAM | 0.095–0.165 | 0.767–0.885 | 0.709–0.847 | 1.838–3.434 | 0.942–1.955 |
| PACERGEN02432 | STATE_MOTION | Delta_INT | 0.117–0.180 | 0.562–0.779 | 0.207–0.726 | 1.104–2.379 | 0.775–1.762 |
| PACERGEN02432 | SIGNED_DRIFT | Delta_AGO | 0.146–0.209 | 0.283–0.423 | -0.017–0.424 | 0.967–1.312 | 0.608–1.094 |
| PACERGEN02432 | SIGNED_DRIFT | Delta_PAM | 0.137–0.204 | 0.288–0.347 | -0.215–0.081 | 0.921–1.050 | 0.631–0.810 |
| PACERGEN02432 | SIGNED_DRIFT | Delta_INT | 0.192–0.308 | 0.323–0.355 | -0.023–0.291 | 0.950–1.132 | 0.701–0.943 |
| PACERGEN00123 | STATE_MOTION | Delta_AGO | 0.102–0.140 | 0.771–0.821 | 0.668–0.829 | 1.696–2.735 | 0.928–1.348 |
| PACERGEN00123 | STATE_MOTION | Delta_PAM | 0.090–0.138 | 0.797–0.840 | 0.770–0.848 | 2.074–2.879 | 0.937–1.845 |
| PACERGEN00123 | STATE_MOTION | Delta_INT | 0.122–0.181 | 0.637–0.721 | 0.523–0.689 | 1.500–2.070 | 0.484–1.453 |
| PACERGEN00123 | SIGNED_DRIFT | Delta_AGO | 0.134–0.176 | 0.331–0.362 | -0.179–0.209 | 0.821–1.338 | 0.730–0.799 |
| PACERGEN00123 | SIGNED_DRIFT | Delta_PAM | 0.131–0.173 | 0.295–0.336 | -0.054–0.172 | 0.983–0.996 | 0.577–0.768 |
| PACERGEN00123 | SIGNED_DRIFT | Delta_INT | 0.189–0.258 | 0.311–0.350 | -0.005–0.073 | 0.930–1.165 | 0.677–0.771 |
| PACERGEN00462 | STATE_MOTION | Delta_AGO | 0.096–0.162 | 0.819–0.890 | 0.828–0.908 | 2.238–3.133 | 1.020–1.649 |
| PACERGEN00462 | STATE_MOTION | Delta_PAM | 0.100–0.157 | 0.810–0.891 | 0.709–0.893 | 1.914–3.361 | 1.482–2.726 |
| PACERGEN00462 | STATE_MOTION | Delta_INT | 0.102–0.167 | 0.640–0.752 | 0.536–0.715 | 1.179–1.514 | 0.697–1.054 |
| PACERGEN00462 | SIGNED_DRIFT | Delta_AGO | 0.135–0.205 | 0.275–0.391 | -0.209–0.293 | 0.787–1.018 | 0.592–0.833 |
| PACERGEN00462 | SIGNED_DRIFT | Delta_PAM | 0.136–0.196 | 0.328–0.334 | 0.158–0.245 | 0.973–1.049 | 0.671–0.839 |
| PACERGEN00462 | SIGNED_DRIFT | Delta_INT | 0.196–0.279 | 0.268–0.374 | -0.041–0.355 | 0.733–1.078 | 0.639–0.833 |

STATE_MOTION 的局部接触信号相对stable控制增强，并不代表相对distal控制始终增强。SIGNED_DRIFT 的接触区范数可与stable控制相近、低于distal控制，且方向常混合；控制区背景变化不能忽略。时间相干度与前/后半段余弦仅描述同一10 ns轨迹内部变化，不解决独立重复冲突，也不能证明达到稳态。

## 6. 科研展示图与可追溯性

每个图有同名PNG（220 dpi）和矢量PDF。英文图例避免跨平台缺字；图内数值对应CSV的原始精度字段，标注仅作四舍五入。全部来源见 `FIGURE_PROVENANCE_v01.csv`。

| 图文件前缀 | 内容与用途 |
| --- | --- |
| STATE_MOTION_replica_direction | 3候选×3差分；全部9区域、全部3重复对，固定−1至1余弦色标 |
| SIGNED_DRIFT_replica_direction | 3候选×3差分；全部9区域、全部3重复对，固定−1至1余弦色标 |
| STATE_MOTION_magnitude_distribution | 三个重复的块q10/median/q90与时间平均向量范数；仅同分支统一幅度轴 |
| SIGNED_DRIFT_magnitude_distribution | 三个重复的块q10/median/q90与时间平均向量范数；仅同分支统一幅度轴 |
| STATE_MOTION_region_specificity | 各重复相对stable-control的范数比，并保留distal-control；1表示相等 |
| SIGNED_DRIFT_region_specificity | 各重复相对stable-control的范数比，并保留distal-control；1表示相等 |
| region_dependence | 冻结区域残基Jaccard与结构扩散共享质量；评估信息依赖 |

`BLOCK_MAGNITUDES_v01.csv` 每行含候选/分支/差分/区域/重复/块时间及NPZ key，可追溯每个分位数。`INPUT_SHA256_v01.csv` 保存输入绝对路径、字节数、哈希及来源预期哈希；分析结束后已复核所读取输入哈希不变。`ANALYSIS_VALIDATION_v01.json` 给出检查计数及数值容差。

分析代码版本 v01，`analyze_v01.py` SHA256 `8141380a2dbc9db68ffc7ad525ba1864a72c7d5b0c88c369dc0e1112ca3a41b5`，仓库基线 HEAD `97e863dd4a0d1f256cde900579dd28e0b0b47d7b`。代码与环境完整记录在 `ANALYSIS_CODE_VERSION_v01.json`，报告渲染脚本为 `write_report_v01.py`。运行命令见该记录，全部输出仅写入本独立目录；没有导入或执行上游模型流水线，没有改动历史结果，也没有提交GitHub。

## 7. 最终回答：动态证据、冲突及实验优先事项

**02432：** 最可重复的是 STATE_MOTION 的 PAM 接触区 ΔPAM，以及局部 ΔAGO；正交 ΔPAM 较弱但三对同向。主要冲突是 PAM ΔINT、胞内 ΔPAM/ΔINT 和 SIGNED_DRIFT 的缺乏一致性。优先检验候选在有/无 probe 背景下是否产生可重复功能变化，尤其检验“有背景差异但未形成可重复交互”是否只是动态表示层面的现象。

**00123：** 最可重复的是 STATE_MOTION cooperativity ΔAGO，PAM consensus ΔAGO/ΔPAM 亦重复。主要冲突是正交 ΔAGO、PAM/extension ΔINT 与远端R2异常突出。优先区分候选自身作用、probe背景调制与非特异响应；胞内弱一致方向只作为次级探索，不应据此提前判定ago-PAM。

**00462：** 最可重复的是 STATE_MOTION PAM/cooperativity 的 ΔAGO和ΔPAM；同区域ΔAGO方向一致性比其余两者更整齐。正交 ΔAGO/ΔINT 有较弱一致支持。最明显冲突是 PAM consensus ΔINT 的R1–R3=−0.606、胞内所有主要差分方向不一致，且stable-control也有变化。优先验证候选自身功能和probe背景调制是否存在，并检验局部接触变化能否对应可重复受体输出；这些优先事项不是功能排序。

**实验设计建议（均为待验证假说）：** 对三个候选采用相同、预先定义的功能实验矩阵：vehicle、probe单独、候选单独、候选+probe；测候选浓度系列和probe完整浓度–反应曲线，分别评估基线效应、probe曲线位移和最大响应变化，不仅测单一浓度。保留跨日独立重复及同批对照，用受体相关功能读数，并加入受体缺失/阻断、细胞状态与检测干扰对照以区分非特异效应。实验是否支持PAM或ago-PAM应由预先约定的实验判据决定，不从本报告学习阈值或权重。

若基础功能矩阵重复，再将PAM接触/协同性区域与正交/胞内区域的结构假说做定点扰动或其他正交验证；选点应遵循冻结结构先验并预先注册，结合表达/基线功能对照。不要从区域RFF向量推出某残基是因果位点，不用单次突变效应直接宣称变构传播路径。三个候选都值得在一致实验条件下检验；当前结果支持的是不同程度的局部动态差异重复性，而不是已验证的PAM/ago-PAM药理结论。

## 附表：全部9区域、两个分支、三主要差分

每行范数为R1/R2/R3的 ||mean(v)||，余弦顺序为R1–R2/R1–R3/R2–R3。表中数值三位小数；CSV保留全精度，接近零的符号以CSV为准。

### PACERGEN02432 / STATE_MOTION

| 差分 | 区域 | 范数 R1/R2/R3 | 余弦 12/13/23 | 幅度CV |
| --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.108/0.094/0.103 | 0.562/0.424/0.744 | 0.071 |
| Delta_AGO | pam_contact_union | 0.093/0.082/0.093 | 0.412/0.327/0.668 | 0.073 |
| Delta_AGO | cooperativity_mutagenesis | 0.078/0.086/0.084 | 0.324/0.274/0.767 | 0.049 |
| Delta_AGO | compound110_extension | 0.082/0.102/0.099 | 0.148/-0.240/-0.127 | 0.114 |
| Delta_AGO | orthosteric_activation_core | 0.050/0.056/0.056 | 0.246/-0.289/0.479 | 0.058 |
| Delta_AGO | orthosteric_contact_union | 0.050/0.054/0.053 | 0.269/-0.314/0.428 | 0.040 |
| Delta_AGO | intracellular_microswitches | 0.051/0.058/0.029 | 0.233/-0.291/0.028 | 0.335 |
| Delta_AGO | stable_core_control | 0.042/0.047/0.044 | 0.024/-0.064/0.228 | 0.062 |
| Delta_AGO | distal_control | 0.076/0.083/0.059 | 0.457/0.289/0.286 | 0.167 |
| Delta_PAM | pam_contact_consensus | 0.132/0.092/0.087 | 0.742/0.729/0.583 | 0.238 |
| Delta_PAM | pam_contact_union | 0.127/0.079/0.086 | 0.727/0.747/0.534 | 0.265 |
| Delta_PAM | cooperativity_mutagenesis | 0.103/0.081/0.068 | 0.666/0.655/0.438 | 0.216 |
| Delta_PAM | compound110_extension | 0.115/0.069/0.071 | 0.137/-0.147/0.055 | 0.306 |
| Delta_PAM | orthosteric_activation_core | 0.044/0.043/0.060 | 0.225/0.349/0.274 | 0.199 |
| Delta_PAM | orthosteric_contact_union | 0.042/0.042/0.058 | 0.229/0.359/0.258 | 0.199 |
| Delta_PAM | intracellular_microswitches | 0.036/0.044/0.032 | 0.322/-0.158/-0.072 | 0.153 |
| Delta_PAM | stable_core_control | 0.039/0.036/0.048 | -0.052/0.267/-0.359 | 0.147 |
| Delta_PAM | distal_control | 0.068/0.078/0.093 | 0.289/0.445/0.140 | 0.159 |
| Delta_INT | pam_contact_consensus | 0.124/0.074/0.079 | 0.287/-0.182/-0.003 | 0.294 |
| Delta_INT | pam_contact_union | 0.133/0.067/0.070 | 0.221/-0.035/0.022 | 0.417 |
| Delta_INT | cooperativity_mutagenesis | 0.107/0.063/0.072 | 0.038/-0.272/0.004 | 0.286 |
| Delta_INT | compound110_extension | 0.147/0.100/0.134 | 0.445/-0.635/-0.340 | 0.192 |
| Delta_INT | orthosteric_activation_core | 0.075/0.053/0.050 | 0.053/-0.065/-0.128 | 0.234 |
| Delta_INT | orthosteric_contact_union | 0.074/0.051/0.050 | 0.113/-0.048/-0.189 | 0.230 |
| Delta_INT | intracellular_microswitches | 0.068/0.070/0.046 | 0.426/-0.445/-0.052 | 0.216 |
| Delta_INT | stable_core_control | 0.052/0.050/0.072 | 0.051/-0.164/-0.218 | 0.209 |
| Delta_INT | distal_control | 0.070/0.092/0.102 | 0.291/-0.117/-0.173 | 0.187 |

### PACERGEN02432 / SIGNED_DRIFT

| 差分 | 区域 | 范数 R1/R2/R3 | 余弦 12/13/23 | 幅度CV |
| --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.070/0.058/0.045 | 0.026/0.032/0.137 | 0.216 |
| Delta_AGO | pam_contact_union | 0.055/0.057/0.047 | -0.145/-0.115/0.205 | 0.095 |
| Delta_AGO | cooperativity_mutagenesis | 0.048/0.054/0.049 | -0.235/-0.118/0.215 | 0.068 |
| Delta_AGO | compound110_extension | 0.066/0.076/0.065 | -0.171/-0.155/0.215 | 0.090 |
| Delta_AGO | orthosteric_activation_core | 0.044/0.053/0.047 | -0.047/-0.232/0.428 | 0.101 |
| Delta_AGO | orthosteric_contact_union | 0.044/0.055/0.045 | -0.067/-0.229/0.347 | 0.126 |
| Delta_AGO | intracellular_microswitches | 0.043/0.050/0.040 | 0.238/0.036/-0.149 | 0.113 |
| Delta_AGO | stable_core_control | 0.054/0.051/0.047 | 0.175/-0.068/-0.183 | 0.067 |
| Delta_AGO | distal_control | 0.064/0.069/0.075 | -0.058/0.019/0.308 | 0.075 |
| Delta_PAM | pam_contact_consensus | 0.044/0.055/0.054 | -0.191/0.121/-0.329 | 0.117 |
| Delta_PAM | pam_contact_union | 0.041/0.051/0.051 | -0.059/-0.008/-0.191 | 0.116 |
| Delta_PAM | cooperativity_mutagenesis | 0.042/0.053/0.057 | -0.073/0.005/-0.438 | 0.153 |
| Delta_PAM | compound110_extension | 0.061/0.066/0.060 | -0.035/-0.209/-0.066 | 0.048 |
| Delta_PAM | orthosteric_activation_core | 0.039/0.037/0.045 | -0.067/-0.002/-0.213 | 0.103 |
| Delta_PAM | orthosteric_contact_union | 0.040/0.037/0.047 | -0.063/0.034/-0.175 | 0.134 |
| Delta_PAM | intracellular_microswitches | 0.046/0.044/0.051 | -0.041/-0.126/0.202 | 0.078 |
| Delta_PAM | stable_core_control | 0.048/0.054/0.051 | -0.220/-0.011/-0.021 | 0.056 |
| Delta_PAM | distal_control | 0.070/0.068/0.069 | 0.004/-0.019/-0.009 | 0.015 |
| Delta_INT | pam_contact_consensus | 0.069/0.090/0.072 | -0.005/0.054/-0.259 | 0.149 |
| Delta_INT | pam_contact_union | 0.065/0.081/0.067 | -0.021/0.038/-0.149 | 0.125 |
| Delta_INT | cooperativity_mutagenesis | 0.057/0.079/0.077 | -0.185/0.033/-0.329 | 0.168 |
| Delta_INT | compound110_extension | 0.090/0.101/0.097 | 0.094/-0.273/0.034 | 0.062 |
| Delta_INT | orthosteric_activation_core | 0.060/0.073/0.048 | -0.151/-0.084/0.134 | 0.204 |
| Delta_INT | orthosteric_contact_union | 0.061/0.072/0.054 | -0.157/-0.064/0.011 | 0.146 |
| Delta_INT | intracellular_microswitches | 0.066/0.066/0.066 | 0.306/-0.106/-0.043 | 0.003 |
| Delta_INT | stable_core_control | 0.070/0.080/0.075 | -0.069/-0.210/-0.235 | 0.066 |
| Delta_INT | distal_control | 0.099/0.096/0.094 | 0.024/-0.006/0.043 | 0.025 |

### PACERGEN00123 / STATE_MOTION

| 差分 | 区域 | 范数 R1/R2/R3 | 余弦 12/13/23 | 幅度CV |
| --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.090/0.105/0.096 | 0.704/0.608/0.430 | 0.082 |
| Delta_AGO | pam_contact_union | 0.090/0.096/0.085 | 0.642/0.441/0.333 | 0.066 |
| Delta_AGO | cooperativity_mutagenesis | 0.096/0.086/0.090 | 0.749/0.761/0.715 | 0.058 |
| Delta_AGO | compound110_extension | 0.137/0.098/0.087 | 0.510/-0.437/-0.100 | 0.248 |
| Delta_AGO | orthosteric_activation_core | 0.048/0.048/0.042 | -0.443/-0.032/-0.046 | 0.083 |
| Delta_AGO | orthosteric_contact_union | 0.049/0.047/0.040 | -0.392/-0.006/-0.040 | 0.106 |
| Delta_AGO | intracellular_microswitches | 0.035/0.051/0.041 | -0.120/-0.194/0.338 | 0.190 |
| Delta_AGO | stable_core_control | 0.053/0.039/0.037 | -0.232/-0.014/-0.148 | 0.201 |
| Delta_AGO | distal_control | 0.067/0.114/0.093 | -0.519/0.394/-0.481 | 0.259 |
| Delta_PAM | pam_contact_consensus | 0.097/0.084/0.101 | 0.624/0.315/0.692 | 0.094 |
| Delta_PAM | pam_contact_union | 0.088/0.078/0.108 | 0.493/0.148/0.724 | 0.168 |
| Delta_PAM | cooperativity_mutagenesis | 0.074/0.083/0.096 | 0.286/0.148/0.811 | 0.129 |
| Delta_PAM | compound110_extension | 0.130/0.071/0.083 | 0.625/-0.290/-0.056 | 0.329 |
| Delta_PAM | orthosteric_activation_core | 0.042/0.057/0.038 | 0.259/0.211/-0.144 | 0.217 |
| Delta_PAM | orthosteric_contact_union | 0.040/0.055/0.036 | 0.236/0.192/-0.216 | 0.228 |
| Delta_PAM | intracellular_microswitches | 0.035/0.036/0.031 | 0.122/0.189/0.413 | 0.086 |
| Delta_PAM | stable_core_control | 0.034/0.040/0.037 | 0.279/-0.094/-0.050 | 0.092 |
| Delta_PAM | distal_control | 0.069/0.089/0.054 | 0.089/0.362/0.045 | 0.248 |
| Delta_INT | pam_contact_consensus | 0.097/0.087/0.099 | 0.377/-0.407/-0.292 | 0.067 |
| Delta_INT | pam_contact_union | 0.114/0.082/0.106 | 0.316/-0.559/-0.177 | 0.162 |
| Delta_INT | cooperativity_mutagenesis | 0.111/0.075/0.067 | 0.063/-0.261/0.226 | 0.281 |
| Delta_INT | compound110_extension | 0.230/0.124/0.131 | 0.668/-0.689/-0.525 | 0.365 |
| Delta_INT | orthosteric_activation_core | 0.060/0.054/0.057 | -0.029/-0.050/-0.012 | 0.055 |
| Delta_INT | orthosteric_contact_union | 0.058/0.052/0.055 | -0.071/-0.016/-0.035 | 0.056 |
| Delta_INT | intracellular_microswitches | 0.044/0.064/0.047 | 0.185/0.003/0.365 | 0.209 |
| Delta_INT | stable_core_control | 0.065/0.052/0.048 | 0.187/-0.115/-0.368 | 0.158 |
| Delta_INT | distal_control | 0.067/0.180/0.099 | -0.332/0.249/-0.637 | 0.507 |

### PACERGEN00123 / SIGNED_DRIFT

| 差分 | 区域 | 范数 R1/R2/R3 | 余弦 12/13/23 | 幅度CV |
| --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.051/0.059/0.051 | -0.088/-0.085/0.238 | 0.089 |
| Delta_AGO | pam_contact_union | 0.051/0.052/0.047 | -0.177/-0.131/0.171 | 0.051 |
| Delta_AGO | cooperativity_mutagenesis | 0.060/0.047/0.055 | -0.121/0.163/-0.037 | 0.121 |
| Delta_AGO | compound110_extension | 0.096/0.066/0.070 | 0.042/-0.291/0.026 | 0.214 |
| Delta_AGO | orthosteric_activation_core | 0.037/0.036/0.037 | -0.171/-0.186/0.129 | 0.011 |
| Delta_AGO | orthosteric_contact_union | 0.036/0.037/0.037 | -0.137/-0.098/0.160 | 0.011 |
| Delta_AGO | intracellular_microswitches | 0.039/0.039/0.047 | 0.116/0.217/0.095 | 0.109 |
| Delta_AGO | stable_core_control | 0.059/0.044/0.062 | -0.273/-0.054/-0.211 | 0.173 |
| Delta_AGO | distal_control | 0.070/0.074/0.069 | -0.309/0.138/0.091 | 0.040 |
| Delta_PAM | pam_contact_consensus | 0.045/0.044/0.048 | 0.083/-0.085/0.035 | 0.047 |
| Delta_PAM | pam_contact_union | 0.043/0.045/0.055 | 0.025/-0.126/0.165 | 0.136 |
| Delta_PAM | cooperativity_mutagenesis | 0.037/0.048/0.054 | -0.026/-0.049/0.287 | 0.183 |
| Delta_PAM | compound110_extension | 0.093/0.071/0.061 | 0.459/0.080/0.018 | 0.219 |
| Delta_PAM | orthosteric_activation_core | 0.037/0.036/0.041 | 0.139/-0.114/-0.163 | 0.070 |
| Delta_PAM | orthosteric_contact_union | 0.038/0.037/0.040 | 0.165/-0.181/-0.087 | 0.035 |
| Delta_PAM | intracellular_microswitches | 0.040/0.043/0.047 | -0.033/0.122/-0.017 | 0.081 |
| Delta_PAM | stable_core_control | 0.046/0.044/0.048 | 0.009/-0.140/0.071 | 0.042 |
| Delta_PAM | distal_control | 0.074/0.076/0.062 | -0.112/-0.047/0.012 | 0.103 |
| Delta_INT | pam_contact_consensus | 0.066/0.074/0.075 | -0.121/-0.227/0.217 | 0.067 |
| Delta_INT | pam_contact_union | 0.066/0.073/0.077 | -0.170/-0.359/0.279 | 0.081 |
| Delta_INT | cooperativity_mutagenesis | 0.075/0.071/0.069 | -0.251/-0.193/0.060 | 0.037 |
| Delta_INT | compound110_extension | 0.163/0.107/0.090 | 0.372/-0.120/-0.046 | 0.318 |
| Delta_INT | orthosteric_activation_core | 0.056/0.053/0.061 | -0.043/-0.186/-0.040 | 0.069 |
| Delta_INT | orthosteric_contact_union | 0.055/0.054/0.058 | -0.042/-0.190/0.039 | 0.044 |
| Delta_INT | intracellular_microswitches | 0.059/0.056/0.074 | 0.087/0.263/-0.010 | 0.149 |
| Delta_INT | stable_core_control | 0.071/0.064/0.064 | -0.146/-0.108/-0.028 | 0.062 |
| Delta_INT | distal_control | 0.096/0.110/0.097 | -0.174/0.125/-0.029 | 0.077 |

### PACERGEN00462 / STATE_MOTION

| 差分 | 区域 | 范数 R1/R2/R3 | 余弦 12/13/23 | 幅度CV |
| --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.127/0.120/0.095 | 0.820/0.640/0.814 | 0.146 |
| Delta_AGO | pam_contact_union | 0.114/0.119/0.090 | 0.799/0.573/0.773 | 0.142 |
| Delta_AGO | cooperativity_mutagenesis | 0.106/0.109/0.087 | 0.868/0.749/0.819 | 0.120 |
| Delta_AGO | compound110_extension | 0.078/0.084/0.104 | 0.165/-0.077/0.780 | 0.150 |
| Delta_AGO | orthosteric_activation_core | 0.041/0.039/0.033 | 0.520/0.322/0.308 | 0.105 |
| Delta_AGO | orthosteric_contact_union | 0.040/0.037/0.033 | 0.452/0.314/0.178 | 0.087 |
| Delta_AGO | intracellular_microswitches | 0.055/0.047/0.034 | -0.287/-0.141/0.197 | 0.233 |
| Delta_AGO | stable_core_control | 0.041/0.054/0.031 | 0.216/0.117/-0.249 | 0.265 |
| Delta_AGO | distal_control | 0.078/0.073/0.094 | 0.088/0.488/0.166 | 0.135 |
| Delta_PAM | pam_contact_consensus | 0.093/0.117/0.127 | 0.746/0.572/0.739 | 0.157 |
| Delta_PAM | pam_contact_union | 0.081/0.109/0.128 | 0.647/0.455/0.716 | 0.225 |
| Delta_PAM | cooperativity_mutagenesis | 0.075/0.109/0.100 | 0.732/0.571/0.784 | 0.186 |
| Delta_PAM | compound110_extension | 0.069/0.089/0.090 | 0.046/-0.011/0.440 | 0.146 |
| Delta_PAM | orthosteric_activation_core | 0.040/0.049/0.038 | -0.009/0.221/0.081 | 0.132 |
| Delta_PAM | orthosteric_contact_union | 0.039/0.049/0.038 | -0.076/0.189/0.108 | 0.151 |
| Delta_PAM | intracellular_microswitches | 0.040/0.069/0.045 | 0.206/-0.258/-0.197 | 0.311 |
| Delta_PAM | stable_core_control | 0.033/0.035/0.067 | 0.085/-0.018/0.112 | 0.423 |
| Delta_PAM | distal_control | 0.063/0.063/0.047 | -0.003/-0.093/-0.060 | 0.163 |
| Delta_INT | pam_contact_consensus | 0.083/0.081/0.107 | 0.052/-0.606/0.077 | 0.163 |
| Delta_INT | pam_contact_union | 0.085/0.074/0.107 | 0.141/-0.575/-0.069 | 0.189 |
| Delta_INT | cooperativity_mutagenesis | 0.064/0.067/0.078 | 0.040/-0.430/0.121 | 0.106 |
| Delta_INT | compound110_extension | 0.097/0.114/0.153 | -0.412/-0.516/0.603 | 0.232 |
| Delta_INT | orthosteric_activation_core | 0.063/0.059/0.049 | 0.260/0.444/0.164 | 0.129 |
| Delta_INT | orthosteric_contact_union | 0.061/0.061/0.049 | 0.156/0.432/0.095 | 0.120 |
| Delta_INT | intracellular_microswitches | 0.058/0.093/0.051 | -0.278/-0.153/0.070 | 0.333 |
| Delta_INT | stable_core_control | 0.055/0.068/0.072 | 0.281/0.183/0.208 | 0.139 |
| Delta_INT | distal_control | 0.119/0.102/0.102 | -0.140/0.196/0.099 | 0.093 |

### PACERGEN00462 / SIGNED_DRIFT

| 差分 | 区域 | 范数 R1/R2/R3 | 余弦 12/13/23 | 幅度CV |
| --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.049/0.060/0.045 | 0.046/0.039/0.208 | 0.151 |
| Delta_AGO | pam_contact_union | 0.049/0.068/0.049 | 0.201/0.093/0.378 | 0.200 |
| Delta_AGO | cooperativity_mutagenesis | 0.045/0.057/0.046 | 0.201/0.260/0.220 | 0.130 |
| Delta_AGO | compound110_extension | 0.063/0.076/0.075 | 0.059/-0.130/0.182 | 0.097 |
| Delta_AGO | orthosteric_activation_core | 0.039/0.038/0.037 | 0.134/0.177/0.237 | 0.017 |
| Delta_AGO | orthosteric_contact_union | 0.039/0.038/0.038 | 0.137/0.247/0.220 | 0.017 |
| Delta_AGO | intracellular_microswitches | 0.044/0.057/0.047 | -0.401/-0.038/0.081 | 0.144 |
| Delta_AGO | stable_core_control | 0.058/0.059/0.058 | -0.113/0.036/-0.293 | 0.016 |
| Delta_AGO | distal_control | 0.075/0.072/0.076 | 0.046/-0.007/-0.223 | 0.029 |
| Delta_PAM | pam_contact_consensus | 0.052/0.052/0.053 | 0.080/0.205/-0.100 | 0.016 |
| Delta_PAM | pam_contact_union | 0.050/0.049/0.056 | 0.094/0.299/-0.045 | 0.068 |
| Delta_PAM | cooperativity_mutagenesis | 0.050/0.055/0.049 | 0.170/0.147/-0.066 | 0.065 |
| Delta_PAM | compound110_extension | 0.066/0.068/0.063 | -0.060/-0.035/-0.066 | 0.033 |
| Delta_PAM | orthosteric_activation_core | 0.038/0.038/0.036 | -0.011/0.003/0.106 | 0.024 |
| Delta_PAM | orthosteric_contact_union | 0.037/0.041/0.037 | -0.033/-0.051/0.128 | 0.056 |
| Delta_PAM | intracellular_microswitches | 0.044/0.048/0.049 | -0.021/0.166/-0.026 | 0.054 |
| Delta_PAM | stable_core_control | 0.052/0.050/0.055 | 0.016/-0.121/0.041 | 0.049 |
| Delta_PAM | distal_control | 0.077/0.062/0.064 | -0.126/0.058/0.128 | 0.122 |
| Delta_INT | pam_contact_consensus | 0.071/0.082/0.064 | -0.123/0.081/-0.074 | 0.131 |
| Delta_INT | pam_contact_union | 0.069/0.079/0.061 | -0.053/0.086/-0.251 | 0.124 |
| Delta_INT | cooperativity_mutagenesis | 0.064/0.074/0.062 | -0.035/0.060/-0.206 | 0.095 |
| Delta_INT | compound110_extension | 0.093/0.105/0.093 | -0.148/-0.103/0.093 | 0.070 |
| Delta_INT | orthosteric_activation_core | 0.059/0.048/0.049 | 0.191/0.172/0.000 | 0.112 |
| Delta_INT | orthosteric_contact_union | 0.058/0.051/0.050 | 0.156/0.159/-0.011 | 0.076 |
| Delta_INT | intracellular_microswitches | 0.062/0.080/0.069 | -0.171/0.186/0.025 | 0.129 |
| Delta_INT | stable_core_control | 0.065/0.078/0.087 | -0.100/0.033/-0.103 | 0.140 |
| Delta_INT | distal_control | 0.093/0.099/0.100 | -0.069/0.082/-0.108 | 0.039 |

