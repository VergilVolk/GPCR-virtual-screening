# Stage4 Prospective PACER-FKG v02 Scientific Evaluation

独立只读科学复核与解释 · 中文正式报告 · 复核日期：2026-10-03。来源为本地冻结产物；未使用此前聊天中的候选排序。新写入内容仅为本报告及配套科学摘要JSON。

## 1. Executive Summary

**完整性核验通过，存在一处需确认的源元数据差异。** 最终产物状态为`STAGE4_FULL_FROZEN_EVALUATION_COMPLETE`。已核验36/36条轨迹、3个候选×4个各自cluster-matched上下文×3个replica，36个Phase1完整缓存、72份Phase2a manifest及18份Phase2b区域向量文件。每条轨迹为10 ns、1000个10 ps原始帧，固定`[::5]`后为200个50 ps分析帧，构成10个名义1 ns的相关时间block。R1/R2/R3均EVALUATION_ONLY；未发现Stage4 R2进入历史拟合、重训练、refit、新阈值或outcome-driven tuning的产物/代码证据。

**必须优先说明：** `C:/projects/PACER_STAGE4_MD_backup/master_manifest.json`的`trajectory_ps=50`与实际认证的原始10 ps间隔不同。实际DCD、state.csv和冻结Stage4契约相互一致；现有特征使用的是10 ps原始数据→stride-5→50 ps分析数据。该差异不构成已发现的冻结数值状态污染，但执行元数据未完全自洽，需要项目组确认manifest是否为早期计划快照、或实际启动时存在未记录覆盖。未核验到证明该覆盖的服务器启动记录，不能自行假定来源；本报告不修正原文件，以下解释以认证的实际时间契约为前提。

**Delta_INT的核心结论：没有候选在多数关键区域、跨两branch具有稳健的三replica一致interaction方向。** PACER0010的STATE_MOTION/PAM-consensus出现R1–R3反向；PACER0027的consensus两branch全部配对为负或近零；PACER0073的SIGNED_DRIFT局部PAM区域同号，但STATE_MOTION对应不一致。局部同号不等于机制确定，不能合理声称已证实cooperative dynamic interaction；该强结论是 **not supported by the current frozen prospective evaluation**。

PACER0010在PAM-related Delta_PAM比Delta_AGO方向更一致，activation-core intrinsic方向不稳定；支持有限模式分离，不能证明无intrinsic效应。

PACER0027在STATE_MOTION的PAM/正交相关Delta_AGO具有较均衡局部复现；Delta_PAM在PAM相关区域有支持，SIGNED_DRIFT支持不全面；应重视intrinsic效应跟进。

PACER0073在PAM/mutagenesis区域的candidate-alone方向比probe-background更稳定，尤其SIGNED_DRIFT；不能据此贴inactive或agonist标签。

控制区域幅度同样明显，限制了严格PAM区域局部化主张；候选之间未构建综合分数、预测概率或药效排序。

## 2. Evaluation Scope and Claim Boundary

冻结结果明确声明：`Prospective frozen PACER-FKG dynamic differential evaluation only. No PAM/ago-PAM labels, probability, efficacy or independent-block inference.`

本轮READ / ANALYZE / REPORT ONLY：执行既有只读`--verify`，读取文件及保存数组，独立复算摘要与余弦；未调用MD、encoder推理、descriptor构建、normalize/RFF应用或图扩散运行路径，也未修改原始数据、历史冻结状态、适配器或已发表结果。表示空间方向的一致性不是激活方向；Delta_AGO、Delta_PAM、Delta_INT是历史contrast名称，不是候选的功能标签。未生成药理分类、预测概率、综合评分、统计显著性或候选药效排名。

## 3. Frozen Analysis Provenance

历史根目录：`C:/projects/GPCR-virtual-screening-fkg-v02/`。冻结manifest为`project/results/pacer_fkg_v02_longmd_v01/calibration/V02_FREEZE_MANIFEST.json`，SHA-256：`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`。与adapter强制pin一致，实际消费的定义、校准inventory、normalization、bandwidth、RFF权重/偏置和相关冻结实现均通过只读哈希认证。

校准源为历史apo、probe_only、compound110__candidate_no_probe、compound110__candidate_probe各自R2，每条50个block，共200个校准block；历史R1/R3不拟合。Stage4三个候选全部replica仅应用既有状态，不能将历史代码中的R2 CALIBRATION_REPLICA标记移植给Stage4 R2。冻结日期为2026-09-29；本轮验证序列化状态，不重建历史校准总体。

| Branch | 原始descriptor宽度 | 冻结bandwidth | RFF宽度 / seed |
| --- | --- | --- | --- |
| STATE_MOTION | 512 | 32.30188361260893 | 512 / 272340 |
| SIGNED_DRIFT | 256 | 29.29934899925964 | 512 / 272084 |

Normalization为逐channel历史median与`max(1.4826×MAD, 1e-6)`尺度；既有Gaussian-RBF RFF权重、偏置直接消费，不重新抽样。STATE_MOTION = `concat(mean(Z), sqrt(mean(diff(Z)^2)))`；SIGNED_DRIFT = `mean(diff(Z)) = (Z_last − Z_first)/19`，差分严格在block内部。两branch映射到不同随机特征坐标，不直接比较跨branch幅度或向量方向。

固定图为`M4_MULTISTRUCTURE_GRAPH_v01`，270节点、1304条序列化边；backbone边权1，contact边按既有support_fraction，去重后构造行归一化转移矩阵。restart diffusion为`F_(t+1)=(1−0.65)F_0+0.65 T F_t`，20步，随后对区域节点做均值。接触/图结构截止分别4.5 Å/8.0 Å。**本次实际九区域来自图文件的`regions`；已认证的历史G2_REGION_MAP_v02还包含旧面板名，不应把它的ACh_pocket等名称或旧边界替换实际Phase2b区域。**

冻结residue mapping为270个embedding_index 0–269，7TRS/R的M32至H466有序位点，匹配图节点。Stage4 OpenMM两条封端片段按唯一冻结序列投影，在原拓扑顺序保留A段195＋B段75残基、2139个ATOM14重原子；resindex 196、274、275、276是排除的构建体标准残基。chainID/segid或重新起算的PDB resid不是embedding身份。

主要认证资产（SHA-256来自冻结回执，当前字节已验证）：

| 资产路径（相对历史根目录） | SHA-256 |
| --- | --- |
| project/results/pacer_fkg_v02_longmd_v01/calibration/V02_FREEZE_MANIFEST.json | `b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd` |
| project/results/pacer_dc_four_context_v01/compound110/G2_RESIDUE_MAPPING_v01.csv | `3cb47de7fa555e635f66e6b1ff16846e6eaa95046c522fea6154d8e1f046b61d` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/V02_BLOCK_DEFINITIONS.json | `dfe5472b2b7ea9f93befa222cfdb15cf3538550fab59f1aeaeed9884e5916e0f` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/V02_CONTRAST_DEFINITIONS.json | `7b897db161e5fb62f23670ca187c889a1f9c998ef661b9839cd05c7c20a5ae2e` |
| project/results/pacer_dc_geom2vec_pilot_v01/M4_MULTISTRUCTURE_GRAPH_v01.json | `25593c11dac0f70ad6811c8c6aa7a4f06c1755f14397f1162eaa1157026be122` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/V02_METHOD_SPEC.json | `ccae0585576332e2ad19eaf8a813ccbc56a3c1fe9e8f0e6c226ec72ebebe9e59` |
| project/results/pacer_dc_four_context_v01/compound110/G2_REGION_MAP_v02.json | `a40fb6ca8131a2f6dfc18caf6f48b9930f1341c072b368e3cb3b7c3c21e593ec` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/V02_SEEDS.json | `1bb280da1cf38884c7b34203fe18a96358e77ffa35ab2e82fdaf8de185ef0bb6` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/V02_CALIBRATION_INVENTORY.json | `2338f26a1dcf10e7760a0cb0293b21a5603906bad39a27f7eb4527cabd8db96a` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/STATE_MOTION/normalization_center.npy | `2e1950951d57ca4775a6dd9e6e27f711c4d62f4c0dc352ff82c36eb83471b1ec` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/STATE_MOTION/normalization_scale.npy | `f9cadd72ac4b68e0a810d9c434c1094e167903f1683a46946226de5cbb2fb430` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/STATE_MOTION/bandwidth.json | `336c4202d0c6b715ddafd8fedfdef48a41e257e2dda5e505f21f1c9dab9597c1` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/STATE_MOTION/rff_weights.npy | `34c6979797ca7d5ca34742f9b9e008bdb6cc9fa1c7e82c4aec2e2406ea1d3619` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/STATE_MOTION/rff_bias.npy | `ae8dc92ed47e9edd49a56498ad78158c796477fa9d408e49b007b082173bf860` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/SIGNED_DRIFT/normalization_center.npy | `ac7d18780e6b096e95d923d4be9f08cb0956065a0fa2f37fae1b27c4bdb5e1bd` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/SIGNED_DRIFT/normalization_scale.npy | `a9caf251cc40d84c31fb0c7c469dc9fb1f9330b5040c572d8357cdfe45a50dd9` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/SIGNED_DRIFT/bandwidth.json | `882173a46651c36f4de1b6252b0f6b5cade813d2fa60c7c5111f6582d7e7ba72` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/SIGNED_DRIFT/rff_weights.npy | `fb1d5a41c2a60e6464b92d57b1be5ee2b318a6f838796a04ee803508712925c0` |
| project/results/pacer_fkg_v02_longmd_v01/calibration/SIGNED_DRIFT/rff_bias.npy | `a4b5bfa9a00a80c57f448981b37a3fbd2ead1d5553e199521cfa544300168069` |

完整只读`--verify`返回`STAGE4_FROZEN_EVALUATION_VERIFIED`。从18份既有NPZ的1458个数组核对810组replica摘要及810个配对余弦，与JSON最大绝对差均为0；保存context与contrast的线性关系最大差为1.0710×10⁻⁸，处于原float32计算舍入量级，不构成结果变化。未重新执行Phase2算法。

## 4. Stage4 Dataset and Temporal Design

| Candidate | Cluster | A / P | C / CP |
| --- | --- | --- | --- |
| PACER0010 | 9 | cluster9__apo / cluster9__probe_only | PACER0010__candidate_no_probe / PACER0010__candidate_probe |
| PACER0027 | 4 | cluster4__apo / cluster4__probe_only | PACER0027__candidate_no_probe / PACER0027__candidate_probe |
| PACER0073 | 0 | cluster0__apo / cluster0__probe_only | PACER0073__candidate_no_probe / PACER0073__candidate_probe |

R1=27101、R2=38201、R3=49301，全部EVALUATION_ONLY。总生产量36×10 ns=360 ns；每个上下文3条轨迹，禁止跨cluster借用apo/probe baseline。种子相同用于标识各上下文replica系列，不保证微观轨迹严格配对。

实际production为5,000,000步×2 fs。state.csv绝对时间记录14–10004 ps，固定10 ps间隔；应按生产步数和相邻间隔解释，不把绝对末时钟误写成生产时长10.004 ns。分析原始帧索引0,5,…,995；无插值。20分析帧组成一个block，19个内部lag×50 ps=950 ps，名义覆盖1 ns。10个block相关，不能当成10个独立样本。

35条progress与物理完成均一致；PACER0010/CP/R3 progress仍为running/9.9 ns，但按既有特定例外由1000帧、最终步数、有限坐标和兼容拓扑认证为物理完整。该例外未泛化，原progress未修改。早期IMPLEMENTATION_VALIDATION的“下载阻断”和builder_input_validation的dry_run/production_started=false均为阶段快照，不是最终缺失；拓扑修复验证和完整冻结回执提供后续状态。源manifest的cadence差异仍需另行确认。

## 5. Contrast Definitions

公式从最终JSON与历史冻结实现双重确认。先在各上下文的冻结RFF表示上做带符号线性组合，再使用固定图扩散和区域均值。以下A/P/C/CP均属于该候选自己的cluster。

### Delta_AGO

`C − A`：candidate-alone相对apo的动态表示扰动；用于检查intrinsic candidate-associated effect，不能直接称激动活性。

### Delta_PAM

`CP − P`：probe背景下加入candidate的动态表示变化；允许描述probe-background modulation-like signature。幅度大不证明PAM，亦未证明“只有probe存在才起作用”。

### Delta_INT

`CP − C − P + A = Delta_PAM − Delta_AGO`：表示空间二阶interaction/nonadditivity contrast。RFF非线性、采样差异、初始状态、各轨迹波动都可能影响该contrast；非零幅度不能单独确证物理协同、因果机制或药理cooperativity。

另有`Probe_effect=P−A`和`Probe_background_effect=CP−C`。两者也完整审阅，关键区域摘要置于附录及JSON；它们分别回答probe在无/有candidate时的背景变化，不能与Delta_PAM混称。

## 6. Region Definitions and Biological Interpretation

区域为冻结结构先验。PAM相关、activation相关和control是解释用途，不是经Stage4训练得到的分类目标。区域重叠（例如Y439跨PAM与正交面板，Y416跨extension与正交面板），经过图扩散后也会共享邻域信息；不能用多个区域当成独立证据票数。

| Region | 节点数 | 冻结sites | 解释边界 |
| --- | --- | --- | --- |
| pam_contact_consensus | 12 | Y89, Y92, I93, G96, F186, L190, N423, Q427, D432, W435, S436, Y439 | 冻结PAM结构接触共识面板，作为allosteric结构先验；接触归属不是功能标签。0010的Delta_PAM两branch有局部支持；0027的STATE_MOTION intrinsic方向最均衡；0073的intrinsic方向跨branch较一致。三者Delta_INT均不跨branch稳健。 |
| pam_contact_union | 16 | Y89, Y92, I93, G96, Q184, F186, L190, S191, N423, C426, Q427, D432, T433, W435, S436, Y439 | 更宽的PAM接触并集，不能当成独立于consensus的重复实验。0027的STATE_MOTION Delta_AGO明显同向，0010 Delta_PAM在两branch弱到中等同向；0073 Delta_PAM及STATE_MOTION Delta_INT出现replica分歧。 |
| cooperativity_mutagenesis | 13 | Y89, Y92, I93, K95, G96, Q184, F186, N423, D432, T433, W435, Y439, Y443 | 历史机制相关位点面板，名称不意味着本项目已做突变实验。0010 Delta_PAM两branch同向；0027/0073 Delta_AGO两branch同向。三者Delta_INT都未形成跨branch一致结论。 |
| compound110_extension | 8 | C185, Y416, M419, V420, S428, C429, I430, P431 | 历史allosteric agonist control结构及compound110背景延伸面板，解释为相邻结构区域，不能泛称PAM专属口袋。两branch的AGO/PAM/INT多数存在分歧，未给出稳定、普适的probe依赖或interaction定位。 |
| orthosteric_activation_core | 13 | D112, Y113, S116, N117, W164, A203, F204, W413, Y416, N417, Y439, C442, Y443 | 预定义activation-related结构位点面板。0027 Delta_AGO STATE_MOTION有局部一致intrinsic perturbation；0073较弱；0010不一致。三者对应SIGNED_DRIFT均不一致，不能标注已激活。 |
| orthosteric_contact_union | 14 | D112, Y113, S116, N117, V120, W164, A203, F204, W413, Y416, N417, Y439, C442, Y443 | 正交配体接触并集，与activation core大量重叠。0027 Delta_INT两branch配对余弦均正，但SIGNED_DRIFT最小0.0160，不能称稳健；0073的STATE_MOTION局部interaction未被SIGNED_DRIFT支持。 |
| intracellular_microswitches | 9 | D129, R130, Y131, F409, W413, P415, N445, P450, Y453 | D129/R130/Y131、F409/W413/P415及N445/P450/Y453面板。0010 Delta_INT STATE_MOTION局部同向，SIGNED_DRIFT混合；0027 intrinsic跨branch不统一；0073 R1的STATE_MOTION intrinsic方向与R2/R3反向。无统一胞内激活链路证据。 |
| stable_core_control | 12 | Y69, L71, F72, I125, I126, I210, M211, T212, V213, L214, H217, I218 | 按跨结构低CA方差、距所有配体≥12 Å且不在声明机制面板中确定的12个位点。0010 Delta_INT STATE_MOTION三配对0.4182/0.5069/0.4303，甚至比PAM-region interaction更一致；显示非机制专属变化。 |
| distal_control | 12 | A221, S222, R223, S224, R225, V226, A392, A393, R394, E395, R396, K397 | 历史远端对照，跨越构建体两段的边界邻近区域。各contrast幅度常高于PAM/正交面板；短轨迹和片段边缘动力学可影响解释，但本次未以数据证明具体原因。不能当作静默对照或用它建立新阈值。 |

数值阅读约定：下列B=block_magnitude_mean，N=mean_vector_norm；B/N三元组按R1/R2/R3，余弦按列R1–R2/R1–R3/R2–R3。B是10个相关block向量范数的均值，N是block平均向量的范数；前者可大而后者因方向抵消较小。q10/median/q90仅为时间block描述，不是置信区间；全精度摘要保存在配套JSON，不复制512维signed_direction_vector。余弦>0只描述符号，不引入“稳健”新阈值。

## 7. PACER0010

### STATE_MOTION

在PAM相关区域，probe背景效应的方向一致性优于candidate-alone：consensus的Delta_PAM三配对为0.4758/0.1001/0.3592，而Delta_AGO为0.0434/−0.4107/0.3363。cooperativity_mutagenesis的Delta_PAM为0.4619/0.2976/0.5283，支持局部probe-background signature；但extension的R1–R3仅0.0589。

| Contrast | Region | B (R1/R2/R3) | N (R1/R2/R3) | cos12 | cos13 | cos23 |
| --- | --- | --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.10804 / 0.09205 / 0.08999 | 0.08214 / 0.05904 / 0.05789 | 0.0434 | -0.4107 | 0.3363 |
| Delta_AGO | pam_contact_union | 0.10382 / 0.08864 / 0.08229 | 0.07750 / 0.05693 / 0.04924 | 0.1987 | -0.3409 | 0.2546 |
| Delta_AGO | cooperativity_mutagenesis | 0.11053 / 0.08598 / 0.08509 | 0.08527 / 0.05515 / 0.05799 | 0.0150 | -0.3654 | 0.4051 |
| Delta_AGO | compound110_extension | 0.10672 / 0.11285 / 0.09877 | 0.05384 / 0.07982 / 0.04548 | 0.4408 | 0.0422 | 0.0235 |
| Delta_AGO | orthosteric_activation_core | 0.08017 / 0.07619 / 0.06951 | 0.05492 / 0.05513 / 0.03418 | -0.0970 | -0.0510 | 0.0736 |
| Delta_AGO | orthosteric_contact_union | 0.07898 / 0.07530 / 0.07010 | 0.05308 / 0.05307 / 0.03393 | -0.1390 | -0.0370 | 0.0382 |
| Delta_AGO | intracellular_microswitches | 0.07391 / 0.09830 / 0.08329 | 0.03906 / 0.08111 / 0.04046 | 0.1660 | 0.1872 | 0.1227 |
| Delta_AGO | stable_core_control | 0.08644 / 0.10361 / 0.07909 | 0.05132 / 0.07622 / 0.03645 | 0.2716 | 0.2938 | 0.3977 |
| Delta_AGO | distal_control | 0.11882 / 0.13177 / 0.12807 | 0.06156 / 0.08128 / 0.08490 | -0.2098 | -0.1758 | -0.2909 |
| Delta_PAM | pam_contact_consensus | 0.08572 / 0.08256 / 0.08738 | 0.05587 / 0.05269 / 0.04874 | 0.4758 | 0.1001 | 0.3592 |
| Delta_PAM | pam_contact_union | 0.08281 / 0.07891 / 0.08338 | 0.05431 / 0.04814 / 0.04265 | 0.3517 | 0.0391 | 0.2561 |
| Delta_PAM | cooperativity_mutagenesis | 0.08160 / 0.07568 / 0.08248 | 0.05859 / 0.04766 / 0.04886 | 0.4619 | 0.2976 | 0.5283 |
| Delta_PAM | compound110_extension | 0.10470 / 0.11728 / 0.11857 | 0.05876 / 0.06263 / 0.05914 | 0.4794 | 0.0589 | 0.1669 |
| Delta_PAM | orthosteric_activation_core | 0.05924 / 0.06278 / 0.05362 | 0.03470 / 0.03933 / 0.02967 | 0.1102 | -0.0671 | 0.5004 |
| Delta_PAM | orthosteric_contact_union | 0.05842 / 0.06404 / 0.05391 | 0.03324 / 0.04096 / 0.02924 | 0.1056 | -0.0746 | 0.5179 |
| Delta_PAM | intracellular_microswitches | 0.08574 / 0.08405 / 0.07859 | 0.05536 / 0.05851 / 0.04436 | -0.0812 | 0.0613 | 0.1784 |
| Delta_PAM | stable_core_control | 0.08408 / 0.07515 / 0.08783 | 0.05016 / 0.03907 / 0.05223 | 0.0063 | 0.2661 | 0.1321 |
| Delta_PAM | distal_control | 0.12168 / 0.15650 / 0.15547 | 0.05012 / 0.11275 / 0.10971 | -0.2827 | 0.2609 | -0.5550 |
| Delta_INT | pam_contact_consensus | 0.14097 / 0.12047 / 0.12643 | 0.09995 / 0.07578 / 0.07392 | 0.1043 | -0.4621 | 0.1329 |
| Delta_INT | pam_contact_union | 0.13282 / 0.11333 / 0.11875 | 0.09169 / 0.06838 / 0.06486 | 0.1783 | -0.3277 | 0.0410 |
| Delta_INT | cooperativity_mutagenesis | 0.14184 / 0.11018 / 0.11973 | 0.10422 / 0.06487 / 0.07555 | -0.0842 | -0.4181 | 0.3368 |
| Delta_INT | compound110_extension | 0.15152 / 0.16715 / 0.15420 | 0.07911 / 0.10611 / 0.07550 | 0.4612 | -0.1251 | -0.1603 |
| Delta_INT | orthosteric_activation_core | 0.10381 / 0.09579 / 0.08675 | 0.06953 / 0.06718 / 0.03918 | 0.0114 | -0.0937 | 0.1777 |
| Delta_INT | orthosteric_contact_union | 0.10203 / 0.09739 / 0.08832 | 0.06667 / 0.06831 / 0.03848 | -0.0356 | -0.1122 | 0.1686 |
| Delta_INT | intracellular_microswitches | 0.11405 / 0.13142 / 0.11261 | 0.06860 / 0.10107 / 0.05551 | 0.4363 | 0.2066 | 0.1793 |
| Delta_INT | stable_core_control | 0.11652 / 0.13543 / 0.12257 | 0.05940 / 0.09788 / 0.07462 | 0.4182 | 0.5069 | 0.4303 |
| Delta_INT | distal_control | 0.16928 / 0.21400 / 0.22947 | 0.07074 / 0.14470 / 0.17118 | 0.0389 | -0.0297 | -0.5922 |

### SIGNED_DRIFT

Delta_PAM在cooperativity_mutagenesis为0.2533/0.3041/0.2773，consensus为0.1246/0.1280/0.1979，方向支持有限但均同号。Delta_AGO在activation core三配对均负（−0.2626/−0.0955/−0.2035），不支持稳定的正交口袋intrinsic方向。

| Contrast | Region | B (R1/R2/R3) | N (R1/R2/R3) | cos12 | cos13 | cos23 |
| --- | --- | --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.15555 / 0.14947 / 0.14939 | 0.04768 / 0.04292 / 0.04540 | 0.0823 | 0.1050 | 0.1141 |
| Delta_AGO | pam_contact_union | 0.15028 / 0.14596 / 0.14314 | 0.04509 / 0.04383 / 0.04120 | -0.0358 | -0.0831 | 0.2041 |
| Delta_AGO | cooperativity_mutagenesis | 0.15132 / 0.14211 / 0.13929 | 0.05065 / 0.04455 / 0.04125 | 0.3446 | 0.0850 | 0.1195 |
| Delta_AGO | compound110_extension | 0.20475 / 0.19334 / 0.20669 | 0.07400 / 0.07058 / 0.06395 | -0.2298 | -0.1122 | -0.0164 |
| Delta_AGO | orthosteric_activation_core | 0.12754 / 0.12525 / 0.13040 | 0.04444 / 0.03723 / 0.04463 | -0.2626 | -0.0955 | -0.2035 |
| Delta_AGO | orthosteric_contact_union | 0.12923 / 0.12634 / 0.13185 | 0.04381 / 0.04064 / 0.04634 | -0.1954 | -0.1136 | -0.3142 |
| Delta_AGO | intracellular_microswitches | 0.13942 / 0.13519 / 0.14287 | 0.03988 / 0.04194 / 0.04867 | 0.0365 | -0.0972 | -0.2728 |
| Delta_AGO | stable_core_control | 0.16685 / 0.16529 / 0.19025 | 0.05830 / 0.06070 / 0.05916 | -0.0645 | 0.2906 | -0.1053 |
| Delta_AGO | distal_control | 0.21934 / 0.23797 / 0.24291 | 0.06560 / 0.08609 / 0.06811 | -0.2684 | 0.2740 | -0.2668 |
| Delta_PAM | pam_contact_consensus | 0.16116 / 0.14402 / 0.15090 | 0.04785 / 0.04427 / 0.04473 | 0.1246 | 0.1280 | 0.1979 |
| Delta_PAM | pam_contact_union | 0.15233 / 0.14183 / 0.14301 | 0.04475 / 0.04480 / 0.04318 | 0.1200 | 0.0529 | 0.0714 |
| Delta_PAM | cooperativity_mutagenesis | 0.15362 / 0.13722 / 0.14511 | 0.05075 / 0.04387 / 0.04414 | 0.2533 | 0.3041 | 0.2773 |
| Delta_PAM | compound110_extension | 0.21117 / 0.21937 / 0.21310 | 0.06913 / 0.07329 / 0.07577 | 0.0717 | 0.2980 | 0.2130 |
| Delta_PAM | orthosteric_activation_core | 0.12377 / 0.12442 / 0.12217 | 0.04168 / 0.04487 / 0.03989 | 0.0692 | -0.0284 | 0.2476 |
| Delta_PAM | orthosteric_contact_union | 0.12243 / 0.12641 / 0.12313 | 0.04256 / 0.04986 / 0.04121 | 0.1539 | -0.0212 | 0.2554 |
| Delta_PAM | intracellular_microswitches | 0.15737 / 0.14675 / 0.14981 | 0.05794 / 0.04703 / 0.04445 | -0.0513 | -0.0081 | 0.0064 |
| Delta_PAM | stable_core_control | 0.16179 / 0.16510 / 0.17142 | 0.04177 / 0.05025 / 0.04971 | 0.2278 | 0.0621 | -0.0138 |
| Delta_PAM | distal_control | 0.23903 / 0.23128 / 0.23658 | 0.07320 / 0.08554 / 0.07752 | -0.1205 | 0.2451 | -0.2660 |
| Delta_INT | pam_contact_consensus | 0.22480 / 0.21502 / 0.21844 | 0.06717 / 0.06100 / 0.06633 | 0.0189 | -0.0658 | 0.0718 |
| Delta_INT | pam_contact_union | 0.21832 / 0.20717 / 0.20620 | 0.06666 / 0.05674 / 0.06167 | -0.0062 | -0.1467 | 0.0579 |
| Delta_INT | cooperativity_mutagenesis | 0.21336 / 0.20362 / 0.20869 | 0.06323 / 0.05585 / 0.06497 | 0.0162 | 0.0370 | 0.0811 |
| Delta_INT | compound110_extension | 0.30014 / 0.29469 / 0.29935 | 0.11360 / 0.09876 / 0.09181 | -0.0383 | 0.1681 | -0.0584 |
| Delta_INT | orthosteric_activation_core | 0.17884 / 0.17449 / 0.18225 | 0.06495 / 0.06144 / 0.05507 | 0.0137 | -0.0060 | -0.1517 |
| Delta_INT | orthosteric_contact_union | 0.18007 / 0.17989 / 0.18234 | 0.06430 / 0.07261 / 0.05542 | 0.0926 | -0.0706 | -0.2284 |
| Delta_INT | intracellular_microswitches | 0.20557 / 0.19720 / 0.20040 | 0.06605 / 0.05398 / 0.06249 | 0.0973 | -0.0490 | -0.0331 |
| Delta_INT | stable_core_control | 0.22562 / 0.23855 / 0.24847 | 0.06129 / 0.08387 / 0.07801 | -0.0316 | 0.1599 | -0.0947 |
| Delta_INT | distal_control | 0.32491 / 0.33676 / 0.34211 | 0.10156 / 0.14377 / 0.12015 | -0.3357 | 0.3065 | -0.4652 |

### Delta_AGO

跨replica稳定的intrinsic效应未获整体支持。STATE_MOTION activation core余弦为−0.0970/−0.0510/0.0736；SIGNED_DRIFT对应均负。局部SIGNED_DRIFT consensus虽三配对均小正值（0.0823/0.1050/0.1141），不能据此否定其他candidate-alone扰动。

### Delta_PAM

在consensus、union和mutagenesis，两branch均出现同号配对，mutagenesis尤为一致。STATE_MOTION consensus幅度反而低于Delta_AGO（0.08572/0.08256/0.08738，相比0.10804/0.09205/0.08999），显示幅度与复现性不可混同。activation core及microswitches仍有跨replica方向分歧。

### Delta_INT

consensus在STATE_MOTION为0.1043/−0.4621/0.1329，SIGNED_DRIFT为0.0189/−0.0658/0.0718。R1与R3在前者明显反向；后者接近正交。mutagenesis也不跨branch一致。microswitches和stable control仅STATE_MOTION三配对同号，不能作为跨branch协同证据。

### Region-level interpretation

PAM相关区域的Delta_PAM方向复现优于activation core的Delta_AGO，支持有限的动态模式分离；但extension、activation及control不呈整齐分区。distal control的Delta_PAM STATE_MOTION幅度0.12168–0.15650，高于consensus，且方向不一致；区域特异性受到限制。

### Replica consistency

consensus Delta_AGO和Delta_INT主要受R1–R3分歧影响，但R1–R2也不强，不能把单一replica判为异常后剔除。Delta_PAM的局部方向较一致，并不等于全受体或两branch全面稳健。

### Candidate-specific interpretation

可作为检验probe-dependent effect与intrinsic effect能否分离的实验候选。证据集中在PAM-related Delta_PAM及activation-core Delta_AGO的不一致，不是“无intrinsic活性”证明。稳健cooperative dynamic interaction：not supported by the current frozen prospective evaluation。

## 8. PACER0027

### STATE_MOTION

candidate-alone方向在consensus为0.7957/0.7167/0.7912，union为0.7938/0.7647/0.8055，mutagenesis为0.7012/0.7783/0.8075；均比其Delta_PAM更均衡。Delta_PAM consensus为0.4716/0.8288/0.4351，仍显示局部复现，但extension的R2与R1、R3均负。

| Contrast | Region | B (R1/R2/R3) | N (R1/R2/R3) | cos12 | cos13 | cos23 |
| --- | --- | --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.11155 / 0.12398 / 0.12505 | 0.08744 / 0.10653 / 0.10796 | 0.7957 | 0.7167 | 0.7912 |
| Delta_AGO | pam_contact_union | 0.11209 / 0.11563 / 0.12509 | 0.08877 / 0.09700 / 0.10807 | 0.7938 | 0.7647 | 0.8055 |
| Delta_AGO | cooperativity_mutagenesis | 0.10313 / 0.10473 / 0.11395 | 0.08172 / 0.08568 / 0.09816 | 0.7012 | 0.7783 | 0.8075 |
| Delta_AGO | compound110_extension | 0.11090 / 0.14848 / 0.14510 | 0.07111 / 0.11224 / 0.11072 | 0.1895 | 0.4202 | 0.7894 |
| Delta_AGO | orthosteric_activation_core | 0.05586 / 0.06931 / 0.08378 | 0.03138 / 0.05146 / 0.06743 | 0.4508 | 0.5122 | 0.6488 |
| Delta_AGO | orthosteric_contact_union | 0.05508 / 0.06921 / 0.08265 | 0.03004 / 0.05002 / 0.06557 | 0.4552 | 0.5091 | 0.6300 |
| Delta_AGO | intracellular_microswitches | 0.08324 / 0.06818 / 0.07529 | 0.06035 / 0.03587 / 0.04441 | 0.1890 | -0.1058 | -0.3338 |
| Delta_AGO | stable_core_control | 0.08027 / 0.07185 / 0.07049 | 0.04515 / 0.03264 / 0.03274 | 0.1445 | -0.0214 | 0.3286 |
| Delta_AGO | distal_control | 0.12312 / 0.11700 / 0.11934 | 0.06566 / 0.05793 / 0.06440 | 0.1501 | 0.0490 | 0.1204 |
| Delta_PAM | pam_contact_consensus | 0.12370 / 0.10061 / 0.10880 | 0.10630 / 0.06712 / 0.08734 | 0.4716 | 0.8288 | 0.4351 |
| Delta_PAM | pam_contact_union | 0.11361 / 0.09218 / 0.10729 | 0.09392 / 0.05771 / 0.08653 | 0.4125 | 0.7351 | 0.3633 |
| Delta_PAM | cooperativity_mutagenesis | 0.11582 / 0.08793 / 0.10550 | 0.10037 / 0.05443 / 0.08588 | 0.5224 | 0.8655 | 0.5331 |
| Delta_PAM | compound110_extension | 0.11662 / 0.10912 / 0.13531 | 0.07276 / 0.05308 / 0.07736 | -0.1129 | 0.6592 | -0.0824 |
| Delta_PAM | orthosteric_activation_core | 0.07336 / 0.06747 / 0.07248 | 0.05356 / 0.03858 / 0.04748 | 0.2311 | 0.1425 | -0.0571 |
| Delta_PAM | orthosteric_contact_union | 0.07375 / 0.06644 / 0.07155 | 0.05321 / 0.03744 / 0.04512 | 0.1620 | 0.1937 | -0.0397 |
| Delta_PAM | intracellular_microswitches | 0.06414 / 0.07681 / 0.07060 | 0.03134 / 0.04180 / 0.03865 | 0.1218 | 0.1094 | 0.1622 |
| Delta_PAM | stable_core_control | 0.07814 / 0.08007 / 0.08959 | 0.04198 / 0.04442 / 0.06404 | -0.4345 | 0.0330 | -0.1861 |
| Delta_PAM | distal_control | 0.11501 / 0.13401 / 0.11812 | 0.06478 / 0.09602 / 0.07242 | -0.0298 | -0.0117 | 0.3586 |
| Delta_INT | pam_contact_consensus | 0.11480 / 0.13833 / 0.10829 | 0.06853 / 0.09781 / 0.06196 | -0.1851 | -0.0912 | -0.0019 |
| Delta_INT | pam_contact_union | 0.11889 / 0.13379 / 0.10561 | 0.07551 / 0.09356 / 0.06042 | 0.1080 | 0.1064 | 0.1209 |
| Delta_INT | cooperativity_mutagenesis | 0.10655 / 0.11657 / 0.09917 | 0.06339 / 0.07569 / 0.05114 | -0.3102 | -0.1790 | 0.2060 |
| Delta_INT | compound110_extension | 0.14292 / 0.18490 / 0.16399 | 0.05950 / 0.13151 / 0.07951 | -0.0647 | -0.0699 | 0.5140 |
| Delta_INT | orthosteric_activation_core | 0.09229 / 0.09815 / 0.09527 | 0.06241 / 0.06845 / 0.06608 | 0.1841 | 0.0778 | 0.5458 |
| Delta_INT | orthosteric_contact_union | 0.09272 / 0.09596 / 0.09760 | 0.06202 / 0.06548 / 0.06811 | 0.1990 | 0.1699 | 0.5439 |
| Delta_INT | intracellular_microswitches | 0.10192 / 0.10365 / 0.10407 | 0.06340 / 0.05853 / 0.05863 | -0.1222 | 0.1459 | 0.2458 |
| Delta_INT | stable_core_control | 0.11421 / 0.11069 / 0.11420 | 0.06036 / 0.06212 / 0.07474 | -0.2307 | -0.2318 | -0.0438 |
| Delta_INT | distal_control | 0.16136 / 0.19252 / 0.16752 | 0.08367 / 0.12802 / 0.09824 | 0.0086 | -0.1031 | 0.3750 |

### SIGNED_DRIFT

Delta_AGO在consensus、union、mutagenesis有三配对正值，程度远低于STATE_MOTION。Delta_PAM consensus为−0.0153/−0.0086/0.2165，union为−0.1426/−0.1628/0.3152；主要由R1相对R2/R3分离。mutagenesis Delta_PAM局部三配对仍均正。

| Contrast | Region | B (R1/R2/R3) | N (R1/R2/R3) | cos12 | cos13 | cos23 |
| --- | --- | --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.15241 / 0.16139 / 0.15698 | 0.04484 / 0.05623 / 0.06168 | 0.1399 | 0.1625 | 0.1339 |
| Delta_AGO | pam_contact_union | 0.14849 / 0.15546 / 0.15519 | 0.04650 / 0.05326 / 0.06348 | 0.2161 | 0.3708 | 0.2711 |
| Delta_AGO | cooperativity_mutagenesis | 0.15042 / 0.15490 / 0.15323 | 0.05122 / 0.05563 / 0.05956 | 0.3053 | 0.3217 | 0.2344 |
| Delta_AGO | compound110_extension | 0.21401 / 0.21971 / 0.21538 | 0.05862 / 0.06776 / 0.07477 | 0.0128 | 0.1954 | 0.2210 |
| Delta_AGO | orthosteric_activation_core | 0.12461 / 0.12734 / 0.12587 | 0.04184 / 0.03606 / 0.03869 | 0.1263 | -0.0503 | -0.0375 |
| Delta_AGO | orthosteric_contact_union | 0.12540 / 0.12777 / 0.12614 | 0.04172 / 0.03703 / 0.03919 | 0.1277 | 0.0276 | 0.0476 |
| Delta_AGO | intracellular_microswitches | 0.14029 / 0.14062 / 0.14861 | 0.04405 / 0.04612 / 0.04815 | 0.1040 | -0.0496 | 0.1276 |
| Delta_AGO | stable_core_control | 0.18482 / 0.16174 / 0.15539 | 0.06466 / 0.05328 / 0.04298 | -0.1373 | -0.1865 | -0.0088 |
| Delta_AGO | distal_control | 0.24513 / 0.23962 / 0.24700 | 0.06773 / 0.07732 / 0.07888 | 0.0711 | -0.0152 | 0.1493 |
| Delta_PAM | pam_contact_consensus | 0.16141 / 0.15067 / 0.15524 | 0.04777 / 0.04815 / 0.05788 | -0.0153 | -0.0086 | 0.2165 |
| Delta_PAM | pam_contact_union | 0.15300 / 0.14188 / 0.15442 | 0.04874 / 0.04802 / 0.06191 | -0.1426 | -0.1628 | 0.3152 |
| Delta_PAM | cooperativity_mutagenesis | 0.15638 / 0.14835 / 0.14477 | 0.04673 / 0.04458 / 0.05784 | 0.2939 | 0.2566 | 0.2439 |
| Delta_PAM | compound110_extension | 0.20947 / 0.21153 / 0.23638 | 0.06671 / 0.07555 / 0.08797 | -0.2297 | -0.1947 | 0.2712 |
| Delta_PAM | orthosteric_activation_core | 0.12404 / 0.12928 / 0.12413 | 0.04856 / 0.04920 / 0.03251 | 0.1122 | 0.0954 | -0.2143 |
| Delta_PAM | orthosteric_contact_union | 0.12302 / 0.12949 / 0.12514 | 0.04639 / 0.04946 / 0.03365 | 0.1410 | 0.1139 | -0.1008 |
| Delta_PAM | intracellular_microswitches | 0.15166 / 0.14971 / 0.14793 | 0.05466 / 0.05189 / 0.04711 | -0.1068 | -0.1445 | -0.1153 |
| Delta_PAM | stable_core_control | 0.16219 / 0.15613 / 0.15984 | 0.05231 / 0.04956 / 0.05031 | 0.0330 | -0.1865 | -0.1938 |
| Delta_PAM | distal_control | 0.23455 / 0.24259 / 0.23042 | 0.06539 / 0.08363 / 0.06344 | 0.0101 | 0.1516 | 0.0156 |
| Delta_INT | pam_contact_consensus | 0.21900 / 0.22549 / 0.21774 | 0.05977 / 0.07148 / 0.07565 | -0.1520 | -0.0079 | -0.1322 |
| Delta_INT | pam_contact_union | 0.20883 / 0.21606 / 0.21353 | 0.06601 / 0.06591 / 0.07160 | -0.1757 | 0.0407 | -0.0614 |
| Delta_INT | cooperativity_mutagenesis | 0.20901 / 0.21341 / 0.20880 | 0.05410 / 0.06135 / 0.07494 | -0.0200 | -0.0764 | -0.0582 |
| Delta_INT | compound110_extension | 0.28677 / 0.30115 / 0.31419 | 0.08786 / 0.10162 / 0.08938 | -0.1217 | 0.0793 | -0.0220 |
| Delta_INT | orthosteric_activation_core | 0.17188 / 0.19014 / 0.18140 | 0.06296 / 0.06459 / 0.05115 | 0.2707 | 0.1932 | -0.1135 |
| Delta_INT | orthosteric_contact_union | 0.17136 / 0.18970 / 0.18263 | 0.05937 / 0.06487 / 0.05361 | 0.3091 | 0.2355 | 0.0160 |
| Delta_INT | intracellular_microswitches | 0.19417 / 0.20602 / 0.20688 | 0.06851 / 0.05977 / 0.06642 | 0.0426 | 0.1136 | -0.0847 |
| Delta_INT | stable_core_control | 0.24015 / 0.22832 / 0.22007 | 0.07771 / 0.07022 / 0.06293 | 0.2948 | -0.1136 | -0.1521 |
| Delta_INT | distal_control | 0.33920 / 0.34386 / 0.32525 | 0.09424 / 0.12826 / 0.09256 | 0.0652 | 0.1062 | 0.0990 |

### Delta_AGO

局部intrinsic candidate-associated perturbation较清楚。STATE_MOTION activation core为0.4508/0.5122/0.6488，且R1/R2/R3均值向量范数0.03138/0.05146/0.06743，幅度与方向有支持。SIGNED_DRIFT activation core为0.1263/−0.0503/−0.0375；microswitches也不一致，因此不是完整稳定的激活通路证据。

### Delta_PAM

STATE_MOTION PAM相关区域的probe-background effect可复现；SIGNED_DRIFT只有mutagenesis较一致，consensus/union不支持广泛跨branch方向复现。candidate-alone consensus幅度0.11155/0.12398/0.12505，与Delta_PAM 0.12370/0.10061/0.10880相近，R2/R3甚至更高；不能把后者强度当作优于其他候选的依据。

### Delta_INT

consensus在STATE_MOTION三配对为−0.1851/−0.0912/−0.0019，SIGNED_DRIFT为−0.1520/−0.0079/−0.1322；不支持一致interaction方向。union仅STATE_MOTION三小正值。orthosteric_contact_union两branch均同号，但SIGNED_DRIFT R2–R3仅0.0160；activation core的SIGNED_DRIFT R2–R3为−0.1135。

### Region-level interpretation

STATE_MOTION intrinsic效应在PAM-related与正交核心均有一致性，stable/distal control方向则较弱或混合，提供部分方向上的区域区分。SIGNED_DRIFT控制幅度可同样高，microswitches不支持统一激活链路；“intrinsic activation-like perturbation risk”可作为跟进理由，不能等同intrinsic agonism。

### Replica consistency

STATE_MOTION Delta_AGO是此候选最明确的局部复现模式；Delta_PAM extension及SIGNED_DRIFT PAM区域具有replica-specific response。consensus Delta_INT没有任何一个配对构成有力正向支持，不能用其他区域单一漂亮配对救济整体interaction结论。

### Candidate-specific interpretation

适合优先安排candidate-alone与组合功能实验，以排查intrinsic效应并检验probe背景变化是否叠加。其强signal带来解释风险与实验需求，不能自动推出更优的PAM候选。稳健cooperative dynamic interaction：not supported by the current frozen prospective evaluation。

## 9. PACER0073

### STATE_MOTION

Delta_AGO在consensus为0.2757/0.5543/0.2464、mutagenesis为0.4705/0.6001/0.2388；Delta_PAM consensus较弱（0.1854/0.0969/0.2787），union有负配对（−0.0486/−0.0274/0.4514）。candidate-alone在这些局部区域更稳定。

| Contrast | Region | B (R1/R2/R3) | N (R1/R2/R3) | cos12 | cos13 | cos23 |
| --- | --- | --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.09005 / 0.10551 / 0.10137 | 0.05784 / 0.08550 / 0.06585 | 0.2757 | 0.5543 | 0.2464 |
| Delta_AGO | pam_contact_union | 0.08886 / 0.09849 / 0.09297 | 0.05528 / 0.07665 / 0.05595 | 0.4225 | 0.4209 | 0.1629 |
| Delta_AGO | cooperativity_mutagenesis | 0.08334 / 0.10136 / 0.08410 | 0.05036 / 0.08389 / 0.04320 | 0.4705 | 0.6001 | 0.2388 |
| Delta_AGO | compound110_extension | 0.11110 / 0.11499 / 0.12247 | 0.05252 / 0.07123 / 0.06910 | 0.3138 | 0.3092 | -0.1631 |
| Delta_AGO | orthosteric_activation_core | 0.06991 / 0.08367 / 0.08569 | 0.04294 / 0.06631 / 0.05453 | 0.3025 | 0.3171 | 0.1354 |
| Delta_AGO | orthosteric_contact_union | 0.06918 / 0.08375 / 0.08403 | 0.04312 / 0.06645 / 0.05303 | 0.2549 | 0.2724 | 0.1065 |
| Delta_AGO | intracellular_microswitches | 0.09366 / 0.08368 / 0.09967 | 0.06378 / 0.05602 / 0.07360 | -0.5368 | -0.3937 | 0.4742 |
| Delta_AGO | stable_core_control | 0.07355 / 0.08134 / 0.08081 | 0.04254 / 0.04760 / 0.04510 | -0.1570 | 0.0852 | -0.2608 |
| Delta_AGO | distal_control | 0.13374 / 0.13758 / 0.12267 | 0.08849 / 0.08266 / 0.08100 | 0.4342 | -0.2051 | 0.0317 |
| Delta_PAM | pam_contact_consensus | 0.08938 / 0.09818 / 0.08691 | 0.05851 / 0.06872 / 0.05942 | 0.1854 | 0.0969 | 0.2787 |
| Delta_PAM | pam_contact_union | 0.09167 / 0.09847 / 0.08415 | 0.06287 / 0.07068 / 0.05717 | -0.0486 | -0.0274 | 0.4514 |
| Delta_PAM | cooperativity_mutagenesis | 0.08142 / 0.08720 / 0.07367 | 0.05272 / 0.05768 / 0.04420 | 0.3929 | -0.0177 | 0.2448 |
| Delta_PAM | compound110_extension | 0.11425 / 0.14495 / 0.13414 | 0.06241 / 0.10135 / 0.09906 | 0.0009 | 0.1032 | 0.6509 |
| Delta_PAM | orthosteric_activation_core | 0.07888 / 0.06757 / 0.07333 | 0.05695 / 0.04547 / 0.05424 | 0.4438 | -0.1143 | 0.1862 |
| Delta_PAM | orthosteric_contact_union | 0.08132 / 0.06791 / 0.07368 | 0.06103 / 0.04585 / 0.05457 | 0.5072 | -0.2122 | 0.0554 |
| Delta_PAM | intracellular_microswitches | 0.09458 / 0.08270 / 0.08459 | 0.07208 / 0.05323 / 0.05643 | -0.3593 | -0.3191 | 0.3503 |
| Delta_PAM | stable_core_control | 0.08206 / 0.07953 / 0.07662 | 0.04987 / 0.04669 / 0.04461 | -0.1053 | 0.0401 | 0.2409 |
| Delta_PAM | distal_control | 0.12031 / 0.13276 / 0.10131 | 0.06708 / 0.08706 / 0.05867 | -0.1394 | -0.2288 | 0.4533 |
| Delta_INT | pam_contact_consensus | 0.12010 / 0.12612 / 0.11633 | 0.07057 / 0.08820 / 0.05726 | -0.1664 | -0.0140 | 0.4230 |
| Delta_INT | pam_contact_union | 0.11641 / 0.12901 / 0.11197 | 0.06720 / 0.09419 / 0.05130 | -0.1719 | -0.0365 | 0.3438 |
| Delta_INT | cooperativity_mutagenesis | 0.10205 / 0.11296 / 0.10648 | 0.05436 / 0.07406 / 0.04750 | -0.1128 | -0.0582 | 0.3410 |
| Delta_INT | compound110_extension | 0.15600 / 0.18656 / 0.17199 | 0.07637 / 0.12902 / 0.10130 | -0.1380 | 0.1008 | 0.2941 |
| Delta_INT | orthosteric_activation_core | 0.11313 / 0.11730 / 0.11473 | 0.07890 / 0.09419 / 0.07893 | 0.6030 | 0.2110 | 0.2615 |
| Delta_INT | orthosteric_contact_union | 0.11327 / 0.11858 / 0.11401 | 0.08007 / 0.09578 / 0.07858 | 0.6148 | 0.0983 | 0.1606 |
| Delta_INT | intracellular_microswitches | 0.11213 / 0.11338 / 0.11782 | 0.06522 / 0.07331 / 0.07276 | -0.1181 | 0.1831 | 0.1915 |
| Delta_INT | stable_core_control | 0.10394 / 0.11988 / 0.10739 | 0.05718 / 0.07342 / 0.05542 | 0.0004 | -0.0502 | -0.1015 |
| Delta_INT | distal_control | 0.19362 / 0.18082 / 0.15724 | 0.13176 / 0.09986 / 0.10114 | 0.2860 | -0.2468 | 0.1389 |

### SIGNED_DRIFT

Delta_AGO consensus为0.3692/0.4680/0.4553，union为0.3456/0.5072/0.4478，mutagenesis为0.3683/0.4691/0.5193；与STATE_MOTION的局部intrinsic模式相互支持。Delta_PAM consensus仅−0.0003/0.0489/0.2457，缺少同样复现。

| Contrast | Region | B (R1/R2/R3) | N (R1/R2/R3) | cos12 | cos13 | cos23 |
| --- | --- | --- | --- | --- | --- | --- |
| Delta_AGO | pam_contact_consensus | 0.16907 / 0.16706 / 0.15736 | 0.06087 / 0.05471 / 0.06588 | 0.3692 | 0.4680 | 0.4553 |
| Delta_AGO | pam_contact_union | 0.15911 / 0.16228 / 0.15247 | 0.06131 / 0.05204 / 0.06451 | 0.3456 | 0.5072 | 0.4478 |
| Delta_AGO | cooperativity_mutagenesis | 0.15789 / 0.16010 / 0.14853 | 0.05213 / 0.05149 / 0.05948 | 0.3683 | 0.4691 | 0.5193 |
| Delta_AGO | compound110_extension | 0.22791 / 0.22108 / 0.22221 | 0.06869 / 0.07688 / 0.06657 | 0.0502 | -0.0348 | 0.1743 |
| Delta_AGO | orthosteric_activation_core | 0.13290 / 0.12931 / 0.14181 | 0.04504 / 0.04117 / 0.04218 | 0.0070 | -0.1132 | 0.0626 |
| Delta_AGO | orthosteric_contact_union | 0.13683 / 0.13115 / 0.14251 | 0.04566 / 0.04020 / 0.04292 | -0.0308 | -0.1182 | -0.0303 |
| Delta_AGO | intracellular_microswitches | 0.14674 / 0.15457 / 0.15958 | 0.05050 / 0.05497 / 0.04873 | -0.2771 | -0.0076 | 0.1312 |
| Delta_AGO | stable_core_control | 0.16426 / 0.16367 / 0.16970 | 0.04876 / 0.05803 / 0.04914 | 0.2060 | 0.1632 | 0.1074 |
| Delta_AGO | distal_control | 0.22798 / 0.23907 / 0.23259 | 0.07366 / 0.08180 / 0.06500 | 0.3866 | 0.2189 | 0.2331 |
| Delta_PAM | pam_contact_consensus | 0.15663 / 0.14760 / 0.15013 | 0.05200 / 0.04878 / 0.05613 | -0.0003 | 0.0489 | 0.2457 |
| Delta_PAM | pam_contact_union | 0.14890 / 0.14559 / 0.14374 | 0.05119 / 0.04653 / 0.05039 | -0.0633 | 0.0310 | 0.1595 |
| Delta_PAM | cooperativity_mutagenesis | 0.15173 / 0.14766 / 0.14305 | 0.04987 / 0.05034 / 0.04915 | 0.1383 | -0.0334 | 0.2894 |
| Delta_PAM | compound110_extension | 0.22819 / 0.21278 / 0.21126 | 0.07755 / 0.06544 / 0.06196 | -0.0186 | 0.1436 | -0.0662 |
| Delta_PAM | orthosteric_activation_core | 0.12777 / 0.12484 / 0.12155 | 0.04944 / 0.04186 / 0.03911 | 0.2542 | 0.3142 | 0.2457 |
| Delta_PAM | orthosteric_contact_union | 0.12929 / 0.12571 / 0.12309 | 0.04999 / 0.04183 / 0.03876 | 0.1822 | 0.2669 | 0.1370 |
| Delta_PAM | intracellular_microswitches | 0.14919 / 0.15771 / 0.14615 | 0.04715 / 0.05878 / 0.04464 | -0.2358 | 0.1094 | -0.1145 |
| Delta_PAM | stable_core_control | 0.16344 / 0.15782 / 0.16201 | 0.04628 / 0.04680 / 0.04638 | 0.0219 | 0.0508 | 0.0751 |
| Delta_PAM | distal_control | 0.22753 / 0.24599 / 0.22502 | 0.06041 / 0.08239 / 0.09589 | 0.0978 | 0.2620 | 0.3133 |
| Delta_INT | pam_contact_consensus | 0.22293 / 0.21518 / 0.21786 | 0.06724 / 0.07690 / 0.07486 | 0.0786 | 0.1597 | 0.3178 |
| Delta_INT | pam_contact_union | 0.21121 / 0.21179 / 0.20717 | 0.06515 / 0.07232 / 0.07526 | 0.0967 | 0.2473 | 0.3159 |
| Delta_INT | cooperativity_mutagenesis | 0.21660 / 0.21023 / 0.21136 | 0.06088 / 0.06933 / 0.07133 | 0.0274 | 0.0958 | 0.2762 |
| Delta_INT | compound110_extension | 0.31452 / 0.29805 / 0.30323 | 0.09382 / 0.10496 / 0.07907 | -0.0883 | 0.1004 | 0.0920 |
| Delta_INT | orthosteric_activation_core | 0.18598 / 0.17860 / 0.18450 | 0.07663 / 0.05776 / 0.05152 | 0.1749 | -0.0343 | -0.0326 |
| Delta_INT | orthosteric_contact_union | 0.19040 / 0.17863 / 0.18739 | 0.07722 / 0.05852 / 0.05370 | 0.1969 | -0.0677 | -0.0760 |
| Delta_INT | intracellular_microswitches | 0.20793 / 0.22217 / 0.21230 | 0.05986 / 0.07080 / 0.06401 | 0.0099 | 0.1318 | 0.0564 |
| Delta_INT | stable_core_control | 0.23938 / 0.22841 / 0.22122 | 0.06640 / 0.07333 / 0.06210 | 0.1325 | 0.2191 | 0.0391 |
| Delta_INT | distal_control | 0.31981 / 0.35134 / 0.32284 | 0.08727 / 0.10260 / 0.10107 | 0.0748 | -0.0374 | -0.1254 |

### Delta_AGO

PAM-related区域的candidate-alone效应比probe背景效应稳定。STATE_MOTION activation core为0.3025/0.3171/0.1354，但SIGNED_DRIFT对应0.0070/−0.1132/0.0626；microswitches的STATE_MOTION为−0.5368/−0.3937/0.4742，R1与R2/R3反向。不能据此推导一致激动方向。

### Delta_PAM

consensus两branch总体弱；mutagenesis不均衡。SIGNED_DRIFT activation core三配对0.2542/0.3142/0.2457全部为正，而STATE_MOTION为0.4438/−0.1143/0.1862，提示局部probe-background扰动，不能泛化成全面无效或全面有效。

### Delta_INT

STATE_MOTION consensus为−0.1664/−0.0140/0.4230，union为−0.1719/−0.0365/0.3438；R2/R3较一致而R1分离。SIGNED_DRIFT consensus为0.0786/0.1597/0.3178，union为0.0967/0.2473/0.3159，为局部同号interaction pattern。正交核心则STATE_MOTION同号、SIGNED_DRIFT两个负配对，跨branch不支持稳健机制。

### Region-level interpretation

PAM-related区域的intrinsic方向是该候选的主要复现特征；orthosteric与intracellular响应并不等价且存在明显分歧。distal control Delta_AGO的SIGNED_DRIFT也有0.3866/0.2189/0.2331同号配对；不能称为严格局部化。

### Replica consistency

R1在STATE_MOTION microswitches的Delta_AGO/Delta_PAM，以及PAM-related Delta_INT中与R2/R3呈不同方向；SIGNED_DRIFT并未统一复现该分组。所有replica应保留，不允许据结果剔除R1。

### Candidate-specific interpretation

在PAM相关面板，更符合candidate-alone driven perturbation的描述，probe效应仍具有区域与branch依赖性。适合检验intrinsic效应及其与probe是否可分离；不能贴inactive、PAM或ago-PAM标签。稳健cooperative dynamic interaction：not supported by the current frozen prospective evaluation。

## 10. Cross-candidate Comparison

| Candidate | Probe-dependent consistency | Intrinsic-effect consistency | Delta_INT consistency | Orthosteric intrinsic perturbation | PAM-region localization | Interpretation |
| --- | --- | --- | --- | --- | --- | --- |
| PACER0010 | PAM/mutagenesis局部两branch同号，activation/microswitches不全面 | 正交核心两branch不稳定 | consensus两branch混合，STATE_MOTION R1–R3负 | 不支持跨replica稳定核心方向；不等于无intrinsic效应 | 有方向分离，控制幅度同样高 | 用于检验probe/intrinsic能否分离 |
| PACER0027 | STATE_MOTION PAM面板较一致，SIGNED_DRIFT有限 | STATE_MOTION PAM/正交核心最均衡，SIGNED_DRIFT局部支持 | consensus两branch不一致；正交contact仅弱局部同号 | STATE_MOTION有局部intrinsic风险，SIGNED_DRIFT/胞内不统一 | 方向上部分区分，非严格局部幅度特异 | 需要candidate-alone功能跟进 |
| PACER0073 | PAM面板弱于intrinsic；正交SIGNED_DRIFT局部同号 | PAM/mutagenesis两branch较一致，activation不统一 | SIGNED_DRIFT PAM局部同号，STATE_MOTION分歧 | STATE_MOTION较弱同号，SIGNED_DRIFT不一致 | intrinsic方向在PAM面板较稳定，distal也有变化 | 局部更像candidate-alone驱动；机制未定 |

以上只比较各自matched contrast的复现模式。没有跨cluster原始baseline合并、没有总体药效优胜者、没有综合分数或outcome-driven rank。对“higher experimental priority”的合理表达应是某个待检验问题优先，而不是PAM成功率更高。

## 11. Replica Robustness

Delta_INT的完整九区域×两branch×三候选余弦如下，避免只报告R1–R3。所有replica保留，未以结果好坏剔除。

| Candidate | Branch | Delta_INT region | cos12 | cos13 | cos23 |
| --- | --- | --- | --- | --- | --- |
| PACER0010 | STATE_MOTION | pam_contact_consensus | 0.1043 | -0.4621 | 0.1329 |
| PACER0010 | STATE_MOTION | pam_contact_union | 0.1783 | -0.3277 | 0.0410 |
| PACER0010 | STATE_MOTION | cooperativity_mutagenesis | -0.0842 | -0.4181 | 0.3368 |
| PACER0010 | STATE_MOTION | compound110_extension | 0.4612 | -0.1251 | -0.1603 |
| PACER0010 | STATE_MOTION | orthosteric_activation_core | 0.0114 | -0.0937 | 0.1777 |
| PACER0010 | STATE_MOTION | orthosteric_contact_union | -0.0356 | -0.1122 | 0.1686 |
| PACER0010 | STATE_MOTION | intracellular_microswitches | 0.4363 | 0.2066 | 0.1793 |
| PACER0010 | STATE_MOTION | stable_core_control | 0.4182 | 0.5069 | 0.4303 |
| PACER0010 | STATE_MOTION | distal_control | 0.0389 | -0.0297 | -0.5922 |
| PACER0010 | SIGNED_DRIFT | pam_contact_consensus | 0.0189 | -0.0658 | 0.0718 |
| PACER0010 | SIGNED_DRIFT | pam_contact_union | -0.0062 | -0.1467 | 0.0579 |
| PACER0010 | SIGNED_DRIFT | cooperativity_mutagenesis | 0.0162 | 0.0370 | 0.0811 |
| PACER0010 | SIGNED_DRIFT | compound110_extension | -0.0383 | 0.1681 | -0.0584 |
| PACER0010 | SIGNED_DRIFT | orthosteric_activation_core | 0.0137 | -0.0060 | -0.1517 |
| PACER0010 | SIGNED_DRIFT | orthosteric_contact_union | 0.0926 | -0.0706 | -0.2284 |
| PACER0010 | SIGNED_DRIFT | intracellular_microswitches | 0.0973 | -0.0490 | -0.0331 |
| PACER0010 | SIGNED_DRIFT | stable_core_control | -0.0316 | 0.1599 | -0.0947 |
| PACER0010 | SIGNED_DRIFT | distal_control | -0.3357 | 0.3065 | -0.4652 |
| PACER0027 | STATE_MOTION | pam_contact_consensus | -0.1851 | -0.0912 | -0.0019 |
| PACER0027 | STATE_MOTION | pam_contact_union | 0.1080 | 0.1064 | 0.1209 |
| PACER0027 | STATE_MOTION | cooperativity_mutagenesis | -0.3102 | -0.1790 | 0.2060 |
| PACER0027 | STATE_MOTION | compound110_extension | -0.0647 | -0.0699 | 0.5140 |
| PACER0027 | STATE_MOTION | orthosteric_activation_core | 0.1841 | 0.0778 | 0.5458 |
| PACER0027 | STATE_MOTION | orthosteric_contact_union | 0.1990 | 0.1699 | 0.5439 |
| PACER0027 | STATE_MOTION | intracellular_microswitches | -0.1222 | 0.1459 | 0.2458 |
| PACER0027 | STATE_MOTION | stable_core_control | -0.2307 | -0.2318 | -0.0438 |
| PACER0027 | STATE_MOTION | distal_control | 0.0086 | -0.1031 | 0.3750 |
| PACER0027 | SIGNED_DRIFT | pam_contact_consensus | -0.1520 | -0.0079 | -0.1322 |
| PACER0027 | SIGNED_DRIFT | pam_contact_union | -0.1757 | 0.0407 | -0.0614 |
| PACER0027 | SIGNED_DRIFT | cooperativity_mutagenesis | -0.0200 | -0.0764 | -0.0582 |
| PACER0027 | SIGNED_DRIFT | compound110_extension | -0.1217 | 0.0793 | -0.0220 |
| PACER0027 | SIGNED_DRIFT | orthosteric_activation_core | 0.2707 | 0.1932 | -0.1135 |
| PACER0027 | SIGNED_DRIFT | orthosteric_contact_union | 0.3091 | 0.2355 | 0.0160 |
| PACER0027 | SIGNED_DRIFT | intracellular_microswitches | 0.0426 | 0.1136 | -0.0847 |
| PACER0027 | SIGNED_DRIFT | stable_core_control | 0.2948 | -0.1136 | -0.1521 |
| PACER0027 | SIGNED_DRIFT | distal_control | 0.0652 | 0.1062 | 0.0990 |
| PACER0073 | STATE_MOTION | pam_contact_consensus | -0.1664 | -0.0140 | 0.4230 |
| PACER0073 | STATE_MOTION | pam_contact_union | -0.1719 | -0.0365 | 0.3438 |
| PACER0073 | STATE_MOTION | cooperativity_mutagenesis | -0.1128 | -0.0582 | 0.3410 |
| PACER0073 | STATE_MOTION | compound110_extension | -0.1380 | 0.1008 | 0.2941 |
| PACER0073 | STATE_MOTION | orthosteric_activation_core | 0.6030 | 0.2110 | 0.2615 |
| PACER0073 | STATE_MOTION | orthosteric_contact_union | 0.6148 | 0.0983 | 0.1606 |
| PACER0073 | STATE_MOTION | intracellular_microswitches | -0.1181 | 0.1831 | 0.1915 |
| PACER0073 | STATE_MOTION | stable_core_control | 0.0004 | -0.0502 | -0.1015 |
| PACER0073 | STATE_MOTION | distal_control | 0.2860 | -0.2468 | 0.1389 |
| PACER0073 | SIGNED_DRIFT | pam_contact_consensus | 0.0786 | 0.1597 | 0.3178 |
| PACER0073 | SIGNED_DRIFT | pam_contact_union | 0.0967 | 0.2473 | 0.3159 |
| PACER0073 | SIGNED_DRIFT | cooperativity_mutagenesis | 0.0274 | 0.0958 | 0.2762 |
| PACER0073 | SIGNED_DRIFT | compound110_extension | -0.0883 | 0.1004 | 0.0920 |
| PACER0073 | SIGNED_DRIFT | orthosteric_activation_core | 0.1749 | -0.0343 | -0.0326 |
| PACER0073 | SIGNED_DRIFT | orthosteric_contact_union | 0.1969 | -0.0677 | -0.0760 |
| PACER0073 | SIGNED_DRIFT | intracellular_microswitches | 0.0099 | 0.1318 | 0.0564 |
| PACER0073 | SIGNED_DRIFT | stable_core_control | 0.1325 | 0.2191 | 0.0391 |
| PACER0073 | SIGNED_DRIFT | distal_control | 0.0748 | -0.0374 | -0.1254 |

符号层面的局部发现：0010仅STATE_MOTION microswitches/stable control、SIGNED_DRIFT mutagenesis三配对皆正，两个branch没有同一面板同时皆正；0027的orthosteric_contact_union在两branch皆正，但SIGNED_DRIFT最低0.0160，其余重要PAM面板不统一；0073的SIGNED_DRIFT PAM/mutagenesis/microswitches/control局部皆正，而STATE_MOTION这些面板并不全正。不是任何候选的跨多数关键区域稳健interaction证据。

显著方向分离的例子：0010 STATE_MOTION consensus Delta_INT cos13=−0.4621；0073 STATE_MOTION microswitches Delta_AGO cos12=−0.5368、cos13=−0.3937、cos23=0.4742；0027 SIGNED_DRIFT union Delta_PAM cos12=−0.1426、cos13=−0.1628、cos23=0.3152。负余弦不能用于指定哪个replica“错了”，也不是拮抗活性证据。

## 12. Region-specific Findings

**pam_contact_consensus**：冻结PAM结构接触共识面板，作为allosteric结构先验；接触归属不是功能标签。0010的Delta_PAM两branch有局部支持；0027的STATE_MOTION intrinsic方向最均衡；0073的intrinsic方向跨branch较一致。三者Delta_INT均不跨branch稳健。

**pam_contact_union**：更宽的PAM接触并集，不能当成独立于consensus的重复实验。0027的STATE_MOTION Delta_AGO明显同向，0010 Delta_PAM在两branch弱到中等同向；0073 Delta_PAM及STATE_MOTION Delta_INT出现replica分歧。

**cooperativity_mutagenesis**：历史机制相关位点面板，名称不意味着本项目已做突变实验。0010 Delta_PAM两branch同向；0027/0073 Delta_AGO两branch同向。三者Delta_INT都未形成跨branch一致结论。

**compound110_extension**：历史allosteric agonist control结构及compound110背景延伸面板，解释为相邻结构区域，不能泛称PAM专属口袋。两branch的AGO/PAM/INT多数存在分歧，未给出稳定、普适的probe依赖或interaction定位。

**orthosteric_activation_core**：预定义activation-related结构位点面板。0027 Delta_AGO STATE_MOTION有局部一致intrinsic perturbation；0073较弱；0010不一致。三者对应SIGNED_DRIFT均不一致，不能标注已激活。

**orthosteric_contact_union**：正交配体接触并集，与activation core大量重叠。0027 Delta_INT两branch配对余弦均正，但SIGNED_DRIFT最小0.0160，不能称稳健；0073的STATE_MOTION局部interaction未被SIGNED_DRIFT支持。

**intracellular_microswitches**：D129/R130/Y131、F409/W413/P415及N445/P450/Y453面板。0010 Delta_INT STATE_MOTION局部同向，SIGNED_DRIFT混合；0027 intrinsic跨branch不统一；0073 R1的STATE_MOTION intrinsic方向与R2/R3反向。无统一胞内激活链路证据。

**stable_core_control**：按跨结构低CA方差、距所有配体≥12 Å且不在声明机制面板中确定的12个位点。0010 Delta_INT STATE_MOTION三配对0.4182/0.5069/0.4303，甚至比PAM-region interaction更一致；显示非机制专属变化。

**distal_control**：历史远端对照，跨越构建体两段的边界邻近区域。各contrast幅度常高于PAM/正交面板；短轨迹和片段边缘动力学可影响解释，但本次未以数据证明具体原因。不能当作静默对照或用它建立新阈值。

不能从这些九区域均值直接判定“全受体完全同向变化”：最终产物没有保存全270节点global-region汇总，也不宜为本次解释新设global score。可以据control与多个面板的明显幅度判断变化不严格局限于机制面板；是否为完全统一的全局模式，则not supported by the current frozen prospective evaluation。

## 13. Controls and Specificity

以下采用Delta_PAM直接列出控制幅度，不用它建立新归一化比率或选择阈值。

| Candidate | Branch | Region | Delta_PAM B | Delta_PAM N | cos12 | cos13 | cos23 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PACER0010 | STATE_MOTION | pam_contact_consensus | 0.08572 / 0.08256 / 0.08738 | 0.05587 / 0.05269 / 0.04874 | 0.4758 | 0.1001 | 0.3592 |
| PACER0010 | STATE_MOTION | orthosteric_activation_core | 0.05924 / 0.06278 / 0.05362 | 0.03470 / 0.03933 / 0.02967 | 0.1102 | -0.0671 | 0.5004 |
| PACER0010 | STATE_MOTION | stable_core_control | 0.08408 / 0.07515 / 0.08783 | 0.05016 / 0.03907 / 0.05223 | 0.0063 | 0.2661 | 0.1321 |
| PACER0010 | STATE_MOTION | distal_control | 0.12168 / 0.15650 / 0.15547 | 0.05012 / 0.11275 / 0.10971 | -0.2827 | 0.2609 | -0.5550 |
| PACER0010 | SIGNED_DRIFT | pam_contact_consensus | 0.16116 / 0.14402 / 0.15090 | 0.04785 / 0.04427 / 0.04473 | 0.1246 | 0.1280 | 0.1979 |
| PACER0010 | SIGNED_DRIFT | orthosteric_activation_core | 0.12377 / 0.12442 / 0.12217 | 0.04168 / 0.04487 / 0.03989 | 0.0692 | -0.0284 | 0.2476 |
| PACER0010 | SIGNED_DRIFT | stable_core_control | 0.16179 / 0.16510 / 0.17142 | 0.04177 / 0.05025 / 0.04971 | 0.2278 | 0.0621 | -0.0138 |
| PACER0010 | SIGNED_DRIFT | distal_control | 0.23903 / 0.23128 / 0.23658 | 0.07320 / 0.08554 / 0.07752 | -0.1205 | 0.2451 | -0.2660 |
| PACER0027 | STATE_MOTION | pam_contact_consensus | 0.12370 / 0.10061 / 0.10880 | 0.10630 / 0.06712 / 0.08734 | 0.4716 | 0.8288 | 0.4351 |
| PACER0027 | STATE_MOTION | orthosteric_activation_core | 0.07336 / 0.06747 / 0.07248 | 0.05356 / 0.03858 / 0.04748 | 0.2311 | 0.1425 | -0.0571 |
| PACER0027 | STATE_MOTION | stable_core_control | 0.07814 / 0.08007 / 0.08959 | 0.04198 / 0.04442 / 0.06404 | -0.4345 | 0.0330 | -0.1861 |
| PACER0027 | STATE_MOTION | distal_control | 0.11501 / 0.13401 / 0.11812 | 0.06478 / 0.09602 / 0.07242 | -0.0298 | -0.0117 | 0.3586 |
| PACER0027 | SIGNED_DRIFT | pam_contact_consensus | 0.16141 / 0.15067 / 0.15524 | 0.04777 / 0.04815 / 0.05788 | -0.0153 | -0.0086 | 0.2165 |
| PACER0027 | SIGNED_DRIFT | orthosteric_activation_core | 0.12404 / 0.12928 / 0.12413 | 0.04856 / 0.04920 / 0.03251 | 0.1122 | 0.0954 | -0.2143 |
| PACER0027 | SIGNED_DRIFT | stable_core_control | 0.16219 / 0.15613 / 0.15984 | 0.05231 / 0.04956 / 0.05031 | 0.0330 | -0.1865 | -0.1938 |
| PACER0027 | SIGNED_DRIFT | distal_control | 0.23455 / 0.24259 / 0.23042 | 0.06539 / 0.08363 / 0.06344 | 0.0101 | 0.1516 | 0.0156 |
| PACER0073 | STATE_MOTION | pam_contact_consensus | 0.08938 / 0.09818 / 0.08691 | 0.05851 / 0.06872 / 0.05942 | 0.1854 | 0.0969 | 0.2787 |
| PACER0073 | STATE_MOTION | orthosteric_activation_core | 0.07888 / 0.06757 / 0.07333 | 0.05695 / 0.04547 / 0.05424 | 0.4438 | -0.1143 | 0.1862 |
| PACER0073 | STATE_MOTION | stable_core_control | 0.08206 / 0.07953 / 0.07662 | 0.04987 / 0.04669 / 0.04461 | -0.1053 | 0.0401 | 0.2409 |
| PACER0073 | STATE_MOTION | distal_control | 0.12031 / 0.13276 / 0.10131 | 0.06708 / 0.08706 / 0.05867 | -0.1394 | -0.2288 | 0.4533 |
| PACER0073 | SIGNED_DRIFT | pam_contact_consensus | 0.15663 / 0.14760 / 0.15013 | 0.05200 / 0.04878 / 0.05613 | -0.0003 | 0.0489 | 0.2457 |
| PACER0073 | SIGNED_DRIFT | orthosteric_activation_core | 0.12777 / 0.12484 / 0.12155 | 0.04944 / 0.04186 / 0.03911 | 0.2542 | 0.3142 | 0.2457 |
| PACER0073 | SIGNED_DRIFT | stable_core_control | 0.16344 / 0.15782 / 0.16201 | 0.04628 / 0.04680 / 0.04638 | 0.0219 | 0.0508 | 0.0751 |
| PACER0073 | SIGNED_DRIFT | distal_control | 0.22753 / 0.24599 / 0.22502 | 0.06041 / 0.08239 / 0.09589 | 0.0978 | 0.2620 | 0.3133 |

控制区域并非静默。distal的Delta_PAM幅度在每个候选/branch均可与PAM区域相当或更大，stable core的SIGNED_DRIFT幅度也常高于consensus。另一方面，0027 STATE_MOTION intrinsic方向在PAM面板明显一致，而stable/distal方向较弱，说明方向上的区域区分仍有信息。准确说法是“局部方向pattern与广泛幅度变化并存”，不能宣称严格PAM特异，也不能断言全受体统一同向变化。

## 14. What the Results Support

- PACER0010在PAM-related Delta_PAM比Delta_AGO方向更一致，activation-core intrinsic方向不稳定；支持有限模式分离，不能证明无intrinsic效应。

- PACER0027在STATE_MOTION的PAM/正交相关Delta_AGO具有较均衡局部复现；Delta_PAM在PAM相关区域有支持，SIGNED_DRIFT支持不全面；应重视intrinsic效应跟进。

- PACER0073在PAM/mutagenesis区域的candidate-alone方向比probe-background更稳定，尤其SIGNED_DRIFT；不能据此贴inactive或agonist标签。

- 固定历史表示下有非零、区域和replica依赖的动态contrast；允许描述局部probe-background/intrinsic signature、interaction contrast及replica一致或不一致响应。

- 允许以这些模式安排功能实验要回答的问题；不得把实验优先问题转化成冻结模型新评分。

## 15. What the Results Do NOT Support

- 稳健跨branch cooperative mechanism、确定的allosteric-to-orthosteric因果传递或激活链路：not supported by the current frozen prospective evaluation。

- 候选的PAM/ago-PAM/inactive药理标签、成功率、效能、效价、拮抗或“无intrinsic活性”的判断。

- 将Delta_PAM大视为PAM，将Delta_AGO大视为agonist，将Delta_INT大视为synergy，或将负余弦视为抑制。

- 把10个block当独立样本，推p值/独立block置信区间，或以只选某个replica/region形成机制结论。

- 认定候选全受体统一变化、已收敛，或跨cluster原始baseline可互换。

## 16. Experimental Implications

- 所有候选采用同一M4功能读出和匹配vehicle对照，同时测candidate alone、ACh/probe alone、candidate+ACh/probe；以浓度系列而非单点回答intrinsic与probe-background效应。

- PACER0010用于检验局部probe-dependent signature是否转化为组合功能调节，并单独检查candidate-alone效应；不是预先认定PAM。

- PACER0027优先排查candidate-alone功能响应并与组合响应并列，鉴别intrinsic agonism、ago-PAM可能性和仅表示扰动；不是给该候选功能分类。

- PACER0073检验candidate-alone driven perturbation能否对应功能以及probe是否改变该响应；不把MD方向分歧当成inactive证据。

- PAM-like modulation、intrinsic agonism、ago-PAM possibility与inactive binding均是待实验区分的目标；功能无响应本身不证明发生结合，需要独立结合证据。

- 本轮不预设实验阈值、成功率或综合优先级；如果资源分期，优先级应按待回答问题而非宣称某候选药效更好。

实验建议只规定比较逻辑，不擅自规定项目尚未实施的读出平台、剂量、阈值或实验成功标准。独立功能数据产生之前，所有功能标签均保持未定。

## 17. Limitations

- 每条10 ns和3个随机种子不足以证明慢转变收敛、完整构象覆盖或稳态功能机制。

- 10个相关时间block不等于10个独立样本；未计算p值、独立block置信区间或新判定阈值。

- 不同上下文是分别模拟的轨迹；同一replica种子和block编号不意味着同一微观事件或严格配对样本。

- 候选来自不同cluster；只能讨论各自matched contrasts，原始A/P baseline不可跨cluster拼接，候选差异也可能包含初始构象背景影响。

- 校准源为历史compound110四上下文R2，而非Stage4；冻结保证没有Stage4拟合污染，但不保证表示在新候选和cluster上功能泛化。

- STATE_MOTION同时含状态均值与波动幅度；SIGNED_DRIFT由端点差决定，可能受短时波动/抵消影响；RFF方向不是物理位移或激活方向。

- 两branch的归一化、带宽和RFF种子不同，不能直接比较跨branch幅度，也不计算两branch向量余弦。

- 区域重叠并经固定图扩散，区域信号不是独立验证，无法从区域均值确认因果传递、直接接触变化或单残基机制。

- 负余弦仅表示表示空间方向夹角大于90度；不是物理“反向运动”、拮抗或抑制活性。小正余弦只是弱方向共享，不构成预注册的稳健阈值。

- 现有MD/PACER-FKG输出不能替代功能实验、证明效能/效价、给出PAM/ago-PAM/inactive分类或预测概率。

- 源master_manifest的50 ps原始输出声明与认证的10 ps原始间隔不同；需确认执行记录。旧progress与早期dry-run/provenance快照不代表最终状态。

- 本轮只核验原始文件哈希与既有物理QC回执，未重跑逐帧坐标QC；对冻结历史状态做序列化哈希认证，未重新拟合或回放calibration。

## 18. Final Scientific Conclusion

冻结输出和实际消费的输入字节完整性通过；源manifest原始输出间隔仍有需确认的元数据差异。以实际认证时间契约为前提，PACER0010支持有限的probe/intrinsic动态模式分离，PACER0027支持较一致的局部intrinsic扰动及相应功能跟进，PACER0073在PAM-related区域更表现为candidate-alone方向复现。三者均缺少跨多数关键区域、跨两branch稳健一致的Delta_INT，因而cooperative dynamic mechanism is not supported by the current frozen prospective evaluation。后续应以四条件功能比较验证这些可检验假说，不提前赋予药理标签。

## Appendix A. Additional Frozen Contrasts

五种contrast全部审阅。以下补充Probe_effect/Probe_background_effect关键区域，其余全精度描述统计可见配套JSON。

| Candidate | Branch | Contrast | Region | B | N | cos12 | cos13 | cos23 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PACER0010 | STATE_MOTION | Probe_effect | pam_contact_consensus | 0.10807 / 0.08800 / 0.09741 | 0.08167 / 0.05623 / 0.05884 | 0.2040 | -0.4341 | 0.1702 |
| PACER0010 | STATE_MOTION | Probe_effect | orthosteric_activation_core | 0.08180 / 0.06765 / 0.06599 | 0.06443 / 0.04481 / 0.03698 | 0.0639 | 0.3175 | 0.5550 |
| PACER0010 | STATE_MOTION | Probe_effect | stable_core_control | 0.08556 / 0.08269 / 0.08104 | 0.05559 / 0.05301 / 0.04552 | 0.3513 | 0.3704 | 0.6018 |
| PACER0010 | STATE_MOTION | Probe_effect | distal_control | 0.12154 / 0.13909 / 0.14870 | 0.06187 / 0.08463 / 0.10235 | 0.3506 | -0.2477 | -0.3304 |
| PACER0010 | STATE_MOTION | Probe_background_effect | pam_contact_consensus | 0.08662 / 0.09873 / 0.08150 | 0.04747 / 0.07586 / 0.05244 | -0.0079 | -0.1253 | 0.6449 |
| PACER0010 | STATE_MOTION | Probe_background_effect | orthosteric_activation_core | 0.07554 / 0.06949 / 0.06694 | 0.04372 / 0.04707 / 0.04043 | 0.1255 | 0.0215 | 0.3069 |
| PACER0010 | STATE_MOTION | Probe_background_effect | stable_core_control | 0.08181 / 0.09652 / 0.08058 | 0.04394 / 0.06627 / 0.04230 | -0.2343 | 0.1688 | -0.0174 |
| PACER0010 | STATE_MOTION | Probe_background_effect | distal_control | 0.11919 / 0.14282 / 0.13647 | 0.05919 / 0.09206 / 0.09182 | -0.1173 | 0.2176 | -0.4677 |
| PACER0010 | SIGNED_DRIFT | Probe_effect | pam_contact_consensus | 0.15429 / 0.15554 / 0.15006 | 0.04726 / 0.04869 / 0.04742 | -0.2573 | -0.2381 | 0.2656 |
| PACER0010 | SIGNED_DRIFT | Probe_effect | orthosteric_activation_core | 0.12323 / 0.12264 / 0.12922 | 0.03750 / 0.04284 / 0.04476 | -0.0927 | 0.2215 | -0.1930 |
| PACER0010 | SIGNED_DRIFT | Probe_effect | stable_core_control | 0.16379 / 0.16239 / 0.17266 | 0.05434 / 0.05289 / 0.05599 | -0.0345 | 0.3081 | -0.2283 |
| PACER0010 | SIGNED_DRIFT | Probe_effect | distal_control | 0.22207 / 0.22948 / 0.24149 | 0.06206 / 0.07893 / 0.07908 | -0.0967 | 0.1091 | -0.1785 |
| PACER0010 | SIGNED_DRIFT | Probe_background_effect | pam_contact_consensus | 0.15793 / 0.15202 / 0.15092 | 0.04555 / 0.05524 / 0.04701 | -0.1800 | -0.0565 | 0.1393 |
| PACER0010 | SIGNED_DRIFT | Probe_background_effect | orthosteric_activation_core | 0.13350 / 0.12153 / 0.12412 | 0.04937 / 0.03851 / 0.03753 | 0.1925 | 0.0443 | 0.1583 |
| PACER0010 | SIGNED_DRIFT | Probe_background_effect | stable_core_control | 0.16296 / 0.16732 / 0.18029 | 0.04444 / 0.05219 / 0.04745 | -0.1249 | 0.0929 | -0.0200 |
| PACER0010 | SIGNED_DRIFT | Probe_background_effect | distal_control | 0.24744 / 0.23172 / 0.24143 | 0.07644 / 0.09401 / 0.07891 | -0.2758 | 0.1255 | -0.4236 |
| PACER0027 | STATE_MOTION | Probe_effect | pam_contact_consensus | 0.08961 / 0.10058 / 0.08120 | 0.05981 / 0.06793 / 0.04664 | -0.1768 | -0.0408 | -0.0338 |
| PACER0027 | STATE_MOTION | Probe_effect | orthosteric_activation_core | 0.06970 / 0.08644 / 0.06890 | 0.05032 / 0.07029 / 0.04809 | -0.0378 | 0.0413 | 0.7184 |
| PACER0027 | STATE_MOTION | Probe_effect | stable_core_control | 0.08081 / 0.07456 / 0.07602 | 0.04383 / 0.03657 / 0.04047 | -0.3587 | -0.1006 | 0.2849 |
| PACER0027 | STATE_MOTION | Probe_effect | distal_control | 0.11784 / 0.13940 / 0.12767 | 0.06249 / 0.09202 / 0.08767 | 0.2393 | 0.0500 | 0.3527 |
| PACER0027 | STATE_MOTION | Probe_background_effect | pam_contact_consensus | 0.07555 / 0.08364 / 0.07589 | 0.03099 / 0.05332 / 0.04495 | 0.1318 | 0.1808 | -0.0071 |
| PACER0027 | STATE_MOTION | Probe_background_effect | orthosteric_activation_core | 0.06787 / 0.07070 / 0.07260 | 0.04601 / 0.04776 / 0.05122 | 0.4710 | 0.2644 | 0.5543 |
| PACER0027 | STATE_MOTION | Probe_background_effect | stable_core_control | 0.07302 / 0.08244 / 0.08560 | 0.03333 / 0.04717 / 0.05817 | 0.0669 | -0.0792 | -0.3111 |
| PACER0027 | STATE_MOTION | Probe_background_effect | distal_control | 0.11951 / 0.13310 / 0.12097 | 0.05973 / 0.09038 / 0.07020 | 0.1577 | 0.2159 | 0.2024 |
| PACER0027 | SIGNED_DRIFT | Probe_effect | pam_contact_consensus | 0.14988 / 0.15991 / 0.16418 | 0.04616 / 0.04600 / 0.04553 | -0.1108 | -0.0286 | 0.0031 |
| PACER0027 | SIGNED_DRIFT | Probe_effect | orthosteric_activation_core | 0.11995 / 0.13337 / 0.12461 | 0.04125 / 0.04384 / 0.03976 | 0.1443 | 0.1146 | -0.0031 |
| PACER0027 | SIGNED_DRIFT | Probe_effect | stable_core_control | 0.18015 / 0.16091 / 0.16059 | 0.05168 / 0.04787 / 0.04457 | -0.2462 | -0.1318 | 0.0706 |
| PACER0027 | SIGNED_DRIFT | Probe_effect | distal_control | 0.23852 / 0.23944 / 0.21757 | 0.06623 / 0.07188 / 0.06413 | 0.0772 | 0.1876 | 0.0557 |
| PACER0027 | SIGNED_DRIFT | Probe_background_effect | pam_contact_consensus | 0.16291 / 0.14618 / 0.14575 | 0.04164 / 0.04834 / 0.05032 | -0.1832 | 0.0235 | -0.2063 |
| PACER0027 | SIGNED_DRIFT | Probe_background_effect | orthosteric_activation_core | 0.12770 / 0.12574 / 0.12839 | 0.04057 / 0.03537 / 0.03513 | 0.2987 | 0.1524 | 0.0014 |
| PACER0027 | SIGNED_DRIFT | Probe_background_effect | stable_core_control | 0.16649 / 0.16234 / 0.16119 | 0.05240 / 0.06407 / 0.04825 | 0.1732 | -0.1117 | 0.0309 |
| PACER0027 | SIGNED_DRIFT | Probe_background_effect | distal_control | 0.22938 / 0.23803 / 0.24115 | 0.06830 / 0.09058 / 0.07246 | -0.0580 | 0.0894 | 0.1052 |
| PACER0073 | STATE_MOTION | Probe_effect | pam_contact_consensus | 0.08124 / 0.10146 / 0.09428 | 0.03668 / 0.07637 / 0.05388 | 0.0045 | 0.2100 | 0.4403 |
| PACER0073 | STATE_MOTION | Probe_effect | orthosteric_activation_core | 0.07392 / 0.07388 / 0.08443 | 0.05050 / 0.05647 / 0.06185 | 0.6978 | -0.0517 | 0.1292 |
| PACER0073 | STATE_MOTION | Probe_effect | stable_core_control | 0.09245 / 0.08161 / 0.07401 | 0.06673 / 0.04154 / 0.03625 | 0.0262 | 0.1567 | -0.4156 |
| PACER0073 | STATE_MOTION | Probe_effect | distal_control | 0.13269 / 0.10705 / 0.11027 | 0.09157 / 0.05291 / 0.06665 | -0.3397 | 0.1003 | -0.3793 |
| PACER0073 | STATE_MOTION | Probe_background_effect | pam_contact_consensus | 0.08662 / 0.08322 / 0.07706 | 0.05825 / 0.05092 / 0.04216 | 0.0715 | -0.2494 | 0.2754 |
| PACER0073 | STATE_MOTION | Probe_background_effect | orthosteric_activation_core | 0.08208 / 0.08653 / 0.09036 | 0.05882 / 0.06971 / 0.07150 | 0.3345 | 0.4347 | 0.7330 |
| PACER0073 | STATE_MOTION | Probe_background_effect | stable_core_control | 0.08097 / 0.08336 / 0.07518 | 0.05136 / 0.05246 / 0.03598 | -0.5092 | -0.1226 | 0.1485 |
| PACER0073 | STATE_MOTION | Probe_background_effect | distal_control | 0.13393 / 0.14994 / 0.11497 | 0.09157 / 0.09661 / 0.07014 | -0.0964 | -0.0496 | -0.0948 |
| PACER0073 | SIGNED_DRIFT | Probe_effect | pam_contact_consensus | 0.16931 / 0.15751 / 0.16777 | 0.04632 / 0.05673 / 0.06425 | 0.0466 | 0.0701 | 0.4744 |
| PACER0073 | SIGNED_DRIFT | Probe_effect | orthosteric_activation_core | 0.13110 / 0.12478 / 0.12968 | 0.05358 / 0.04236 / 0.04965 | 0.3808 | 0.4132 | 0.2313 |
| PACER0073 | SIGNED_DRIFT | Probe_effect | stable_core_control | 0.16628 / 0.16196 / 0.15081 | 0.04814 / 0.04947 / 0.04453 | 0.2393 | 0.0616 | -0.0615 |
| PACER0073 | SIGNED_DRIFT | Probe_effect | distal_control | 0.22394 / 0.21703 / 0.23005 | 0.06803 / 0.07843 / 0.07901 | 0.0041 | 0.2655 | -0.3149 |
| PACER0073 | SIGNED_DRIFT | Probe_background_effect | pam_contact_consensus | 0.15228 / 0.15429 / 0.14395 | 0.05083 / 0.04873 / 0.04304 | 0.1357 | 0.1802 | 0.2270 |
| PACER0073 | SIGNED_DRIFT | Probe_background_effect | orthosteric_activation_core | 0.12818 / 0.12935 / 0.12953 | 0.04330 / 0.04175 / 0.05074 | 0.0557 | -0.3225 | 0.0929 |
| PACER0073 | SIGNED_DRIFT | Probe_background_effect | stable_core_control | 0.17238 / 0.16601 / 0.16566 | 0.05078 / 0.04901 / 0.04496 | -0.1659 | 0.0740 | -0.0581 |
| PACER0073 | SIGNED_DRIFT | Probe_background_effect | distal_control | 0.23682 / 0.25783 / 0.22952 | 0.07340 / 0.07219 / 0.06490 | -0.1555 | 0.1179 | -0.0520 |

补充contrast没有恢复普遍probe基准复现：各cluster的Probe_effect在PAM区域同样有分歧；如0010 STATE_MOTION consensus为0.2040/−0.4341/0.1702，0027为−0.1768/−0.0408/−0.0338，0073为0.0045/0.2100/0.4403。这增加了背景依赖和采样不确定性，不能解释为probe无功能。0073 Probe_background_effect的STATE_MOTION正交核心为0.3345/0.4347/0.7330，但SIGNED_DRIFT对应0.0557/−0.3225/0.0929；仍不形成跨branch机制闭环。

## Appendix B. Read-only Audit and Source Index

已清点Stage4结果目录原有147个文件；完整读取JSON和NPZ，并核验拓扑文件哈希及上游回执关系：input_audit及12个sanitized topology；Phase1 full/smoke回执和37个manifest；Phase2a full回执和72个manifest；Phase2b final结果、冻结回执、18份BLOCK_REGIONAL_VECTORS.npz；provenance两份阶段验证记录。原始36个DCD/state/progress、12个原始PDB及job_matrix的当前字节哈希与输入认证回执一致。Phase1与Phase2a数组形状、float32、有限性及manifest身份由只读验证复核；本轮另以DCDReader独立读取36条原始DCD头部，均为1000帧、约10 ps间隔；不重跑逐帧坐标QC。

相关代码已审阅：stage4_prospective_common_v01.py及三份run_stage4_prospective_phase1/phase2a/phase2b脚本，历史run_phase1_bs256.py、run_phase2_calibration.py、run_phase3_fkg.py的契约/公式与冻结源码身份；历史V02_FREEZE_MANIFEST、V02_METHOD_SPEC、V02_BLOCK_DEFINITIONS、V02_CONTRAST_DEFINITIONS、V02_SEEDS、V02_CALIBRATION_INVENTORY、两branch normalization/bandwidth/RFF资产、M4_MULTISTRUCTURE_GRAPH、G2_REGION_MAP_v02及G2_RESIDUE_MAPPING。还读取Stage4源master_manifest、builder_input_validation与production progress，识别源cadence差异和历史快照。

最终结果SHA-256：`ab224fa812588ccd89e5f988fd21d47a24f133a8c6cd25945564256f21be9a70`；Phase2a回执SHA-256：`04e95f7fa9a2e85f95ce25bcc01daa59ad84cdfaf016cecfd02fe3b0dd940f60`。复核命令为`python -B -u project/pacer_fkg_v02/run_stage4_prospective_phase2b_graph_region_v01.py --verify --source-root C:/projects/PACER_STAGE4_MD_backup`，退出码0。摘要与余弦独立核算仅消费既有NPZ，不执行graph/region重新生成。

工作区进入本轮时已存在4个脏的历史管理文件：BASELINE_TRACKED_TREE_SHA256.json、PHASE0_GATE.json、REPORT_PHASE0.md、V01_COMPONENT_INVENTORY.json，以及Stage4适配器/测试未跟踪文件。本轮未改变这些文件；未用git restore或重写冻结资产。报告目录由git忽略，git status不显示两个报告本身，因此只凭git status不能证明原始/忽略资产不变，应同时参考消费资产哈希和写入前后文件元数据检查。

输出：`C:/projects/GPCR-virtual-screening/project/results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md`；`C:/projects/GPCR-virtual-screening/project/results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_SUMMARY_v01.json`。

写入范围复核：498个既有相关文件的大小/mtime均未改变，仅新增上述两个报告；git status与本轮起始状态一致。既有消费资产另已通过SHA-256认证。结果目录为忽略路径，原有4个管理文件的脏状态与5个Stage4脚本/测试未跟踪状态均保留。
