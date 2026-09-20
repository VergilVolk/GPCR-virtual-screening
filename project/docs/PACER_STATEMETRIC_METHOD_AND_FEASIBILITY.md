# PACER-StateMetric / StateMIL 方法定义与可行性评估

> 状态：概念冻结、尚未训练；2026-08-30。本文记录下一代动态主线，防止其退化为 RMSD 聚类或多构象 Vina 加权。

## 1. 核心假设

传统 ensemble docking 在严格 M4 PAM/inactive 真值上接近随机，原因是它只检验候选对既有构象的结合相容性。PACER-StateMetric 改为学习：

> 哪些受体动态状态代表由 PAM 诱导的正构配体稳定、变构协同性和 G 蛋白耦合，而不是哪些 snapshot 在原子坐标上接近。

它不在固定 embedding 上事后聚类，而是用功能正负关系重塑 embedding，使功能相近状态形成可学习的软状态原型。

## 2. 模型架构

### 2.1 无标签动态编码器

每个 snapshot `x_t` 经基础编码器得到 `h_t = E0(x_t)`。输入优先采用可解释的受体状态特征：

- PAM/正构配体残基接触与水桥；
- ACh/iperoxo 位姿和正构口袋稳定性；
- DRY、NPxxY、TM3-TM6、TM7 等微开关；
- G 蛋白 alpha5 螺旋和受体胞内界面；
- 口袋体积、水网络和残基动态相关；
- 可选的局部几何图表示。

E0 先使用时间邻近、跨 replica 一致性、遮蔽重建或慢模态目标预训练；功能标签不进入此阶段。

### 2.2 功能度量适配器

小型 adapter 将基础表示与实验上下文映射为 `z_t = A(h_t, context)`。上下文至少包含 species、orthosteric probe、probe 浓度/ECxx、readout、G 蛋白/嵌合系统与 endpoint。

基础编码器默认冻结，只训练 adapter；完全解冻必须作为独立消融，防止小样本过拟合。

### 2.3 可学习状态原型与 attention

设置 K 个可学习 prototype `c_k`，snapshot 的软归属为：

`a_tk = softmax(-||z_t-c_k||^2 / temperature)`

prototype 不预先命名为 PAM/inactive。模型训练后再根据训练体系富集、特征和对照轨迹解释。轨迹输出包括 prototype occupancy、状态转移、跨 replica 方差和 attention 加权表示。

### 2.4 Multiple-Instance Learning

功能标签属于整条轨迹/分子，而非单个 snapshot。每条轨迹是 bag，snapshot 是 instance。attention 从 bag 中识别与 `log alpha`、`log alpha beta`、`log tauB` 或 PAM/inactive 标签相关的稀有状态；禁止把 PAM 轨迹中的每个 snapshot 全部标成 PAM 状态。

## 3. Triplet 的保留方式

Triplet 不是独立大模型，而是两级训练约束。

### 3.1 分子级 Functional Triplet

- Anchor：分子 + assay context；
- Positive：同 endpoint/context 且功能接近；
- Hard negative：ECFP4/理化性质接近但 potency、PAM/inactive、intrinsic agonism、M2 selectivity、species 或 readout 表型显著不同。

它替代 PACER-FS 中较普通的分子表示，但保留“系列内相对 SAR + 新系列功能锚点标定”。

### 3.2 轨迹级 State Triplet

- Anchor：已知 PAM 的一条 trajectory bag；
- Positive：另一 replica 或不同 chemotype、相同功能耦合表型的 trajectory bag；
- Hard negative：ACh-only、inactive 近邻、弱 PAM 或高 intrinsic-agonism 体系，并优先选择基础几何 embedding 接近的对照。

snapshot 级只使用时间连续、跨 replica prototype 一致性和弱对照来源，不赋予伪造的逐帧功能标签。

## 4. 必需对照体系

每个训练 chemotype 最低需要：

1. `M4 + orthosteric agonist + PAM + G protein`；
2. `M4 + orthosteric agonist + G protein`；
3. `M4 + PAM + G protein`；
4. `M4 + orthosteric agonist + inactive matched analog + G protein`；
5. 强/弱 PAM 梯度及至少 3 个独立 replica。

所有配对体系使用一致起始结构、膜、质子化、离子、力场、平衡和采样协议。

## 5. 防塌缩和防泄漏

- 先按完整 chemotype/source/document 划分，再在训练部分构造 pair/triplet；
- 测试配体及其 snapshot 不参与预训练后功能 adapter、prototype 数量或阈值选择；
- prototype 利用率和 attention 熵约束防止单簇/单帧塌缩；
- replica consistency 检查同体系状态分布可重复性；
- chemotype-adversarial 或系列留出检验模型不是在识别骨架；
- K、margin、temperature 和损失权重全部在内层训练系列选择；
- 外层按 ligand/chemotype 同时留出所有 replicas 和 endpoints。

## 6. 冻结 baseline 与消融

1. RMSD average-linkage k=10；
2. TICA/VAMP + k-means/HDBSCAN；
3. PACER-FS / ECFP-LightGBM；
4. 普通 graph encoder，无 triplet；
5. 随机 negative triplet；
6. functional hard-negative triplet；
7. prototype attention，无 MIL；
8. StateMIL，无化学分支；
9. StateMetric chemistry-only；
10. StateMetric + MD；
11. 去除 context、去除 G protein、去除 orthosteric-control 的消融。

任何 MD 增量必须比较相同查询集上的 chemistry-only 模型，并报告 paired bootstrap CI。

## 7. 成功标准

### 分子级 Go/No-Go

- leave-source/series-out macro Spearman 高于冻结 PACER-FS 0.263；
- paired series-bootstrap 增益 CI 下界 > 0；
- activity-cliff ranking 明显改善；
- worst-series 不低于 PACER-FS；
- 至少一套未用于开发的外部 campaign 支持。

### 动态级 Go/No-Go

- 在 leave-chemotype-out 中 PAM/inactive 或强/弱 PAM 的状态分布可分；
- PAM-state occupancy 与 `log alpha beta` 的方向和相关性跨 replica 稳定；
- intrinsic-agonism 状态与 `log tauB` 可分，不把两者混成单一活性；
- MD 分支相对 chemistry-only 的主要指标 paired CI 下界 > 0；
- 随机化轨迹标签和随机 prototype 对照回到机会水平。

未达到任一主条件时，不将该模块用于候选效力排序。

## 8. 当前资源审计

### 已具备

- MK-97–M4R holo-GaMD：6 条 500 ns replica；公开 1 ns stride 共 3000 帧；
- 完整 Amber `sys.parm7`、膜/水/离子体系和输入文件；
- 7TRQ 正构/别构配体结构；
- Monash LY2033298：16 分子，13 个完整 `pKB/log alpha beta/log tauB`，3 个明确 inactive；
- mechanistic multidomain：27 分子，均有 cAMP `log alpha beta/log tauB/pKB`，部分有 GoB、arrestin 和 binding endpoints；
- 430 分子 potency、66 strict inactive、93 activity cliffs，以及 Acadia/US 多终点外部数据；
- MDAnalysis 可读取现有轨迹，已完成 3000 帧几何重聚类。

### 缺失

- 与公开 PAM 轨迹严格匹配的 ACh-only、PAM-only、inactive 和弱 PAM 轨迹；
- 多 chemotype 配对轨迹和足够 replica；
- 逐分子已参数化膜蛋白体系；
- 当前 Python 3.13 环境无可用 OpenMM/ParmEd；ParmEd 构建需要 MSVC，生产环境应使用独立 Python 3.11/Conda；
- CPU 工作站不适合多体系 100–500 ns 生产模拟。

## 9. 可行性结论

| 部分 | 科学可行性 | 当前执行可行性 | 结论 |
|---|---|---|---|
| 分子 Functional Triplet | 中高 | 高，CPU可做 | 立即实现严格baseline |
| 无标签公开GaMD表示预训练 | 中高 | 高，已有3000帧 | 可立即做特征/MIL原型pilot |
| 单PAM replica一致性 | 中 | 高 | 只能验证技术稳定性，不能验证PAM特异性 |
| PAM vs ACh-only StateShift | 高 | 中低 | 需新对照轨迹和模拟环境 |
| 多chemotype leave-out StateMIL | 高 | 低 | 需GPU与批量参数化 |
| 候选功能预测 | 潜力高 | 当前不足 | 通过上述Go/No-Go后才能应用 |

总体判断：方法在科学上成立，且比 ensemble docking 加权更直接针对 PAM 功能；现有数据足以启动分子级 triplet、无标签动态表示和单体系 replica pilot，但不足以训练或验证可泛化的 PAM 状态模型。下一步必须先获得配对对照轨迹，不能用六条同一 PAM replica 冒充多 chemotype 训练集。

### 量化审计结果

按“先同系列约束、后构造triplet”的冻结定义，430个potency分子提供1084个功能近邻正对、405个同系列效力cliff困难负对；229个anchor同时具备正负伙伴，覆盖10个系列，可形成4182个原始triplet组合。463 PAM/66 experimental inactive之间另有411个Tanimoto>=0.55的困难对。故分子Functional Triplet为Go，但这些组合是相关增强而非4182个独立样本。

动态侧仍为No-Go：公开3000帧来自六条replica，但只有一个PAM chemotype，且没有匹配的ACh-only、PAM-only、inactive和弱PAM轨迹。完整机器可读审计见 `results/pacer_statemetric_feasibility_v01/audit.json`。

## 10. 结论边界

PACER-StateMetric 当前是预注册开发假设，不是完成的算法、SOTA 或候选验证证据。模型生成的是功能 embedding、soft prototype membership 和轨迹状态分布，不生成新的物理 snapshot。只有在严格留出、对照轨迹和外部/前瞻功能实验中胜出后，才能晋级为候选排序模块。

## 11. 分子 Functional Triplet v0.1 冻结结果（2026-08-30）

预注册的11个整 `source_component` 留出已完成。三锚点标签仅用于MAE偏移，不参与系列内排序；所有triplet均在外层划分后从训练系列重新构造。

| 方法 | macro Spearman | worst-series | cliff方向准确率 |
|---|---:|---:|---:|
| PACER-FS | 0.263 | -0.005 | 0.693 |
| Absolute-LightGBM | 0.199 | -0.157 | 0.725 |
| MLP-Reg | 0.138 | -0.297 | 0.586 |
| MLP-RandomTriplet | 0.142 | -0.220 | 0.601 |
| PACER-FM-HardTriplet | 0.141 | -0.341 | 0.574 |

HardTriplet相对PACER-FS的macro差值为 -0.122，series-bootstrap 95% CI `[-0.245,+0.002]`，优于PACER-FS的系列仅3/11；相对MLP-Reg的差值仅+0.003，CI跨零。内部Go条件全部未通过。

失败不是单一实现错误可以解释：普通MLP已显著弱于树模型；hard-triplet没有改善cliff方向；原始4182个triplet中44.9%来自一个系列、79.3%来自三个系列，相关增强不能替代独立chemotype和多终点功能标签。因此v0.1冻结为No-Go，不在相同外层结果上调margin、网络宽度或损失权重追分。

这不否定轨迹级StateMetric，但收窄了路线：分子potency triplet不再作为主创新；动态路线必须以匹配对照轨迹、轨迹袋级标签、prototype occupancy和chemotype留出为核心。在获得ACh-only/inactive/弱PAM/多chemotype对照前，仅允许无标签snapshot表示与replica一致性pilot，不允许声称PAM功能预测。

## 12. 无标签状态 pilot v0.1 冻结结果（2026-08-30）

在公开3000帧上比较PCA、tICA与无标签TemporalTriplet。PCA/tICA/TemporalTriplet的mean replica occupancy JSD分别为0.536/0.743/约0.737；tICA的1 ns状态持续率为0.930，但只有1个状态在六条replica均达到1%占比。TemporalTriplet三种子的平均AMI为0.623、有效状态数约9.56，没有数值塌缩，但六replica共同覆盖仅1–2个状态，技术Go失败。

该结果揭示了两个边界。第一，六条500 ns replica的状态占比差异很大，单体系的MD采样尚不能被当作已收敛的功能状态分布。第二，v0.1将其他replica随机帧设为negative，在目标上与“跨replica共同功能状态应接近”相矛盾，会鼓励replica shortcut。该版本保持冻结，不通过调margin追分。若继续，必须预注册为：同replica时间远且几何远的帧作negative、跨replica几何互近邻作弱positive、replica-invariance约束，并使用leave-one-replica-out评价，防止直接优化评估集占比分布。

## 13. Replica-invariant StateMetric v0.2 冻结结果（2026-09-01）

按上述纠正完成六折leave-one-replica-out与三个种子。ReplicaInvariantTriplet将held-out occupancy JSD由PCA的0.313降至0.285，并将有效状态数由3.79提高至4.94，说明replica shortcut有所减弱；但时间持续率由0.653降至0.530，共同状态覆盖由3.17降至2.78，三种子AMI仅0.234，远低于0.60门槛。技术Go失败，功能StateMIL继续No-Go。

该失败说明“去除replica身份”与“得到稳定功能簇”不是同一个问题。无标签triplet可重塑状态占比，却没有足够跨体系监督确定唯一、可复现的离散边界。

## 14. 三静态锚点与第三PAM轨迹检验（2026-09-01）

使用7TRQ/VU0467154与7TRP/LY2033298作为两个独立PAM正锚，7TRS/ACh-only作为负锚，盲注释第三PAM MK-97的六条公开GaMD轨迹。126个共同口袋原子的静态RMSD满足 `7TRQ-7TRP=0.785 A < 7TRQ-7TRS=0.874 A < 7TRP-7TRS=1.026 A`，证明两个PAM结构之间存在弱但方向一致的几何共性；但整体RMSD PAM-likeness仅4/6 replica中位数为正，冻结premise未通过。

为排除整体对齐与构建体差异，预注册了21个C-alpha、210个内部距离的Consensus-PAM signature。主阈值0.10 A保留141个在两个PAM相对ACh-only中同向变化的距离。MK-97六条replica中5/6中位数为正，0/0.05/0.10/0.20 A四个阈值均保持5/6；replica 2为唯一反向，故严格Go仍失败。replica 2前100 ns为正，100 ns后转入并长期停留于反向盆地，主要差异涉及Y97-D432、Y97-N423、I187-Y439、W108-L190、W98-D432等跨区距离。它们是结构投影贡献，不是因果残基或效力决定因子。

## 15. 冻结triplet cluster的外部注释（2026-09-01）

Consensus-PAM signature在ReplicaInvariantTriplet训练后才构造，因此可作回顾性外部注释。cluster对该连续签名的mean eta-squared为：PCA 0.410、tICA 0.208、ReplicaInvariantTriplet 0.321；三者相对保持时间自相关的循环移位null均为p=0.001，但triplet相对PCA差值为-0.089，且仅1/6折更好。结论是原无标签triplet捕获了真实轨迹结构，却没有比线性PCA更好地组织PAM共同方向。

## 16. Anchor-Guided Trajectory Triplet v0.3（2026-09-01）

进一步严格留一replica，用训练五条轨迹内“跨replica、相同静态锚投影”构造positive，用“整体几何相近但锚投影差至少0.50”构造hard negative；held-out replica的signature从未参与训练。AnchorGuidedTriplet的held-out eta-squared为0.695，高于同特征PCA的0.503，6/6折均提高，平均增量+0.192；有效状态数也由2.65提高至3.32。

然而三种子AMI仅0.449，时间持续率由PCA的0.635降至0.497，未通过冻结Go。最准确的结论是：**静态跨PAM共识可以引导embedding学到可跨replica迁移的连续结构轴，但现有单体系数据不足以把该轴离散成稳定、时间连贯、可称为功能状态的cluster。** 这保留了anchor-guided metric作为候选创新模块，同时否定了当前版本已经完成PAM效力预测或达到SOTA。

## 17. 外部静态特异性与probe-matched差分（2026-09-01）

冻结141距离签名在独立第三PAM 7V68/LY2119620上为1.015，高于iperoxo-only 7TRK的0.377与inactive 5DSG的0.094，三个敏感性阈值顺序一致。但同研究allosteric agonist 7V6A/compound-110为1.177，高于真实PAM；故旧轴是变构占位/激活耦合轴，不是PAM-specific。

改用相同iperoxo背景的`7TRQ-7TRK`与`7TRP-7TRK`建立paired-delta轴后，独立`7V68-7V69`仍稳定复现+0.54左右的PAM结构增量，但compound-110仍更高。结论是probe-matched结构差分可识别PAM引起的结构增量，却不能从单一上下文区分PAM与allosteric agonist；必须加入无正构探针上下文。

## 18. 公开多副本动态端点（2026-09-01）

从Wang et al. 2022官方Source Data恢复M4R-iperoxo与M4R-iperoxo-LY2119620各3条独立1 us轨迹的iperoxo RMSD。以trajectory为独立单位，LY2119620使iperoxo RMSD平均降低0.651 A，三个副本差值为+0.706/+0.467/+0.780 A；20种精确标签置换单侧p=0.05，50 ns分层block-bootstrap 95% CI `[0.527,0.762]` A，去除首100 ns后效应0.656 A。该结果支持orthosteric stabilization作为机制相关动态端点，但只有一个PAM体系。

在独立MK-97六副本中，冻结蛋白coupling coordinate与`-ACh pose RMSD`在6/6副本同向，mean within-replica Spearman 0.196，保持时间自相关的循环移位p<0.0001；encoder未读取配体坐标。连续结构轴因此获得独立动态注释支持。相反，AnchorGuidedTriplet hard cluster对ACh RMSD的eta-squared比PCA低0.027，仅2/6折改善，再次否定硬聚类作为主输出。

## 19. Soft/sticky prototype head（2026-09-01）

有序soft prototype在6/6折保持ACh稳定性单调，mean occupancy JSD 0.127，但瞬时连续性仅0.762。使用训练副本自相关时间（4--5 ns）确定sticky transition并前向-后向推断后，连续性升至0.928、JSD 0.216、max occupancy 0.542，仍保持6/6单调；effective prototype count为2.979，略低于冻结3.0门槛，因此开发Go保持false。该模块保留为确定性trajectory attention head，不称功能状态预测器。
