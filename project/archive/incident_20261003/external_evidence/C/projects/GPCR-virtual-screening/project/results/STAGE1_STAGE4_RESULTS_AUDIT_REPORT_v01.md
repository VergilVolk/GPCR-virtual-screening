# Stage1–Stage4 实际结果只读审计报告 v01

审计日期：2026-10-03（Asia/Shanghai）。目标仓库：`C:/projects/GPCR-virtual-screening`。

本报告只汇总现存候选推理产物；未执行生成、过滤、训练、DrugCLIP 推理、docking、MD、PACER-FKG、阶段验证脚本或科学指标重算。审计操作限于读取代码、表格、JSON、Markdown、provenance、文本/结果文件 SHA-256 和文件元数据，并新增本报告。文中的既有 QC、测试、坐标认证和独立科学复核均为读取已有记录，不表示本次重新执行。

## 0. 阶段定义、证据范围与总结果

阶段定义来自实际运行规格，而不是用一般虚拟筛选流程推断：上游 `STAGE1_RUN_SPEC_v01.json` / `README.md` 明确 Stage1 为分子生成、Stage2 为理化/类药性过滤、Stage3 为冻结 DrugCLIP 筛选；当前仓库的 Stage3D 代码与冻结 manifest 将 Stage3 扩展至十构象结构门控、骨架多样性和 MD pose 冻结；Stage4 最终回执定义为候选四上下文 prospective MD 及冻结 PACER-FKG 评价。Stage4 内部 Phase1/2a/2b 不是项目 Stage1/2/3，历史 Stage A/B 也不是本报告的 Stage1/2。

**重要来源位置说明：**当前仓库仅保留本轮上游的 `drugclip_corrected_v01` 目录，Stage1/Stage2 的回执及原始生成 CSV 实际位于 provenance 引用的另一工作副本 `C:/projects/GPCR-virtual-screening-final`。本次只读追溯该来源，未复制或修补源文件。该副本的 Stage2 输出与当前仓库 `predock_portfolio.csv` 实测字节 SHA-256 完全一致，因而可以接入当前最终候选链。不能把只含当前 Git checkout 的环境视为已含所有上游证据。

以下路径前缀均为实际绝对目录，后文括号中的目录别名仅为减少重复：

| 别名 | 绝对目录 |
| --- | --- |
| R | `C:/projects/GPCR-virtual-screening` |
| U | `C:/projects/GPCR-virtual-screening-final/project/results/pacer_prospective_run_20261002_v01` |
| G | `C:/projects/GPCR-virtual-screening-final/project/results/generated` |
| B | `C:/projects/GPCR-virtual-screening/project/results/pacer_candidates_v01` |
| D | `C:/projects/GPCR-virtual-screening/project/results/pacer_prospective_run_20261002_v01/drugclip_corrected_v01` |
| K | `C:/projects/GPCR-virtual-screening/project/results/m4_gamd_ensemble/pacer_stage3_200_docking_stage3_v01` |
| S | `C:/projects/GPCR-virtual-screening/project/results/pacer_stage3d_corrected_20261002_v01` |
| M | `C:/projects/GPCR-virtual-screening/project/results/pacer_stage4_prospective_10ns_v01` |
| F | `C:/projects/GPCR-virtual-screening/project/results/pacer_stage4_prospective_fkg_v02_v01` |
| RAW | `C:/projects/PACER_STAGE4_MD_backup` |

| 阶段 | 实际输入 | 实际完成/留下 | 淘汰或未晋级 | 实际交接 |
| --- | --- | --- | --- | --- |
| Stage1 | 12 个冻结已知分子种子；17,146 次 join 尝试 | 5,000 个去重生成结构中 2,605 个通过生成器规则 | 2,395 个生成器规则不通过 | 2,605 行生成 CSV 进入 Stage2 |
| Stage2 | 2,605 个本轮生成分子 | 1,646 个满足资格；配额留下 200 个 | 合计 2,405；其中配额未入选 1,446 | `PACER0001`–`PACER0200`，160 local + 40 exploratory |
| Stage3（最终修正链） | 同一 200 个 canonical 分子 | 两模型各 200 个有效分数；2,000/2,000 docking 任务成功；200 个结构门通过；140 个不同骨架代表 | 结构淘汰 0；骨架代表压缩 60；140 代表中 137 未入三分子 MD 配额 | PACER0010、PACER0073、PACER0027 |
| Stage4 | 冻结三分子、各自 cluster-matched A/P/C/CP | 12 systems、36/36 条 10 ns 轨迹被认证；360 ns；三分子两 branch 评价全部冻结 | 未据动态结果剔除分子或 replica；不存在冻结功能合格率 | 三分子计算证据包；功能实验尚未开展 |

证据优先级为本轮冻结结果、manifest、receipt 与 provenance；日期更早的总览、准备状态、历史排名不得覆盖最终状态。整个链不是“2,605 → 五个双优 → 最终三分子”：早期五分子交集与最终 M4-safe 路由是不同结果版本，详见第 6 节。

## 1. Stage1 — 本轮 prospective 分子生成

### 输入

运行标识 `PACER-PROSPECTIVE-20261002`。基线为 `integration/final-freeze-v01`、commit `ddc11f970f508d89a75dac8ca30c41cfcd361731`、tag `pacer-code-freeze-v01`。

种子来自冻结 `CHRM4_PAM_Modeling_Handoff_v1.0/data/modeling_potency_exact_calcium.csv`（665 条测量、430 个不同分子），按分子取 max(pEC50)、按效力降序及 ID 排序选 top12，再映射 `compounds.csv` 的结构。实际 12 个种子为 CM00482、CM00681、CM00366、CM00825、CM00353、CM00123、CM00240、CM00478、CM00202、CM00651、CM00243、CM00758；没有使用六分子 fallback。这些已知分子是生成输入，不是本轮发现的候选。

### 处理过程

实际实现为 `R/project/scripts/fragment_generation.py`：BRICS 拆分、受控重组和生成器内药效团/理化过滤。冻结上限为 5 fragments、5,000 个生成结构、200,000 次 join 尝试；实际权威运行固定 `PYTHONHASHSEED=0`。生成器先用 canonical SMILES set 汇集结构，随后筛选。

### 关键数量

| 项目 | 结果 |
| --- | ---: |
| 输入种子 / 成功重建 | 12 / 12 |
| join 尝试 | 17,146 |
| 到达生成上限的去重 raw structures | 5,000 |
| 通过生成器规则、实际输出 | 2,605 |
| 未通过生成器药效团/类药性规则 | 2,395 |
| 输出 CSV 无效 SMILES | 0 |
| 输出 raw / canonical 重复 | 0 / 0 |
| 输出规则违反 | 0 |

回执字段 `invalid=2395` 指生成器规则不通过，并非交付 CSV 中 2,395 个 SMILES 解析失败。生成过程重复 join 已由 set 吸收，原脚本没有输出它们的精确次数；不能把交付重复数 0 外推为尝试过程中从未重复。

### 主要结果

状态 `STAGE1_GENERATION_COMPLETE`，QC PASS。实际权威库是 2,605 个 fragment 分子，文件 SHA-256 为 `a832c9cb2cfbed9b2e3db7902b48db5c384980cc70c551961c87f57db9e3be5d`。两次 pinned 运行字节一致是已有回执记录，本审计没有生成新库。

### 晋级到下一阶段的结果

完整 2,605 行 `G/generated_pam_analogs.csv` 进入 Stage2。此时未分配 PACER ID；三个最终候选的 exact SMILES 均在该 CSV 中出现一次，数据行分别为 0010→1756、0027→124、0073→95（不计表头）。该行号仅为现存表格定位，不是生成时间或母体编号。

### 对应核心文件路径

- `U/STAGE1_RUN_SPEC_v01.json`
- `U/STAGE1_SEED_DERIVATION_RECEIPT_v01.json`
- `U/STAGE1_GENERATION_RECEIPT_v01.json`
- `U/STAGE1_QC_v01.json`
- `G/generated_pam_analogs.csv`
- `C:/projects/GPCR-virtual-screening-final/project/data/generated/actives.smi`
- `R/project/scripts/fragment_generation.py`

## 2. Stage2 — 理化、类药性、适用域与配额过滤

### 输入

仅使用上述本次 Stage1 的 2,605 行 CSV，输入 SHA 与 Stage1 输出一致。不是历史 3,271 个混合生成分子，也没有加入 LSTM/GPT/diffusion 的历史库。

### 处理过程

实际实现 `R/project/scripts/build_candidate_portfolio.py`：标准化及 exact-known 排除、PAINS/类药性、与既有结构的相似度/近邻适用域、strict-inactive risk 门控，然后 local/exploratory 配额。

冻结资格：MW 250–550、logP −1–5.5、TPSA≤140、rotatable bonds≤10、PAINS=0；local 为最大 Tanimoto 0.55–0.85 且至少 3 个 Tanimoto≥0.45 近邻；exploratory 为最大相似度 0.35–<0.55；risk<0.45。local 配额 160、exploratory 配额 40。局部 kNN 效力和全数据 LGBM 字段是参考证据，不是实验效力；risk 也不是经本轮验证的非活性概率。

代码含既有参考模型拟合路径，属于当时生产步骤；本审计只读该代码及其已有输出，没有训练任何模型，也不将这些参考模型的训练/benchmark 分数列为本轮阶段成果。

### 关键数量

| 项目 | 数量 |
| --- | ---: |
| 输入 | 2,605 |
| exact known 排除 | 13 |
| 去除 exact known 后进入属性审计 | 2,592 |
| 审计域分类：local / exploratory / reject_domain | 415 / 1,246 / 931 |
| 配额前 eligible | 1,646 |
| 最终保留 / 唯一分子 | 200 / 200 |
| 最终 local / exploratory | 160 / 40 |
| 总未晋级 | 2,405 |
| ID 重复 / 独立规则违反 | 0 / 0 |

互斥首要原因按回执顺序统计，不能把重叠规则命中相加：

| 首要淘汰原因 | 数量 |
| --- | ---: |
| exact known | 13 |
| outside AD | 931 |
| druglike fail | 0 |
| strict-inactive risk | 15 |
| portfolio quota cap | 1,446 |
| 合计 | 2,405 |

单规则重叠命中另为 risk≥0.45 共 55，PAINS 0、druglike composite 0；55 与互斥 risk 15 口径不同，并非文件冲突。1,446 个是符合资格但未入配额，不能称为化学失败。

### 主要结果

状态 `STAGE2_FILTERING_COMPLETE`，QC PASS；留下 200 个 fragment 候选。`predock_portfolio.csv` SHA 为 `0dd7689860879eaae8abe806699a05a1c261a439e7355b99628e982abe46f52b`，上游工作副本与当前仓库文件实测相同。

### 晋级到下一阶段的结果

`PACER0001`–`PACER0200` 全量 200 个进入最终 Stage3 的评分与结构门控，完整 ID—SMILES 清单以 `B/predock_portfolio.csv` 为准。PACER0010、PACER0027、PACER0073 均属于 local，均为 PAINS=0、druglike=1、近邻数 3。

### 对应核心文件路径

- `U/STAGE2_FILTERING_RECEIPT_v01.json`
- `U/STAGE2_QC_v01.json`
- `C:/projects/GPCR-virtual-screening-final/project/results/pacer_candidates_v01/all_generated_audit.csv`（2,592 行）
- `C:/projects/GPCR-virtual-screening-final/project/results/pacer_candidates_v01/audit.json`
- `B/predock_portfolio.csv`（当前最终链 canonical input）
- `R/project/scripts/build_candidate_portfolio.py`

## 3. Stage3 — 修正 DrugCLIP、十构象结构门控与最终冻结选择

### 输入

同一 Stage2 的 200 个 canonical 分子。修正 DrugCLIP、docking 与 Stage3D 冻结 manifest 共同引用上述 `0dd768…e52b` 输入。冻结身份检查采用完整 isomeric graph；docking 微状态经显式质子化 neutralization 后映射 parent，不使用截断 InChIKey 或只按复用 PACER ID 拼接。

### 处理过程

实际最终链包含两个证据分支，再在 Stage3D 合并：

1. 修正后的 DrugCLIP 2023 backbone + M4 held-out GPCR-LOTO 三种子投影，使用本次分子表示和冻结口袋表示，得到 200 个分数；三种子 projected cosine 均值作为 2023 分数。Science2026 则独立使用对应的 Science2026 backbone 和 family-aug 三种子 head，对本次 7TRQ/7TRP/7TRS 三口袋先逐种子 max cosine，再取三种子均值。
2. M4-safe router 在 M4 上只使用 2023 percentile 排序；2026 仅为第二意见。配置文件的 0.5/0.5 权重只用于非 M4，不构成本轮 M4 融合分数。
3. 200 个候选在十个 M4R GaMD cluster（0–9）产生结构证据：Vina 1.2.7、exhaustiveness=8、num_modes=9、30 Å box；pH7 状态由 molscrub 0.1.1 提供。实际每个候选只有一个有效微状态，因此共 2,000 个状态×cluster 任务；未形成可报告的多微状态敏感性结果。
4. 冻结结构门：cluster coverage≥0.8、median native pocket coverage≥0.5、median contacts≥8、max centroid outlier≤20 Å。门控只能剔除，不重排 router。
5. 按 M4-safe final_rank 每个 Murcko scaffold 取一个代表，再取前 3 个。pose 选择为距候选 ensemble-median centroid 最近，其次口袋覆盖、再 Vina；并非直接取最低 Vina。

### 关键数量

| 项目 | 结果 |
| --- | ---: |
| 2023 / 2026 成功计分 | 200 / 200 |
| 每份评分及路由表 missing / extra / duplicates | 0 / 0 / 0 |
| canonical、DrugCLIP、docking 分子图交集 | 三组两两均 200 |
| 2023 分数范围 | −0.32773045 至 0.09303563 |
| 2026 分数范围 | −0.2368780822 至 0.1724239141 |
| docking expected / completed / successful / failed | 2,000 / 2,000 / 2,000 / 0 |
| ledger / statewise / framewise 行数 | 各 2,000 |
| docking pose 文件 / ligand 文件 | 2,000 / 200 |
| 10/10 cluster 全覆盖候选 | 200 |
| 结构门 PASS / FAIL | 200 / 0 |
| eligible distinct scaffolds / diverse shortlist | 140 / 140 |
| 骨架代表压缩未保留 | 60 |
| 最终 MD 配额 | 3 |
| 200 个中未进入 MD | 197（没有因此获得 inactive 标签） |

### 主要结果

`D/COMPLETION_RECEIPT.json` 为 COMPLETE、`D/FINAL_VALIDATION.json` 为 PASS；修正 Stage3D validation 为 PASS，无 historical DrugCLIP contamination。冻结选择保留原始排序和阈值，不生成新功能分数。

| 分子 | 2023 分数 | 2026 第二意见分数 | router final_rank | diverse_rank | 冻结 cluster / state | median pocket coverage | median contacts |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| PACER0010 | 0.09303563 | −0.019286387 | 1 | 1 | 9 / 0 | 0.775210 | 14 |
| PACER0073 | 0.07242324 | −0.12625457 | 2 | 2 | 0 / 0 | 0.723810 | 13 |
| PACER0027 | 0.05368191 | 0.0563544 | 4 | 3 | 4 / 0 | 0.678571 | 14 |

它们是不同骨架的结合优先候选，不是经功能验证的 PAM。2026 分数较低不会在冻结 M4 路由中自动淘汰，尤其 PACER0073 的晋级不能解释成两模型同时高排名。`stage3d_selected_candidates.csv` 的 `stage3d_rank` 对 PACER0027 保存为 4，而 diverse shortlist/selection manifest 的 `diverse_rank` 为 3：前者沿用 router 排位，后者是不同骨架代表序号，不能混用。

### 晋级到下一阶段的结果

按实际选择顺序：**PACER0010 → PACER0073 → PACER0027**。进入 Stage4 的唯一三分子清单由 `S/stage3d_selection_manifest.json`、`S/stage3d_selected_candidates.csv`、`S/md_shortlist_pose_manifest.csv` 共同冻结。140 个代表的完整清单保存在 `S/pacer_stage3d_diverse_shortlist.csv`，其余 137 个未进入本轮 MD。

### 对应核心文件路径

- `D/COMPLETION_RECEIPT.json`、`D/FINAL_VALIDATION.json`
- `D/2023_provenance_audit.json`、`D/2026_provenance_audit.json`、`D/router_provenance_audit.json`
- `D/prospective200_2023_gpcr_loto_scores.csv`、`D/prospective200_2026_family_aug_scores.csv`
- `D/prospective200_two_model_scores.csv`、`D/prospective200_m4_safe_routed_ranking.csv`
- `K/metadata.json`、`K/ledger.jsonl`、`K/statewise_scores.csv`、`K/framewise_scores.csv`、`K/poses/`
- `R/project/results/pacer_stage3_v01/stage3a_docking_qc.json`（其 docking SHA 与修正 manifest 一致，仅复用该 docking QC）
- `S/pacer200_structural_gate.csv`、`S/pacer200_stage3d_merged.csv`、`S/pacer_stage3d_diverse_shortlist.csv`
- `S/stage3d_selection_manifest.json`、`S/stage3d_validation.json`、`S/stage3d_provenance_audit.json`、`S/stage3d_stdout.log`
- `S/run_corrected_stage3d.py`、`S/identity_contract.py`、`S/stage3d_c_gate.py`、`S/stage3d_def_merge_diversity_pose.py`
- `R/project/scripts/pacer_drugclip_router.py`

## 4. Stage4 — 三候选 prospective MD 与冻结 PACER-FKG

### 输入

输入为 Stage3D 三个冻结分子及 concrete pose。各自使用自身受体 cluster 的 A/P/C/CP 四上下文，A=apo、P=ACh-only、C=candidate alone、CP=candidate+ACh；没有借用其他 cluster 或历史 7TRS 原始 A/P 对照。源 master manifest 明确 apo/probe-only reuse 为 NO。

| 分子 | cluster | A | P | C | CP | 轨迹 / 总时长 |
| --- | ---: | --- | --- | --- | --- | --- |
| PACER0010 | 9 | cluster9__apo | cluster9__probe_only | PACER0010__candidate_no_probe | PACER0010__candidate_probe | 12 / 120 ns |
| PACER0027 | 4 | cluster4__apo | cluster4__probe_only | PACER0027__candidate_no_probe | PACER0027__candidate_probe | 12 / 120 ns |
| PACER0073 | 0 | cluster0__apo | cluster0__probe_only | PACER0073__candidate_no_probe | PACER0073__candidate_probe | 12 / 120 ns |

### 处理过程

既有生产方案为 POPC 膜、300 K、1 bar、0.15 M、2 fs、每轨迹 5,000,000 steps / 10 ns，配对 seeds 27101/38201/49301；生产约束强度为 0。原始数据以 `RAW/systems/<system>/minimized.pdb` 与 `RAW/production/<system>/replica_01..03/{trajectory.dcd,state.csv,progress.json}` 保存，精确路径和 SHA 以输入认证回执为准。

已完成的冻结应用顺序：Phase1 C1-BS256 缓存 → Phase2a 固定 block/normalization/RFF → Phase2b 固定 graph diffusion/region/contrast → 科学解释。Stage4 R1/R2/R3 全部 EVALUATION_ONLY；历史 long-MD compound110 四上下文 R2 的 200 个 calibration blocks 是既有数值资产来源，不是本轮训练成果；Stage4 R2 未参与重新校准。

实际原始间隔 10 ps、每条 1,000 帧；`frames[::5]` 留 200 个 50 ps 分析帧；20 帧/block，10 个相关 block/trajectory，名义 1 ns、内部跨度 950 ps。冻结受体顺序投影保留 270 残基、2,139 ATOM14 重原子，处理 OpenMM A/B 两片段。STATE_MOTION 为状态均值和内部差分 RMS 的 512 维拼接；SIGNED_DRIFT 为内部平均差分 256 维；分别应用历史冻结 RFF 到 512 维，graph alpha=0.65、20 步 restart diffusion、九区域均值。

冻结 contrasts 为 `Delta_AGO=C−A`、`Delta_PAM=CP−P`、`Delta_INT=CP−C−P+A`，另有 Probe_effect=P−A 和 Probe_background_effect=CP−C。名称是方法术语：分别描述 candidate-alone 扰动、probe 背景变化、表示空间二阶非加性交互，不直接对应药理标签。

### 关键数量

| 层级 | 实际完成结果 |
| --- | --- |
| MD 输入认证 | `36_OF_36_ACCEPTED`；12 systems × 3 replicas |
| MD 总量 | 36 × 10 ns = 360 ns；每条 1,000 raw frames |
| progress | 35 complete；1 个已认证陈旧 progress 特例 |
| Phase1 | `PHASE1_FULL_FROZEN`；36 full caches，float32[1000,270,256] |
| Phase2a | `PHASE2A_FULL_FROZEN`；36 jobs，72 manifests（两 branch） |
| Phase2b | `STAGE4_FULL_FROZEN_EVALUATION_COMPLETE`；3 candidates × 2 branches；18 regional-vector NPZ |
| 动态淘汰/剔除 | 未按结果剔除任何候选或 R1/R2/R3 |
| 功能实验 | NOT YET PERFORMED；没有已确认 PAM/ago-PAM 候选 |

特殊 progress 是 PACER0010 / CP / R3，显示 running、9.9 ns。认证回执理由为 `KNOWN_STALE_PROGRESS_EXCEPTION_PACER0010_CP_R3_PHYSICAL_COMPLETE`，依据既有完整 1,000 帧、最终 5,000,000 steps、有限坐标和拓扑兼容性认定完成；不据 progress 将它记作失败，也不推广此特例。

### 主要结果

三候选均已得到完整四上下文、多 replica、两 branch 的 prospective 动态差分证据包。最终科学解释如下：

| 分子 | 已有冻结科学结果 | 结论边界 |
| --- | --- | --- |
| PACER0010 | PAM-related 区域中 probe-background 方向较 intrinsic 更一致；activation-core intrinsic 不稳定，支持有限的 probe/intrinsic 动态模式分离 | 不能证明无 intrinsic 功能效应，不能定为 PAM |
| PACER0027 | STATE_MOTION 在 PAM/正交相关区域有较一致 candidate-alone 扰动；probe-background 局部有支持，SIGNED_DRIFT 支持不全面 | candidate-alone 功能响应是实验问题，不能定为 ago-PAM/agonist |
| PACER0073 | PAM/mutagenesis 区域的 candidate-alone 方向比 probe-background 更稳定，尤其 SIGNED_DRIFT；probe 影响依赖区域和 branch | 不能据此定为 inactive 或 agonist |

下表直接转录现存科学摘要的 `pam_contact_consensus` 配对 direction cosine，顺序固定 R1/R2、R1/R3、R2/R3，四舍五入至四位；没有重新计算向量、余弦或设阈值。该区域仅作共同数值例示，最终解释以九区域、所有 replica、两 branch 的完整结果为准。

| 分子 | Branch | Delta_AGO 三配对余弦 | Delta_PAM 三配对余弦 | Delta_INT 三配对余弦 |
| --- | --- | --- | --- | --- |
| PACER0010 | STATE_MOTION | 0.0434 / −0.4107 / 0.3363 | 0.4758 / 0.1001 / 0.3592 | 0.1043 / −0.4621 / 0.1329 |
| PACER0010 | SIGNED_DRIFT | 0.0823 / 0.1050 / 0.1141 | 0.1246 / 0.1280 / 0.1979 | 0.0189 / −0.0658 / 0.0718 |
| PACER0027 | STATE_MOTION | 0.7957 / 0.7167 / 0.7912 | 0.4716 / 0.8288 / 0.4351 | −0.1851 / −0.0912 / −0.0019 |
| PACER0027 | SIGNED_DRIFT | 0.1399 / 0.1625 / 0.1339 | −0.0153 / −0.0086 / 0.2165 | −0.1520 / −0.0079 / −0.1322 |
| PACER0073 | STATE_MOTION | 0.2757 / 0.5543 / 0.2464 | 0.1854 / 0.0969 / 0.2787 | −0.1664 / −0.0140 / 0.4230 |
| PACER0073 | SIGNED_DRIFT | 0.3692 / 0.4680 / 0.4553 | −0.0003 / 0.0489 / 0.2457 | 0.0786 / 0.1597 / 0.3178 |

**总体实际结论：**没有任何候选在多数关键区域且跨 STATE_MOTION/SIGNED_DRIFT 两 branch 呈现稳健一致的 Delta_INT；当前冻结 prospective 评价不支持确定 cooperative PAM mechanism。控制区域也有明显幅度变化，区域相互重叠且经图扩散，不能声称严格 PAM-region 特异或统一全受体激活。

严格保留原 claim boundary：**Prospective frozen PACER-FKG dynamic differential evaluation only. No PAM/ago-PAM labels, probability, efficacy or independent-block inference.** 非零 Delta_INT 不是协同证明；负余弦不是物理反向运动或抑制证据；两 branch 不直接比较幅度/跨 branch 向量余弦；10 个相关 block 不当独立样本；10 ns × 3 replicas 不足以证明慢转变收敛。

### 晋级到下一阶段的结果

仓库没有在 Stage4 后冻结新的功能合格名单或药效排名。三者均保留为带完整计算证据的待验证候选，下一科学步骤是功能实验：比较 candidate alone、ACh alone、candidate+ACh 与匹配 vehicle。0010 检验模式分离是否转为功能调节；0027 检查 candidate-alone 响应；0073 检验 candidate-alone 及 probe 背景影响。没有实验结果，不能把“计算链完成”写成“PAM 发现完成”。

### 对应核心文件路径

- `M/master_manifest.json`、`M/job_matrix.csv`；`R/project/scripts/pacer_stage4_prospective.py`、`R/project/scripts/build_pacer_stage4_prospective_membrane.py`
- `RAW/master_manifest.json`、`RAW/job_matrix.csv`、`RAW/systems/`、`RAW/production/`
- `F/input_audit/INPUT_AUTHENTICATION_v01.json`
- `F/phase1_bs256/PHASE1_FREEZE_RECEIPT_v01.json`、`F/phase1_bs256/manifests/full/`
- `F/phase2a_frozen_apply/PHASE2A_FREEZE_RECEIPT_v01.json`、`F/phase2a_frozen_apply/manifests/full/`
- `F/phase2b_graph_region/PHASE2B_FREEZE_RECEIPT_v01.json`、`F/phase2b_graph_region/STAGE4_PROSPECTIVE_RESULTS_v01.json`
- `F/STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md`、`F/STAGE4_PROSPECTIVE_SCIENTIFIC_SUMMARY_v01.json`
- `F/provenance/STAGE4_TEMPORAL_METADATA_NOTE_v01.md`
- `R/project/docs/PACER_STAGE4_FINAL_HANDOFF_v01.md`
- `R/project/pacer_fkg_v02/stage4_prospective_common_v01.py` 及 `run_stage4_prospective_phase1_bs256_v01.py`、`run_stage4_prospective_phase2a_frozen_apply_v01.py`、`run_stage4_prospective_phase2b_graph_region_v01.py`
- 历史冻结数值资产：`C:/projects/GPCR-virtual-screening-fkg-v02/project/results/pacer_fkg_v02_longmd_v01/calibration/V02_FREEZE_MANIFEST.json`

## 5. 三个最终候选的完整可证 lineage

共同链：12 个已知种子 → pinned fragment 库（2,605）→ 本轮药化/AD/risk/配额（200）→ 修正两模型评分与 M4-safe 2023-only 路由 → 原有十 cluster docking 结构证据 → 结构 PASS → 不同 Murcko 代表 → frozen pose/state/cluster → 自身 cluster A/P/C/CP × R1/R2/R3 → C1-BS256 → frozen PACER-FKG 两 branch/九区域 → prospective 计算解释。

### 5.1 结构身份与药化证据

| 字段 | PACER0010 | PACER0027 | PACER0073 |
| --- | --- | --- | --- |
| Stage1 CSV 数据行 | 1756 | 124 | 95 |
| 生成来源 | fragment | fragment | fragment |
| Stage2 域 | local | local | local |
| nearest_known_id（近邻，不是已证母体） | CM00482 | CM00123 | CM00594 |
| max Tanimoto | 0.695652 | 0.605634 | 0.582090 |
| n_neighbors_045 | 3 | 3 | 3 |
| MW | 477.981 | 480.041 | 463.901 |
| logP / TPSA | 4.58514 / 103.13 | 3.11782 / 74.72 | 3.08824 / 101.44 |
| QED | 0.395438 | 0.468204 | 0.452665 |
| PAINS / druglike | 0 / 1 | 0 / 1 | 0 / 1 |
| strict_inactive_risk_ref | 0.051188 | 0.091177 | 0.087244 |
| router rank / diverse rank | 1 / 1 | 4 / 3 | 2 / 2 |
| 最终 source cluster / state | 9 / 0 | 4 / 0 | 0 / 0 |

canonical SMILES（取冻结 selection manifest，不依赖 ID 单独对齐）：

```text
PACER0010  Cc1c(Cl)c2nnc(C)n2c2sc(C(=O)Nc3cnn(C4Cc5ccccc5C4)c3)c(N)c12
PACER0027  Cc1c(Cl)c2nncn2c2sc(N3CC(N4CC(NC5Cc6ccccc6C5)C4)C3)c(N)c12
PACER0073  Cc1c(Cl)c2nnc(C)n2c2sc(N3CC(NC(=O)c4cc(F)nc(F)c4)C3)c(N)c12
```

PACER0027 的冻结 docking/MD source state 不是上述中性 parent 文本，而是：

```text
Cc1c(Cl)c2nncn2c2sc(N3CC([NH+]4CC([NH2+]C5Cc6ccccc6C5)C4)C3)c(N)c12
```

0010/0073 的 source_state_smiles 与各自 canonical 文本相同。parent/state 差异已经由冻结身份契约处理，不能据 charged 文本差异认定换分子，也不能将中性 parent 误报成 0027 实际模拟电荷状态。

**lineage 能力边界：**Stage1 CSV 只保存最终结构及描述符，没有逐分子 BRICS fragment/具体母体组合事件 provenance。因此可证明三者来自该 12-seed 生成库，不能恢复每个候选唯一母体、精确拼接序列或随机事件。nearest_known_id 是 Stage2 相似度近邻，不等于生成母体；尤其 CM00594 不是此轮 12 seeds 之一。报告不补造这部分记录。

### 5.2 冻结 pose 与 MD 交接

| 分子 | 原 docking pose（K/poses 下） | MD 冻结 pose（S/md_shortlist_poses 下） | 两文件冻结 SHA-256 |
| --- | --- | --- | --- |
| PACER0010 | `126806d1a73a_s00_c09.pdbqt` | `PACER0010_c09_s00.pdbqt` | `836a68d0cced23dc6e345a55d43e66119c27f26b9d49803635f97aed56e1e9de` |
| PACER0027 | `b55b59821a4b_s00_c04.pdbqt` | `PACER0027_c04_s00.pdbqt` | `88cd9f7b228a25f65e82c9d82a49f65653cd67eb2abe53c3280146cc001d737e` |
| PACER0073 | `699298bb49d2_s00_c00.pdbqt` | `PACER0073_c00_s00.pdbqt` | `3c503fc44a107cca2612dc71d9b0d529577bb5447eb2ea329dbe3adc32bcc37a` |

Stage4 master manifest 引用 `S/stage3d_selected_candidates.csv` SHA `5672be147ee2d6b1cc4ece4d0b95b0e90512fe6b54dfacb88de976353a507dd3`，pose manifest SHA `8ef3c4be8ff3800543d09e3b37c11331671e961c7bd3f29cc70850b787e311a6`；本次实测均相符。每分子均沿第 4 节自己的 cluster/control namespaces 进入 12 条轨迹，没有最终动态重排或替换。

## 6. 冲突、版本替代与历史结果隔离

| 现存差异 | 审计采用的依据与处理 |
| --- | --- |
| 历史报告 fragment=2,954、mixed generated=3,271、160+40、早期24候选，与本轮生成数不同 | 本轮 Stage1 receipt/QC 为权威：5,000→2,605；Stage2 2,605→200。历史2,954与未固定 hash seed 的生成有关，不能复制为本次产出。旧200与新200的 PACER ID 可复用，必须按完整图及批次 SHA 区分。 |
| 上游最初 Stage3 dual-model intersection 五分子 | `U/STAGE3_DRUGCLIP_RECEIPT_v01.json`、`U/STAGE3_OUTPUTS/STAGE3_DUAL_MODEL_DECISION_v01.json` 记录 old_rank≤50 AND new_rank≤25，200计分、5入选：PACER0108、PACER0157、PACER0080、PACER0149、PACER0089，195未入交集。这是本轮早期已实际发生的结果版本，不是最终三候选的上游 shortlist。其 Model B 复用了 Model A 表示；后续 D 的独立 Science2026 backbone provenance、M4-safe router 与 S 冻结选择明确取代该晋级链。不能把五分子记作最终进入 MD。 |
| 历史双优 PACER0040/0154/0125/0057/0053；旧 `pacer_stage3_v01` 排名与修正链并存 | 只作历史推理/政策来源。S manifest 明确 `historical_drugclip_ranking_used=false`，旧 Stage3D 文件仅用于冻结阈值和实现规则；最终 ranking 来自 D。 |
| Stage3 的 3A 名称复用 | 上游 FINAL_RUN_REPORT 的 3A 是 asset restoration；当前 docking QC 的 Stage3A 是200候选十簇 docking。按各自文件实际定义解释，不把二者任务数合并。 |
| 上游 FINAL_RUN_REPORT 前段 Stage3 环境阻断 / 3B–3D 未开始，后段 COMPLETE | 原文是追加的历史状态；优先最终 receipt、后段 resumed result；且最终候选链进一步采用 D/S 修正结果。不得以早期 BLOCKED 判定最终失败。 |
| 当前主线2026-10-02文档称候选推理未完成，旧 PIPELINE 称 MD 未启动 | 日期早的开发状态。最终 Stage4 输入认证、Phase1/2a/2b freeze receipts、2026-10-03 final handoff 证明 prospective computational evaluation 已完成；湿实验仍未完成。 |
| master_manifest `trajectory_ps=50`，认证 raw=10 ps | 以 INPUT_AUTHENTICATION、冻结 temporal contract 和实际既有 DCD 头部检查记录为准：10 ps raw×stride5=50 ps analysis。最新 `F/provenance/STAGE4_TEMPORAL_METADATA_NOTE_v01.md` 定性为源 raw-frame 描述错误，不影响既有分析；原文件、早期科学 caveat 保持不改。没有证据证明服务器覆盖原因，不自行推测。 |
| PACER0010 CP R3 progress running/9.9 ns | 以已冻结单一特例物理完整认证为准，36/36 accepted。原 progress 保持原字节。 |
| selected CSV 的 rank4 与 manifest diverse_rank3 | PACER0027 router/final rank4、不同骨架代表序号3。以实际冻结规则和明确字段含义区分。 |
| 上游 README 的 freeze SHA 拼写与权威 receipt 不一致 | README 写成包含 `...dad698e346e11bd` 的错误串；本报告使用 Stage1 run spec/receipt 及 Stage4 freeze receipt 的完整 SHA：`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`。 |

以下成果不计入本轮 Stage1–Stage4 新候选结果：DrugCLIP 13/20靶点 AUC、LOTO/LOSO、EF/benchmark；PACER-FS/LGBM 历史效力验证；LY2119620/compound110/CM00734 的 Stage A/B、long-MD 回顾性验证；OneProt/G0–G4/encoder 消融；历史24候选 docking/15 Pareto 集。这些只能解释算法/参考资产来源，不能替代三候选本次 prospective 输出，不能作为候选功能确认。

## 7. 本次只读核对与交付范围

本次通过直接读取现存表格定位三候选生成行及 Stage2/3身份，并核对以下关键文件的当前字节 SHA-256；未调用分子标准化、模型 forward、轨迹分析或阶段脚本来再现结果。

| 关键文件 | 当前 SHA-256（与回执相符） |
| --- | --- |
| G/generated_pam_analogs.csv | `a832c9cb2cfbed9b2e3db7902b48db5c384980cc70c551961c87f57db9e3be5d` |
| 上游及 B/predock_portfolio.csv | `0dd7689860879eaae8abe806699a05a1c261a439e7355b99628e982abe46f52b` |
| D/prospective200_m4_safe_routed_ranking.csv | `bda22c9354a9ebf4a263ce92ccf1257c8b27e8e3e747404b0f4fa97c3fb34d89` |
| S/stage3d_selected_candidates.csv | `5672be147ee2d6b1cc4ece4d0b95b0e90512fe6b54dfacb88de976353a507dd3` |
| S/md_shortlist_pose_manifest.csv | `8ef3c4be8ff3800543d09e3b37c11331671e961c7bd3f29cc70850b787e311a6` |
| F/phase2b_graph_region/STAGE4_PROSPECTIVE_RESULTS_v01.json | `ab224fa812588ccd89e5f988fd21d47a24f133a8c6cd25945564256f21be9a70` |

MD 原始完整性、Phase1/2数组有限性、810组摘要/810个余弦独立复核和27个测试通过均引用已有认证/科学报告，不声称本次重新逐帧检查、重测或独立复算。未核验所有外部大文件当前 SHA；本报告的 MD/动态完整性结论依赖现有冻结认证及 provenance。这不妨碍梳理实际产物，但不等同于重新做整条计算复现。

新增文件仅为 `C:/projects/GPCR-virtual-screening/project/results/STAGE1_STAGE4_RESULTS_AUDIT_REPORT_v01.md`。未修改已有结果、代码、配置、manifest、receipt、provenance、原始 MD 或训练资产；未 commit/push。报告目录属于 Git 忽略路径，前后 `git status --short` 均没有文件变更输出，不能仅据此证明忽略产物未变。

最终写入范围核对：对 `R/project/results`、`R/project/config`、`R/project/scripts`、`R/project/pacer_fkg_v02` 的 4,106 个既有文件执行写入前后 path/length/LastWriteTimeUtc 比对，变化 0、删除 0；仅新增本报告。该元数据检查与关键消费文件 SHA 核对共同支持本次只读范围声明，未将元数据相同宣称为所有文件全量字节认证。
