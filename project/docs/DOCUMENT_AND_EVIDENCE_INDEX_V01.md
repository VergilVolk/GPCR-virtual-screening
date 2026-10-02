# PACER-M4 文档与证据索引 v01

更新时间：2026-10-02

## 使用说明

仓库积累了大量开发记录。为避免把计划、失败支线或旧结论误当成当前结果，本索引采用四类标记：

- **主线**：当前对外叙事和代码联调应采用；
- **结果**：已经运行并有输出，但须遵守文件中的数据集与主张边界；
- **协议**：预注册、运行手册或交接说明，不代表实验已经完成；
- **历史/消融**：保留用于解释方法演化和失败原因，不作为当前默认方法。

## 一、项目全貌与当前状态

| 文档 | 类型 | 用途 |
|---|---|---|
| `START_HERE_FOR_COLLABORATORS_V01.md` | 主线 | 合作者的第一阅读入口 |
| `PACER_M4_PROJECT_MAINLINE_DEVELOPMENT_AND_INFERENCE_V01.md` | 主线 | 区分算法开发与候选推理，说明完整数据链 |
| `PACER_M4_INTEGRATED_ALGORITHM_AND_BENCHMARK_V02.md` | 主线/结果 | 当前核心结果、基线与严格主张边界 |
| `PROJECT_INTEGRATION_AND_BENCHMARK_V02.md` | 主线 | 模块联调、benchmark 矩阵与缺口 |
| `RELEASE_VALIDATION_V01.md` | 主线/审计 | 安装、测试、wheel 构建和 PACER-200 端到端复算 |
| `LOCAL_API_AND_DOCKER_V01.md` | 接口 | 统一 HTTP 接口、安全执行开关与容器边界 |
| `../PROJECT_STATUS.md` | 主线/日志 | 按日期记录最新冻结状态与结果 |
| `PACER_M4_RESEARCH_CHARTER.md` | 历史/约束 | 初始科学问题、验证原则和禁止事项 |
| `PACER_M4_ABSTRACT_EVIDENCE_BANK.md` | 写作材料 | 摘要可用数字与禁止表述；部分内容需用最新结果替换 |
| `../paper/PAPER_SKELETON.md` | 写作材料 | 论文结构草案 |

## 二、数据、效力建模与 PACER-FS

| 文档 | 类型 | 用途 |
|---|---|---|
| `PACER_FS_METHOD.md` | 结果/方法 | 系列中心化、三锚点适配、拒绝规则和复现协议 |
| `PACER_FS_V02_PREREGISTRATION.md` | 协议 | PACER-FS 后续版预注册 |
| `PACER_M4_ENRICHMENT_ROUTE_AUDIT.md` | 结果 | 大型 active/decoy 与真实功能标签的任务差异审计 |
| `PACER_M4_FUNCTIONAL_VALIDATION_PLAN.md` | 协议 | 功能性 PAM 验证设计 |
| `M4_RESIDUE_POTENCY_ASSOCIATION_PROTOCOL.md` | 协议 | 残基特征与效力关联的防泄漏设计 |
| `PACER_SEQUENTIAL_SAR_PREREGISTRATION.md` | 协议 | 序贯 SAR 与新系列验证 |
| `PACER_METRIC_PREREGISTRATION.md` | 协议 | 指标、统计与判定规则 |

## 三、DrugCLIP 结合检索与路由

### 当前应读

| 文档 | 类型 | 用途 |
|---|---|---|
| `PACER_DRUGCLIP_ROUTER_INTERFACE_V01.md` | 主线/接口 | 2023 GPCR-LOTO 与 2026 family-aug 的冻结路由和代码接口 |
| `DRUGCLIP_GENERATION_BAKEOFF_PROTOCOL_V01.md` | 协议 | 多 checkpoint、同数据同口袋的统一盲测协议 |
| `DRUGCLIP_GENERATION_BAKEOFF_RESULT_V01.md` | 结果 | 4-GPCR、28,618 对统一跑分与配对置信区间 |
| `DRUGCLIP_FINETUNING_PIPELINE_INTEGRATION.md` | 主线/交接 | 微调权重、推理数据和全流程接入方式 |
| `DRUGCLIP_SCIENCE2026_FINETUNING_HEADLINE.md` | 结果 | Science-2026 表征上的主要微调结果 |
| `DRUGCLIP_SCIENCE2026_13TARGET_MIGRATION_RESULT.md` | 结果 | 13-target family-aug 结果 |
| `DRUGCLIP_PROVENANCE_AUDIT_20260926.md` | 审计 | checkpoint、数据和输出来源核验 |

### 补充与历史支线

| 文档 | 类型 | 用途 |
|---|---|---|
| `DRUGCLIP_GPCR_SYSTEMATIC_FINETUNING_RESULT.md` | 结果 | GPCR 系统微调早期结果 |
| `DRUGCLIP_TARGET_TRANSFER_ABLATION.md` | 结果/消融 | 靶点迁移消融 |
| `DRUGCLIP_MUSCARINIC_TRIPLET_RESULT.md` | 结果/消融 | muscarinic triplet 试验 |
| `DRUGCLIP_FUNCTIONAL_PROBE_RESULT.md` | 结果/消融 | DrugCLIP 对功能标签的能力边界 |
| `DRUGCLIP_CGDA_RESULT.md` | 历史/关闭支线 | CGDA 结果，已不作为结合线头条 |
| `PACER_CGM_SCIENCE2026_PILOT_20260927.md` | 历史/试验 | 早期 CGM 试验 |
| `DRUGCLIP_CPU_HANDOFF.md`、`DRUGCLIP_LARGE_FINETUNE_RUNBOOK.md` | 协议 | 环境与运行交接 |
| `DRUGCLIP_LOCKBOX_ACQUISITION_V01.md` | 协议 | 外部 lockbox 获取边界 |

## 四、静态、多状态、GaMD 与六路 docking

| 文档 | 类型 | 用途 |
|---|---|---|
| `THOMPSON_MIAO_2026_EXACT_REPRODUCTION.md` | 主线/结果 | 论文六路分数、公共 GaMD 聚类和本地 Vina 的精确复现审计 |
| `PACER_XR_CASCADE_METHOD.md` | 主线/结果 | 六路 rank 共识与 1% 级联算法 |
| `M4_GAMD_ENSEMBLE_REPLICATION_PROTOCOL.md` | 协议/复现 | GaMD ensemble docking 复现步骤 |
| `DOCKING_METHODS_SURVEY.md` | 调研/历史 | docking 工具与受限环境可行性；部分早期结论已被后续复现更新 |
| `PACER_STATE_MOE_METHOD.md` | 结果/消融 | 三状态结构 mixture-of-experts |
| `PACER_STATEMETRIC_METHOD_AND_FEASIBILITY.md` | 结果/消融 | 受体状态特征可行性 |
| `PACER_M4_DYNAMIC_VALIDATION_PROTOCOL.md` | 协议 | 动态结构验证规则 |
| `PACER_GENERATION4_STEREO_ENSEMBLE_PROTOCOL.md` | 结果/压力测试 | 立体异构体与 GaMD ensemble 的机制压力测试 |

六路固定为：Glide PDB、Glide BEmin、Glide BEavg、Vina PDB、Vina BEmin、Vina BEavg。任何缺路结果不得用其他分数填补。

## 五、四上下文 MD 与 PACER-FKG

### 当前主方法

| 文档 | 类型 | 用途 |
|---|---|---|
| `PACER_FACTORIAL_KERNEL_GRAPH_METHOD.md` | 主线/方法 | A、P、C、CP 四上下文，核分布交互与残基图定义 |
| `PACER_FKG_IMPLEMENTATION_AUDIT_20260927.md` | 主线/审计 | 实现正确性、修正和测试证据 |
| `PACER_FKG_LONG_MD_FINAL_REPORT_v01.md` | 结果 | 600 ns 描述性 long-MD 结果与边界 |
| `PACER_DC_CLOSE_LOOP_CONTRACT_v01.md` | 协议 | 20 ns 闭环的预冻结规则 |
| `PACER_DC_CLOSE_LOOP_20NS_FINAL_REPORT_v01.md` | 主线/结果 | LY2119620 与 CM00734 hard-negative 闭环结果 |
| `PACER_DC_CLOSE_LOOP_20NS_HANDOFF_v01.md` | 交接 | 闭环数据和运行接口 |
| `PACER_DC_PUBLICATION_PROTOCOL.md` | 协议 | 统计单位、baseline、消融和发表门槛 |

### 编码器探索与失败消融

| 文档 | 类型 | 用途 |
|---|---|---|
| `PACER_DC_DYNAMIC_ENCODER_DECISION.md` | 结果/决策 | PCA、tICA、VAMP、OneProt 等动态表示选择 |
| `ONEPROT_PACER_DC_BRANCH_REVIEW_20260921.md` | 历史/审计 | OneProt 分支代码与证据审核 |
| `PACER_DC_ENCODER_DIAGNOSTIC_AND_V2.md` | 历史/诊断 | OneProt 表示不稳定的诊断 |
| `PACER_DC_GEOM2VEC_HANDOFF.md` | 历史/交接 | Geom2Vec 试验接口；最终采用其 C1-BS256 中间 readout 思路，而非原 final embedding |
| `PACER_MCV_METHOD_AND_HANDOFF.md` | 基线 | 不依赖预训练编码器的可解释物理特征基线 |
| `PACER_STATIC_ANCHORED_TRAJECTORY_TEST.md` | 消融 | 静态锚定轨迹测试 |

### 机制与后续验证

| 文档 | 类型 | 用途 |
|---|---|---|
| `M4_ALLOSTERIC_COUPLING_HYPOTHESES.md` | 假设 | 变构耦合区域和可检验机制 |
| `PACER_M4_FUNCTIONAL_VALIDATION_PLAN.md` | 协议 | 如何用药理实验确认 PAM、ago-PAM 与 inactive |
| `PACER_PUBLIC_PAM_DYNAMIC_ENDPOINT_PREREGISTRATION.md` | 协议 | 公共 PAM 动态端点预注册 |
| `PACER_PAIRED_DELTA_PAM_SIGNATURE_PREREGISTRATION.md` | 协议 | 配对 Delta_PAM 验证 |
| `PACER_REPLICA_INVARIANT_STATE_PREREGISTRATION.md` | 协议 | 跨 replica 稳定性验证 |

## 六、分子生成、候选选择与实验假设

| 文档 | 类型 | 用途 |
|---|---|---|
| `PACER_GENERATION4_STEREO_ENSEMBLE_PROTOCOL.md` | 协议/结果 | 立体化学候选与 ensemble 压力测试 |
| `PACER_M4_LEAD_PAM_HYPOTHESES.md` | 主线/候选 | 首轮候选证据卡和预设湿实验成功标准 |
| `PACER_M4_ABSTRACT_EVIDENCE_BANK.md` | 写作材料 | 生成路线数量、适用域结果和摘要素材 |
| `../PACER_M4_FINAL_REPORT.md` | 历史结果 | 2026-08-30 前的数据、生成、PACER-FS 和候选闭环；不是最新总报告 |

## 七、外部验证与预注册

以下文件用于冻结外部评价规则，文件存在不等于评价已经完成：

- `PACER_ACADIA_2025_EXTERNAL_PREREGISTRATION.md`
- `PACER_SUVEN_PATENT_EXTERNAL_PREREGISTRATION.md`
- `PACER_US20260055116_EXTERNAL_PREREGISTRATION.md`
- `PACER_EXTERNAL_2026_PREREGISTRATION.md`
- `PACER_EXTERNAL_STATIC_PAM_SIGNATURE_PREREGISTRATION.md`
- `PACER_LY2033298_ALLOSTERY_PREREGISTRATION.md`
- `PACER_ACM_V2_PREREGISTRATION.md`
- `CGDA_EXTERNAL_EVALUATION_PROTOCOL.md`

对外汇报前必须到 `PROJECT_STATUS.md` 或相应 `RESULT`/`FINAL_REPORT` 文件确认是否已经产生结果。

## 八、当前不应作为入口的文件

- `../PIPELINE.md`：仍保留“OneProt-MD 是唯一主线”等历史表述，不能代表 2026-10-02 状态。
- `../PACER_M4_FINAL_REPORT.md`：是 8 月阶段性冻结报告，不含最新 DrugCLIP 路由和 PACER-FKG hard-negative 闭环。
- 所有 `PLAN`、`PREREGISTRATION`、`CONTRACT`、`RUNBOOK`、`HANDOFF`：只能说明预设或执行方式，不能自动视为结果。
- 旧版 `PACER_M4_INTEGRATED_ALGORITHM_AND_BENCHMARK_V01.md`：由 v02 取代。
- CGDA、OneProt final embedding、原始 Geom2Vec final embedding、State-MoE 等：保留为历史或消融，不再作为默认主方法。

## 九、证据追溯原则

1. 总结文件只负责解释；定量主张必须能追溯到 `results/` 中的机器可读输出。
2. 计划文件和运行成功日志不能替代 benchmark 结果。
3. docking、DrugCLIP 和 PACER-FKG 分属不同任务，不得把各自最优数字拼成单一端到端 AUC。
4. 动态轨迹窗口不是独立分子样本；统计单位优先为分子、药化系列、靶点或独立 replica。
5. 未经湿实验确认的新分子统一称为“计算优先候选”。
