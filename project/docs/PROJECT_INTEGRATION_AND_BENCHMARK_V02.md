# PACER 项目联调与基准评估 v02

日期：2026-10-01  
分支：`codex/project-integration-benchmark-v02`

## 1. 项目到底在做什么

项目不是用一个总分直接宣布“发现了 PAM”，而是把问题拆成三道可审计的门：

1. **结合候选检索**：分子是否具有与 M4 口袋匹配的跨靶点结合表征；
2. **结构可行性检查**：候选能否形成合理姿势、关键接触及多构象稳定性；
3. **功能机制检查**：在 ACh 条件下，候选引起的受体动态变化是否具有可重复方向，并与已知 PAM、无功能近邻相区分。

因此，当前闭环输出的是“优先验证的 PAM 假设”，不是湿实验确认的 PAM。

```text
带靶点标签的配体数据
  -> DrugCLIP / 2D 基线与严格 target-LOSO
  -> 两种检索模型共同保留候选
  -> docking / IFP / PACER-FS 结构门控
  -> 四上下文 MD: R0, RA, RC, RCA
  -> DeltaPAM=(RCA-RA), DeltaAGO=(RC-R0), DeltaINT=RCA-RA-RC+R0
  -> C1-BS256 表征 + PACER-FKG 区域方向审计
  -> 已知 PAM / hard negative / 候选并列比较
  -> 候选、证据、拒绝原因和主张边界
```

## 2. 各模块的输入、算法、输出和验收条件

| 模块 | 输入 | 当前算法 | 输出 | 必须通过的基准或门槛 |
|---|---|---|---|---|
| 数据与协议 | 分子、靶点、活性标签、来源、骨架 | 规范化、去冲突、target/scaffold/source/series 切分 | 冻结 CSV、audit、哈希 | 训练/测试靶点和同系列泄漏检查 |
| 结合检索 | pocket 表征、分子表征 | 官方 DrugCLIP；ECFP4 logistic；ep80；family-augmented ensemble | 每个靶点的分数与排名 | 同一 13T 或 20T 协议内比较；ROC、PR、BEDROC20、EF1%、EF5% |
| 候选交接 | 同一批 PACER-200 的旧、新分数 | ID+canonical SMILES 双重核对；双模型门控 | 200 分子审计表、交集短名单 | 必须 200/200 身份一致；禁止跨库按行拼接 |
| 静态结构 | 受体构象、候选姿势 | Vina、IFP、2D/结构融合、PACER-FS | 姿势、相互作用和静态优先级 | 与 2D 基线、LOSO 和外部富集比较；静态分数不冒充 PAM 功能 |
| MD 表征 | 匹配的四上下文轨迹 | PCA/tICA/VAMP；OneProt-MD；Geom2Vec；冻结 C1-BS256 | 每个窗口、区域和 replica 的动态表征 | replica 留出、匹配 seed、公共原子核、方向一致性 |
| 功能差分 | R0/RA/RC/RCA 表征 | DeltaPAM、DeltaAGO、DeltaINT；PACER-FKG 区域方向 | 功能相关区域的方向和量级描述符 | 已知 PAM 正控、近邻无功能 hard negative、远端对照；禁止事后调参 |
| 最终决策 | 结合、结构、动态证据 | 分层门控，不做不可解释加权总分 | shortlist + evidence card + reject reason | 缺任一上下文就拒绝功能结论；无湿实验时只称候选假设 |

## 3. v02 的严格 benchmarking

### 3.1 结合检索：统一从原始预测重算

指标实现统一为 sklearn ROC/AP、RDKit `CalcBEDROC(alpha=20)` 和 ceil-based EF。13T 与 20T 是不同协议，只能分别比较。

| 协议/方法 | ROC-AUC | PR-AUC | BEDROC20 | EF1% | EF5% |
|---|---:|---:|---:|---:|---:|
| 13T 官方 DrugCLIP | 0.544 | 0.178 | 0.188 | 1.788 | 1.127 |
| 13T ECFP4 logistic | 0.568 | 0.201 | 0.239 | **2.377** | 1.708 |
| 13T ep80 ensemble | 0.608 | 0.215 | 0.228 | 1.207 | 1.323 |
| 13T family-aug ensemble | **0.641** | **0.248** | **0.294** | 2.105 | **2.058** |
| 20T 官方 DrugCLIP | 0.582 | 0.145 | 0.198 | **3.990** | **2.161** |
| 20T ECFP4 logistic | 0.533 | 0.144 | 0.184 | 1.944 | 1.470 |
| 20T ep80 ensemble | **0.628** | 0.159 | **0.214** | 1.728 | 2.109 |
| 20T family-aug ensemble | 0.625 | **0.163** | 0.213 | 2.528 | 1.999 |

成对靶点 bootstrap 的关键结论：

- 13T family-aug 相对官方 DrugCLIP：ROC `+0.097`、PR `+0.070`、BEDROC20 `+0.107`、EF5% `+0.931` 的 95% CI 均不跨零；EF1% 不确定。
- 13T family-aug 相对 ECFP4：ROC 与 PR 改善的 CI 不跨零，但 BEDROC20、EF1%、EF5% 不能宣称稳定优于 ECFP4。
- 20T family augmentation 相对 ep80 没有稳定总体增益；它是范围受限的 13T 改进，不是普遍 SOTA。

旧脚本的 ensemble BEDROC 公式有误；v02 的 RDKit 重算值取代旧值。

### 3.2 候选交接：修复了错误数据表联结

旧检查把 PACER-200 与另一套 28,519 分子库比较，得到零重叠，没有科学意义。v02 在同一 200 分子上完成：

- candidate ID 与 canonical SMILES：`200/200` 一致；
- 旧 DrugCLIP 排名与 family-aug 排名 Spearman：`-0.059`，说明两者提供互补而非重复排序；
- 规则“旧排名 <= 50 且新排名 <= 25”只得到 **5 个**候选：`PACER0040`、`PACER0154`、`PACER0125`、`PACER0057`、`PACER0053`；
- 这 5 个是双模型交集假设，不是已确认 PAM。

### 3.3 功能动态线：目前能说和不能说的

冻结端点为 `STATE_MOTION / DeltaINT / R1-R3 direction cosine`：

| 分子 | compound110 extension | PAM contact consensus | 解释 |
|---|---:|---:|---|
| compound110 | +0.584 | -0.308 | 局部区域有重复信号，但并非所有机制区一致 |
| LY2119620 已知 PAM | +0.634 | +0.375 | 多个机制区方向为正，支持正控机制一致性 |
| CM00734 无功能近邻 | -0.327 | -0.251 | 与正控方向相反，构成部分特异性证据 |

该结果支持“方向描述符可能区分已知 PAM 与近邻无功能分子”的机制假设；样本仍只有三类代表，不能称为通用 PAM 分类器、效力预测器或临床候选确认。

### 3.4 已覆盖的基线

方法注册表共 47 条，分为五层：

- 2D/QSAR：TanimotoKNN、RandomForest、ExtraTrees；
- 静态结构：Vina、IFP、结构/2D 融合、PACER-FS；
- 结合检索：官方 DrugCLIP、ECFP4、ep80、family-aug、随机标签、CGDA 外部迁移；
- MD 表征：PCA、tICA、VAMP、OneProt-MD、Geom2Vec、C1-BS256、PACER-MCV；
- 相关公开先进方法：GNINA、RTMScore、DeepRLI、EquiScore、VAMPnet、SPIB。

公开先进方法若没有同一输入、同一切分和原始预测，只在覆盖表中标为 `NOT_RUN` 或 `PILOT`，不出现在胜负柱状图中。

## 4. 图表与机器可读产物

- `fig_v02_01_binding_fair_bars.png`：13T/20T 分协议柱状图；
- `fig_v02_02_binding_delta_ci_heatmap.png`：方法差值与 target-bootstrap CI；
- `fig_v02_03_binding_target_heatmap.png`：13 个靶点逐靶 ROC 热图；
- `fig_v02_04_functional_three_class_heatmap.png`：功能正控、候选、hard negative 动态方向；
- `fig_v02_05_baseline_coverage_grid.png`：全部基线的运行状态格子图；
- `fig_v02_06_candidate_handoff_scatter.png`：同一 PACER-200 双模型候选交接；
- `binding_recomputed_v02.json/csv`：统一指标和 bootstrap；
- `benchmark_method_registry_v02.csv`：基线协议、状态、来源；
- `candidate_handoff_audit_v02.json`：候选身份和哈希审计；
- `integration_run_receipt_v02.json`：完整命令、退出码和产物 SHA256。

## 5. 分支联调状态

以下远端分支均已成为本整合分支的 Git 祖先：OneProt G0、encoder optimization、FKG long-MD、Geom2Vec pilot、20 ns close-loop、CM00734 Stage B、DrugCLIP handoff、project benchmark v01。历史实现和冻结证据被保留，但默认主线仍是：

`family-aug binding -> corrected PACER-200 handoff -> static gate -> matched four-context MD -> frozen C1-BS256/PACER-FKG audit`。

## 6. 一条命令复现

```powershell
python project/scripts/run_project_benchmark_v02.py `
  --repo-root D:\CLC `
  --artifact-root D:\CLC `
  --old-ranking D:\CLC\project\results\drugclip_m4_candidate_screen_v01\ranked_candidates.csv `
  --famaug-scores D:\CLC\project\data\drugclip_muscarinic_family_aug_v01\cand200\cands200_famaug_scores.csv `
  --bootstrap 20000
```

原始预测和大轨迹不进入 Git；它们通过路径、行身份与 SHA256 纳入审计。若输入缺失、候选身份不一致或任一步退出非零，编排器立即失败，不生成成功 receipt。

## 7. 目前最诚实的结论

已经完成的是：可复现的多模块联调、13T 结合检索的实质改善、同库候选交接修复，以及三类分子的功能动态方向 pilot。尚未完成的是：足够规模的前瞻功能标签验证、同协议运行全部结构 SOTA、稳定 EF1% 优势和湿实验 PAM 确认。因此目前可写成“有严格基准和机制线索的候选优选框架”，不能写成“已经发现有效 PAM”或“全流程达到 SOTA”。

## 8. 测试状态与环境缺口

- 联调、指标、路径保护、Phase-2 校准和 FKG 数学测试：`13 passed + 20 subtests passed`。
- `test_phase1_bs256.py` 在当前 CPU 环境缺少 `torch_geometric`，测试收集失败；它需要 GPU/图网络环境。
- `test_phase3_fkg.py` 因 git-ignored 的 `V02_FREEZE_MANIFEST.json` 不在本机而 fail-closed；预期 SHA256 为 `b48bc74a...e11bd`。这证明冻结门在生效，但也说明 GPU 机器必须补交该工件后才能完成完整 Phase-3 复验。
- 以上两项不影响 v02 的结合指标重算和候选交接审计，但阻止我们宣称“任意新机器上完整动态训练链全部通过”。
