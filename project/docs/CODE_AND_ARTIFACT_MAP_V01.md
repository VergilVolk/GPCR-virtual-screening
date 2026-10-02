# PACER-M4 代码与产物地图 v01

更新时间：2026-10-02  
目的：区分可复用主代码、复现实验脚本、历史探索和大型外部依赖，便于合作者运行与审计。

## 1. 推荐代码入口

| 任务 | 入口 | 输入 | 输出 |
|---|---|---|---|
| DrugCLIP M4-safe 路由 | `project/scripts/pacer_drugclip_router.py` | 同一 target-molecule 行上的 2023/2026 两列分数 | 靶点内百分位、路由分数、最终排名、审计 JSON |
| PACER-200 离线重评分 | `project/scripts/score_pacer200_drugclip2023_m4_loto.py` | 打包的冻结表示、M4 口袋表示和三 seed 适配器 | 三 seed 分数、均值、标准差、排名 |
| DrugCLIP 代际统一跑分 | `project/scripts/run_drugclip_generation_bakeoff.py` | 冻结 benchmark 与各模型预测 | ROC、PR、BEDROC、EF 和配对 bootstrap |
| 六路 PACER-XR | `project/scripts/apply_pacer_xr_candidates.py` | Glide/Vina × PDB/BEmin/BEavg 六列原始分数 | 六路秩共识、1% cascade 排名、审计 JSON |
| Thompson–Miao 复现 | `project/scripts/run_thompson_miao_reproduction.py` | 作者公开分数、公共轨迹与配置 | 官方指标复算、聚类和复现审计 |
| GaMD ensemble docking | `project/scripts/run_m4_gamd_complete.py` | 候选、10 个代表构象、Vina 环境 | 多构象 poses、BEmin/BEavg 与任务台账 |
| 功能效力基准 | `project/scripts/run_pacer_baselines.py`、`run_lightgbm_loso.py` | 清洗的 M4 功能数据与冻结划分 | 系列留出 baseline |
| PACER-FS | `project/scripts/run_pacer_fs_complete.py`、`apply_pacer_fs_candidates.py` | 训练系列、三锚点或候选 | 系列内效力排序、校准与拒绝标记 |
| 候选 Pareto 选择 | `project/scripts/select_pareto_candidates.py` | 结构、药化、适用域和多样性证据 | 非加权 Pareto 候选表 |
| 四上下文体系构建 | `project/scripts/build_pacer_dc_reference_complexes.py` 至 `run_pacer_dc_production_md.py` | A/P/C/CP 体系定义、配体和膜环境 | 匹配 MD 轨迹与审计记录 |
| 四上下文端点提取 | `project/scripts/extract_pacer_dc_production_endpoints.py` | 匹配轨迹 | 轨迹 QC 与端点表 |
| PACER-FKG | `project/pacer_fkg_v02/` 与相应分析脚本 | 四上下文残基级轨迹表示 | Delta_PAM、Delta_AGO、Delta_INT、区域方向一致性 |
| 传统 PACER-DC 计分 | `project/scripts/pacer_dc_score.py` | 多 replica 证据表 | 动态证据与 Pareto 标记；保留作兼容/基线 |
| 项目 benchmark 汇总 | `project/scripts/run_project_benchmark_v02.py` | 各模块冻结指标 | 统一模块级 benchmark 与图表 |

## 2. DrugCLIP 可直接使用的推理包

目录：`project/artifacts/drugclip2023_m4_loto_pacer200_v01/`

包含：

- 三个 M4-held-out 投影适配器；
- PACER-200 的 512 维冻结分子表示；
- M4 cluster0 的 512 维冻结口袋表示；
- 三 seed 分数、两模型输入、M4-safe 最终排名；
- manifest 与 SHA256 审计。

该包不需要下载基础 DrugCLIP checkpoint，也不需要 GPU。它只覆盖已打包的 PACER-200 表示；对任意新分子进行完整编码仍需要相应 DrugCLIP backbone 或预先提取的冻结表示。

## 3. 六路 docking 的固定语义

PACER-XR 需要以下六列全部存在：

1. `Glide_PDB_raw`
2. `Glide_BEmin_raw`
3. `Glide_BEavg_raw`
4. `Vina_PDB_raw`
5. `Vina_BEmin_raw`
6. `Vina_BEavg_raw`

每一路先相对冻结 M4 参考分布转换成经验秩百分位，`PACER_XR` 为六路等权平均。Cascade 将 Glide-BEmin 前 1% 保持原顺序，其余按 PACER-XR 排序。缺少 Glide 许可或任一路分数时必须失败闭锁，不允许插值或用 Vina 复制填充。

与 docking 相关的其他入口：

- `dock_candidate_portfolio.py`：单一 M4 结构上的候选 docking；
- `dock_potency_coupling_features.py`：7TRQ/7TRP/7TRS 状态特征；
- `dock_m4_gamd_ensemble.py`：10 个 GaMD 代表构象；
- `analyze_m4_gamd_ensemble.py`：ensemble 分析；
- `rank_m4_gamd_candidates.py`：候选多构象排序；
- `audit_pacer_xr_functional_semantics.py`：检查 broad-AM 与 functional-PAM 标签不可混用。

## 4. 四上下文与 PACER-FKG 代码链

四上下文使用同一受体底座、膜、参数和匹配 seed：

```text
A  = apo
P  = ACh-only
C  = candidate-only
CP = candidate+ACh
```

主要步骤：

1. `prepare_pacer_dc_common_protein.py`：共同蛋白准备；
2. `build_pacer_dc_reference_complexes.py`：四上下文参考复合物；
3. `build_pacer_dc_openmm_reference_systems.py`：无溶剂 QC；
4. `build_pacer_dc_membrane_reference.py`：膜体系；
5. `run_pacer_dc_short_equilibration.py` 与 `run_pacer_dc_restraint_release.py`：平衡与分级去约束；
6. `run_pacer_dc_production_md.py`：生产轨迹；
7. `extract_pacer_dc_production_endpoints.py`：端点与 QC；
8. `project/pacer_fkg_v02/`：冻结 C1-BS256 表示、核交互、图区域与跨 replica 审计。

运行前必须阅读 `PACER_DC_CLOSE_LOOP_CONTRACT_v01.md`。轨迹窗口只能作为时间采样，不能冒充独立分子或生物学重复。

## 5. 结果和证据位置

Git 默认忽略 `project/results/` 中的大多数运行产物；已经纳入版本控制的轻量结果用于审计。关键目录包括：

- `project/results/drugclip_generation_bakeoff_v01/`
- `project/results/project_wide_integration_benchmark_v02/`
- `project/results/pacer_dc_close_loop_20ns_v01/`
- `project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01/`
- `project/results/pacer_xr_candidate_transfer_v01/`（若本地存在完整六路候选分数）

大型轨迹、第三方 checkpoint、外部源码和商业软件不应提交。文档中的 SHA256、commit、tag 和 manifest 用于恢复其来源。

## 6. 哪些脚本不是当前生产入口

`project/scripts/` 保留了大量模型探索和消融。以下命名通常表示研究过程，而非默认部署：

- `evaluate_*`、`compare_*`、`analyze_*`：评价或诊断；
- `finetune_*`、`train_*`：开发期训练；
- `*_preregistration` 对应脚本：预注册验证；
- OneProt、CGDA、State-MoE、早期 triplet/soft-coupling：历史或消融；
- `run_pacer_complete.py`：早期简化流水线，不代表当前完整主线。

是否属于当前主方法，以 `START_HERE_FOR_COLLABORATORS_V01.md` 和 `PROJECT_STATUS.md` 为准，而不是以文件名是否含 `final` 判断。

## 7. 最小提交检查

提交给评审或合作者前至少执行：

```bash
python -m pytest project/tests/test_pacer_drugclip_router.py
python -m pytest project/tests/test_drugclip_generation_bakeoff.py
python -m pytest project/tests/test_pacer_factorial_kernel_graph.py
python -m pytest project/tests/test_project_benchmark_v02.py
```

同时检查：

- 主结果是否能追溯到机器可读表；
- 是否误把 protocol/plan 写成 completed result；
- 是否把 broad-AM、binding、potency 和 functional PAM 指标混为一谈；
- 是否提交了轨迹、商业软件或不允许再分发的第三方权重；
- 新候选是否明确标注为未经过湿实验确认。

