# PACER-M4 最小复现指南 v01

更新时间：2026-10-02

本指南用于验证仓库核心代码和轻量产物。它不下载大型轨迹，不重新运行数百纳秒 MD，也不要求 Schrödinger 许可证。

## 1. 环境

提交级核心推理环境冻结为 CPython 3.9.23：

```bash
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu
pip install -r requirements.txt
python project/scripts/audit_core_environment.py
```

若机器只使用 CPU，建议先从 PyTorch 官方 CPU index 安装 `torch==2.8.0`，再安装 requirements；其余固定版本均存在 CPython 3.9 Windows wheel。开发测试另行安装：

```bash
pip install pytest
```

完整模块环境分别记录在：

- `project/environment_drugclip_cpu.yml`
- `project/environment_vs.yml`
- `project/environment_sbdd.yml`
- `project/environment_pacer_dc_md.yml`
- `project/environment_pacer_dc_geom2vec.yml`

## 2. 验证提交内容

```bash
python project/scripts/validate_pacer_m4_release.py \
  --output project/release_audit_v01.json
```

该命令验证正式入口、核心脚本、测试和轻量产物是否存在，并写出 SHA256。大型 checkpoint 和 MD 轨迹属于可选外部资产；缺失时会报告，但不会把源码发布判为失败。

也可以先查看统一阶段注册表：

```bash
python -m pacer_m4 stages
python -m pacer_m4 capabilities
```

需要调用原脚本时使用：

```bash
python -m pacer_m4 run <stage-id> -- <原脚本参数>
```

加 `--dry-run` 可只生成命令和审计回执而不执行计算。

## 3. 运行核心测试

```bash
python -m pytest project/tests -q
```

未安装可选大型 checkpoint 时，对应 checkpoint 完整性测试应显示 `skipped`，而不是失败。当前冻结仓库预期为 `25 passed, 2 skipped`。

## 4. 复算 PACER-200 M4 DrugCLIP 排名

无需基础 DrugCLIP checkpoint 或 GPU：

```bash
python project/scripts/score_pacer200_drugclip2023_m4_loto.py \
  --bundle-dir project/artifacts/drugclip2023_m4_loto_pacer200_v01 \
  --output project/artifacts/drugclip2023_m4_loto_pacer200_v01/reproduced_scores.csv
```

随后执行冻结 M4-safe 路由：

```bash
python project/scripts/pacer_drugclip_router.py \
  --input project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_two_model_scores.csv \
  --output project/artifacts/drugclip2023_m4_loto_pacer200_v01/reproduced_routed_ranking.csv
```

检查输出中的：

- `route = m4_2023_gpcr_loto`；
- `eligible_for_pam_claim = false`；
- 生成的 `.audit.json` 包含输入输出 SHA256。

## 5. 复核六路 PACER-XR 合同

先生成候选模板：

```bash
python project/scripts/apply_pacer_xr_candidates.py --create-template
```

填入六路 protocol-matched 原始分数后运行：

```bash
python project/scripts/apply_pacer_xr_candidates.py \
  --input project/results/pacer_xr_candidate_transfer_v01/candidate_six_channel_raw_scores.csv
```

若缺少任一路，程序必须失败并指出缺失列。三条 Glide 分数需要外部合法 Schrödinger 环境；仓库不提供该商业软件，也不允许用 Vina 或均值填充。

## 6. 复核 PACER-FKG 数学核心

```bash
python -m pytest project/tests/test_pacer_factorial_kernel_graph.py -q
```

测试覆盖：

- 非加和四上下文交互可被 `Delta_INT` 捕获；
- candidate-only 自激活轴与交互轴可分离；
- 分布统计不依赖轨迹帧的排列顺序。

这只验证估计器实现，不等于 M4 药理性能。真实闭环证据见：

- `project/docs/PACER_DC_CLOSE_LOOP_20NS_FINAL_REPORT_v01.md`
- `project/results/pacer_dc_close_loop_20ns_v01/`
- `project/results/pacer_dc_cm00734_stage_b_20ns_analysis_v01/`

## 7. 完整四上下文 MD

完整 MD 不是轻量复现的一部分。运行者必须使用相同受体底座、膜环境、参数化协议和匹配 seed，依次完成体系构建、平衡、去约束、生产模拟和端点提取。开始前阅读：

- `project/docs/PACER_DC_CLOSE_LOOP_CONTRACT_v01.md`
- `project/docs/PACER_DC_GPU_HANDOFF.md`
- `project/docs/PACER_FKG_IMPLEMENTATION_AUDIT_20260927.md`

短平衡轨迹只能标记为 `equilibration_pilot`，不得用于 PAM 功能结论。独立轨迹的相同帧号没有物理配对关系，不得逐帧相减。

## 8. 结果解释

| 输出 | 可以说明 | 不能说明 |
|---|---|---|
| DrugCLIP 路由分数 | 结合候选的相对排序 | PAM 功能、效力、协同性 |
| docking / PACER-XR | pose 与口袋相容性、广义变构配体富集 | 功能 PAM 或 ACh 依赖增强 |
| PACER-FS | 有少量锚点时的系列内效力排序 | 零样本 PAM 判定 |
| PACER-FKG | 四上下文动态交互的跨 replica 重现性 | 湿实验确认或通用 PAM 概率 |

所有新分子在功能实验确认前统一称为“计算优先候选”。
