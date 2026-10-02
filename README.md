# PACER-M4

PACER-M4 是面向人源 CHRM4/M4 正性别构调节剂（PAM）的分层计算筛选项目。项目将“可能结合 M4 别构口袋”和“可能在 ACh 条件下产生功能性正变构调节”作为两个不同任务，分别用结合检索、结构门控和四上下文动态分析处理。

## 从这里开始

- 项目全貌与阅读顺序：[`project/docs/START_HERE_FOR_COLLABORATORS_V01.md`](project/docs/START_HERE_FOR_COLLABORATORS_V01.md)
- 全部文档与证据索引：[`project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md`](project/docs/DOCUMENT_AND_EVIDENCE_INDEX_V01.md)
- 开发与候选推理主线：[`project/docs/PACER_M4_PROJECT_MAINLINE_DEVELOPMENT_AND_INFERENCE_V01.md`](project/docs/PACER_M4_PROJECT_MAINLINE_DEVELOPMENT_AND_INFERENCE_V01.md)
- 当前整体结果与边界：[`project/docs/PACER_M4_INTEGRATED_ALGORITHM_AND_BENCHMARK_V02.md`](project/docs/PACER_M4_INTEGRATED_ALGORITHM_AND_BENCHMARK_V02.md)
- 代码、数据和产物地图：[`project/docs/CODE_AND_ARTIFACT_MAP_V01.md`](project/docs/CODE_AND_ARTIFACT_MAP_V01.md)
- 最小复现指南：[`project/docs/REPRODUCTION_QUICKSTART_V01.md`](project/docs/REPRODUCTION_QUICKSTART_V01.md)
- 本地 API 与 Docker：[`project/docs/LOCAL_API_AND_DOCKER_V01.md`](project/docs/LOCAL_API_AND_DOCKER_V01.md)

## 方法概览

```text
分子生成或独立分子库
  -> 结构标准化与药化过滤
  -> DrugCLIP 路由结合检索
  -> 单结构 / 多状态 / GaMD ensemble docking
  -> 六路 PACER-XR 结构级联与 IFP 门控
  -> 适用域、化学多样性与不确定性筛选
  -> 四上下文 MD（A、P、C、CP）
  -> PACER-FKG 动态功能复核
  -> 候选证据卡与实验建议
```

四上下文定义：

- `A`：apo；
- `P`：ACh-only；
- `C`：candidate-only；
- `CP`：candidate + ACh；
- `Delta_INT = CP - P - C + A`：候选与 ACh 的非加和动态交互。

## 已冻结的主要结果

- DrugCLIP M4-safe 路由在统一 4-GPCR 开发基准上达到 ROC-AUC `0.735`、PR-AUC `0.323`、BEDROC20 `0.398`、EF5% `4.474`；EF1% 增益仍不确定。
- PACER-XR Cascade 在公开 M4 broad allosteric-modulator benchmark 上达到 AUC `0.7775`、EF0.5% `20.40`、EF1% `13.73`。该任务不等于功能 PAM 预测。
- PACER-FS 在共同 11 系列三锚点协议中将 macro Spearman 从 `0.199` 提高到 `0.263`，增益 95% CI 为 `[+0.017,+0.113]`。
- 冻结 PACER-FKG 回顾性检验中，主区域跨 replica 方向一致性为 LY2119620 `+0.634`、compound110 `+0.584`、实验 inactive 近邻 CM00734 `-0.327`。该结果是初步功能特异性证据，不是通用 PAM 分类器。

## 快速复核

轻量测试不需要 GPU：

```bash
python -m pytest project/tests/test_pacer_drugclip_router.py \
  project/tests/test_drugclip_generation_bakeoff.py \
  project/tests/test_pacer_factorial_kernel_graph.py
```

统一命令入口不会改变各模块算法，只负责发现、校验和调用已经审计的脚本：

```bash
python -m pacer_m4 stages
python -m pacer_m4 capabilities --output project/capabilities_v01.json
python -m pacer_m4 run drugclip-route --dry-run -- \
  --input scores.csv --output routed.csv
```

每个阶段保留自己的输入、输出和科学主张边界；命令入口不会把不同任务的分数合成虚假的端到端概率。

可选本地 API：

```bash
pip install -e ".[api,core]"
uvicorn pacer_m4.api:app --host 127.0.0.1 --port 8000
```

API 默认只允许查询和 dry-run。真正执行必须在可信本地环境显式设置 `PACER_M4_API_ALLOW_EXECUTION=1`。

复算已经打包的 PACER-200 M4 DrugCLIP 排名：

```bash
python project/scripts/score_pacer200_drugclip2023_m4_loto.py \
  --bundle-dir project/artifacts/drugclip2023_m4_loto_pacer200_v01 \
  --output project/artifacts/drugclip2023_m4_loto_pacer200_v01/reproduced_scores.csv
```

重新执行冻结 DrugCLIP 路由：

```bash
python project/scripts/pacer_drugclip_router.py \
  --input project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_two_model_scores.csv \
  --output project/artifacts/drugclip2023_m4_loto_pacer200_v01/reproduced_routed_ranking.csv
```

环境说明见：

- `project/environment_drugclip_cpu.yml`
- `project/environment_vs.yml`
- `project/environment_sbdd.yml`
- `project/environment_pacer_dc_md.yml`
- `project/environment_pacer_dc_geom2vec.yml`

## 仓库边界

- DrugCLIP 与 docking 输出衡量结合检索或结构相容性，不能单独证明 PAM 功能。
- Glide 三路依赖 Schrödinger 商业许可；仓库不会分发该软件，也不会在缺失三路时填补 PACER-XR 分数。
- MD 轨迹和大型第三方 checkpoint 不进入 Git；仓库保留协议、校验值、轻量产物和可追溯结果。
- 新生成分子未经湿实验确认时，只称“计算优先候选”。

## 许可证

本仓库代码遵循 [`LICENSE`](LICENSE)。第三方数据、模型和商业软件仍受各自许可证约束。
