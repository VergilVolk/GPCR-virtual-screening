# PACER-M4 代码发布验证 v01

验证日期：2026-10-02  
分支：`codex/pacer-m4-service-v01`

## 1. 源码与文档完整性

执行：

```bash
python project/scripts/validate_pacer_m4_release.py \
  --output project/release_audit_v01.json
```

结果：`valid = true`。manifest 中的项目入口、核心文档、DrugCLIP 路由、六路 docking、PACER-FS、Pareto、PACER-FKG 和测试文件均存在。四个可选轻量证据目录也全部存在。

核心 `requirements.txt` 的 12 个固定直接依赖均已确认存在 CPython 3.9/Windows wheel；报告生成器所需的 `reportlab==4.2.5` 已显式加入。项目元数据已由 Python `>=3.10` 修正为 `>=3.9`，7 个包入口文件通过 Python 3.9 grammar parse，生成 wheel 的 `Requires-Python` 为 `>=3.9`。当前本机只有 Python 3.11/3.13，因此完整精确版本导入闸门仍应由队友在 `chrm4_vs` Python 3.9.23 环境运行 `audit_core_environment.py` 最终确认。

## 2. 自动测试

执行：

```bash
python -m pytest project/tests -q
```

结果：`37 passed, 2 skipped`。

两项 skip 是未分发的 Science-2026 大型开发 checkpoint 完整性测试。测试现在能区分“整套外部 checkpoint 未安装”和“checkpoint 只复制了一部分”：前者跳过，后者仍然失败。

## 3. Python 包与命令入口

以下操作均成功：

```bash
python -m pip install -e . --no-deps
pacer-m4 stages
python -m pip wheel . --no-deps --wheel-dir tmp/wheels
```

生成标准 wheel：`pacer_m4-0.1.0-py3-none-any.whl`。阶段注册表包含 19 个入口，覆盖仓库校验、DrugCLIP、四种 docking 层级、PACER-FS、Pareto、四上下文 MD、动态证据计分和项目 benchmark。

## 4. PACER-200 真实推理复算

使用打包的 512 维分子/口袋表示和三个 M4-held-out 适配器，重新计算 200 个候选：

- 行数：200；
- 三 seed 权重和表示均通过 SHA256 读取；
- 分数全部有限；
- 分数范围：`[-0.384211, 0.152515]`。

随后通过统一入口重新运行冻结 M4-safe 路由。与仓库冻结排名比较：

- `pacer_binding_score` 最大绝对差：`0.0`；
- 排名不一致：`0/200`；
- 路由不一致：`0/200`；
- 所有分子均使用 `m4_2023_gpcr_loto` 路由。

这证明打包的 M4 结合检索结果可在 CPU 环境逐位复现。它不证明候选具有 PAM 功能。

## 5. 候选选择代码修复

`select_pareto_candidates.py` 已从固定本地路径和单行实现改为：

- 显式 `--predock`、`--docking-features`、`--known-potency`、`--output-dir`；
- 可测试的 Pareto front、ECFP/Butina 多样性和 portfolio 函数；
- 缺列、空连接和非法 SMILES 的明确错误；
- 保持原冻结科学规则：Vina affinity 与未验证效力不进入 Pareto 目标。

新增 3 个测试验证非支配关系、front 分层和多样性选择。

## 6. 本地 API 烟雾测试

实际启动 Uvicorn 后完成：

- `GET /health`：200，`execution_allowed=false`；
- `GET /v1/capabilities`：返回 19 个注册阶段；
- `POST /v1/stages/drugclip-route/run`：返回 `dry_run` 和精确命令；
- 默认真实执行请求：403；
- API 自动测试：5 项通过。

## 7. Docker 一致性

容器入口已从 Python 3.11 改为与冻结核心环境一致的 `python:3.9.23-slim-bookworm`。Docker 构建按以下顺序执行：

1. 从 PyTorch CPU 索引安装 `torch==2.8.0`；
2. 安装 `requirements.txt` 的冻结核心依赖；
3. 安装 `requirements-api.txt` 中仍支持 Python 3.9 的 FastAPI/Uvicorn；
4. 以 `--no-deps` 安装 PACER-M4 包；
5. 执行核心环境审计和阶段能力 smoke test。

Docker Hub 已确认基础镜像 tag 存在。本机没有 Docker 命令，因此这里不声称已经完成镜像实构建；最终提交前须在装有 Docker 的机器运行：

```bash
docker compose build --no-cache
docker compose up -d
curl http://127.0.0.1:8000/health
docker compose down
```

## 8. 仍需外部资源的部分

- 任意新分子的完整 DrugCLIP backbone 编码；
- Glide PDB/BEmin/BEavg 三路商业计算；
- 四上下文生产 MD 的 GPU 时间和大型轨迹；
- 新候选的功能湿实验。

统一命令入口会报告这些依赖，但不会伪造、插补或绕过它们。
