# PACER-M4 本地 API 与 Docker v01

更新时间：2026-10-02

## 1. 设计边界

API 复用 `pacer_m4/stages.py` 中同一份阶段注册表，不重新实现 DrugCLIP、docking、PACER-FS 或 PACER-FKG。每个阶段保留自己的参数、输出和科学主张边界。

服务默认禁止真正执行脚本，只允许：

- 查看健康状态；
- 查询 19 个阶段及依赖；
- 获取某一阶段的说明；
- 生成 dry-run 命令和审计信息。

只有在可信本地机器显式设置 `PACER_M4_API_ALLOW_EXECUTION=1` 后，API 才接受执行任务。不要把允许执行的服务暴露到公网。

## 2. 本地启动

```bash
pip install -e ".[api,core]"
uvicorn pacer_m4.api:app --host 127.0.0.1 --port 8000
```

浏览器接口：

- OpenAPI UI：`http://127.0.0.1:8000/docs`
- 健康检查：`http://127.0.0.1:8000/health`
- 能力列表：`http://127.0.0.1:8000/v1/capabilities`

## 3. API

### 查询能力

```bash
curl http://127.0.0.1:8000/v1/capabilities
```

返回每个阶段的代码状态、Python 模块状态、硬件说明、外部依赖和主张边界。`python_ready=true` 不表示商业软件、轨迹或受体文件一定存在；应同时检查 `execution_readiness` 与 `external_requirement`。

### 预演任务

```bash
curl -X POST http://127.0.0.1:8000/v1/stages/drugclip-route/run \
  -H "Content-Type: application/json" \
  -d '{"args":["--input","scores.csv","--output","routed.csv"]}'
```

默认返回 `status=dry_run`，不会运行脚本。

### 允许本地执行

PowerShell：

```powershell
$env:PACER_M4_API_ALLOW_EXECUTION = "1"
uvicorn pacer_m4.api:app --host 127.0.0.1 --port 8000
```

提交时在 JSON 中增加 `"execute": true`。服务返回 job id；状态通过：

```text
GET /v1/jobs/{job_id}
```

日志和回执写入 `runs/pacer_m4_api/`。API 使用后台任务避免请求一直等待，但它不是分布式调度系统；服务器重启后应以磁盘回执为准。

## 4. Docker

```bash
docker compose up --build
```

默认容器同样禁止执行。需要挂载外部数据、合法软件和结果目录后，才可在可信环境打开执行开关。Docker 镜像不包含：

- Schrödinger Glide；
- MD 轨迹；
- 基础 DrugCLIP 大 checkpoint；
- OpenMM GPU/CUDA 环境。

因此默认镜像适合文档、能力查询、轻量路由和 CPU 分析，不是假装“一镜像运行全部重计算”。GPU MD 应继续使用已冻结的 WSL/conda 交接环境。

## 5. 已完成验证

- FastAPI 健康检查返回 200；
- 19 个阶段可查询；
- DrugCLIP 路由 dry-run 返回精确命令；
- 默认执行请求返回 403；
- API 自动测试 5 项通过；
- 本机未安装 Docker，因此 Dockerfile 和 compose 已完成静态配置，但尚未在本机实际构建镜像。

