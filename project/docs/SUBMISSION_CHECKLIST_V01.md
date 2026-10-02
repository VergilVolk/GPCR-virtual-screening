# PACER-M4 比赛提交检查表 v01

更新时间：2026-10-02

## 必交材料

- [ ] 项目评审表 PDF；
- [ ] 带讲解录音的项目介绍 PPT；
- [ ] 项目原始代码；
- [ ] 文件名符合赛事系统要求。

Docker 是可复现性补充，不替代原始代码。提交包必须保留 Python 源码、环境文件、运行说明和证据索引。

## 匿名与合规

- [ ] 评审表、PPT 和代码中不出现学校名称、校徽、成员或指导老师姓名、照片、邮箱和本机用户名；
- [ ] 不提交口令、访问令牌、许可证文件或私人路径；
- [ ] 第三方模型、数据和商业软件写明来源与许可证边界；
- [ ] 未经湿实验确认的分子只称“计算优先候选”；
- [ ] 团队成员能够解释提交代码、数据来源、失败结果和主要结论。

赛事公开通知强调原创、独立完成、模型可复现和结论可解释，未列出所谓“AIGC率”阈值。提交前应人工重读文字和代码，删除空泛、重复或夸张表述，但不以规避检测器为目标。

## 代码验收

```bash
python project/scripts/validate_pacer_m4_release.py \
  --output project/release_audit_v01.json
python -m pytest project/tests -q
python -m pacer_m4 capabilities \
  --output project/capabilities_v01.json
```

若有 Docker 主机，再执行：

```bash
docker compose build --no-cache
docker compose up -d
curl http://127.0.0.1:8000/health
docker compose down
```

## 提交前人工确认

- [ ] README 中的指标均能追溯到仓库证据；
- [ ] 没有把结合检索、docking 或回顾性动态结果写成功能 PAM 实验证明；
- [ ] `git status` 仅包含计划提交的文件；
- [ ] 从干净目录解压提交包后，按 README 复核一次；
- [ ] 最终提交包不包含 `.git`、缓存、轨迹、大型 checkpoint 和临时日志。
