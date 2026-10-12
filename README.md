# 面向 CHRM4 正性别构调节剂发现的动态—药理协同智能筛选

[![Frozen replay](https://github.com/VergilVolk/GPCR-virtual-screening/actions/workflows/replay.yml/badge.svg)](https://github.com/VergilVolk/GPCR-virtual-screening/actions/workflows/replay.yml)
[![License: MIT](https://img.shields.io/badge/code-MIT-6f4cff.svg)](LICENSE)
[![Release](https://img.shields.io/github/v/release/VergilVolk/GPCR-virtual-screening?color=6f4cff)](https://github.com/VergilVolk/GPCR-virtual-screening/releases)

这是华北赛区本科生组赛道三参赛项目的公开复现仓库。公开包聚焦一件可核验的事：把提交材料中的候选筛选结果，按照冻结输入、冻结模型、冻结规则重新计算，并由独立校验脚本检查模型哈希、模块输出和最终结果表。

> **结论边界：** 本仓库提供计算筛选与动态证据，不把模型分数解释为 PAM 活性。3 个候选的四条件受体动力学状态均为 `complete_unresolved`，湿实验状态均为 `not_performed`；PAM 功能仍需实验验证。

## 提交材料中的主要结果

| 评估项 | 提交结果 |
| --- | ---: |
| M4 测试集 ROC-AUC | 0.540 → **0.714** |
| M4 测试集 PR-AUC | 0.107 → **0.155** |
| M4 测试集 BEDROC20 | 0.125 → **0.180** |
| 13 靶点严格留出 ROC-AUC | 0.544 → **0.641** |
| 13 靶点严格留出 PR-AUC | 0.178 → **0.248** |
| 六路共识 ROC-AUC | **0.771** |
| 六路级联 ROC-AUC | **0.777** |

冻结筛选链路为 `2,605 生成分子 → 2,592 去参考泄漏分子 → 200 个 DrugCLIP 候选 → 8 个六通道记录 → 3 个四条件动力学候选`。最后 8 个记录由 5 个新生成候选和 3 个历史锚点组成，不能把“8”误写成 8 个全新候选。

## 一分钟复现

```bash
cd submission
python -m pip install -r requirements.txt
python predict.py
python verify_submission.py
```

成功时最后一行应为：

```text
PASS: frozen four-module replay, model hashes and final table
```

程序仅重放仓库内可再分发的冻结数据，不重新运行 Glide、分子动力学或训练流程。完整输入、输出和证据范围见 [`submission/docs/REPRODUCIBILITY_MATRIX.md`](submission/docs/REPRODUCIBILITY_MATRIX.md)。

## 三个进入四条件动力学分析的候选

| 答辩编号 | 生成编号 | 六路综合分数 | 计算状态 | 湿实验 |
| --- | --- | ---: | --- | --- |
| PACER0010 | PACERGEN01755 | 0.7304 | complete_unresolved | not_performed |
| PACER0027 | PACERGEN00123 | 0.7307 | complete_unresolved | not_performed |
| PACER0073 | PACERGEN00094 | 0.7041 | complete_unresolved | not_performed |

PACER0073 在正构激活核心中出现重复同向变化（方向余弦 `+0.603 / +0.211 / +0.261`），但这只是需要后续验证的动力学证据，不构成 PAM 功能证明。

## 仓库导航

- [`submission/`](submission/)：答辩与评审优先入口，可在普通 CPU 环境复现。
- [`submission/docs/SUBMISSION_CLAIMS.md`](submission/docs/SUBMISSION_CLAIMS.md)：申报书表述与公开证据的对应关系。
- [`submission/docs/REPRODUCIBILITY.md`](submission/docs/REPRODUCIBILITY.md)：逐步复现说明。
- [`docs/`](docs/)：静态展示面板源文件；不收集数据、不提供在线预测。
- [`docs/DOMAIN_DEPLOYMENT.md`](docs/DOMAIN_DEPLOYMENT.md)：独立域名选择与安全绑定方案。
- [`project/`](project/)：完整研究过程、算法开发和审计材料。
- [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)：第三方模型、代码和不可再分发组件说明。

## 开源与匿名边界

本项目代码按 [`LICENSE`](LICENSE) 中的 MIT 条款发布。DrugCLIP 模型权重及其输出受 CC BY-NC 4.0 约束；Geom2Vec 相关代码还包含 MIT/GPL-3.0 来源，详见第三方声明。Glide 为商业软件，本仓库不分发其程序或许可证。

为遵守盲审要求，公开材料不写学校、指导教师和成员姓名；引用时可使用 [`CITATION.cff`](CITATION.cff) 中的匿名团队信息。
