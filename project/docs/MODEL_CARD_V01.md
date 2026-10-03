# Model Card — PACER-M4 提交模型说明 v01

依据附件5第三节"模型说明文件"要求编写。适用日期：2026-10-02。

## 1. 模型清单与来源披露

| 模型 | 类型 | 来源 | 本项目工作 |
|---|---|---|---|
| DrugCLIP（2026, "Science-2026" checkpoint） | 开源预训练口袋-分子对比模型 | 官方发布（未训练时 AUC 见仓库基准） | 参数高效微调（双投影），3 随机种子，权重=投影适配器 |
| DrugCLIP 2023 checkpoint | 开源预训练 | 官方发布 | GPCR-LOTO 适配器训练（3 种子） |
| family-augmented 微调 | 上项的微调变体 | 本项目 | ChEMBL M2/M4 家族增强数据 |
| Geom2Vec / C1-BS256 编码器 | 开源 ViSNet 基座 | arXiv:2409.19838 + vendored checkpoint | 冻结编码器优化链（C0→C1→BS256），仓库含冻结 tag |
| PACER-FKG 判别规则 | 本项目原创 | — | 四上下文差分 + 区域方向一致性，全程冻结 |

**权重提交方式**：投影适配器（.projection.pt ×9）与 200 候选冻结表示（.npz）已入仓库
`project/artifacts/drugclip2023_m4_loto_pacer200_v01/`；base checkpoint 因体积与许可不入库，
复算不需要它（表示已冻结，见 `upstream_frozen_inputs/README.md`）。

## 2. 训练与微调说明

- 训练入口：`project/scripts/`（finetune 系列，见仓库 CODE_AND_ARTIFACT_MAP_V01.md 索引）
- 数据划分：严格整靶点留出（LOSO/LOTO）+ Murcko 骨架隔离 + 3 种子（20260925/26/27）
- 关键冻结结果：13 靶 famaug 集成 ROC 0.6414（CI [+0.067,+0.102] vs 官方）；20 靶 [+0.028,+0.062]
- 泄漏防控：随机标签对照（0.496）、预注册协议、冻结 SHA 审计、独立复算 6/6 对账
- 随机种子：全部记录于结果 JSON 与审计文件

## 3. 适用范围与已知局限

**适用**：M4 别构口袋结合候选的排序与检索；已冻结四上下文动态复核对"方向可复现性"的判别
（回顾性三类：LY2119620 +0.634 / compound110 +0.584 / CM00734 −0.327）。

**不适用/局限**：
1. 不输出"PAM 概率"，无校准阈值；EF1% 早期富集未超过二维基线
2. 功能判别验证样本 n=3，无通用分类器主张
3. 静态结合表示在外部功能终点上接近随机（外部基准 v02 第一层实测）
4. 候选均为计算假设，确认 PAM 数量 = 0

## 4. 输入输出

- 主入口：`python -m pacer_m4`（CLI/API，含 dry-run 与 claim boundary 强制随行）
- 最终结果：`project/scripts/generate_submission_results_v01.py` → `results.csv`（UTF-8，
  字段符合附件5第四节；Stage3D manifest 存在时输出 3 候选冻结短名单，否则 200 全量 pre-MD）

## 5. 第三方组件版本

requirements.txt（Python 3.9.23 冻结）+ Dockerfile + docker-compose.yml；DrugCLIP /
Geom2Vec / Vina / Schrödinger 的名称、版本、许可与调用参数见
`project/docs/CODE_AND_ARTIFACT_MAP_V01.md` 与各阶段审计 JSON。
