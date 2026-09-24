# DrugCLIP 功能适配器：CPU 初步结果与 GPU 分工

## 结论

DrugCLIP 可以纳入 PACER，但现阶段不应直接全量微调。官方 DrugCLIP 静态结合分数无法直接区分 M4 functional PAM 与 inactive；冻结其预训练表征后训练低容量分类头，在 scaffold 留出上表现较高，但在更严格的 source/series 留出上失效，说明当前主要瓶颈是化学系列与文献来源混杂，而非算力不足。

## 已完成实验

- 数据：529 个 M4 分子（463 PAM、66 inactive）。
- 结构上下文：7TRQ、7TRP、7TRS 三个 M4 别构口袋。
- DrugCLIP：官方 checkpoint，编码器冻结，128 维分子表征。
- CPU 适配器：标准化 + L2 logistic regression，class weight balanced。
- 主设定：`C=0.1`，测试折不参与调参。

| 方法 | scaffold holdout AUC | source holdout AUC | series holdout AUC |
|---|---:|---:|---:|
| DrugCLIP 三状态 cosine 适配器 | 0.681 | 0.610 | 0.620 |
| DrugCLIP 分子 embedding 适配器 | 0.855 | 0.396 | 0.346 |
| embedding + 三状态 cosine | 0.856 | 0.400 | 0.348 |
| 既有 2D 最佳基线 | 0.919 | 0.630 | 0.624 |

零样本 DrugCLIP 的三个单结构 AUC 为 0.497、0.446、0.540。高 AP 不能单独解释为性能好，因为阳性率为 87.5%。

正则化敏感性（`C=0.01/0.1/1.0`）显示，source holdout 的分子 embedding AUC 在 0.306–0.564 之间，结论不稳定且未超过既有 2D baseline。不得将 scaffold AUC=0.855 宣称为可迁移性能。

## 科学解释

1. DrugCLIP 预训练表征含有较强的化学结构信息，因此在 scaffold 被相对均匀分散的划分中能够拟合标签。
2. 当整篇文献或完整化学系列被留出时，标签分布和化学空间同时改变，模型无法迁移。
3. DrugCLIP 原任务是蛋白口袋—配体结合匹配，不是 ACh 条件下的 PAM 功能预测。此结果支持“结合兼容性不等于 PAM 功能”，但尚不证明 PACER-DC 有效。

## 设备分工

- GPU 游戏本：优先连续完成缺失的四上下文 MD、OneProt-MD embedding；不要与 DrugCLIP 全模型训练并行抢占同一块 GPU。
- CPU 设备：数据审计、冻结 embedding、linear/低秩 adapter、source/series holdout、bootstrap 和消融。
- 只有当冻结 adapter 在严格 source/series 外推上稳定超过 2D baseline 后，才安排 GPU 做 LoRA 或末层微调。

## 复现命令

```bash
python project/scripts/train_drugclip_m4_functional_probe.py \
  --embeddings project/results/drugclip_m4_benchmark_v01/embeddings.npz \
  --benchmark project/data/benchmarks/m4_pam_v1/pam_vs_inactive.csv \
  --output project/results/drugclip_m4_benchmark_v01/functional_probe.json \
  --c 0.1 --bootstrap 2000
```

证据等级：`retrospective_frozen_drugclip_functional_probe`。该实验不是湿实验 PAM 验证，也不是 DrugCLIP backbone 全量微调。
