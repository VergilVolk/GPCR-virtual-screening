# PACER-200 DrugCLIP 2023 M4-LOTO 离线推理包

本目录已经包含完整推理所需材料，不需要下载 1.18 GB DrugCLIP 基础 checkpoint，也不需要 DrugCLIP、Uni-Core、LMDB、GaMD 或GPU。

## 已包含

- 三个固定 seed 的 M4-held-out 投影权重；
- PACER-200 的 2023 DrugCLIP 512维冻结分子表示；
- M4 cluster0 的512维冻结口袋表示；
- 已计算的三seed分数及均值；
- 与2026 family-aug分数合并后的输入；
- M4-safe路由最终排名；
- SHA256审计文件。

三个适配器仅使用 B2AR、CCR2、M2R 训练，M4 标签未进入训练。适配器集成对原冻结LOTO预测的最大复现误差为 `4.97e-08`。

## 最快使用方式

最终结果已经算好：

```text
pacer200_m4_safe_routed_ranking.csv
```

关键字段：

```text
pair_id,target,canonical_smiles,pacer_binding_score,final_rank,route
```

## 在队友电脑重新推理

最低依赖：

```bash
pip install torch numpy pandas
```

在仓库根目录运行：

```powershell
python project/scripts/score_pacer200_drugclip2023_m4_loto.py `
  --bundle-dir project/artifacts/drugclip2023_m4_loto_pacer200_v01 `
  --output project/artifacts/drugclip2023_m4_loto_pacer200_v01/reproduced_scores.csv
```

输出包含三个seed分数、平均分、seed标准差和最终排名。

随后可重新运行路由：

```powershell
python project/scripts/pacer_drugclip_router.py `
  --input project/artifacts/drugclip2023_m4_loto_pacer200_v01/pacer200_two_model_scores.csv `
  --output project/artifacts/drugclip2023_m4_loto_pacer200_v01/reproduced_routed_ranking.csv
```

## 科学边界

这些分数衡量 M4 结合候选检索，不是PAM功能、协同性或效力预测。最终PAM判断仍需PACER-FKG四上下文分析或实验验证。
