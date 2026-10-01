# DrugCLIP 未见靶点迁移：适配器与损失消融

更新日期：2026-09-26

## 目的

本实验不再把多个模块合并后只报告一个总分，而是回答三个具体问题：

1. 未见靶点迁移主要来自分子端还是口袋端适配？
2. target retrieval 与表示保持约束是否真正贡献性能？
3. 针对最差训练靶点加权的 Target-DRO 是否提高外推稳健性？

## 严格验证协议

- 4 个 GPCR 与 9 个 LIT-PCBA 靶点，共 37,524 个配对；
- 每折完整留出一个靶点；
- 从其余 12 个靶点删除所有与留出靶点重叠的 Murcko scaffold；
- 留出口袋不进入 BCE、retrieval、正则化或选模；
- 三个独立随机种子；
- 不确定性按靶点内 Murcko scaffold 成簇配对 bootstrap 计算。

## 第一轮完整消融

以下为 seed 20260926 的宏平均，用于筛除无效方向：

| 模式 | ROC-AUC | PR-AUC | BEDROC20 | EF5% |
|---|---:|---:|---:|---:|
| 官方 DrugCLIP | 0.562 | 0.201 | 0.241 | 1.78 |
| 仅分子端 LoRA | 0.580 | 0.200 | 0.241 | 1.59 |
| 仅口袋端 LoRA | 0.603 | 0.228 | 0.293 | 2.13 |
| 双侧 BCE | 0.609 | 0.239 | 0.286 | 2.16 |
| 双侧 BCE + retrieval + preservation | 0.617 | 0.244 | 0.300 | 2.25 |
| 双侧 BCE + retrieval，无 preservation | **0.626** | **0.251** | **0.307** | 2.34 |
| Target-DRO | 0.559 | 0.192 | 0.228 | 1.66 |
| 完整末层投影微调 | 0.591 | 0.221 | 0.266 | 2.00 |

由此可得：口袋端适配是主要来源，双侧适配进一步提高性能，target retrieval 有贡献；当前 Target-DRO 明显有害，已否决。分子端单独适配不能解释主结果。

完整末层投影微调更新 131,328 个参数，约为 rank-8 LoRA（10,240 个参数）的 12.8 倍，但同一种子下仍明显低于 LoRA。该结果支持低秩更新的结构正则化价值；由于 full-last 仅使用一个预设学习率，现阶段将其作为受控参数效率基线，不声称已穷尽完整微调的最优超参数。

## 三种子确认：去除表示保持约束

| 模式 | ROC-AUC | PR-AUC | BEDROC20 | EF1% | EF5% |
|---|---:|---:|---:|---:|---:|
| 标准双侧适配 | 0.623 | 0.247 | 0.300 | **2.39** | 2.28 |
| 无 preservation 双侧适配 | **0.637** | **0.258** | **0.311** | 2.24 | **2.35** |

无 preservation 相对标准版的 scaffold-cluster bootstrap 95% CI：

- ROC-AUC：`[+0.0079, +0.0197]`；
- PR-AUC：`[+0.0062, +0.0160]`；
- BEDROC20：`[+0.0010, +0.0210]`；
- EF1% 与 EF5% 区间跨零。

因此可以确认 preservation 在当前未见靶点迁移任务上限制了排序性能，但不能声称它改善了最前端 EF1%。

## 分子端与口袋端锚定

三种子平均：

| 锚定方式 | ROC-AUC | PR-AUC | BEDROC20 | EF5% |
|---|---:|---:|---:|---:|
| 保留口袋空间、放开分子端 | 0.626 | 0.249 | 0.302 | 2.29 |
| 保留分子空间、放开口袋端 | **0.635** | **0.256** | **0.310** | **2.33** |

分子锚定相对口袋锚定的 ROC-AUC、PR-AUC、BEDROC20 95% CI 均严格为正。结果支持以下机制解释：

> DrugCLIP 的通用分子语义可以适度保留，但口袋投影必须允许针对新的靶点分布发生重排。

完全不锚定的总体 ROC-AUC 仍略高，因此当前部署模型采用 `dual_no_preserve`；非对称锚定保留为后续跨域稳定性方案，而不是当前最佳模型。

## 当前方法边界

本结果证明的是 active/decoy 虚拟筛选和未见靶点早期富集迁移，不证明：

- M4 配体一定是 PAM；
- PAM 效力、协同性或内在激动；
- 通用虚拟筛选 SOTA；
- prospective 湿实验有效性。

## 复现文件

- `scripts/evaluate_drugclip_13target_ablation.py`
- `scripts/summarize_drugclip_13target_ablation.py`
- `results/litpcba_drugclip_external_v01/combined_13target_adapter_ablation_seed2026092*.json`
- `results/litpcba_drugclip_external_v01/combined_13target_no_preserve_vs_standard_3seed_bootstrap.json`
- `results/litpcba_drugclip_external_v01/combined_13target_mol_vs_pocket_anchor_3seed_bootstrap.json`

## 部署权重

推理时建议平均三个 `dual_no_preserve` 投影：

- `drugclip_13target_dual_no_preserve.seed20260925.projection.pt`  
  SHA256 `E80429516EE37717229225F024EC6C20F3149CA639F1FFEBBEE297EAAFB7BB1F`
- `drugclip_13target_dual_no_preserve.seed20260926.projection.pt`  
  SHA256 `32D71BBA40DF7D19DD6AF2464FD5752F985C5D6D3D70CD3D9B8AA50DA372AA18`
- `drugclip_13target_dual_no_preserve.seed20260927.projection.pt`  
  SHA256 `A5C2994F3D79AC1BEFFA3152BAA097C798B405BA20FDF4930E5AB2C0DC783912`

部署权重由全部13靶点训练，仅用于后续筛选；论文性能必须引用 LOSO 结果，不能引用全数据训练分数。
