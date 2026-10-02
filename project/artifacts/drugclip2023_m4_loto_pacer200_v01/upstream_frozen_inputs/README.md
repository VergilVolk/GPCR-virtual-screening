# Upstream frozen inputs（2026 侧补齐包）

本目录补齐 router 上游的全部"本机孤儿"文件，使 Stage3 及其复算在任何克隆上闭环，**不需要 base checkpoint、不需要 GPU、不需要 molecules.lmdb**。

## 文件清单

| 文件 | 作用 |
|---|---|
| `cands200_famaug_scores.csv` | 200 分子的 famaug 三种子 + 集成分数（`pacer200_two_model_scores.csv` 的 2026 列的原始来源，SHA 19de17ae… 与 handoff 审计一致） |
| `cands_science2026.npz` | 200 分子的 Science2026 冻结表示（forward inference 已完成的产物） |
| `m4r_pocket_science2026_subsets.npz` | 10 个 M4R 口袋的 Science2026 冻结表示（从 69MB 全口袋 npz 中抽取的 M4R 子集） |
| `predock_portfolio.csv` | 候选推理第一棒的 200 分子组合（与 docking 2000 结果的分子集对齐基准） |

## 为什么不需要 base checkpoint

分数链 = 冻结表示（npz，本目录）→ 仓库内 3-seed 投影（`../drugclip2023_m4_loto_seed*.projection.pt` 与 `project/results/drugclip_science2026/family_aug_v01/*.projection.pt`）→ 点积。lmdb → checkpoint 的 forward 是一次性历史步骤，产物即 npz，不需要也无法在克隆上重演（base checkpoint 按大小政策不入库）。

## 最小复算路径（任意机器，CPU，秒级）

```python
import numpy as np, torch
mol = np.load('cands_science2026.npz', allow_pickle=False)
pkt = np.load('m4r_pocket_science2026_subsets.npz', allow_pickle=False)
m = torch.as_tensor(mol['molecule_representations'].astype(np.float32))
p = torch.as_tensor(pkt['pocket_representations'].astype(np.float32))
# 对每个 seed 加载 family_aug 投影后：proj_mol(m) @ proj_pkt(p).T 取 max——
# 与 cands200_famaug_scores.csv 的 seed 列逐一对照（脚本见
# D:\CLC 侧 rescore_cands200_famaug.py，本目录文件路径可替换）
```

## Router 使用（零复算）

`pacer_drugclip_router.py` 直接吃 `../pacer200_two_model_scores.csv`（列契约
`pair_id,target,score_2023_gpcr_loto,score_2026_13t_famaug` 已满足）；最终排名
`../pacer200_m4_safe_routed_ranking.csv` 及审计 JSON 已在本 bundle。
