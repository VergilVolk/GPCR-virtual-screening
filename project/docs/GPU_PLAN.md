# GPU 资源执行计划（2026-08-22 定稿）

> 背景：调研确认"变构位点是 2026 年文献公认的 DL/通用打分盲点"，我们的几何接触证据是正面解法。GPU 到位后按优先级执行；免费算力备选：Colab（T4/A100）、Kaggle 30h/周、AlphaFold Server 网页版、Boltz-2 自托管。

## P0 结构锚点（GPU 前可做）
- 下载 M4 PAM 结合 cryo-EM 结构：**MK-97 结合 M4**（Kaoullas 2026 bioRxiv 10.64898/2026.05.06.723386，PDB ID 待查证）、VU6016235 相关结构、ModelArchive ma-vohv7；
- 与 7TRQ（已有，M4+Gi1+iperoxo+VU0467154）对比口袋差异，建立多结构口袋锚点；
- 解决 M4R GaMD 集合簇（编号 629–908，280 残基截短构建）与 7TRQ 编号的映射——用结构比对（OBAlign）传递口袋中心/残基。

## P1 大样本几何基准（GPU 前可跑，pose 策略已备）
- `scripts/benchmark_geometry_m4r.py`：M4R 2314 PAM + 11.5 万 decoy，抽样 200+200 → Vina 出 pose（seed 42）→ 几何特征 AUC/EF 大样本立证；
- 产出：几何证据（0.87）在大样本上的正式统计 + 与文献 Vina PDB/BEmin/BEavg 对比表。

## P2 变构反向筛选（纯 CPU 可先行，第二创新点）
四步：① ligand 快筛（候选 vs 已知 M1–M5 PAM 的 Tanimoto 相似性，CPU 秒级）；② 五亚型逆对接（M1–M5 结构下载 → 候选对接各亚型变构口袋）；③ **接触证据替代打分**（各亚型变构 hub 对应残基，GPCRdb 注释）；④ 选择性归因（cryptic pocket 机制，Hollingsworth 2019；eLife 2023 M4 变构标志残基）。

## P3 GPU 工具链部署
| 工具 | 用途 | 优先级 |
|---|---|---|
| Uni-Dock（GPU） | 亿级库初筛（物理打分≈Vina，GPU 100–300×） | 高 |
| GNINA 1.1 | 重打分 + pose 生成（CNN 打分稳健） | 高 |
| AutoDock-GPU | 柔性侧链对接，配 GaMD 集合 | 中 |
| OpenMM/AMBER GPU + GaMD | 受体集合生成（M4 活性态弛豫） | 高 |
| SurfDock / DiffDock 2.0 | DL 姿态生成（初筛命中/生成分子重打分前） | 中 |
| RTMScore | 终筛重打分（图 Transformer） | 中 |
| Boltz-2 / AF3 server | 共折叠结构补全（须做口袋采样修复+弛豫，JCIM 2026 教训） | 低 |

## P4 3D 口袋条件生成闭环（核心创新落地）
- TargetDiff / FlowDock / TD3B 类：以 M4 变构口袋为条件生成 PAM 骨架；
- **几何证据条件化**（文献空白）：生成分子必须命中 hub 接触（TYR439 对应残基）+ 口袋接触原子对阈值；
- 生成 → 几何打分 → 再生成迭代闭环；候选经 BEmin 二次排序。

## P5 跨位点泛化验证（防过拟合）
- 几何接触证据在 M2R / M4R / A1R 至少三靶点重测（M2R/A1R 的 actives+decoys 同仓库可取）；
- 证明 0.87 不是单靶点过拟合。

## 里程碑
- GPU 到位当天：P0 + P3 部署
- +1 天：P1 大样本基准出数
- +2 天：P2 反向筛选出数
- +3~5 天：P4 生成闭环初版
- +7 天：P5 泛化验证 + 论文方法学定稿
