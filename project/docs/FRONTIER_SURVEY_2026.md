# 前沿调研：GPCR 反向虚拟筛选与变构位点 VS 算法（2023–2026）

> 目的：回答"最新颖的 GPCR 反向虚拟筛选怎么做算法？与我们之前的双线思路（算法线 + 机制线）差多少？"并为即将到位的 GPU 资源制定执行计划。
> 来源：3 路并行 web 调研（2026-08-22），全部结论附文献来源；与 docs 里 4 篇论文的分析（见 LITERATURE_REUSE_ANALYSIS.md）互补。

---

## 一、三个关键结论（先看这个）

1. **变构位点是 2026 年文献公认的"DL 对接与通用打分系统性盲点"**。至少 3 篇独立工作（Allosteric Blind Spot、Decoding the Allosteric Paradox、GPCR AM 对接基准）一致证明：DiffDock/SurfDock/Boltz-2 类模型对变构位点早期富集弱，且传统打分同样弱。**我们的几何接触证据（AUC 0.87 vs Vina 0.33）恰好正面回应了这个盲点**——这不是自嗨，是文献留出的位置。
2. **两个可占领的空白**：
   - ① 位点特异"几何感知评分 + 生成式条件化闭环 + 自建 M4 变构基准"——2024–2026 无公开工作（生成式变构设计集中在激酶，TD3B 等未绑定 GPCR PAM）；
   - ② 变构反向虚拟筛选（给定 PAM 候选预测 GPCR 亚型选择性 M4 vs M1/M2/M3/M5）——无成熟专用工具。
3. **GPU 工具链已明确**：Uni-Dock（亿级/天）、GNINA、AutoDock-GPU、GaMD(OpenMM/AMBER)、SurfDock、DiffDock 2.0、Boltz-2/AF3 共折叠。VSDS-VD 基准（Nat. Mach. Intell. 2025）的教训：**AI 对接在 VS 富集上仍不及传统打分，推荐"AI 出姿态 + 物理/几何打分验证"混合**——与我们"打分失效→几何证据"的思路天然一致。

---

## 二、与"我们之前的双线思路"差距分析

### 算法线（QSAR + 生成 + docking）
| 维度 | 我们现状 | 2023–2026 前沿 | 差距 |
|---|---|---|---|
| 对接 | Vina 1.2.5（CPU，打分已被证伪） | Uni-Dock / Vina-GPU（亿级/天）、DiffDock 2.0 / SurfDock（DL 姿态）、RTMScore（重打分） | 大：缺 GPU 工具；但"打分不可信"反而是我们的论据 |
| 受体柔性 | 单结构 7TRQ | GaMD 集合对接（Miao lab，M4R 数据已开源） | 中：已拿到 M4R 10 构象簇 + PMF，可零成本接入 |
| 生成 | SMILES 级 LSTM/GPT/BRICS | 3D 口袋条件生成（TargetDiff/DecompDiff/FlowDock）、TD3B 变构生成 | 大：生成维度（2D→3D）差距，GPU 后可补 |
| 基准 | n=20 手工集 | M4R 2314 PAM + 11.5 万 decoy（已复现 Table 4） | 已追上：大样本基准在手 |

### 机制线（变构-正位偶联，8.03 Å / TYR439 hub）
- 文献对照：Miao 2016 的"前庭双低能位点 + 几何反应坐标"与我们同构；Nguyen 2025 的"浅袋 + 膜环境"解释了为什么打分失效。
- **差距不在方向上，在叙事闭环上**：文献证明"通用打分对变构失效"（盲点），我们只有"单个靶点的几何证据 AUC 0.87"。**升级为"跨位点泛化"（至少 M2/M4/A1 三靶）才能证明不是过拟合**——这是机制线的首要补强点。

### 反向虚拟筛选（全新空白）
- 定义确认：IVS/RVS = 以配体为查询，反推靶点 / 亚型选择性 / off-target 谱。
- 三路线：ligand 相似性反查（SEA/SwissTargetPrediction/PASS，CPU 秒级）、逆对接（TarFisDock/INVDOCK 类）、DL 靶点预测（DeepPurpose 等）。
- **未发现"给定 PAM 预测 GPCR 亚型/变构位点"的成熟工具** → 我们的四步落地：ligand 快筛 → 五亚型（M1–M5）变构口袋逆对接 → **接触证据替代打分**（与我们 QSAR 失败→几何证据的观察同构！）→ GPCRdb/ASD2023 选择性归因（cryptic pocket 机制，Hollingsworth 2019；eLife 2023 M4 变构标志残基）。

---

## 三、推荐组合创新（文献无人做过）

**"位点特异几何感知的变构筛选闭环"**：
1. **几何接触证据**（已验证 AUC 0.87）→ 升级为位点特异评分函数，在大样本 M4R 基准（2314 PAM/11.5 万 decoy）上正式立证（已备好 pipeline，待 pose 策略确认）；
2. **GaMD 集合 + BEmin 重排序**（Miao 2026 协议，数据已复现）作为受体柔性层；
3. **几何证据做生成条件化**：约束 3D 生成（TargetDiff/FlowDock/TD3B 类）必须命中 hub 接触，形成"生成→打分→再生成"闭环（文献空白）；
4. **自建 M4 变构富集基准**（文献尚无，Miao 2026 只发了打分数据）；
5. **变构反向筛选**：候选 → M1–M5 亚型选择性预测（第二创新点，纯 CPU 可先行）。

风险提示：AUC 0.87 需跨位点（M2/M4/A1）泛化验证；生成候选需类药性/合成可行性过滤。

---

## 四、GPU 到位后执行计划（优先级排序）

| 优先级 | 任务 | 工具 | 产出 |
|---|---|---|---|
| P0 | 下载 M4 PAM 结合结构（MK-97、VU6016235，ModelArchive ma-vohv7）做口袋锚点 | 网络 + PDB | 变构口袋定义（比 7TRQ 更权威） |
| P1 | 大样本几何基准正式跑（400–800 pose） | 现有 CPU 管线（Vina 仅出 pose） | 几何证据 AUC/EF 大样本立证 |
| P2 | 变构反向筛选四步（CPU 先行，不依赖 GPU） | SEA/SwissTargetPrediction + 逆对接 + 接触证据 | 亚型选择性创新点 |
| P3 | GPU 工具链部署：Uni-Dock / GNINA / AutoDock-GPU / OpenMM-GaMD | GPU | 亿级库初筛 + GaMD 集合 |
| P4 | 3D 口袋条件生成（TargetDiff 类）微调/条件化到 hub 接触 | GPU | 生成闭环 |
| P5 | 跨位点泛化验证（M2/M4/A1） | 全链 | 防过拟合证据 |

免费算力备选：Colab（T4/A100）、Kaggle 30h/周、AlphaFold Server 网页版、Boltz-2 自托管。

---

## 五、关键文献速查

- Nguyen 2025 PNAS：A1R PAM 膜感知对接（[链接](https://www.pnas.org/doi/10.1073/pnas.2421687122)）
- Thompson & Miao 2026 bioRxiv：GPCR AM 对接基准（[链接](https://www.biorxiv.org/content/10.64898/2026.08.12.744492v1)），数据仓库 [tylerdt1/gpcr-am-ensemble-docking](https://github.com/tylerdt1/gpcr-am-ensemble-docking)
- Allosteric Blind Spot 2026（[ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2666386426003541)）
- Allosteric Paradox 2026（[bioRxiv](https://www.biorxiv.org/content/10.64898/2026.02.24.707829v1)）
- DiffDock-L 变构微调 2026（[JCIS](https://pubs.acs.org/jcisd8/article/66/6/3036/5080410/)）
- Hollingsworth 2019：cryptic pocket 与 M1–M5 选择性（[PMC6650467](https://pmc.ncbi.nlm.nih.gov/articles/PMC6650467/)）
- eLife 2023：M4 变构药理学标志残基（[PMID 37248726](https://pubmed.ncbi.nlm.nih.gov/37248726/)）
- VSDS-VD：AI 对接 VS 富集评测（[zenodo](https://zenodo.org/records/14874127)）
- TD3B：转换引导离散扩散变构生成（[arXiv:2605.09810](https://huggingface.co/papers/2605.09810)）
