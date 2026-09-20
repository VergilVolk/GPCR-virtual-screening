# M4R GaMD ensemble docking：复现边界与执行协议

## 为什么做这一层

静态 docking 只能在一个冻结的受体微状态上搜索配体姿势。M4 PAM 口袋处于胞外前庭，形状和侧链构象会随正构配体、G 蛋白和 PAM 改变。因此本层检验的是：**对已采样的多个受体微状态进行 docking，并用微状态自由能重加权，能否提高 PAM 与 inactive 的结构富集。**

它不能直接证明候选物具有正向协同性、功能效力或特定正构探针依赖性；这些是后续功能模型和湿实验的责任。

## 对方仓库实际做法

来源：`tylerdt1/gpcr-am-ensemble-docking` 和配套 Miao et al. 2026 预印本。

冻结仓库提交：`44798c841ee77230b1f89fe41074e41b458c7070`（2026-08-19）。M4R ensemble PDB SHA256：`14F61D066E328E4D4E525B0E22BE240A7D830B163ECA2D6145DE43FB49EE77D4`；PMF SHA256：`F8F3294A809775F865B92A127340B28320869EF2296BCD676D6191317AC9278E`。

1. 对 ACh–M4R–G protein–MK-97 holo 复合物做 GaMD。
2. 用 PAM 口袋、ACh（907）和 MK-97（908）的重原子联合聚类，得到 10 个代表构象。
3. 对每个构象的别构口袋分别 docking。
4. 以 GaMD 重加权 PMF 校正每个构象的 docking score：
   `BE_i = docking_score_i + PMF_i`。
5. 两种排序：`BEmin = min(BE_i)`；`BEavg = mean(BE_i)`。
6. 用已知变构调节剂和 property-matched decoys 测 ROC-AUC、logAUC、EF 和 EF'。

正式 Vina 参数：Vina 1.2.7；30 Å search box；exhaustiveness 8；最多 9 poses；energy range 3 kcal/mol。蛋白保留 GPCR 与正构配体，水删除；配体按 pH 7 标准化后由 RDKit/Meeko 制备。

## 我们的两套验证，不能混为一谈

### A. 论文原始 active/decoy benchmark（代码复算）

- PDB 单结构 ROC-AUC：0.7253。
- GaMD BEmin ROC-AUC：0.7338。
- GaMD BEavg ROC-AUC：0.6958。
- 但 Vina 的 top-0.5% EF：PDB 0.778、BEmin 0.349、BEavg 0.959，均未形成可靠早期富集。

因此，Vina ensemble 的全局 AUC 有轻微增益，不代表前排候选更可靠；原论文对 M4R 的强结果来自商业 Glide BEmin，而不是 Vina。

### B. 我们的 assay-specific functional benchmark

- 主分析只用 27 个 **A-tier primary-source-confirmed inactive** + 27 个性质匹配 PAM。
- 66 + 66 全集包含 database-text “Not Active” 阴性，只保留作敏感性分析，不进入主结论。
- 不使用“未测过所以当阴性”的人工 decoy。
- 这是更贴近 PAM/inactive 区分的压力测试，但样本更小，必须报告 bootstrap 置信区间。

早期工程 pilot（22 Å、exhaustiveness 4、1 pose）得到 cluster-00 AUC 0.548。该值因参数不符合论文正式协议，已冻结为 **pilot only**，不得用于最终方法比较。正式结果只接受 Vina 1.2.7 / 30 Å / exhaustiveness 8 / 9 poses。

### C. 标签语义审计（关键发现）

529 个 assay-specific 标签中有 300 个与 Miao/ASD 库精确 SMILES 重叠。33 个重叠的 functional negative 全部在 ASD benchmark 中至少被标为一次 active；6 个结构同时出现在官方 active 与 decoy 行。33 个冲突阴性均来自 database-text tier，而非 A-tier primary-confirmed tier。

这说明 ASD 的 broad allosteric-modulator membership 与“指定 ACh 探针下的人 M4 功能 PAM”不是同一 endpoint。前者能评估文献库 enrichment，不能直接作为 PAM efficacy ground truth。完整交叉表见 `results/m4_gamd_ensemble/miao_vs_functional_label_crosswalk.csv`。

## 当前实现

- `prepare_m4_gamd_ensemble.py`：拆分 10 簇、保留 ACh、移除 MK-97、制备 PDBQT，并以每簇共模拟 MK-97 的位置映射 30 Å 口袋。
- `dock_m4_gamd_ensemble.py`：molscrub 0.1.1 在 pH 7 枚举质子化/互变异构状态，逐状态 CPU 并行 docking，每个分子/构象保留最低能状态，并输出 rank-1 姿势和原生口袋覆盖。
- `analyze_m4_gamd_ensemble.py`：按官方公式计算 BEmin/BEavg，并与 cluster-00 做成对分层 bootstrap。

## 复现等级与剩余差异

- **完全一致**：公开 10 簇坐标、PMF、ACh 保留、Vina 1.2.7、30 Å、exhaustiveness 8、9 poses、energy range 3、molscrub 0.1.1 pH 7 状态枚举、BEmin/BEavg 公式。
- **工程等价但非同一软件**：论文用 AutoDockTools 1.5.7 制备受体；本机用 Meeko 制备刚性 PDBQT。ACh 在所得 PDBQT 中净电荷为 +1.000。
- **数据集有意不同**：论文用 ASD actives + property-matched 未知活性 decoys；我们另加 66 PAM + 66 实验 inactive 的严格 benchmark。
- **仓库限制**：对方仓库提供最终 ensemble、逐簇 docking scores、化合物库和分析脚本，但没有完整的原始 receptor-preparation/docking job 文件。因此不能声称逐字节复现，只能做参数级复现并记录软件替代。

## 预注册判据

只有在严格 experimental-inactive benchmark 上，ensemble 相对 cluster-00 的成对 bootstrap `ΔAUC` 95% CI 不跨 0，且 EF10% 不恶化，才能声称“动态集合带来验证增益”。否则将其作为候选构象稳健性/状态选择性证据，不作为 PAM 分类器。

## 最终结果（2026-08-29）

- A-tier 主集：54 分子（27 PAM / 27 primary-confirmed inactive），十簇全部完成，零 docking 失败。
- cluster-00：ROC-AUC 0.499、PR-AUC 0.511、EF10%=1.00。
- PMF-BEmin：ROC-AUC 0.499、PR-AUC 0.511、EF10%=1.00；ΔAUC=0。
- PMF-BEavg：ROC-AUC 0.457、PR-AUC 0.460、EF10%=0.33；ΔAUC=-0.044，95% CI `[-0.147,+0.045]`。
- 预注册晋级条件未满足。GaMD ensemble 降级为姿势稳健性/状态敏感性证据，不作为 PAM 分类器。
- 24 个候选、300 个 pH7 状态/构象任务全部完成且零失败；15 个处于动态四目标 Pareto 前沿。所有候选的 PMF-BEmin 均由 cluster-00 决定。

主统计见 `results/m4_gamd_ensemble/VALIDATION_REPORT.md`；候选证据见 `results/m4_gamd_ensemble/candidate_ensemble_evidence.csv`。
