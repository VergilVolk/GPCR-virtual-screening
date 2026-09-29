# PACER-M4 当前状态

更新时间：2026-09-26

## 2026-09-26 CGDA：Full LIT-PCBA 未见靶点适配

- 已实现 Context-Gated DrugCLIP Adapter（CGDA）：由目标口袋与共晶参考配体共同产生门控权重，组合两个 rank-2 低秩专家；DrugCLIP 编码器保持冻结。
- 使用官方 15 个 LIT-PCBA 靶点、2,807,612 个候选、129 个真实口袋构象，执行严格 leave-one-entire-target-out；留出靶点活性标签只用于最终评价。
- 三种子困难负样本版相对 reference-ligand retrieval：宏 ROC-AUC `0.5693→0.5788`，配对 target-bootstrap 95% CI `[+0.00134,+0.01840]`；BEDROC80.5 `0.07103→0.07378`，CI `[+0.00045,+0.00504]`，15 靶点中 12 个改善。
- EF1% `6.06→6.31`，但仅少数靶点产生离散增益；EF5%下降，不声称所有截断点均改善。PR-AUC区间仍跨零。
- 四专家 rank-4 版本过拟合；二专家 rank-2 更稳。去除 retrieval、直接口袋到配体对齐、共享监督 LoRA、简单堆叠与多口袋聚合均未通过。
- 正式主张限定为：在参考配体辅助、严格未见靶点的回顾性筛选中，CGDA 显著改善整体 ROC-AUC 与 BEDROC；不是通用 SOTA、PAM 身份或药效预测。
- 完整方法、消融和复现入口见 `docs/DRUGCLIP_CGDA_RESULT.md`。
- 已完成冻结后的 DUD-E 外部 GPCR 压力测试（AA2AR、ADRB1、CXCR4、DRD3、M3R；91,311候选/1,342活性物）。相对 reference retrieval，CGDA 宏 ROC-AUC `0.7233→0.7437`、PR-AUC `0.2658→0.2762`、BEDROC80.5 `0.3637→0.3738`，三者的 target-bootstrap 95% CI 均为正且5/5靶点改善。
- 但原始 pocket DrugCLIP 在该 DUD-E 面板更强（ROC-AUC/BEDROC80.5 `0.8105/0.3936`）；CGDA 在 DRD3、M3R 退化。DUD-E 含生成 decoy，只作为迁移压力测试，不作为实验 inactive、PAM 或真实 HTS 命中证据。

## 2026-09-26 DrugCLIP 大样本双侧投影微调

- 已完成13靶点适配器与损失拆解：仅分子端适配基本无增益，口袋端适配贡献主要信号，双侧进一步提高；Target-DRO 退化至 ROC-AUC `0.559`，已否决。
- 三种子中，去除表示保持约束将宏 ROC-AUC/PR-AUC/BEDROC20 从 `0.623/0.247/0.300` 提高到 `0.637/0.258/0.311`；paired scaffold-cluster bootstrap 95% CI 分别为 `[+0.0079,+0.0197]`、`[+0.0062,+0.0160]`、`[+0.0010,+0.0210]`。EF1%/EF5% 差异仍跨零。
- 分子锚定显著优于口袋锚定，支持“分子语义可保留、口袋空间需允许域重排”；当前最高总体迁移模型仍为完全不锚定的 `dual_no_preserve`。
- 已冻结三份 `dual_no_preserve` 全数据部署权重并核验可加载；完整消融与 SHA256 见 `docs/DRUGCLIP_TARGET_TRANSFER_ABLATION.md`。
- 完整末层投影微调使用 131,328 个参数（LoRA 的 12.8 倍），同种子 ROC-AUC/PR-AUC/BEDROC20/EF5% 仅 `0.591/0.221/0.266/2.00`，低于 LoRA 的 `0.626/0.251/0.307/2.34`；支持参数效率，但该 full-last 基线尚未做大范围调参。

- 已直接微调 DrugCLIP 内部 `mol_project` 与 `pocket_project` 的 rank-8 LoRA，共 10,240 个可训练参数；不是 ECFP 外接分类器。
- 使用 28,618 个 B2AR/CCR2/M2R/M4R active/decoy 配对，联合平衡 BCE、活性分子靶点检索、表示保持和蒸馏；全局 Murcko-scaffold 五折 OOF 无骨架重叠。
- 域内完整多任务模型将官方 DrugCLIP 的宏 ROC-AUC/PR-AUC/BEDROC20 从 `0.629/0.195/0.245` 提高到 `0.951/0.864/0.912`；全随机标签对照为 `0.493/0.116/0.128`。
- 真实 muscarinic triplet 与随机 triplet 几乎无差异，因此当前 triplet 分支未证明有增益；推荐模型为 `BCE + retrieval`。
- 三种子严格 target-LOSO 中，微调模型的宏 ROC-AUC/PR-AUC/BEDROC20/EF5% 为 `0.655/0.218/0.313/3.57`，官方 DrugCLIP 为 `0.629/0.195/0.245/2.39`；相对官方的 scaffold-bootstrap 95% CI 均为正（EF1% 除外）。
- pooled ECFP4 的 ROC-AUC 略高（`0.663`），但 PR-AUC/BEDROC20/EF5% 仅 `0.193/0.229/2.32`，低于 DrugCLIP 微调，形成当前最有价值的未见靶点早期富集证据。
- 冻结 M4 外部药理面板未改善：Acadia AUC 持平，Monash AUC 下降，效力排序无可靠信号；不得将该模型称为 PAM 功能或效力预测器。
- GaMD 十构象聚合和严格嵌套 docking 融合均未稳定超过单独 target-LOSO DrugCLIP，已保留为负结果。
- 已冻结 `bce_retrieval` 与 `full_multitask` 两套部署投影权重；详细结果见 `docs/DRUGCLIP_LARGE_MULTITASK_RESULT.md`。
- 9 个独立 LIT-PCBA 靶点、8,906 个配对的冻结外部测试已完成：完整多任务模型宏 ROC-AUC `0.541`，官方为 `0.532`，差值 95% CI `[+0.0007,+0.0168]`；PR-AUC、BEDROC 和早期富集区间均跨零。模型未出现明显灾难性遗忘，但尚无通用早期富集提升。
- 已扩大为 4 个 GPCR + 9 个 LIT-PCBA 靶点、37,524 配对的统一 13-target LOSO。三种子 DrugCLIP BCE+retrieval 的宏 ROC-AUC/PR-AUC/BEDROC20/EF5% 为 `0.623/0.247/0.300/2.28`，官方为 `0.562/0.201/0.241/1.78`，严格 ECFP4 为 `0.568/0.201/0.239/1.71`。
- 13-target 模型相对官方与 ECFP4 的 ROC-AUC、PR-AUC、BEDROC20、EF5% scaffold-bootstrap 95% CI 均为正；EF1% 仍未通过。13 个靶点中 ROC-AUC/PR-AUC/BEDROC20 分别有 `10/11/10` 个改善，CCR2 仍失败。
- 已按 `20260925/26/27` 三个 seed 冻结全数据部署权重，推理时建议三模型平均；这构成目前最强的方法学结果，但仍是回顾性 13 靶点验证，不称通用 SOTA。

## 2026-09-26 DrugCLIP 非对称分子适配器

- 冻结 DrugCLIP 三维口袋 embedding，以同分子跨亚型 triplet 训练 ECFP→口袋空间的轻量分子适配器。
- 五亚型严格 target-LOSO 同时满足待测口袋未参与监督、训练/测试分子重叠为 0；ECFP3 适配器宏准确率 `73.0%`，高于 DrugCLIP triplet 的 `55.4%`。
- 完整 M4 折为 `87.5%` 对 `53.3%`，分子成簇 bootstrap 95% CI 为 `+28.2` 至 `+39.9` 个百分点。
- 指纹形式只在非 M4 校准集选择；一次性冻结 M4 面板上，ECFP3 为 `68.6%`、ECFP2 为 `66.7%`、DrugCLIP triplet 为 `66.0%`。ECFP3 相对 DrugCLIP 的 CI 跨零，仍需更大独立集。
- ECFP+DrugCLIP 残差、拼接适配器均未在完整 LOSO 与冻结 M4 上同时改善，暂不作为主模型。
- ECFP3 直接跨 B2AR/CCR2/M2R/M4R 家族迁移失败；DrugCLIP 路由加局部适配器的宏 Recall@1 增益不显著，两特征门控外迁后有害。家族内适配器不能冒充通用 GPCR 模型。
- 当前主张限定为“毒蕈碱家族内未见分子、未见亚型的活性排序迁移”；不等同于 PAM 身份、效力或通用 GPCR SOTA。
- 已对 200 个 PACER 候选完成 M1–M5 选择性审计；24 个候选为 M4 top-1。原 12 个多证据 shortlist 中仅 `PACER0054` 同时保留较强 M4 选择性支持（选择性 rank 6/200，margin `0.0854`，三种子 SD `0.0064`），现列为优先验证假设而非 PAM 结论。
- 完整结果与负结果见 `docs/DRUGCLIP_ASYMMETRIC_ADAPTER_RESULT.md`。

## 2026-09-25 DrugCLIP 四 GPCR 系统微调

- 已在 B2AR、CCR2、M2R、M4R 的 2,500 个可编码别构调节剂上完成五折 Murcko-scaffold OOF；四个口袋来自公开 GaMD benchmark 的位点定义和代表构象。
- 冻结两个 Uni-Mol encoder，仅在分子、口袋末端投影训练双侧 rank-4 LoRA，共 5,120 个参数。
- 靶点平衡 CE 将官方 DrugCLIP 的宏 Recall@1 从 `0.413` 提高到 `0.925`，pair ROC-AUC 从 `0.743` 提高到 `0.9965`；hardest-pocket triplet 未优于 CE，已否决为主模型。
- 强二维 ECFP4-logistic 的宏 Recall@1 为 `0.963`、pair ROC-AUC 为 `0.9996`，仍高于 DrugCLIP 微调；因此不能声称通用 SOTA。
- 随机标签模型的总体 Recall@1 可虚高到 `0.840`，但宏 Recall@1 仅 `0.235`，证明不平衡条件下必须报告宏平均和逐靶点指标。
- 低相似审计中 Tanimoto `<0.30` 仅 18 个分子，CE 未优于官方或 ECFP；当前大幅提升主要是已知 GPCR 化学域特化。
- 严格留一整靶点时，CE 的宏 Recall@1 为 `0.370`，低于官方模型 `0.413`；专项 adapter 不具备已证明的未见 GPCR 迁移能力。
- 已冻结三种子 CE 部署权重。部署规则：已知四靶点域可使用 CE adapter；未见 GPCR 保留官方 DrugCLIP 或已验证的跨亚型关系模型。
- 详细证据与边界见 `docs/DRUGCLIP_GPCR_SYSTEMATIC_FINETUNING_RESULT.md`。

## 2026-09-24 DrugCLIP GPCR triplet 结论

- 已完成官方 DrugCLIP CPU 复现，并开放 512 维投影前表示；微调对象为 DrugCLIP 内部 `pocket_project` / `mol_project`，不是把冻结 embedding 接一个冒充微调的外部分类器。
- 同分子跨亚型 triplet 在 33 个未见分子、54 个 M1/M2/M3/M5 配对上将官方 DrugCLIP 准确率从 `59.3%` 提到 `79.6%`，相对随机方向对照的分子成簇 95% CI 为 `+5.8` 至 `+32.1` 个百分点。
- 强二维基线 ECFP4-Ridge 在同一 54 对上达到 `98.1%`；因此该集合是化学系列内插问题，不能用于宣称 DrugCLIP SOTA。
- 在训练完全不含 M4 口袋、测试分子完全未见的冻结零样本测试中，官方 DrugCLIP 为 `57.5%`，GPCR-triplet 为 `66.0%`；提升的分子成簇 95% CI 为 `+3.1` 至 `+14.3` 个百分点，分子置换检验单侧 `p=0.018`。这是当前最强阳性证据。
- 五亚型严格 leave-one-target-out 上，官方模型宏平均 `53.6%`，全投影 triplet `55.4%`，rank-4 pocket LoRA `55.6%`；靶点间差异明显，尚不支持通用 GPCR SOTA。
- 外部 M4 PAM/inactive 数据上，GPCR-triplet 在 Acadia/Monash 的 AUC 分别为 `0.930/0.788`，高于官方模型 `0.837/0.576`，但两集合合计仅 5 个阴性，差值置信区间跨零。
- Suven、VU6025733 和 US20260055116 的效力排序未通过；DrugCLIP 结构分数不能作为 PAM EC50 预测器。
- M4 全投影功能微调出现跨来源退化；rank-4 molecule LoRA 系列留出 AUC `0.630`，低于随机配对对照 `0.641`，该路线已否决。
- 当前可成立的创新边界：GPCR triplet 可修复一部分未见 M4 口袋的零样本结构迁移，但 PAM 身份、效力和内在激动必须由独立功能药理/动态模块承担。完整报告见 `docs/DRUGCLIP_MUSCARINIC_TRIPLET_RESULT.md`。

## 2026-09-19 动态表示对照结果

- 已完成公开 M4–iperoxo–MK-97 六 replica 的逐 replica 留出表示 benchmark；
- VAMP-5 的留出 lagged R² 为 `0.750`，相对 PCA-10 提高 `0.334`，6/6 folds 改善，95% CI `[0.287, 0.378]`；
- VAMP-5 与 tICA-5 基本相同（平均差 `0.001`，95% CI 跨零），因此不能宣称 VAMP 优于简单慢模态基线；
- tICA/VAMP 的跨 replica 状态 JSD 高于 PCA，说明慢状态仍明显受 replica 采样影响；
- 该结果仅验证单一 PAM 体系的动力学表示，不构成功能性 PAM 预测证据；
- OneProt-MD 已完成源码可执行性审计，暂列为可选冻结基线，不作为既定主干模型；
- 详细决策见 `docs/PACER_DC_DYNAMIC_ENCODER_DECISION.md`，原始结果见 `results/pacer_dynamic_representation_bakeoff_v01/`。

## 2026-09-19 四上下文体系构建进展

- 六个无溶剂参考体系的拓扑与受限最小化 QC 已全部通过；
- 已使用 OPM 7TRS 膜取向（270 个受体 CA 对齐 RMSD `0.365 Å`）构建统一 POPC/水/0.15 M 离子底座；
- ACh-only、apo、LY2119620 有/无 ACh、compound-110 有/无 ACh 六个周期体系已在同一次共享底座运行中构建完成；
- 每个体系约 `227,000` 原子，六套体系能量均有限且最小化后下降；
- 这只证明体系可构建，不是平衡完成或 PAM 功能证据；下一闸门为短 NPT 数值稳定性测试。
- ACh-only 首轮 1 ps 及分阶段 0.5 fs 升温测试均在早期出现坐标 NaN，动力学闸门当前**未通过**；不得启动生产轨迹，正在按最大受力原子和约束偏差定位几何问题。
- 已定位 NaN 根因：旧蛋白氢中 Gα-THR51 HG1 与 CG2 几乎重叠；改为保留重原子、按当前力场重建全部蛋白氢，并只约束蛋白重原子。最大力由约 `1.4e12` 降至 `2.8e3 kJ mol-1 nm-1`。
- 修复后的 ACh-only 完成 50→100→200→300 K、0.5 fs、0.2 ps 数值稳定性烟雾测试，能量和坐标均有限；这只通过单体系数值闸门，六体系统一重建与更长平衡仍待完成。
- 六体系已按最终“重建蛋白氢、仅约束重原子”的协议完成一次统一构建。ACh-only、apo 和 compound-110 两上下文在初始最小化后正常；LY2119620 两上下文在100步最小化时仍为高能，暂不通过。
- 对 LY2119620 无 ACh 体系追加1000步最小化后，势能降至 `-3.19e6 kJ/mol`，最大力约 `3.12e3 kJ mol-1 nm-1`，最大受力原子为水氧且约束误差约 `1.05e-8 nm`。这表明LY体系并非拓扑失败，但需要更充分的统一预平衡，不能直接沿用100步最小化结果。
- 已增加可复现的蛋白氢几何闸门：固定 seed `1701`，非键 H–重原子最小距离必须不低于 `0.06 nm`；最终统一底座实际为 `0.160 nm`。
- 最终统一底座上的6个参考上下文均已完成深度最小化及50→300 K、0.5 fs、0.2 ps NVT烟雾测试，6/6能量与坐标有限；深度最小化后最大力约 `2.70–2.76e3 kJ mol-1 nm-1`，最大受力原子均为溶剂水氧。
- 该结果只关闭了数值稳定性阻塞；尚未达到平衡、独立 replica 或功能差分证据要求。
- 六体系均完成首轮 `0.25 ps NVT + 0.25 ps NPT` 压力耦合试验，能量/坐标有限、盒体积变化为 `-0.08%` 至 `-0.72%`；但终温仅约 `208–211 K`，因此数值稳定为 `6/6`、热平衡闸门为 `0/6`，不得提取功能结论。
- 升级后的 apo 协议（`2 ps NVT + 1 ps NPT`、`1 fs`）达到 `289.2 K`，盒体积变化 `-0.73%` 并通过热平衡协议闸门；其余五体系正按相同冻结协议运行。
- 已增加短平衡证据隔离：QC 轨迹必须标记为 `equilibration_pilot`，不能满足 PACER-DC 的正式 `trajectory` 证据要求。
- 已实现蛋白重原子约束分级释放模块（`500→100→10→0 kJ mol-1 nm-2`）；只有零约束末态通过温度/体积/有限性检查后，才输出生产轨迹起点。
- 升级热平衡协议现已在六体系全部通过：终温 `289.2–290.3 K`，盒体积变化 `-0.73%` 至 `-1.16%`，数值稳定/热化联合闸门为 `6/6`。
- 六体系均完成 `500→100→10→0 kJ mol-1 nm-2` 分级去约束；零约束末态终温 `297.6–300.1 K`，`production_start_eligible=6/6`。这表示已获得生产 MD 起点，不表示动力学收敛。
- 已补齐并实测 `MDAnalysis 2.10.0` 与 `pandas 3.0.6`；修复了膜体系重编号（原始受体链 R → 膜体系链 E）及十六进制 `CONECT` 记录导致的轨迹解析问题。
- 六体系 QC 端点提取完成，共 13 条 replica 级记录；全部强制标记为 `equilibration_pilot` 且 `eligible_for_functional_scoring=false`，不得用于 PAM 功能、性能提升或生物机制结论。
- 首次端点审计发现周期边界未处理导致 ACh RMSD 约 `200 Å`、口袋接触为零的伪异常；已加入最邻近周期镜像和最小镜像接触距离并重跑。修复后 ACh RMSD 为 `1.75–1.93 Å`、配体口袋接触覆盖为 `0.775–0.908`，仅说明提取链与短时位姿保持表现合理，仍不是功能或收敛证据。
- 已冻结当前最强可汇报机制结果：公开三对 replica 中，LY2119620 使 iperoxo 平均 RMSD 降低 `0.651 Å`，分层 block-bootstrap 95% CI `[0.527, 0.762] Å`，三对方向一致且 100 ns 后效应仍为 `0.656 Å`；外部 probe-matched 静态结构增量为 `0.539`。该结果支持一个已知 PAM 的正构配体稳定化机制，不支持通用 PAM 分类或新候选确认。
- 已实现可断点续跑的无约束生产 MD 程序，并通过 apo CPU `0.2 ps` 端到端测试；已生成六体系×三 replica×`5 ns` 的 18-job GPU Phase-1 清单和 Slurm job array。Phase-1 只用于轨迹 QC，未达到功能推断门槛。
- 生产程序已再次通过同一任务 `0.2→0.4 ps` 检查点续跑测试，并增加跨GPU节点不兼容时的 XML state 便携恢复路径。
- 已实现18任务统一 Phase-1 审计：温度、密度、体积漂移、有限性、帧数和匹配上下文完整性全部通过后，才允许延长到 `20 ns`。
- 已实现生产轨迹端点提取闭环。`5 ns`结果只能标记为 `production_pilot`，评分器会拒绝；只有显式通过正式收敛审核后，才允许标记为 `trajectory` 并进入四上下文计分。

## 2026-09-15 执行更新

- 论文级 baseline、切分、统计检验和晋级标准已冻结于 `docs/PACER_DC_PUBLICATION_PROTOCOL.md`；
- WSL2 原生 Linux 环境已建立，OpenMM、OpenFF、AmberTools、ParmEd 和 MDAnalysis 审计通过；
- 修复了 OpenFF Toolkit 0.16.x 与 Interchange 0.4.5 的真实兼容性错误，环境已固定到通过实测的版本组合；
- LY2119620、compound-110、ACh 和 iperoxo 四个参考配体的 OpenFF 参数化烟雾测试为 4/4 通过；
- 42-run v2 清单的六个配体均已有明确 canonical SMILES；
- 已在统一 7TRS 底座上构建 ACh-only、apo 和 LY2119620 两个上下文的无严重碰撞起始结构；
- compound-110 从 7V6A 直接转移到 7TRS 后出现 0.32 Å 严重碰撞，已按 QC 规则阻止进入 MD，需重新对接或采用受控受体松弛后再构建。
- compound-110 已完成 7TRS 口袋重对接：Vina `-7.838 kcal/mol`、预声明口袋覆盖率 `0.846`、接触 11/13 个口袋残基；该位姿仅作为无碰撞起始构象提案，仍需拓扑映射和最小化 QC。

## 当前主线

开发 PACER-DC：用 M4 在有、无 ACh 条件下的匹配动态差分，区分 ACh 依赖的 PAM 协同作用与候选物自身的内在激动作用。

## 已完成

- 四上下文、三 replica 的 30-run 试验清单已冻结；
- PACER-DC 多维证据计分器已实现；
- 缺失上下文禁止排名，静态证据禁止冒充完整动态证据；
- 第一阶段自动闸门和测试已实现；
- 主流程已从加权总分改为结合门控、动态差分、baseline/消融和可拒判输出。
- 已锁定同来源近邻对 CM00717（PAM）/CM00734（实验非活性，ECFP4 Tanimoto 约 0.955），用于困难负样本验证；其结合能力仍待检验。

## 当前真实闸门

| 闸门 | 状态 | 原因 |
|---|---:|---|
| 试验清单结构 | 通过 | 四上下文和 replica 对齐正确 |
| 生产 MD 环境 | 未通过 | Windows 无 AmberTools；正在建立 WSL2 Linux 隔离环境 |
| PACER-DC 完整计分 | 未通过 | 当前完整四上下文候选数为 0 |
| Triplet/监督轨迹学习 | 未通过 | 仅一个 PAM 动态正对照，且缺少无功能结合对照 |
| 功能性 PAM 结论 | 未通过 | 尚无前瞻功能实验 |

## 下一执行点

1. 完成 WSL2 中 OpenFF/AmberTools/OpenMM 环境；
2. 对已知 PAM、compound-110、CM00717/CM00734 和共享对照做参数化与最小化 QC；
3. 先跑短轨迹检查体系稳定性，QC 通过后才进入 3 × 100 ns；
4. 先完成简单动态差分 baseline，再比较 VAMP/SPIB；
5. 数据闸门通过后才做 Triplet 微调。

自动报告：`project/results/pacer_dc_phase1_audit/REPORT.md`。
# 当前协作任务入口

当前队员交付以仓库根目录 `START_HERE_ONEPROT_PACER_DC.md` 为准：
**OneProt-MD 复现 → M4 四上下文 embedding → PACER-DC Adapter/Triplet 训练与 baseline 比较**。
仅完成重对接、体系构建或短平衡不等于完成当前算法任务。
