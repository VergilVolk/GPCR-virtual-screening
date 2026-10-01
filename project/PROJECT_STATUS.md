# PACER-M4 当前状态

更新时间：2026-10-01

## 2026-10-01 DrugCLIP 结合筛选线闭环

- DrugCLIP 内部双侧投影微调已迁移到官方 Science-2026 表征，并以严格 target-LOSO、Murcko scaffold 隔离和三随机种子评估。
- 13-target family-augmented ep80 三种子集成宏 ROC-AUC 为 `0.6414`；对应官方 DrugCLIP 为 `0.5442`、ECFP4 logistic 为 `0.5681`。原始预测独立复算已通过。
- 20-target 全覆盖时 family augmentation 不增益，冻结路由为：13-target 使用 family augmentation，20-target 使用纯 ep80；不得跨协议择优拼接。
- EF1% 仍未稳定优于强二维基线，因此正式主张限定为整体排序与部分早期识别改善，不称通用虚拟筛选 SOTA。
- 已完成 28,519 分子库和 PACER 200 候选重评分。新旧权重排序近乎独立，当前作为双模型第二意见；双优交集进入 PACER-DC 功能复核，不直接称 PAM。
- 结合线完整交接见 `docs/DRUGCLIP_FINETUNING_PIPELINE_INTEGRATION.md`；功能线闭环与主张边界见下节。

## 2026-10-01 PACER-FKG v02 Stage A + Stage B 闭环

### 当前冻结结论

- PACER-FKG v02 已完成从历史 600 ns long-MD、LY2119620 20 ns matched closure（Stage A）到 CM00734 20 ns hard-negative test（Stage B）的同一冻结表示链验证。
- 三批分析统一使用冻结的优化 encoder：原始 Geom2Vec/ViSNet 预训练 checkpoint 不微调，但使用冻结的 `C1-BS256` 中间层/残基 readout，而不是原版 Geom2Vec final embedding。
- 冻结 PACER-FKG v02 数值锚点：`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`。
- Stage A 的 matched 20 ns 结果中，主要信号位于 `STATE_MOTION / Delta_INT / compound110_extension`：
  - compound110 reference：R1/R3 direction cosine = `+0.583656`
  - LY2119620：R1/R3 direction cosine = `+0.634304`
- Stage B 对实验 inactive hard negative `CM00734` 使用相同冻结 encoder、normalization、RFF、graph、region 和 contrast，无任何 outcome-driven refit。相同主要信号下：
  - CM00734：R1/R3 direction cosine = `-0.326839`
  - CM00734 / LY2119620 的 Delta_INT magnitude 比例为 R1 `0.644`、R2 `0.798`、R3 `0.525`
- 在 `STATE_MOTION / Delta_INT` 中，CM00734 的全部 9 个报告区域均为负的 R1/R3 direction cosine；LY2119620 在多个机制相关区域为正。
- 但 CM00734 在多个区域仍产生与 LY2119620 相当甚至更大的 Delta_INT magnitude；`SIGNED_DRIFT` 分支的 magnitude 也大量重叠。因此 **Delta_INT magnitude 单独不能作为 PAM/cooperativity 判别量**。
- 当前总体科学结论为：`PARTIAL_SPECIFICITY_SUPPORT`。被冻结方法的可区分信息主要来自 `STATE_MOTION` 分支中 Delta_INT **跨 replica 方向一致性**，而不是单纯幅度。
- 这一结论属于 representation-space mechanistic evidence；不构成通用 PAM classifier、效力预测、药理活性证明或前瞻候选确认。
- R2 是历史 calibration/fitting replica；冻结外推解释优先看 R1/R3。20 个时间 block 是相关时间样本，不是独立生物学 replica，不做独立样本 p-value 推断。
- Stage B 结果已冻结：
  - branch：`experiment/pacer-dc-cm00734-stage-b-20ns-v01`
  - commit：`f2f73403c1e965f63cf6fd68d3b6dc60c7dba26a`
  - tag：`pacer-dc-cm00734-stage-b-20ns-v01`
  - final result SHA256：`be663c7d45cc4ade1e053029e6aa09e65ca24c1e2fc22dae5398a26bf49d7d99`

### 方法主张边界

- `Delta_INT = CP - P - C + A`；不得用简单 `CP - C` 替代。
- `Delta_AGO` 仅为历史 representation axis；不得据此直接宣称 intrinsic agonism。
- docking 仅用于 pose proposal / pocket compatibility，不用于 efficacy 或 PAM 判别。
- CM00734 hard-negative 结果已接受，不允许根据其结果重新调 encoder、normalization、bandwidth、RFF、graph、region 或 threshold。
- 旧状态中“Geom2Vec/FKG 已转为失败/消融基线”的表述只代表 2026-09-27 当时阶段，现已被后续 encoder optimization + frozen PACER-FKG v02 结果取代。

## 2026-09-27 PACER-MCV 无编码器物理基线

- OneProt-MD 与 Geom2Vec/FKG 均未在 compound-110 的 R2/R3 上给出稳定的协同方向；它们保留为失败/消融基线，不再继续充当默认主表示。
- 四上下文设计不变。新增 PACER-MCV，直接提取 PAM 接触区、正构核心和胞内微开关的 53 个距离、侧链质心与区域紧致度特征，再计算 `ΔC|A`、`ΔAGO` 与 factorial `ΔINT`。
- 特征清单、按统计量分别进行的 R2-only 标准化和 R3 留出判据已经冻结；单元测试及端到端注入信号测试通过，R1/W0 四条真实 atom14 轨迹的 53 特征提取通过。
- 当前仓库缺少 R2/R3 原始 atom14，仅有 Geom2Vec 派生结果，因此 PACER-MCV 的跨 replica 实值尚未产生。该计算不需要 GPU；有原始坐标的队友应按 `docs/PACER_MCV_METHOD_AND_HANDOFF.md` 运行。
- 在 R2/R3 门禁通过之前，只能称“物理构象基线已实现”，不能声称检测到协同或预测 PAM。

## 2026-09-27 Geom2Vec 替代编码器小试

- OneProt最终表示因局部四上下文差分跨replica不稳定，停止作为默认主编码器；四上下文MD和`dINT/dAGO`定义不变。
- 已将冻结Geom2Vec ViSNet接入现有atom14窗口，无需重跑MD；真实M4输入为100帧、270残基、2,139个重原子，输出为每帧`270 x 128`残基级不变量特征。
- CPU实测约`3.0–4.8 s/frame`；identical repeat与旋转相对误差为0，平移相对误差为`4.42e-7`。
- replica 1五帧描述性小试中，compound110口袋的context separation/temporal variation为`1.18`，global与distal control均约`0.68`，提示局部信号未被全局池化完全抹除。
- 当前只通过工程可用性闸门，未通过跨replica替换闸门；R2/R3必须运行冻结批处理和matched-window审计，随后还需LY2119620药理控制。
- 运行与主张边界见`docs/PACER_DC_GEOM2VEC_HANDOFF.md`。
- 已补充 PACER-FKG 原型：用三套 PAM 结构、一个变构激动剂结构和两套正构结构构建冻结的多结构残基图，并以核均值的四上下文 factorial interaction 取代错误的独立轨迹逐帧相减。
- 三种 PAM 的4.5 Å共识接触为`Y89/Y92/I93/G96/F186/L190/N423/Q427/D432/W435/S436/Y439`；compound-110另有向`Y416/M419/V420/I430/P431`等残基延伸的结构差异。
- 无偏核统计的3个合成测试均通过。真实compound-110单replica五帧烟雾测试同时发现旧ICL3 distal control存在强伪信号；新的跨结构稳定核心对照较低，但尚无跨replica和药理对照证据。
- 100-seed合成基准中，均值保持但分布形状改变时，FKG核交互相对空对照AUC为`0.9966`，线性均值差为`0.5973`；这只证明非线性分布检验的数学增量，不是M4药理性能。
- 方法、baseline和晋级门槛见`docs/PACER_FACTORIAL_KERNEL_GRAPH_METHOD.md`；不得把当前烟雾结果写成PAM判别、性能提升或机制结论。

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

PACER-DC 当前主线已从“继续寻找新 encoder”转为：

**冻结的 C1-BS256 encoder + PACER-FKG v02 → matched four-context Delta_INT → 跨 replica 方向一致性 → hard-negative specificity 审计。**

Stage A 与 Stage B 已完成，不再进行 outcome-driven tuning。当前工作重点是结果归档、外部审计和为下一轮真正前瞻验证冻结 protocol。

## 已完成

- 四上下文、三 replica 的 30-run 试验清单已冻结；
- PACER-DC 多维证据计分器已实现；
- 缺失上下文禁止排名，静态证据禁止冒充完整动态证据；
- 第一阶段自动闸门和测试已实现；
- 主流程已从加权总分改为结合门控、动态差分、baseline/消融和可拒判输出。
- 已锁定同来源近邻对 CM00717（PAM）/CM00734（实验非活性，ECFP4 Tanimoto 约 0.955），用于困难负样本验证；其结合能力仍待检验。

## 当前真实闸门

| 闸门 | 状态 | 说明 |
|---|---:|---|
| 优化 encoder 冻结 | 通过 | C1-BS256 candidate 已冻结；原始 Geom2Vec/ViSNet 权重不微调 |
| 历史 long-MD PACER-FKG v02 | 通过 | 4 systems × 3 replicas × 50 ns 已进入统一冻结表示 |
| Stage A：LY2119620 matched 20 ns | 通过 / 描述性正向 | 主要 `STATE_MOTION / Delta_INT / compound110_extension` R1/R3 cosine `+0.634` |
| Stage B：CM00734 hard negative | 部分特异性支持 | 同一主要信号 R1/R3 cosine `-0.327`，但 magnitude 存在明显重叠 |
| magnitude-only PAM 判别 | 未通过 | inactive CM00734 在多个区域仍有较大 Delta_INT magnitude |
| 通用 PAM classifier | 未建立 | 无冻结分类阈值、无足够独立正负样本、无前瞻验证 |
| 药理效力 / potency 预测 | 未建立 | 当前输出为 representation-space mechanistic descriptors |
| outcome-driven retuning | 禁止 | Stage B 后不得反向调参以改善 hard-negative 结果 |

## 下一执行点

1. 完成 `PROJECT_STATUS.md` 与 Stage A+B 总结报告归档；
2. 保留 Stage A / Stage B 冻结代码、receipt、result SHA 与 tag，不重跑、不调参；
3. 下一轮若继续科学验证，必须在运行前冻结新增正/负样本清单与评价规则；
4. 优先增加新的实验 inactive/binder hard negatives 和独立已知 PAM，而不是继续在 LY2119620 / CM00734 上优化；
5. 只有在独立前瞻样本中复现方向性 discrimination 后，才讨论 classifier、threshold 或候选排名。

自动报告：`project/results/pacer_dc_phase1_audit/REPORT.md`。
# 当前协作任务入口

PACER-FKG v02 的 Stage A + Stage B 已闭环。当前不再以 PACER-MCV 或原版 Geom2Vec 作为默认主线，也不继续对现有 LY2119620 / CM00734 结果做参数优化。

当前入口为：
**冻结结果归档 → 独立审计 → 下一轮预注册样本扩展 / 前瞻验证。**

总结合并报告：
`project/results/pacer_dc_stage_a_b_final_synthesis_v01/PACER_DC_STAGE_A_B_FINAL_SYNTHESIS_v01.md`。
