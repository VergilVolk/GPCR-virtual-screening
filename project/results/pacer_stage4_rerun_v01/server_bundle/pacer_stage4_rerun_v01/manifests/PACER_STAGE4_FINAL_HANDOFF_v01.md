# PACER Stage4 Final Handoff v01

最终交接日期：2026-10-03。用途：供未参与当前聊天的新成员接手冻结产物与功能实验验证。

## 1. Final status

```text
Stage4 prospective MD: COMPLETE
PACER-FKG v02 frozen evaluation: COMPLETE
Scientific synthesis: COMPLETE
Stage4 prospective computational evaluation: COMPLETE
Wet-lab functional validation: NOT YET PERFORMED
```

已完成计算链：candidate generation / filtering → multi-conformation docking → M4-safe DrugCLIP → candidate selection → cluster-matched four-context prospective MD → PACER-FKG v02 frozen evaluation → scientific interpretation。

本次交接只归档文档、脚本和轻量文本结果，不重新计算、重新测试或改变已完成的科学结果。完整结果状态为`STAGE4_FULL_FROZEN_EVALUATION_COMPLETE`，不是“候选已获药理确认”。

## 2. Frozen candidates

| Candidate | Frozen source cluster | A | P |
|---|---|---|---|
| PACER0010 | 9 | cluster9__apo | cluster9__probe_only |
| PACER0027 | 4 | cluster4__apo | cluster4__probe_only |
| PACER0073 | 0 | cluster0__apo | cluster0__probe_only |

候选身份、cluster和上下文已冻结；本次没有重新筛选、重排或替换候选。不得借用另一个cluster或历史7TRS的原始A/P baseline。

## 3. Four-context design

| Alias | Context | System naming |
|---|---|---|
| A | apo | clusterN__apo |
| P | probe_only / ACh alone | clusterN__probe_only |
| C | candidate_no_probe | PACERxxxx__candidate_no_probe |
| CP | candidate_probe / candidate + ACh | PACERxxxx__candidate_probe |

每个候选只使用自身cluster-matched的A/P/C/CP。contrast在每个replica系列分别计算，保留全部R1/R2/R3，不依据结果剔除replica。

## 4. MD design

- 3 candidates × 4 contexts = 12 systems；每system 3 replicas，共36 trajectories。
- 每条10 ns、5,000,000 production steps、2 fs timestep，共360 ns aggregate。
- Paired seeds：R1=27101、R2=38201、R3=49301。跨上下文相同seed标识同一replica系列，不保证微观状态严格配对。
- 每条1000 raw frames，实际10 ps间隔。数据根目录`C:/projects/PACER_STAGE4_MD_backup`保留原始systems/production及来源记录。
- 35条正常progress完成；PACER0010/CP/R3的running/9.9 ns progress属于既有特定陈旧记录例外，由完整帧数、最终步数、有限坐标与拓扑兼容性认证为物理完整。原progress不改，例外不泛化。
- 10 ns轨迹和3个种子不等于动力学已收敛或药理机制已证实。

## 5. PACER-FKG pipeline

```text
36 authenticated trajectories
→ Phase1 frozen C1-BS256 encoder
→ Phase2a fixed temporal blocks + historical frozen-state apply
→ Phase2b fixed graph / regions / matched contrasts
→ frozen prospective evaluation
→ independent scientific interpretation
```

| Phase | Final receipt status | Inventory |
|---|---|---|
| Phase1 | PHASE1_FULL_FROZEN | cache_count=36；每条float32[1000,270,256] |
| Phase2a | PHASE2A_FULL_FROZEN | job_count=36；manifest_count=72 |
| Phase2b | STAGE4_FULL_FROZEN_EVALUATION_COMPLETE | 3 candidates × 2 branches；18份区域向量NPZ |

实际1000帧通过`frames[::5]`得到200个50 ps分析帧，按20帧/block划分为10个名义1 ns相关block（内部跨度950 ps）。STATE_MOTION为状态均值与内部差分RMS的512维拼接；SIGNED_DRIFT为内部平均差分的256维表示。两branch分别使用既有normalization/bandwidth/RFF资产，映射为512维RFF，再按冻结图执行alpha=0.65、20步restart diffusion和九区域均值。

历史freeze manifest SHA-256：

```text
b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd
```

历史calibration来自long-MD四上下文的R2（compound110/apo/probe背景，共200个历史block）；**Stage4 R1/R2/R3全部EVALUATION_ONLY，Stage4 R2从未重新用于calibration**。没有Stage4-based refit、retraining或threshold selection。

Stage4拓扑选择修复已在common adapter中：OpenMM把构建体分为A/B两段，依据冻结序列在原拓扑顺序上做唯一投影，保留270残基/2139个ATOM14重原子；不要求单一PDB chain，不取第一条链，缺失/歧义均拒绝。历史encoder和270节点身份不变。

已有验证记录：27 tests（原23＋新增4）全部通过，0 failures/errors/skips；真实CUDA smoke通过；后续独立只读科学复核验证完整冻结产物，并从保存NPZ独立核对810组摘要及810个余弦，差值均0。**这些是已有记录，本次最终收尾未重新执行测试或阶段脚本。**

## 6. Frozen contrasts

```text
Delta_AGO = C - A
Delta_PAM = CP - P
Delta_INT = CP - C - P + A
```

这些名称是frozen method terminology，不是功能分类：Delta_AGO表示candidate-alone相对apo的动态扰动；Delta_PAM表示probe背景下加入candidate的变化；Delta_INT表示representation-space二阶nonadditivity/interaction contrast。非零幅度不能单独证明激动、正调节或协同。

另有Probe_effect=P−A及Probe_background_effect=CP−C。两branch的RFF空间不同，不直接比较跨branch幅度或计算跨branch向量余弦。10个block相关，不当作独立统计样本。

## 7. Final scientific findings

以[正式科学报告](../results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md)为完整数值与限制来源：

- PACER0010：支持有限的probe-dependent/intrinsic dynamic pattern separation；不能证明无intrinsic功能效应。
- PACER0027：STATE_MOTION candidate-alone/intrinsic perturbation较一致；应在功能实验中专门检查candidate-alone响应。
- PACER0073：PAM-related regions中candidate-alone directional response较稳定；probe效应仍有区域和branch依赖。
- 没有候选在多数关键区域、跨STATE_MOTION与SIGNED_DRIFT两branch表现稳健一致的Delta_INT。确定cooperative PAM mechanism is **not supported by the current frozen prospective evaluation**。
- 控制区域同样有明显幅度变化，区域重叠且经图扩散；不能宣称严格PAM区域特异或统一全受体激活。

这些结果用于安排功能实验问题与机制假说，不构建新综合评分或候选药效排名。

## 8. Claim boundary

> Prospective frozen PACER-FKG dynamic differential evaluation only.
> **No PAM/ago-PAM labels, probability, efficacy or independent-block inference.**

不得将候选称为confirmed PAM、ago-PAM或inactive，不输出模型成功率/效能，不根据某个漂亮region/replica宣称协同。负余弦是表示空间方向分歧，不是物理反向运动或药理抑制证据。

## 9. Temporal metadata note

```text
source master_manifest: trajectory_ps=50 (raw metadata description error)
actual DCD: 1000 frames / 10 ns, raw spacing 10 ps
PACER-FKG: frames[::5], 10 ps × 5 = 50 ps analysis spacing
```

正式说明：[STAGE4_TEMPORAL_METADATA_NOTE_v01.md](../results/pacer_stage4_prospective_fkg_v02_v01/provenance/STAGE4_TEMPORAL_METADATA_NOTE_v01.md)。实际输入与冻结Stage4分析契约一致；该源字段描述错误不影响已完成结果。原manifest/DCD与早期科学复核报告保持原字节，不静默修改、不重跑。服务器参数覆盖的来源不作无证据推测。

## 10. Important artifacts

本仓库相对路径如下，Windows绝对根目录见下一节。

| Artifact | Path |
|---|---|
| Final result JSON | project/results/pacer_stage4_prospective_fkg_v02_v01/phase2b_graph_region/STAGE4_PROSPECTIVE_RESULTS_v01.json |
| Scientific report | project/results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_REPORT_v01.md |
| Scientific summary | project/results/pacer_stage4_prospective_fkg_v02_v01/STAGE4_PROSPECTIVE_SCIENTIFIC_SUMMARY_v01.json |
| Provenance | project/results/pacer_stage4_prospective_fkg_v02_v01/provenance/ |
| Input authentication | project/results/pacer_stage4_prospective_fkg_v02_v01/input_audit/INPUT_AUTHENTICATION_v01.json |
| Phase1 full / smoke receipts | project/results/pacer_stage4_prospective_fkg_v02_v01/phase1_bs256/PHASE1_FREEZE_RECEIPT_v01.json；SMOKE_RECEIPT_v01.json |
| Phase2a receipt | project/results/pacer_stage4_prospective_fkg_v02_v01/phase2a_frozen_apply/PHASE2A_FREEZE_RECEIPT_v01.json |
| Phase2b receipt | project/results/pacer_stage4_prospective_fkg_v02_v01/phase2b_graph_region/PHASE2B_FREEZE_RECEIPT_v01.json |
| Per-input manifests | phase1_bs256/manifests/full/、manifests/smoke/；phase2a_frozen_apply/manifests/full/（相对上述结果根目录） |
| Frozen manifest（历史资产根目录） | project/results/pacer_fkg_v02_longmd_v01/calibration/V02_FREEZE_MANIFEST.json |
| Adapter / Phase1 | project/pacer_fkg_v02/stage4_prospective_common_v01.py；run_stage4_prospective_phase1_bs256_v01.py |
| Phase2a / Phase2b | project/pacer_fkg_v02/run_stage4_prospective_phase2a_frozen_apply_v01.py；run_stage4_prospective_phase2b_graph_region_v01.py |
| Tests | project/tests/test_stage4_prospective_fkg_v02_v01.py |
| Prior 27-test / topology / smoke evidence | provenance/TOPOLOGY_SELECTION_REPAIR_VALIDATION_v01.json（相对结果根目录） |

Git归档只收脚本/测试、状态/交接文档和本Stage4结果根目录内明确选择的JSON/Markdown。最终结果JSON约19.4 MB，保存为文本结果证据；不纳入DCD、PDB解析副本、NPY/NPZ缓存、checkpoint、临时文件、server copy或PACER_STAGE4_MD_backup。`.gitignore`保持不变，只对选定文本文件使用精确force-add。

Git内的回执仍引用外部原始/缓存资产及其哈希；新成员必须另行取得这些原始存档和历史冻结资产，不能把只有Git checkout视为完整数据环境。禁止用缺失资产触发“自动重跑”补齐。provenance中的早期下载阻断/未开始状态属于历史快照，不能覆盖最终冻结回执。

仓库当前`core.autocrlf=true`；本次新增的`.gitattributes`仅对Stage4结果根目录及这4份adapter脚本/1份测试设置`-text`，避免跨平台checkout改写其换行、破坏回执中精确字节SHA。未修改这些产物或脚本内容，也未改变其他历史文件的换行规则。

## 11. Reproducibility / environment

```text
Windows repository: C:\projects\GPCR-virtual-screening
Stage4 raw MD backup: C:\projects\PACER_STAGE4_MD_backup
Historical frozen PACER-FKG assets: C:\projects\GPCR-virtual-screening-fkg-v02
Conda environment: pacer_dc_geom2vec
Python: C:\Users\barne\anaconda3\envs\pacer_dc_geom2vec\python.exe
Pinned geom2vec source: C:\projects\geom2vec-source
```

历史checkpoint SHA：`b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`；geom2vec source pin：`371d642ec1061664f16e49fcac702d07fc8d0b51`。历史CUDA smoke/strict-load证据见已有provenance，不重新加载模型推理。环境说明文件为`project/environment_pacer_dc_geom2vec.yml`，不能用README的通用CPU核心环境替代该隔离环境。

数值资产、source hash、每个结果的实现身份、seed和输入路径以冻结回执为准。本次收尾只复核已有记录、文本哈希及freeze SHA；不重新调用Phase1/Phase2或测试命令。

任务开始前已存在4个历史管理文件变更：BASELINE_TRACKED_TREE_SHA256.json、PHASE0_GATE.json、REPORT_PHASE0.md、V01_COMPONENT_INVENTORY.json。本次保持其工作区字节不变，排除在收尾commit外。不能为了clean tree而restore/覆盖它们。本次commit留在本地，不push。

## 12. What NOT to redo

- No docking rerun / DrugCLIP rerun / Stage3 rerun。
- No MD rerun / Phase1 rerun / Phase2 rerun or refit。
- No retraining / recalibration / Stage4-based calibration。
- No threshold tuning / outcome-driven tuning / v03 redesign。
- 不改变frozen residue mapping、encoder、normalization、bandwidth、RFF、graph、regions、contrasts或候选。
- 不修改历史冻结资产、Stage4原始PDB/DCD/progress/master_manifest或已完成的科学结果。

## 13. Next scientific step

**Functional experimental validation**。

至少比较candidate alone、ACh alone、candidate + ACh，并使用匹配vehicle对照。intrinsic agonism、positive modulation、ago-PAM possibility与inactive binder均为待实验区分的hypotheses，不是当前预测标签；inactive binding还需独立结合证据，功能无响应本身不能证明结合。

以最终报告的candidate-specific问题安排验证：0010检查probe/intrinsic是否可分离；0027排查candidate-alone功能响应；0073检验candidate-alone driven动态模式及probe影响。这里不增加新的实验阈值、成功概率、评分或药效排名。
