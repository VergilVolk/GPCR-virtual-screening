# PACER-DC 20 ns 三类闭环执行计划 v01

**冻结日期：2026-09-29**

**状态：EXECUTION PLAN — PRE-OUTCOME FROZEN**

本文件定义当前 PACER-DC 四上下文功能判别闭环的实际执行路线。

它不覆盖历史 50 ns / 600 ns 结果，不修改已经冻结的 PACER-FKG
common-kernel v03 或 PACER-FKG v02 数值状态。

由于赛程与 GPU 资源限制，从本闭环开始，所有新增候选体系统一使用
**20 ns matched closure MD**。

---

## 1. 闭环目标

当前目标不是继续优化 encoder 或 PACER-FKG。

目标是在已经冻结的方法下，对已知功能类型不同的配体进行回顾性测试：

1. compound110：已有 ago/interacting control；
2. LY2119620：known functional PAM；
3. CM00734：实验无功能 hard negative，binding compatibility 尚待验证。

最终回答：

> 同一个冻结四上下文框架，能否在已知 PAM、compound110 类
> interaction/agonism control、以及实验无功能 hard negative 中产生
> 可区分且可重复的动态表示模式？

这是 retrospective framework validation。

不是 PAM 药效预测，不是候选确认，也不是 supervised classifier validation。

---

## 2. 当前已经完成的事实

### 2.1 compound110 long-MD

已经完成：

- apo；
- probe_only；
- compound110 candidate_no_probe；
- compound110 candidate_probe；
- R1 / R2 / R3；
- 每条 50 ns；
- 总计 600 ns；
- 12,000 stored frames。

### 2.2 common-kernel PACER-FKG

已完成并冻结：

- historical 128D Geom2Vec representation；
- complete preprocessing provenance；
- frozen G2-B channel median/MAD；
- shared bandwidth；
- common RFF basis；
- dAGO / dPAM / dINT；
- R1/R2/R3 replica consistency；
- stable/distal specificity audit。

当前主要结果：

`dINT / compound110_extension`

是 compound110 最可复现的 interaction-associated signal。

当前结果不支持稳定的三-replica intrinsic-agonism signature。

qualification gate = NOT_DEFINED。

historical kernel_u2 functional gate = REVOKED。

### 2.3 optimized encoder / PACER-FKG v02

已经完成 encoder A-H 优化。

冻结 encoder candidate：

C1 intermediate ViSNet
→ atomic [x, ||v||]
→ backbone / side-chain separated readout
→ BS256。

冻结 anchor：

- encoder candidate commit:
  `0a4b4f238c73add75a09a0875f698250a39bdf33`
- encoder tag:
  `encoder-candidate-ah-frozen-20260929`
- ViSNet checkpoint SHA256:
  `b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`
- Geom2Vec source commit:
  `371d642ec1061664f16e49fcac702d07fc8d0b51`

PACER-FKG v02 Phase 2 frozen numerical state：

- freeze manifest SHA256:
  `b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`
- STATE_MOTION bandwidth:
  `32.30188361260893`
- STATE_MOTION RFF seed:
  `272340`
- SIGNED_DRIFT bandwidth:
  `29.29934899925964`
- SIGNED_DRIFT RFF seed:
  `272084`

R2 是 v02 calibration replica。

R1 / R3 是 frozen application replicas。

这些定义后续不得重新拟合。

---

## 3. 为什么新的闭环统一改为 20 ns

50 ns compound110 long-MD 已经证明：

- 数据链可运行；
- PACER-FKG 可运行；
- common-kernel representation 可审计；
- encoder/FKG v02 可冻结；
- replica consistency 可以正式评估。

后续任务已经不是探索轨迹长度，而是用冻结方法测试新的已知分子。

由于剩余时间和 GPU 资源不足，新体系统一限制为：

**20 ns / trajectory**

为了避免 compound110 50 ns 与 LY/CM00734 20 ns 不公平比较，
闭环的跨分子直接比较统一使用 **前 20 ns**。

历史完整 50 ns compound110 分析继续保留，作为独立 robustness evidence。

不得用 50 ns compound110 数值直接与 20 ns LY/CM00734 magnitude 做主比较。

---

## 4. 20 ns matched closure sampling contract

所有 closure molecule 使用：

- independent replicas: R1 / R2 / R3
- seeds:
  - R1 = 27101
  - R2 = 38201
  - R3 = 49301
- trajectory length: **20 ns**
- frame spacing: **50 ps**
- frames per trajectory: **400**
- block size: **20 stored frames**
- block duration: **1 ns**
- non-overlapping blocks per trajectory: **20**

统计单位仍然是：

**replica / molecule**

20 个 blocks 是 trajectory temporal subsamples。

它们不是 20 个独立 biological replicates。

现有 apo / probe_only / compound110 50 ns trajectories 在 closure layer 中统一使用：

**frames 0–399，即前 20 ns**

不得根据结果选择其它时间窗。

---

## 5. 三条分析轨道

所有闭环分子都尽可能运行三条预先定义的分析轨道。

### Track A — historical 128D common-kernel PACER-FKG v03

复用已经冻结的：

- preprocessing；
- channel median；
- channel scale；
- bandwidth；
- RFF weights/bias；
- graph；
- regions；
- dAGO / dPAM / dINT definitions。

新分子只能 APPLY。

禁止重新 fit。

20 ns closure 时：

400 frames
→ 20 × 20-frame blocks
→ common RFF
→ region vectors
→ per-replica pooled vectors
→ R1/R2/R3 consistency
→ stable/distal specificity。

### Track B — frozen C1-BS256 PACER-FKG v02

新分子使用 frozen C1-BS256 encoder。

复用现有 v02 numerical freeze：

`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`

禁止在 LY 或 CM00734 上运行新的 Phase 2 calibration。

新分子只能作为 frozen application data。

STATE_MOTION 与 SIGNED_DRIFT 分开报告。

不得根据结果选择 branch。

### Track C — PACER-MCV physical baseline

复用：

`project/config/pacer_m4_mechanism_cv_v01.json`

物理 feature panel 不得改。

需要新建适配 long/closure 数据的 runner，使其支持：

- R1/R2/R3；
- 400 frames；
- 20 × 1 ns blocks；
- 与 closure axes 相同的 dAGO / dPAM / dINT。

旧 R2/R3 5-window MCV 结果保留为历史结果，不覆盖。

---

## 6. Stage 0 — execution freeze

正式运行任何新的 biological outcome 前，先完成：

1. 建立本 execution branch；
2. 记录当前 Git HEAD；
3. 冻结本执行计划 SHA256；
4. 建立 `CLOSE_LOOP_EXECUTION_MANIFEST_v01.json`；
5. 核对：
   - encoder anchor；
   - checkpoint；
   - Geom2Vec source；
   - common-kernel contract；
   - common-kernel final freeze；
   - FKG v02 freeze；
   - MCV config；
   - seeds；
   - 20 ns sampling contract。

从这一步以后：

**不得 outcome-driven tuning。**

---

## 7. Stage A — LY2119620 known-PAM closure

### 7.1 先做 preflight，禁止直接开 MD

现有仓库 MD pipeline 已经支持：

- `LY2119620__candidate_no_probe`
- `LY2119620__candidate_probe`

所以第一步不是重建体系。

先核查实际 GPU / result filesystem 中：

- membrane reference；
- minimized.pdb；
- system.xml；
- restraint-release output；
- production_start_state.xml；
- audit.json；
- SHA256；
- `production_start_eligible=true`。

然后做 control-lineage compatibility audit。

至少比较 LY 与已有 apo/probe_only 的：

- receptor origin；
- force field；
- ligand parameterization policy；
- OPM transform；
- POPC/water/ions；
- temperature；
- barostat；
- restraint release；
- final k=0 condition；
- receptor sequence/mapping；
- periodic box setup；
- timestep；
- seeds；
- frame spacing。

### 7.2 control reuse decision

如果 compatibility audit PASS：

复用现有 apo / probe_only closure 20 ns reference。

只新增：

- LY no-probe R1/R2/R3；
- LY + probe R1/R2/R3。

总新增 MD：

**6 × 20 ns = 120 ns**

如果 FAIL：

不得称现有 controls 为 paired controls。

必须单独记录 lineage incompatibility，再决定是否有资源重建 controls。

### 7.3 LY production

必须显式使用：

`--ns 20`

以及：

`--trajectory-ps 50`

注意现有 production runner 默认 trajectory interval 是 10 ps，
因此不能依赖默认值。

新 LY closure MD 必须写入新的 result root。

不得覆盖历史 `pacer_dc_production_v01`。

### 7.4 LY QC

每条 trajectory 完成后检查：

- completed_ns = 20；
- 400 frames；
- finite；
- temperature/density；
- box；
- receptor 270 residues；
- sequence mapping；
- atom14 mapping；
- PBC / CA continuity；
- exact seed provenance。

QC PASS 后才进入 functional analysis。

### 7.5 LY analysis

依次运行：

1. common-kernel v03 apply-only；
2. frozen C1-BS256 extraction；
3. frozen FKG v02 apply-only；
4. MCV physical baseline。

最终比较：

compound110 20 ns
vs
LY2119620 20 ns。

Stage A 的完成条件是：

**所有预注册轨道有结果。**

不是要求 LY 一定表现出预期 PAM pattern。

支持或否定预期都必须报告。

---

## 8. Stage B — CM00734 inactive hard-negative closure

CM00734 当前不是 runnable production system。

它现在只存在于：

- benchmark label；
- manifest；
- hard-negative design。

其当前合法称呼是：

**experimentally inactive hard negative, binding compatibility to be tested**

不能提前称为 confirmed inactive binder。

### 8.1 必须先建立系统

在 LY GPU 运行期间并行完成：

pose/source freeze
→ ligand parameterization
→ static pocket compatibility
→ minimization
→ membrane reference
→ short equilibration
→ restraint release
→ zero-restraint production-start audit
→ short binding/pose-stability QC。

只有 production eligibility PASS 才能进入 20 ns closure MD。

### 8.2 CM00734 production

如果 binding/system gate PASS：

运行：

- candidate_no_probe R1/R2/R3；
- candidate_probe R1/R2/R3；

每条 20 ns。

总新增：

**120 ns**

协议与 LY 完全一致。

### 8.3 CM00734 analysis

使用和 LY 完全相同的 frozen pipeline：

- Track A common-kernel v03；
- Track B frozen FKG v02；
- Track C MCV。

禁止因 CM00734 是 negative 而调整 region、threshold、kernel 或 preprocessing。

### 8.4 若 CM00734 gate FAIL

必须如实记录：

`CM00734 NOT ELIGIBLE FOR FUNCTIONAL MD CLOSURE`

可以把 GPU 队列切换到 CM00717 作为第二 functional PAM evidence。

但 CM00717 **不能替代 negative class**。

因此 CM00734 FAIL 时不得宣布 three-class closure 成立。

---

## 9. 最终三类矩阵

最终目标矩阵：

| role | molecule |
|---|---|
| ago/interacting control | compound110 |
| known functional PAM | LY2119620 |
| experimental inactive hard negative | CM00734 |

所有体系使用统一 20 ns closure horizon。

统一报告：

- dAGO；
- dPAM；
- dINT；
- replica pooled norms；
- pairwise direction cosines；
- coherence；
- stable-core specificity；
- distal-control specificity；
- temporal block behavior；
- STATE_MOTION；
- SIGNED_DRIFT；
- MCV physical response。

不得生成一个临时 PAM 总分。

不得为了实现三类分离重新定义指标。

---

## 10. MCV 的作用

MCV 不是第三个竞争模型。

它回答：

> learned representation 提供的信息，是不是超出了简单受体物理构象特征？

无论结果：

- FKG > MCV；
- FKG ≈ MCV；
- MCV > FKG；

都必须报告。

它属于 framework interpretation，而不是为了证明 FKG 更好。

---

## 11. Training gate

当前闭环不启动 supervised training。

即使 compound110 / LY2119620 / CM00734 三类结果完整，
也不自动打开 training gate。

必须另行满足当前仓库 training audit 的 molecule-level 数据要求、
class coverage 和 leak-free train/validation/test 条件。

20 ns closure 的目的：

**framework validation**

不是 classifier training。

---

## 12. GPU / CPU 并行路线

GPU critical path：

LY2119620 120 ns
→ CM00734 120 ns（若 eligible）。

CPU / engineering parallel path：

- execution manifest；
- LY preflight；
- 20 ns compound110 reference；
- common-kernel apply-only runner；
- FKG v02 apply-only runner；
- 20 ns MCV runner；
- CM00734 system preparation；
- provenance / audit / reporting。

LY GPU 启动后，不应等待 CPU 工作。

所有分析工具必须在 LY outcome inspection 前尽可能冻结。

---

## 13. 明确移出 critical path 的任务

当前不要执行：

- 新 encoder optimization；
- PACER-FKG recalibration；
- v01 exact replay；
- v1-v2 superiority audit；
- PACER0076；
- PACER0057；
- dual-context supervised heads；
- training-gate model fitting；
- outcome-driven region selection；
- outcome-driven cosine threshold selection。

这些任务都不能阻塞三类 closure。

---

## 14. Claim boundary

20 ns matched closure 可以支持：

- frozen-framework retrospective discrimination evidence；
- replica consistency comparison；
- learned-vs-physical representation comparison。

20 ns 不能自动支持：

- conformational convergence；
- equilibrium free-energy claim；
- PAM efficacy prediction；
- potency prediction；
- candidate confirmation；
- prospective success claim。

blocks 不是独立 biological N。

已知标签用于 retrospective testing，
不能被描述为 prospective discovery。

---

## 15. Stage A / Stage B 完成定义

### Stage A COMPLETE

必须具有：

- compound110 matched 20 ns reference；
- LY2119620 3-replica × 2-context 20 ns data；
- Track A；
- Track B；
- MCV baseline；
- provenance + QC。

结果方向不影响 COMPLETE。

### Stage B COMPLETE

必须具有：

Stage A

加：

- eligible CM00734 3-replica × 2-context 20 ns data；
- 同样三条分析轨道；
- negative-class comparison；
- final three-class matrix。

结果是否符合原假设不影响 COMPLETE。

---

## 16. 最终归档

最终至少保存：

- execution contract / plan；
- execution manifest；
- trajectory manifests；
- SHA256；
- QC；
- feature manifests；
- common-kernel outputs；
- FKG v02 outputs；
- MCV outputs；
- replica audits；
- specificity audits；
- final comparison table；
- final report；
- final freeze receipt。

大型 DCD / cache 不进入 Git。

Git 保存：

代码、配置、manifest、small results、reports、SHA256 provenance。

最终建立独立 immutable tag。

---

## 17. 第一条实际执行动作

新对话开始后，不直接启动 MD。

第一项任务必须是：

**完整读取本 execution plan + 当前仓库 + 本地/GPU LY assets，
生成并运行 LY2119620 preflight / control-lineage compatibility audit。**

只有明确得到：

`LY2119620_CLOSE_LOOP_PREFLIGHT_PASS`

才启动 20 ns production MD。
