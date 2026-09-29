# PACER-DC 20 ns 闭环 — 新对话交接 v01

## 新对话第一原则

不要依赖旧聊天记忆。

先读取：

`project/docs/PACER_DC_CLOSE_LOOP_20NS_EXECUTION_PLAN_v01.md`

再读取当前 execution branch、关键 freeze 文件和实际 filesystem assets。

仓库与 result artifacts 是事实源。

---

## 当前任务

执行 PACER-DC 三类 20 ns retrospective closure：

1. compound110：使用已有 50 ns 数据的前 20 ns作为 matched reference；
2. LY2119620：新增 6 × 20 ns；
3. CM00734：系统准备 / binding gate PASS 后新增 6 × 20 ns。

每条：

- 3 replicas；
- seeds 27101 / 38201 / 49301；
- 20 ns；
- 50 ps/frame；
- 400 frames；
- 20 frames/block；
- 20 × 1 ns blocks。

---

## 已冻结、禁止重拟合

### historical/common-kernel line

使用已经完成的 PACER-FKG common-kernel v03 numerical state。

不得重新估计：

- channel median；
- channel scale；
- bandwidth；
- RFF；
- graph；
- regions；
- dAGO/dPAM/dINT definitions。

### optimized FKG v02 line

encoder anchor：

`0a4b4f238c73add75a09a0875f698250a39bdf33`

tag：

`encoder-candidate-ah-frozen-20260929`

checkpoint SHA256：

`b8f1ef9b591c57f7687566bd60d3664800956ac280cf2ae8e64f54196cd8d417`

Geom2Vec source：

`371d642ec1061664f16e49fcac702d07fc8d0b51`

FKG v02 freeze SHA256：

`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`

不得对 LY/CM00734 重新运行 outcome calibration。

---

## compound110 已有事实

600 ns long-MD 已完成并归档。

common-kernel 最可靠 descriptive result：

`dINT / compound110_extension`

为最稳定 interaction-associated shift。

dAGO 三-replica intrinsic agonism signature 不稳健。

historical kernel_u2 gate = REVOKED。

qualification gate = NOT_DEFINED。

50 ns 完整结果继续保留。

closure comparison 只使用前 20 ns。

---

## 当前关键工程事实

LY2119620 已经存在于现有：

- membrane/reference preparation pipeline；
- restraint-release pipeline；
- production runner。

因此 LY 第一任务是资产 + lineage preflight，
不是重新设计体系。

CM00734 当前主要存在于：

- benchmark；
- labels；
- manifest。

它尚未接入现有 runnable membrane/reference/production system。

因此 CM00734 在 20 ns production 之前必须先完成完整 system preparation。

旧 `run_pacer_mechanism_cv.py` 是 R2/R3 + 5-window runner。

它不能直接作为 20 ns three-replica closure runner。

需要新建 long/closure adapter，但 frozen physical feature config 不得修改。

---

## 新对话第一阶段必须执行

1. `git status`
2. `git rev-parse HEAD`
3. 阅读 execution plan
4. 阅读当前 close-loop source-of-truth files
5. 检查本地/GPU LY assets：
   - membrane reference
   - restraint-release state
   - audit
   - hashes
6. 自动比较 LY 与现有 apo/probe controls lineage
7. 输出机器可读 preflight JSON
8. 只有 PASS 才启动 LY MD

预期第一阶段最终状态：

`LY2119620_CLOSE_LOOP_PREFLIGHT_PASS`

随后才生成 6 × 20 ns production command。

---

## 不要做

不要：

- 修改 historical frozen results；
- 调 encoder；
- 重拟合 FKG；
- 按 outcome 选 region；
- 按 outcome 设 threshold；
- 直接假定 CM00734 是 confirmed binder；
- 把 temporal blocks 当独立 biological replicates；
- 启动 supervised training；
- 直接跑 prospective candidates。

---

## 当前无关工作区文件

历史未跟踪文件：

`repair_g2a_pbc_report.py`

它与 20 ns closure 无关。

不要误提交或修改，除非独立审计证明需要它。

---

## 新对话首句建议

“请先完整读取
`project/docs/PACER_DC_CLOSE_LOOP_20NS_EXECUTION_PLAN_v01.md`
和
`project/docs/PACER_DC_CLOSE_LOOP_20NS_HANDOFF_v01.md`，
然后以仓库和实际 result filesystem 为唯一事实源，
不要依赖旧聊天记忆。先完成 Stage 0 execution audit 和
LY2119620 preflight；在 preflight PASS 前不要启动 MD。”
