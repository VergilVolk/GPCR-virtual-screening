# PACER-DC Stage A + Stage B Final Scientific Synthesis v01

更新时间：2026-10-01

## 1. 目的

本报告合并 PACER-FKG v02 的两个 matched-20-ns 验证阶段：

- **Stage A**：已知 M4 PAM `LY2119620` 的 matched closure；
- **Stage B**：实验 inactive hard negative `CM00734` 的冻结 specificity test。

两阶段均使用同一冻结 encoder / PACER-FKG 数值状态。Stage B 不参与任何调参。

## 2. 统一 representation provenance

当前主分析不是原版 Geom2Vec final embedding，而是优化后冻结的 `C1-BS256` candidate：

1. 使用原始 Geom2Vec / ViSNet 预训练 checkpoint，权重不重新训练、不微调；
2. 从 message-passing layers 0/1/2 后的 `C1` accumulated state 读取中间表示；
3. 原子级 invariant descriptor 为 scalar 64 + vector-norm 64 = 128D；
4. residue readout：backbone mean = B128，sidechain mean = S128，concat(B,S) = BS256；
5. 后续 PACER-FKG v02 使用冻结 temporal branches、normalization、RFF、graph 与 region。

冻结 PACER-FKG v02 anchor：

`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`

因此 historical 600 ns long-MD、Stage A LY2119620 和 Stage B CM00734 属于同一 encoder generation 下的连续验证链。

## 3. 冻结 contrast

四上下文：

- A = apo
- P = probe_only
- C = candidate_no_probe
- CP = candidate_probe

冻结 contrasts：

- `Delta_PAM = CP - P`
- `Delta_AGO = C - A`
- `Delta_INT = CP - P - C + A`

本报告的 specificity 讨论以 `Delta_INT` 为主。

## 4. Stage A：LY2119620 matched closure

Stage A 使用 20 ns / trajectory、50 ps / frame、400 frames / trajectory、20 frames / temporal block、20 blocks / trajectory；A/P 为 historical first 20 ns，C/CP 为 LY2119620 no-probe / +probe，3 个 paired replicas。

Stage A result status：

`TRACK_B_MATCHED_20NS_COMPLETE`

Stage A claim boundary：

> Retrospective matched-20ns descriptive mechanism comparison. No threshold tuning, p-values, classifier fitting, efficacy prediction, or independent-block inference.

### 4.1 主要 Stage-A signal

`STATE_MOTION / Delta_INT / compound110_extension`

| Panel | R1 magnitude | R2 magnitude | R3 magnitude | R1/R3 direction cosine |
|---|---:|---:|---:|---:|
| compound110 reference | 0.140523 | 0.084188 | 0.105948 | +0.583656 |
| LY2119620 | 0.104348 | 0.085594 | 0.104993 | +0.634304 |

LY2119620 在该 region 中表现出明显正的 R1/R3 directional agreement，并且 magnitude 与 compound110 reference 同一量级。

### 4.2 LY2119620 的 region-level directional pattern

| Region | LY2119620 R1/R3 cosine |
|---|---:|
| compound110_extension | +0.634 |
| cooperativity_mutagenesis | +0.235 |
| distal_control | -0.529 |
| intracellular_microswitches | +0.180 |
| orthosteric_activation_core | +0.321 |
| orthosteric_contact_union | +0.300 |
| pam_contact_consensus | +0.375 |
| pam_contact_union | +0.306 |
| stable_core_control | -0.319 |

这一结果支持“机制相关区域中的方向性重复性”作为后续 hard-negative 检查对象，但 Stage A 单独不能建立 specificity。

## 5. Stage B：CM00734 frozen hard-negative test

CM00734 是预先指定的实验 inactive hard negative。

Stage-B production：

- 2 candidate contexts
- 3 paired replicas / context
- 20 ns / trajectory
- 120 ns 新 production
- 2400 stored frames
- paired seeds R1/R2/R3 = 27101 / 38201 / 49301

所有 6 条轨迹通过 production QC。

Stage-B final result：

`STAGE_B_MATCHED_20NS_COMPLETE`

Stage-B final SHA256：

`be663c7d45cc4ade1e053029e6aa09e65ca24c1e2fc22dae5398a26bf49d7d99`

冻结点：

- branch: `experiment/pacer-dc-cm00734-stage-b-20ns-v01`
- commit: `f2f73403c1e965f63cf6fd68d3b6dc60c7dba26a`
- tag: `pacer-dc-cm00734-stage-b-20ns-v01`

Stage B 未执行 encoder tuning、normalization refit、bandwidth refit、RFF refit、graph refit、region refit、threshold selection 或 outcome-driven tuning。

## 6. Stage A vs Stage B 核心比较

同一主 signal：

`STATE_MOTION / Delta_INT / compound110_extension`

| Panel | R1 magnitude | R2 magnitude | R3 magnitude | R1/R3 direction cosine |
|---|---:|---:|---:|---:|
| compound110 reference | 0.140523 | 0.084188 | 0.105948 | +0.583656 |
| LY2119620 | 0.104348 | 0.085594 | 0.104993 | +0.634304 |
| CM00734 | 0.067244 | 0.068286 | 0.055084 | -0.326839 |

CM00734 / LY2119620 magnitude ratio：

- R1 = 0.644
- R2 = 0.798
- R3 = 0.525

最重要的区别不是 magnitude 是否为零，而是跨 replica direction：compound110 reference 与 LY2119620 为正向一致；CM00734 为负。

## 7. Region-level specificity evidence

`STATE_MOTION / Delta_INT`：

| Region | CM00734 | LY2119620 | compound110 |
|---|---:|---:|---:|
| compound110_extension | -0.327 | +0.634 | +0.584 |
| cooperativity_mutagenesis | -0.095 | +0.235 | -0.246 |
| intracellular_microswitches | -0.006 | +0.180 | -0.111 |
| orthosteric_activation_core | -0.385 | +0.321 | -0.074 |
| orthosteric_contact_union | -0.376 | +0.300 | -0.052 |
| pam_contact_consensus | -0.251 | +0.375 | -0.308 |
| pam_contact_union | -0.241 | +0.306 | -0.183 |

CM00734 在全部 9 个报告 region 中的 R1/R3 cosine 均为负；LY2119620 则在多个机制相关 region 中为正。

注意：这些数值是**每个 panel 内部的 R1-vs-R3 directional reproducibility**，不是 CM00734 与 LY2119620 两个 mean vector 之间的直接 cosine。

## 8. 为什么不能把 magnitude 当 classifier

CM00734 并不是“没有 Delta_INT”。

在 `STATE_MOTION` 中：

- cooperativity_mutagenesis：CM/LY = 1.148 / 0.984 / 1.320
- intracellular_microswitches：CM/LY = 1.065 / 1.177 / 1.033
- orthosteric_activation_core：CM/LY = 1.255 / 1.089 / 0.972

在 `SIGNED_DRIFT` 中，CM00734、LY2119620 和 compound110 的 magnitude 也经常处于同一量级，而 R1/R3 directional agreement 整体较弱。

因此：

**Delta_INT magnitude alone is not specific.**

任何仅凭 magnitude 排序、阈值或单次轨迹幅度来宣称 PAM 的做法，都不受当前数据支持。

## 9. 最终科学结论

总体 outcome：

`PARTIAL_SPECIFICITY_SUPPORT`

当前最可辩护的结论是：

> Frozen PACER-FKG v02 shows partial hard-negative specificity primarily through cross-replica consistency of Delta_INT direction in the STATE_MOTION branch, rather than through Delta_INT magnitude alone.

中文解释：

优化后的冻结 C1-BS256 + PACER-FKG v02 能在当前 matched design 中区分“已知 PAM 的可重复方向模式”和“实验 inactive hard negative 的不一致方向模式”，尤其是在 `STATE_MOTION / Delta_INT / compound110_extension` 及若干机制相关区域。

但当前数据**没有**证明：

- 一个通用 PAM classifier 已建立；
- 存在可直接使用的 magnitude threshold；
- `Delta_AGO` 等价于真实 intrinsic agonism；
- docking affinity 能推断 efficacy；
- temporal blocks 是独立 biological replicates；
- 当前 representation 可以预测 potency；
- 任意新候选具有 PAM 活性。

## 10. 统计与验证边界

- R2 是历史 calibration/fitting replica；
- R1/R3 是冻结外推解释的重点；
- 20 个 temporal blocks 是相关时间样本，不是独立 replica；
- 不进行基于 blocks 的独立样本 p-value 解释；
- Stage B 后禁止根据 CM00734 outcome 回头调整 encoder / normalization / RFF / graph / regions / thresholds。

## 11. 对下一阶段的含义

当前模型开发阶段应结束于冻结结果，而不是继续针对这两个化合物优化。

下一步若继续，应采用预注册式扩展：

1. 在运行前冻结新的 positive / hard-negative 清单；
2. 优先加入独立实验 inactive binders / near-neighbor negatives；
3. 加入独立已知 PAM；
4. 继续使用同一冻结 encoder 和 PACER-FKG numerical state；
5. 预先规定 directional reproducibility 的描述规则；
6. 在独立前瞻样本复现前，不建立 classifier threshold 或候选排名。

## 12. 最终项目定位

PACER-DC / PACER-FKG v02 当前已经从“encoder 可行性探索”进入：

**frozen mechanistic representation with preliminary hard-negative specificity support**

而不是：

**validated PAM prediction model**。
