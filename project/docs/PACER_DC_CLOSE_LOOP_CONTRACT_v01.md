# 四上下文闭环合同 CONTRACT v01（2026-09-29 冻结）

取代 `PACER_DC_CLOSE_LOOP_PLAN_20260929.md`（同日 v01 草案）。
执行中不得依据 outcome 回改本合同的判读定义；变更需新版本号并说明理由。

## 0. 目标

把 PACER-DC 四上下文从"compound110 单体系信号"推进到
**形成框架层面的回顾性判别证据**：区分 functional PAM / ago-PAM / 无功能 binder
三类已知配体。本合同不产生药理效力预测、候选确认或 PAM 身份主张。

## 1. 两条并存的冻结基线（判读的出发点，不得择一引用）

同一批 600 ns 长程 MD（4 体系 × 3 replica × 50 ns）上存在两条独立冻结分析线：

| 线 | 冻结日 | 表示 | 关键结论 |
|---|---|---|---|
| common-kernel（CONTRACT_v03 / FINAL_FREEZE_v01） | 09-28 | 128D 历史表示 | dINT/compound110_extension 最强：pooled pairwise cosines R1R2/R1R3/R2R3 = 0.568/0.661/0.172；dAGO 三 replica 不稳；distal 对照 2/3 为正 |
| FKG v02 STATE_MOTION（C1-BS256） | 09-29 | 冻结 encoder candidate 0a4b4f2 | ΔAGO/ΔINT 于 compound110_extension（R1/R3 应用副本）= 0.78 / 0.62；ΔPAM = −0.52 不稳 |

v1-v2 对比审计未执行；两线互不替代。**本合同下的所有新分析必须两条轨道都跑、
都报告，禁止依据结果择优引用。**

## 2. 里程碑定义

- **Stage A（presentation closure，10-04 前）**：compound110（已有）+
  LY2119620 镜像 + MCV 物理基线对照。产出"框架在两类已知配体上呈现相反轴模式
  + 物理基线位置"的回顾性证据。
- **Stage B（full three-class closure，10-10 区赛前）**：Stage A +
  CM00734 硬负样本（无功能 binder）。三条都有数（无论方向）才算三类判别闭环。
- CM00717（第二 PAM，chemotype 泛化）、PACER0076/0057（候选）、training-gate
  复评全部在 Stage B 之后，不进关键路径。

## 3. P0 — LY2119620 镜像测试（预注册判读，先冻结后跑数）

**预注册判读标准**（模式级，非单一数字门）：

1. 判读单位：三 replica 全报告。报告 pooled pairwise cosines（R1R2/R1R3/R2R3）、
   coherence、相对 stable-core 与 distal 对照的 pooled excess、block fraction
   ——完全沿用 FINAL_FREEZE_v01 的口径；v02 轨道另报 R1/R3（R2 为其校准 replica
   的既定设计）。
2. 镜像成立 = LY2119620 满足：**PAM 轴（dPAM/dINT，pam_contact 区域）的
   跨 replica 方向一致性高于其自身 dAGO**，且轴序与 compound110 相反
   （compound110 现状：common-kernel 线 dINT>dAGO、v02 线 ΔAGO≫ΔPAM）。
3. 方向余弦 0.50 仅作描述性参照（源自 ENCODER_DIAGNOSTIC §五 与 MCV 冻结判据
   的项目惯例），**不作为 pass/fail 门**；qualification_gate 维持 NOT_DEFINED。
4. 两条轨道（common-kernel、FKG v02）都计算；结论按"两线一致 / 两线分歧"
   分层表述，分歧时如实报告分歧。

**执行**：
- Preflight（先于 MD）：LY2119620 ±probe restraint-release state SHA 核对；
  **manifest compatibility audit**——LY 生产起点与既有 apo/probe_only controls 在
  force field、膜/水、受体构象来源、去约束协议、seed 定义、box、采样节奏上的
  兼容性。PASS 才允许复用既有对照轨迹并称"paired"；FAIL 则 LY 需自建全套四上下文
  （含自身对照），工期顺延并如实标注 lineage。
- 协议：3 配对 seed（27101/38201/49301）、50 ns/条、50 ps 存帧、同 atom14 映射，
  共 6 × 50 ns；分析复用两条线的冻结数值状态（common-kernel 合同 + v02 freeze
  manifest b48bc74a…），**不得重新校准**。

## 4. P1 — CPU 并行（不占 GPU，不阻塞 P0）

1. **MCV 600 ns 基线**：对既有 12 条轨迹跑 `run_pacer_mechanism_cv.py`。
   注意：DCD 在 GPU 机；本机执行前需先传轨迹或直接在 GPU 机跑（二选一，
   记录执行位置）。输出与两条 FKG 轨道同口径的跨 replica 方向表。
2. v01 exact replay / v1-v2 对比审计：**暂缓**，移出关键路径（Stage B 后再做）。

## 5. P2 — CM00734 硬负样本（GPU 第二队列）

1. 启动前置：binding compatibility 核验（v2 清单遗留项）。
   **PASS → 6 × 50 ns**（协议同 P0）；
   **FAIL → 该体系无法启动，GPU 队列改排 CM00717**（其体系准备与 CM00734
   并行进行，Plan B 常备），Stage B 以 CM00717 完成时点重估。
2. 判读预期（预注册）：CM00734 作为无功能 binder，dPAM/dINT 均不出现
   跨 replica 稳定方向信号；出现稳定信号即为对框架的否定性证据，如实报告。

## 6. 日程映射

- 今日：冻结本合同 → LY preflight audit → PASS 即启动 LY 300 ns；
- 本周内：MCV 基线（CPU）；
- ~10-02 至 10-04：LY 分析（两条轨道）→ Stage A 交付；
- 10-04 后：CM00734（或 fallback）→ Stage B，区赛前；
- 决赛前：CM00717 / 候选 / 训练门。

## 7. 主张边界

- 全部结论限于四上下文框架的回顾性判别证据；
- 统计单位是 replica / 分子；50 个 1-ns blocks（20 frames/block）是时间子样本；
- 镜像成立只主张"框架区分了两类已知配体"；
- 候选只输出证据向量；确认 PAM 数保持 0 直到湿实验。

## 8. 责任分工

- GPU 机（4070 笔记本）：LY / CM00734 的 MD、BS256 提取、（可选）MCV；
- 本机 geom2vec 工作树：MCV（若传轨迹）、合同权威、文档；
- 本机 binding 分支：LIT-PCBA 15 靶点验收门（看门狗护航，独立推进）。
