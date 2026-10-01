# 四上下文闭环计划（2026-09-29 冻结草案）

> **状态：已被同日 `PACER_DC_CLOSE_LOOP_CONTRACT_v01.md` 取代。**
> 修订要点：双轨基线（common-kernel v03 合同 + FKG v02 并存）、三 replica 全报告、
> 0.50 仅作描述性参照、CM00734 提前、Stage A/B 分里程碑、v01 replay 移出关键路径。
> 本文件保留作历史草稿，判读以 CONTRACT v01 为准。

目标：把 PACER-DC 四上下文从"compound110 单体系有方向一致信号"推进到
"框架层面证明能区分 functional PAM / ago-PAM / 无功能 binder"。

依据：PACER-FKG v02 已在 600 ns 长程 MD（apo / ACh / compound110±ACh，R2 校准、
R1/R3 冻结应用）上得到 STATE_MOTION 的 ΔAGO/ΔINT 在 7/9 区域跨 replica 方向一致
（compound110_extension：ΔAGO cosine 0.78、ΔINT 0.62），ΔPAM 尚不稳定。
encoder candidate 已冻结（tag `encoder-candidate-ah-frozen-20260929`，commit 0a4b4f2），
FKG v02 数值状态已冻结（manifest SHA256 b48bc74a…）。

## 闭环判据（本计划的终点定义）

满足以下三条，四上下文闭环成立：

1. **镜像判别**：LY2119620（已知 PAM）在同一冻结管线下呈现与 compound110 相反的
   轴模式——ΔPAM（尤其 pam_contact 区域）方向稳定、ΔAGO 弱；
2. **binder 判别**：CM00734（实验非活性、与 CM00717 Tanimoto≈0.955 的近邻硬负样本）
   在 ΔPAM/ΔINT 上均不出现稳定方向信号；
3. **基线比较**：同批轨迹上 MCV 物理特征基线与 FKG v02 的可复现性对照完成，
   明确 learned 表示相对物理基线的位置（无论正负，都要有数）。

三条全部有数后，无论结果方向如何，都形成可写进报告的完整证据链；
不满足判据也是结果，按边界如实报告。

## 阶段与优先级

### P0 — LY2119620 镜像测试（最高优先级，先于一切新 MD）

理由：这是用最小算力买到最大判别力的实验——compound110 的 ΔAGO 0.78 已经是
"框架能抓到内在激动"的证据，镜像测试直接回答"框架能不能区分 PAM 与 ago-PAM"。

1. **先冻结镜像判读标准，再启动 MD**（防止事后解释）：
   - LY2119620 的 ΔPAM 在 pam_contact_consensus / pam_contact_union 区域
     R1/R3 方向余弦 ≥ 0.50 记为稳定；
   - 同体系 ΔAGO 不高于 compound110 的 ΔAGO 水平（对照既有 0.78）；
   - STATE_MOTION 为主分支，SIGNED_DRIFT 只作记录。
2. 生产起点：沿用已过分级去约束门的 LY2119620 ±ACh 零约束末态
   （restraint_release_v01 产物，运行前核对 state manifest 与 SHA）。
3. 协议与 compound110 长程完全一致：3 配对 seed（27101/38201/49301）、
   50 ns/条、50 ps 存帧、同 atom14 映射；共 6 条 × 50 ns = 300 ns。
4. 分析严格复用冻结 v02 数值状态：Phase 1 用同一 frozen C1-BS256 提取，
   Phase 2 **不得重新校准**（bandwidth/RFF seed/normalization 原样套用），
   Phase 3 以 LY2119620 的 R1/R3 为应用副本。
5. 共享对照 apo / probe_only 不重跑——直接复用 600 ns 现有轨迹。

### P1 — 不需要新 MD 的两个对照（CPU，与 P0 并行）

1. **MCV 物理基线**：对现有 12 条长程轨迹跑 `run_pacer_mechanism_cv.py`
   （600 ns 的 DCD 在 GPU 机上，atom14 提取脚本已在 long-MD 管线中验证），
   输出与 FKG v02 同口径的 ΔPAM/ΔAGO/ΔINT 跨 replica 方向表。
   这一步同时解决 MCV 长久以来"缺 R2/R3 原始坐标"的等待——长程数据里
   R1/R2/R3 全部在场。
2. **v01 exact replay + v1-v2 对比审计**：按 9/29 审计报告的停止点执行。
   若 v01 RFF 状态无法精确重放，按报告规定另立 matched-calibration comparator，
   不冒称历史 v01。

### P2 — 特异性对照矩阵（GPU 排队第二位）

1. CM00717（同源 PAM 正对照）±ACh：检验框架跨 chemotype 的 PAM 模式泛化；
2. CM00734（实验非活性）±ACh：binder-vs-function 硬负样本；
   两者各 6 条 × 50 ns，协议同 P0。
3. CM00734 使用前先完成其结合相容性核验（v2 清单中的遗留项）。

### P3 — 前瞻候选（视算力与赛程再启动）

1. PACER0076 / PACER0057 ±ACh 各 300 ns；
   若 10-04 前算力不足，明确推迟到区赛（10-10）后、决赛（10-16）前完成；
2. 候选判读只输出证据向量（CoupledShift / OrthostericStabilization /
   IntrinsicActivationRisk / BindingCompatibility / Uncertainty），
   不输出单一 PAM 分数，不称确认 PAM。

### P4 — 训练门复评（数据齐后）

对照 TRAINING_GATE 条目盘点：PAM 正例（LY2119620、CM00717）、
ago/agonist 正例（compound110，第二个激动阳性是否纳入 LY2033298 由数据覆盖决定）、
负例（CM00734 + 是否补第二个阴性）。若仍不满足 ≥2/≥2/≥2，记录缺口并保持门禁关闭；
满足则在 chemotype 隔离拆分下运行 `train_dual_context_heads.py` 与七行基线矩阵。

## 算力与日程（单张 RTX 4070 Laptop 现实约束）

- P0：300 ns，约 3–6 天（对照 compound110 600 ns 的实跑速率核实后更新）；
- P2：600 ns；P3：600 ns。三者串行约 12–24 天，不可能全部落在 10-04 前。
- 因此排序为：**P0 立即开始 → 10-04 交付用"框架 + encoder 优化 + compound110 信号 +
  （若赶上）LY 镜像"叙事；P2 力争 10-10 区赛前出数；P3 决赛前。**
- P1 全程 CPU（本机即可），本周内完成。

## 主张边界（沿既有口径，不因闭环成功而放宽）

- 全部结论限于"四上下文动态框架的回顾性判别证据"，不含药理效力预测；
- 50 帧块是时间子样本，统计单位是 replica / 分子；
- LY2119620 镜像成立也只能主张"框架区分了两类已知配体"，不能主张候选预测；
- 候选永远输出为"计算优先假设"，确认 PAM 数保持 0 直到湿实验。

## 资产核对清单（执行前逐项确认）

- [ ] LY2119620 ±ACh 生产起点 state 文件与 SHA（restraint_release_v01）
- [ ] 长程 MD 运行脚本（compound110 600 ns 同款）可直接换体系运行
- [ ] FKG v02 冻结数值状态 manifest（b48bc74a…）在新体系上只读套用
- [ ] 12 条既有轨迹 checksum 复核（Phase 0 报告已列）
- [ ] MCV 特征清单 config（pacer_m4_mechanism_cv_v01.json）未改动
- [ ] GPU 机与本仓库分支同步（experiment/pacer-fkg-v02-longmd-v01 及合并后的主线）

## 责任分工建议

- GPU 机（4070 笔记本）：P0 全部、P2、P3 的 MD 与 BS256 提取；
- 本机（D:\CLC_geom2vec_pilot）：P1 的 MCV 长程分析、v1-v2 审计、文档与仓库权威；
- 本机（D:\CLC）：binding 线 LIT-PCBA 15 靶点验收门（看门狗护航，独立推进，
  不与功能线抢占人力）。
