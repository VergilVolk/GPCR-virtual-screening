# PACER-DC 闭环实际执行计划 v01

## 一、现在项目到底处于什么状态

当前不能把 `main` 当作执行权威。

目前真正包含最新工作的是：

`codex/pacer-dc-geom2vec-pilot`

这个分支已经包含并整合：

- encoder A–H 优化；
- frozen encoder candidate；
- PACER-FKG v02；
- 600 ns common-kernel PACER-FKG v03 最终分析；
- long-MD replica/specificity audit；
- 新的 close-loop contract。

`main` 明显落后于当前研究状态，因此接下来的闭环实验必须从当前活跃分支再建立一个专门 execution branch，而不是切回 main。

当前已经完成的科学资产如下：

| 模块 | 当前状态 |
|---|---|
| compound110 600 ns MD | COMPLETE |
| 原版 128D Geom2Vec long-MD | COMPLETE |
| common-kernel PACER-FKG v03 | COMPLETE / FROZEN |
| common-kernel replica audit | COMPLETE |
| common-kernel specificity audit | COMPLETE |
| encoder optimization A–H | COMPLETE |
| C1-BS256 encoder candidate | FROZEN |
| PACER-FKG v02 STATE_MOTION/SIGNED_DRIFT | COMPLETE / FROZEN |
| LY2119620 long-MD | NOT STARTED |
| long-MD MCV baseline | NOT DONE |
| CM00734 runnable MD system | NOT YET BUILT |
| three-class closure | NOT DONE |
| supervised training | GATE CLOSED |

所以项目已经不再是“继续开发 PACER-FKG”。

现在的任务已经变成：

**冻结方法 → 拿新的已知分子做外部回顾性测试。**

---

## 二、两条 PACER-FKG 线都必须保留

现在事实上存在两条合法、冻结而且含义不同的分析线。

第一条是我们刚完成的：

**common-kernel v03 / historical 128D representation**

它使用三 replica，当前 compound110 的主要结论是：

`dINT / compound110_extension`

最稳定，而 `dAGO` 并不稳健。

第二条是 9 月 29 日完成的：

**PACER-FKG v02 / C1-BS256 / STATE_MOTION**

其冻结 encoder 为：

`C1 → backbone/sidechain → BS256`

数值状态为：

- encoder anchor commit：`0a4b4f2...`
- encoder tag：`encoder-candidate-ah-frozen-20260929`
- v02 freeze SHA256：`b48bc74a757a3d1421acb5c4bc0544ce5ca590e3959b971dad698e32146e11bd`
- STATE_MOTION bandwidth：`32.30188361260893`
- STATE_MOTION RFF seed：`272340`
- SIGNED_DRIFT bandwidth：`29.29934899925964`
- SIGNED_DRIFT RFF seed：`272084`

compound110 在这条线上表现为较强的 `ΔAGO / ΔINT` R1/R3 direction consistency，而 `ΔPAM` 不稳定。

这两条分析**不能选结果好看的那一条**。

后面的 LY2119620、CM00734 必须：

**两条都跑、两条都报告。**

如果两条结论一致，这是加强证据。

如果两条结论不一致，分歧本身就是结果。

---

## 三、实际关键路径

1. **先冻结一个新的闭环执行锚点，不再改方法。**

   从当前最新活跃分支创建专门的：

   `experiment/pacer-dc-close-loop-v01`

   并自动生成一个 `CLOSE_LOOP_EXECUTION_MANIFEST_v01.json`。

   其中至少冻结：

   当前 Git HEAD、encoder candidate commit/tag、ViSNet checkpoint SHA、Geom2Vec source commit、common-kernel v03 contract SHA、common-kernel final freeze SHA、FKG v02 freeze SHA、MCV config SHA、3 个 seed、50 ns、50 ps/frame、20 frames/block。

   从这一刻开始，LY 或 CM 的结果无论好坏，都不能再修改 encoder、normalization、bandwidth、RFF、graph、region 或 contrast definitions。

2. **立即做 LY2119620 production preflight。**

   这一步不是重新构建 LY。

   当前仓库的 MD 管线已经原生支持：

   `LY2119620__candidate_no_probe`

   `LY2119620__candidate_probe`

   现有 production runner、membrane-reference builder、restraint-release 路线都已经包含 LY。

   所以现在真正需要检查的是 GPU/本地结果盘上的实际资产：

   `pacer_dc_membrane_reference_v01/LY...`

   和

   `pacer_dc_restraint_release_v01/LY...`

   必须自动核对：

   `minimized.pdb`、`system.xml`、`production_start_state.xml`、`audit.json`、SHA256、`production_start_eligible=true`。

   更重要的是做 control compatibility audit，比较 LY 与现有 apo/probe_only 的：

   force field、受体来源、OPM transform、膜/水构建、离子条件、去约束协议、box、sequence/mapping、积分步长、barostat、seed 定义和采样间隔。

   如果全部兼容：

   **复用现有 apo/probe_only 600 ns controls。**

   因此只需要新增：

   `LY-no-probe × R1/R2/R3`

   `LY+probe × R1/R2/R3`

   即 **6 × 50 ns = 300 ns**。

   如果 compatibility FAIL：

   不允许偷偷复用现有 controls。

   那时才需要重新建立同 lineage 的 apo/probe controls。

3. **在 LY MD 运行之前，把“新分子应用器”写好并冻结。**

   这是现在代码层最大的实际缺口。

   目前很多分析 runner 都把 `compound110` 硬编码在里面，所以不能直接拿去分析 LY。

   但我们不能去修改已经冻结的历史 runner。

   应新建一个独立目录：

   `project/pacer_close_loop_v01/`

   在这里建立三个 apply-only 工具。

   第一套是：

   `apply_common_kernel_v03_candidate.py`

   它只允许读取 frozen v03 preprocessing、bandwidth、RFF、graph 和现有 apo/probe controls。

   对新候选只能生成 candidate/no-probe 和 candidate/probe 两组特征，然后代入：

   `dAGO`

   `dPAM`

   `dINT`

   禁止 fitting。

   第二套是：

   `extract_frozen_bs256_candidate.py`

   用冻结的 C1-BS256 encoder 提取 LY 的 `[1000,270,256]`。

   第三套是：

   `apply_frozen_fkg_v02_candidate.py`

   直接读取 SHA 为 `b48bc74a...` 的 v02 numerical state。

   必须明确：

   **绝对不能对 LY 再运行 Phase 2 calibration。**

   LY 是 application molecule，不是 calibration molecule。

   这三个脚本必须先通过 synthetic signs、shape、RFF/hash、no-fitting-path 测试，然后在看到 LY biological outcome 前冻结。

4. **随后立即启动 LY 300 ns，同时并行完成 MCV long-MD baseline。**

   LY production 使用已有 runner。

   参数必须显式固定为：

   R1 = `27101`

   R2 = `38201`

   R3 = `49301`

   `50 ns`

   `50 ps/frame`

   注意这是一个真实代码细节：

   现有 `run_pacer_dc_production_md.py` 的默认 trajectory spacing 是 **10 ps**。

   因此这次必须显式传：

   `--trajectory-ps 50`

   否则会和 compound110 600 ns 数据采样合同不一致。

   LY 输出必须写入全新的 close-loop result root，不能复用或覆盖旧 `pacer_dc_production_v01`。

   与此同时 CPU 跑 MCV。

   这里也发现了一个重要问题：

   现有 `run_pacer_mechanism_cv.py` 是旧的 R2/R3、5-window runner，并且代码会强制：

   `replicas == [2,3]`

   所以不能简单把 600 ns DCD 塞进去。

   正确做法是新建：

   `run_pacer_mechanism_cv_long_md_v01.py`

   保留完全相同的 frozen physical feature config，但改成长程合同：

   `R1/R2/R3`

   `1000 frames`

   `50 × 20-frame = 50 个 1-ns blocks`

   R2-only scaling 如需延续旧 MCV calibration policy，则先冻结 scaling，再报告三 replica；不允许看结果以后改 feature panel。

5. **LY 完成后形成 Stage A，而不是再优化方法。**

   每条 LY 轨迹首先做与 compound110 相同的：

   completion、frame count、finite、temperature/density、PBC、270-residue mapping、atom14 QC。

   全部通过后才进入分析。

   对 LY 同时生成：

   historical 128D features

   和

   frozen C1-BS256 features。

   然后得到三套证据：

   `common-kernel v03`

   `FKG v02`

   `MCV physical baseline`

   common-kernel 报 R1/R2/R3。

   FKG v02 保持其原始设计：

   R2 = calibration-role anchor

   R1/R3 = frozen application replicas。

   LY 是否“符合我们预期的 PAM 模式”不是 Stage A 是否完成的条件。

   **只要预注册分析全部跑完，无论结果支持还是反对假设，Stage A 都完成。**

   Stage A 要回答的只有：

   > frozen framework 是否在 compound110 与已知 PAM LY2119620 上产生可区分的动态模式？

6. **LY 占用 GPU 时，CPU/工程线提前准备 CM00734，但不要马上跑长 MD。**

   这是当前合同中还不够详细的地方。

   当前仓库里：

   CM00734 已经存在于 label、benchmark 和 MD manifest 中；

   但它**还没有进入现有 membrane-reference / restraint-release / production runner 的可运行系统列表**。

   因此 CM00734 实际上不能直接执行“6 × 50 ns”。

   必须先建立完整 system-preparation lane：

   pose/source freeze

   → static pocket compatibility

   → ligand parameterization

   → minimization QC

   → membrane reference construction

   → short equilibration

   → staged restraint release

   → zero-restraint production start audit

   → short binding/pose stability QC

   → production eligibility。

   只有这些通过后，它才能成为真正的 hard-negative MD experiment。

   特别注意它现在的合法称呼是：

   **实验无功能 hard negative，binding compatibility 待验证。**

   不能提前称为“已确认无功能 binder”。

   如果 CM00734 binding compatibility PASS：

   启动同样的 **6 × 50 ns = 300 ns**。

   然后使用和 LY 完全一样的两条 frozen FKG pipeline + MCV。

   如果 FAIL：

   这个结果本身必须记录。

   CM00717 可以接替下一张 GPU 队列，增加 PAM positive evidence，但：

   **CM00717 不能代替 CM00734 完成 negative class closure。**

   也就是说，CM00734 FAIL 时不能宣布“三类闭环成功”。

7. **最后才做 Stage B 三类闭环，不启动训练。**

   Stage B 最终矩阵应是：

   | 类型 | 分子 | 目的 |
   |---|---|---|
   | ago/interacting control | compound110 | 已完成 |
   | known functional PAM | LY2119620 | Stage A |
   | inactive hard negative | CM00734 | Stage B |

   每个分子统一报告：

   `dAGO`

   `dPAM`

   `dINT`

   replica direction consistency

   stable/distal specificity

   STATE_MOTION

   SIGNED_DRIFT

   MCV physical baseline。

   不构造一个临时“PAM 总分”。

   不重新选 region。

   不重新拟合 kernel。

   不用 block 当 biological N。

   不根据结果改变成功标准。

   三类无论是否按照我们期待的方向分开，只要都有完整、冻结、可审计的数据，就形成真正的 retrospective closure。

   此时生成：

   `PACER_DC_CLOSE_LOOP_FINAL_REPORT_v01.md`

   `CLOSE_LOOP_FINAL_FREEZE_v01.json`

   全部脚本、manifest、small results、SHA256、commit/tag 一并进入 GitHub。

   **训练仍然不应该启动。**

   当前真正的 `audit_pacer_dc_training_gate.py` 比旧计划中的“≥2/≥2/≥2”更严格：它还要求至少 6 个独立训练分子，并且数据规模足以做 molecule-level train/validation/test。

   所以即使 compound110 + LY + CM00734 三类闭环成功，也只是完成 framework validation，不等于 training gate 打开。

---

## 四、最快的并行方式

GPU 长任务严格串行：

**LY2119620 → CM00734**

CPU/工程任务并行：

**execution freeze + LY preflight + generic application scripts + MCV long-MD + CM00734 system preparation**

也就是说，GPU 一旦开始 LY，就不再等代码开发。

LY 在跑的时候，把后面所有分析器和 CM00734 前处理全部准备好。

这样 LY 一结束即可立即分析，随后 GPU 无缝进入 CM00734。

---

## 五、哪些事情现在不要做

现在不要：

重新优化 encoder；

重新 fit PACER-FKG bandwidth/RFF；

运行 v01 exact replay；

继续做 v1-v2 superiority；

启动 PACER0076/PACER0057；

启动 dual-head training；

为了得到漂亮结果修改 region 或 cosine threshold。

这些都不在当前闭环 critical path。

---

## 六、真正的终点定义

Stage A：

**compound110 + LY2119620 + MCV baseline**

完成后可以说：

> 已用冻结方法对两种已知功能类型进行了回顾性差异测试。

Stage B：

**Stage A + 可用的 CM00734 hard negative**

完成后才能说：

> 同一冻结四上下文框架已经在已知 PAM、ago/interacting ligand 和实验无功能 hard negative 上完成三类回顾性判别测试。

这里的“闭环成功”不是要求结果一定符合我们的期待。

如果 LY 没有出现预期 PAM pattern，或者 CM00734 也出现很强的 functional-like signal，那仍然是完整闭环，只是结果会否定或削弱当前框架。

这比“为了闭环必须得到阳性结果”科学上强得多。