# 全项目联调与基准跑分报告 v01

日期：2026-10-01。执行分支：`integration/project-benchmark-v01`（基于
`experiment/pacer-dc-cm00734-stage-b-20ns-v01` @ 5e86149b）。

本报告把三条线（结合 / 功能 / QSAR-静态）的全部已提交结果汇集为一套可审计的
联调记录与基准图表。所有数字取自仓库内已提交的 JSON/CSV/冻结文档，
每个数字带出处；本机无法执行的项如实标注 SKIP，不做任何虚构覆盖。

## 1. 联调结果

### 1.1 分支与冻结锚点（9/9 PASS）

| 检查 | 期望 | 实测 | 判定 |
|---|---|---|---|
| tag encoder-candidate-ah-frozen-20260929 | 0a4b4f23 | 0a4b4f23 | PASS |
| tag pacer-fkg-v02-phase2-frozen-20260929 | ecaf7a60 | ecaf7a60 | PASS |
| tag pacer-dc-ly2119620-close-loop-20ns-v01 | b719f01e | b719f01e | PASS |
| tag pacer-dc-cm00734-stage-b-20ns-v01 | f2f73403 | f2f73403 | PASS |
| tag pacer-dc-stage-a-b-synthesis-v01 | 5e86149b | 5e86149b | PASS |
| tag pacer-fkg-long-md-v01 | df219d70 | df219d70 | PASS |
| Stage B 引用的 FKG v02 冻结 SHA | b48bc74a… | b48bc74a757a | PASS |
| TRACK_B 上游冻结 SHA | b48bc74a… | b48bc74a757a | PASS |
| TRACK_B 禁止操作全 false | 全 false | 全 false | PASS |

分支清单（按提交时间）：`codex/drugclip-pacer-handoff`（结合线，78d4e6a7）、
`experiment/pacer-dc-cm00734-stage-b-20ns-v01`（Stage B+合成，5e86149b）、
`experiment/pacer-dc-close-loop-20ns-analysis-v01`（853aaa12）、
`experiment/pacer-dc-close-loop-20ns-v01`（6d12bcef）、
`codex/pacer-dc-geom2vec-pilot`（14d86967）、
`experiment/pacer-fkg-v02-longmd-v01`（07ce8ae1）、
`experiment/pacer-encoder-opt-v01`（0a4b4f23）、
`audit/oneprot-pacer-dc-g0`、`main`（c2fb61f3）。

### 1.2 项目自带测试套件

| 套件 | 结果 |
|---|---|
| test_path_guard.py | **4 passed + 20 subtests passed** |
| test_phase2_calibration.py | **7 passed** |
| test_phase1_bs256.py | SKIP——需要 torch_geometric（仅 GPU 机装有 ViSNet 栈） |
| test_phase3_fkg.py | SKIP——冻结锚点验证逻辑正确地拒绝在缺 `V02_FREEZE_MANIFEST` 工件的克隆上运行（工件按设计 git-ignored，在 GPU 机） |

两个 SKIP 都是环境性缺资产，不是逻辑失败；phase3 的拒绝行为本身是
冻结纪律按设计工作的证据。

### 1.3 本机联调执行到的真实入口

- 汇编脚本：`project/scripts/assemble_project_wide_benchmark_v01.py`（本报告全部数据来源）
- 图表脚本：`make_benchmark_figures_v01.py` + `make_benchmark_figures_supp_v01.py`
- 历史分支数据提取：FGK v01 common-kernel 四件套经 `git show` 从
  `codex/pacer-dc-geom2vec-pilot` 提取（UTF-8 校验通过）

## 2. 基准跑分：结合线（DrugCLIP，13/20 靶严格 LOSO）

宏指标全表见 `csv/binding_macros.csv`。要点（ROC-AUC）：

| 方法 | 13 靶 | 20 靶 |
|---|---:|---:|
| 官方 Science-2026 raw | 0.5442 | 0.5816 |
| ECFP4 logistic | 0.5681 | 0.5332 |
| 随机标签对照 | 0.4962 | 0.5144 |
| 微调 ep40（3种子集成） | 0.6121 | — |
| 微调 ep80（3种子均值） | 0.6062 | 0.6259 |
| **family-augmented ep80（3种子）** | **0.6400** | — |
| **family-augmented 集成** | **0.6414** | — |
| 纯 GPCR 域内（4靶，M4 探针） | 0.6003→0.6887 | — |

要点判读：
1. family augmentation 是 13 靶协议新头条（0.6414，三种子 0.6353–0.6434）；
2. 20 靶全覆盖：微调 0.6245/0.6255/0.6278 vs 官方 0.5816，CI [+0.028, +0.062]；
3. EF1% 边界保持：13 靶 ECFP 2.38 vs 微调 1.26；20 靶 ECFP 1.94 vs 1.97/1.75——
   早期富集仍无优势（fig12 显示官方 delta 的 EF1% CI 为负）；
4. 纯 GPCR 域内微调使 M4R 探针面板 0.6003→0.6887，确认"异源共训伤 M4"机制。

图：fig01（宏指标柱状）、fig02（13 靶×方法逐靶热图）、fig03（种子稳定性）、
fig12（delta bootstrap CI）。

## 3. 基准跑分：功能线（四上下文，matched 20 ns 三类闭环）

核心结果（STATE_MOTION / ΔINT 的 R1/R3 方向余弦，fig04）：

| Region | compound110 | LY2119620 | CM00734 |
|---|---:|---:|---:|
| compound110_extension | +0.584 | **+0.634** | **−0.327** |
| cooperativity_mutagenesis | −0.246 | +0.235 | −0.095 |
| intracellular_microswitches | −0.111 | +0.180 | −0.006 |
| orthosteric_activation_core | −0.074 | +0.321 | −0.385 |
| orthosteric_contact_union | −0.052 | +0.300 | −0.376 |
| pam_contact_consensus | −0.308 | +0.375 | −0.251 |
| pam_contact_union | −0.183 | +0.306 | −0.241 |
| distal_control / stable_core | （TRACK_B 中有值） | 有值 | 未列（合成表仅 7 region） |

判读：
1. CM00734 在全部报告 region 方向余弦为负，两个功能分子在机制相关 region
   为正——**方向可复现性（而非量级）区分功能与无功能**（fig05 右图：CM/LY
   量级比 0.53–0.80，无分离）；
2. SIGNED_DRIFT 分支整体弱于 STATE_MOTION（fig06 全景）；
3. LY 与 compound110 的直接方向一致性在 ΔINT/compound110_extension 上
   R1/R3 双正（fig13），两个功能分子在别构延伸区推动受体的方向一致；
4. FKG v01 common-kernel（600ns 历史线）：dINT/compound110_extension
   mean pairwise cosine 0.467（fig07），与 v02 线互为独立佐证；
5. 表示演进全景（fig08）：OneProt 全失败 → C0-M128 0.467 → C1-BS256
   c110 0.584 / LY 0.634——三代表示的跨 replica 方向稳定性单调改善。

数据边界：CM00734 逐 region 原始 JSON 在 GPU 机（git-ignored），本报告用
合成文档表格数字（出处已标注）；ΔPAM 在 20ns 不稳（两分子皆然），不做
分类端点主张。

## 4. 基准跑分：QSAR 与静态结构线（历史，冻结参照）

1. **QSAR 泛化崩塌**（fig09）：PAM-vs-inactive 在 scaffold 切分 ROC 0.87–0.92，
   换到 source/series 切分跌到 0.57–0.63、MCC≈0——8 月立项动机的量化证据；
2. **静态结构失效**（fig10 左）：VinaOnly 宏 Spearman 0.043，结构融合全部
   不敌 2D（LightGBM 0.199）；PACER-FS 系列中心化 0.263 为静态线最优
   （CI [+0.017, +0.113]）；
3. **外部大库筛选**（fig10 右）：Miao-2026 十万分子级，各法 ROC 0.65–0.75
   但 EF1%≈随机——高 AUC ≠ 早期富集，与结合线 EF1% 边界互证。

## 5. 全项目总览

fig11 把 14 个头条数字按线归一着色。一句话结论：

> 结合线（0.54→0.64，CI 为正，边界在 EF1%）与功能线（三类方向判别成立，
> 边界在无阈值无分类器）各自的贡献与边界都已量化；QSAR/静态线的历史失败
> 构成两条新线存在理由的完整证据链。

## 6. SKIP 清单（如实记录）

| 项 | 原因 |
|---|---|
| phase1/phase3 测试套件 | torch_geometric / V02_FREEZE_MANIFEST 工件在 GPU 机 |
| MCV 物理基线长程版 | 未执行（R2 scaling 工件未提前冻结，合同禁止事后补） |
| Track A（128D）LY 复算 | GEOM2VEC_LONG_MD_MANIFEST 缺失，apply-only 纪律禁止重建 |
| CM00734 逐 region 原始 JSON | 在 GPU 机结果目录（git-ignored） |
| 600ns DCD/BS256 缓存 | 同上 |

## 7. 复现

```powershell
python project/scripts/assemble_project_wide_benchmark_v01.py
python project/scripts/make_benchmark_figures_v01.py
python project/scripts/make_benchmark_figures_supp_v01.py
```

产物：`project_wide_benchmark_master_v01.json`（全部数字+出处）、
`csv/*.csv`（4 张平面表）、`figures/*.png`（13 张图）。

## 8. 图表索引

| 图 | 内容 |
|---|---|
| fig01 | 结合线宏指标柱状（ROC/PR/EF1%） |
| fig02 | 13 靶 × 7 方法逐靶 ROC 热图 |
| fig03 | ep80/famaug 三种子稳定性 |
| fig04 | **三类方向余弦热图（主结果）** |
| fig05 | 量级格子图 + CM/LY 量级比 |
| fig06 | 分支 × 对比 × panel 全景格子 |
| fig07 | FKG v01 三轴 × 九区热图（600ns 历史线） |
| fig08 | 表示演进三代对比 |
| fig09 | QSAR 切分崩塌柱状 |
| fig10 | 静态结构 + Miao 外部 |
| fig11 | 全项目 14 头条总览 |
| fig12 | 结合线 delta bootstrap CI |
| fig13 | LY vs c110 直接方向一致性 |
