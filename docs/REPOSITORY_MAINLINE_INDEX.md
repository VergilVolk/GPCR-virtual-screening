# 仓库主线索引（REPOSITORY MAINLINE INDEX）

更新：2026-10-01 ｜ 权威分支：`codex/project-integration-benchmark-v02`
（本地尖端 `77f44dc6`；全部业务分支内容已并入，无独有提交）

## 一、主线是什么

**四模块串行级联（官方术语，2026-10 定名），有拒绝机制的证据链**：

```
模块1  M4 PAM候选分子生成与药化筛选（片段生成+药化过滤+新颖性 → 候选池）
模块2  基于DrugCLIP迁移学习的M4智能靶点筛选算法
       （微调线：13靶 famaug 集成 ROC 0.6414；M4 决胜配置 2023基座+mean池化
        ROC 0.7136/BEDROC 0.1802；提交包 submission/ = 本模块）
模块3  静态-动态构象组合分子对接算法（晶体+GaMD系综 Vina + IFP/PACER-FS，
       pose proposal 不判效力）
模块4  基于四条件受体动力学的PAM功能判别算法（A/P/C/CP 20ns → ΔINT
       跨replica方向一致性 → 证据卡，可拒绝）
```

逐级闸门漏斗：28,519 → ~1,400（模块2 top-5%）→ 几十（模块3 pose 通过）
→ 5 个优先验证候选（模块4 幸存）。模块3 失败退回模块2 取下批切片。

**当前最强冻结主张**：预冻结 hard-negative 检验——已知 PAM（LY2119620 /
compound110）跨 replica 方向一致性为正（R1/R3 cosine +0.63/+0.58），实验
inactive 高相似近邻 CM00734（ECFP4 Tanimoto≈0.955）为负（−0.327，9/9
区域全负）；判定 `PARTIAL_SPECIFICITY_SUPPORT`。零 outcome-driven 调参。

## 二、主入口文件（按用途）

| 用途 | 文件 |
|---|---|
| 项目当前状态（最权威） | `project/PROJECT_STATUS.md` |
| 整体算法+基准（v02 报告源） | `project/docs/PACER_M4_INTEGRATED_ALGORITHM_AND_BENCHMARK_V02.md` |
| 生成 PDF 报告 | `project/scripts/build_pacer_m4_algorithm_report_v01.py` → `output/pdf/` |
| Stage A+B 合成报告 | `project/results/pacer_dc_stage_a_b_final_synthesis_v01/…FINAL_SYNTHESIS_v01.md` |
| 结合线交接 | `project/docs/DRUGCLIP_FINETUNING_PIPELINE_INTEGRATION.md`、`DRUGCLIP_SCIENCE2026_FINETUNING_HEADLINE.md` |
| 模块计分板/级联漏斗图 | `project/results/project_wide_integration_benchmark_v02/figures/` |
| 项目级基准套件 | `project/results/project_wide_integration_benchmark_v02/` |

## 三、工作树与分支地图

| 工作树 | 分支 | 状态 |
|---|---|---|
| `D:\CLC` | `codex/drugclip-pacer-handoff`（已 ff 到 v02 尖端） | **主工作树，主线内容** |
| `D:\CLC_integration_benchmark` | `codex/project-integration-benchmark-v02` | 集成线专用 |
| `D:\CLC_geom2vec_pilot` | `integration/project-benchmark-v01` | **已被 v02 取代**，仅存档 |
| `D:\CLC\.codex-worktrees\oneprot-*` | oneprot 审计/后续 | 内容已在 v02，worktree 存档 |

所有 `experiment/*`、`audit/*`、`codex/*` 业务分支：内容均已并入 v02
（逐一验证 `git log v02..<branch>` 为空）。分支保留作历史线索，不再开发。

## 四、各线冻结结论速查

| 线 | 冻结结论 | 关键文档 |
|---|---|---|
| DrugCLIP 微调（结合） | 13 靶 famaug 0.6414（CI [+0.080,+0.114]）；20 靶不叠加；EF1% 不稳定优于 2D 基线 | `DRUGCLIP_SCIENCE2026_FINETUNING_HEADLINE.md` |
| 管线融入 | 28,519 库 + PACER-200 重评分；新旧权重排序独立；双优交集 25 进功能复核；M4 决策线守旧 2023 权重 | `DRUGCLIP_FINETUNING_PIPELINE_INTEGRATION.md` |
| PACER-FKG v02（Stage A） | matched 20 ns：主信号 `STATE_MOTION/ΔINT/compound110_extension` R1/R3 cosine +0.634 | Stage A 报告（frozen tag） |
| PACER-FKG v02（Stage B） | hard-negative CM00734：−0.327，9/9 区域负；magnitude 重叠 → magnitude 不能单独判别 | tag `pacer-dc-cm00734-stage-b-20ns-v01` |
| encoder 优化 | C1-BS256 冻结（原始 Geom2Vec/ViSNet 权重不微调）；锚点 SHA256 b48bc74a… | encoder-opt 审计报告 |
| OneProt-PaCeR-DC G0-G4 | 局部四上下文差分跨 replica 不稳定 → 退役为消融基线；H 验证已交接 | `ONEPROT_PACER_DC_BRANCH_REVIEW_20260921.md` |
| CGDA | 侧线关闭（7 小时分层子集门控决策）；LIT-PCBA LOTO ROC 0.5788 | `DRUGCLIP_CGDA_RESULT.md` |
| PACER-MCV | 物理构象基线已实现；R2/R3 原始 atom14 缺失，跨 replica 实值未产生 | `PACER_MCV_METHOD_AND_HANDOFF.md` |
| 四上下文 MD 工程 | 六体系×生产起点 6/6；iperoxo 稳定化 +0.651 Å [0.527,0.762]（机制支持，非分类） | `PROJECT_STATUS.md` 2026-09 节 |

## 五、下一执行点（PROJECT_STATUS 冻结）

1. 结果归档 + **独立审计**（不重跑、不调参）；
2. 下一轮**预注册**新增 inactive/binder hard negatives 与独立已知 PAM；
3. 只有独立前瞻样本复现方向性 discrimination 后，才谈 classifier/threshold/排名。

## 六、红线（整理时不得破坏）

- Stage A/B 冻结代码、receipt、result SHA、tag 不动；
- CM00734 结果后禁止任何 outcome-driven 重调 encoder/normalization/graph/threshold；
- 不得把不同协议数字拼成总 AUC；不得称端到端 PAM SOTA；
- 未追踪大数据（`project/data/*` LMDB/npz）不入 git，由各自 audit JSON 登记。
