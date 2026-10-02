# PACER-M4 合作者阅读入口

更新时间：2026-10-02  
用途：帮助新合作者在不翻阅全部开发记录的情况下，理解项目问题、算法主线、已有结果、候选推理流程和证据边界。

## 1. 项目在解决什么问题

PACER-M4 面向人源毒蕈碱型乙酰胆碱 M4 受体正性别构调节剂（PAM）的计算筛选。

项目的核心判断不是“分子能否进入 M4 别构口袋”，而是连续回答三个不同问题：

1. 分子是否具有合理的理化和药化性质；
2. 分子是否可能与 M4 别构口袋形成合理结合；
3. 分子是否可能在乙酰胆碱（ACh）存在时产生正性协同，而不是仅能结合、没有功能，或在无 ACh 时自行推动受体激活。

这三个问题不能由一个 docking score 或一个分类概率替代。因此，PACER-M4 是分阶段证据级联，而不是把所有分数加权成一个总分。

## 2. 完整工作分为开发与推理

### 2.1 算法开发

开发阶段使用公开数据、已知 M4 配体和实验标签，建立并检验各模块：

- M4 PAM 功能数据整理与防泄漏基准；
- ligand-only QSAR、ECFP、Chemprop、ChemBERTa 等强基线；
- DrugCLIP 结合检索的 GPCR 迁移、模型路由和统一跑分；
- 单结构、多状态、GaMD ensemble 和六路 docking 的复现与审计；
- 四上下文 MD 与 PACER-FKG 动态功能表示；
- 已知 PAM、ago-PAM 和实验 inactive 近邻的回顾性功能特异性检验。

这些已知分子用于开发和验证算法，不是项目发现的新候选。

### 2.2 候选推理

推理阶段把冻结方法应用于未参与开发的新分子：

```text
分子生成或独立分子库
    -> 结构标准化与药化性质过滤
    -> DrugCLIP 路由结合检索
    -> docking / IFP / 多构象结构门控
    -> 适用域、化学多样性与不确定性筛选
    -> 少量高优先级候选四上下文 MD
    -> 冻结 PACER-FKG 功能复核
    -> 候选证据卡与实验建议
```

推理输出只能称为“计算优先候选”。没有湿实验时，不能称为已发现或已确认的 PAM。

## 3. 当前算法模块

| 模块 | 回答的问题 | 主要输入 | 主要输出 | 当前证据 |
|---|---|---|---|---|
| 数据与功能效力基准 | 模型是否真正跨系列泛化 | M4 PAM 功能实验记录 | 清洗数据、系列/时间/骨架划分 | 430 个精确效力分子、严格外推基准 |
| 分子生成与候选库 | 从哪里获得可筛选的新分子 | 已知骨架、片段、生成模型或商业库 | 标准化候选 SMILES/立体异构体 | 已完成多条生成路线与适用域审计 |
| 药化门控 | 分子是否具备基本开发合理性 | 候选结构 | PAINS、反应性、MW、logP、TPSA、QED、SA 等 | 可运行；不预测 PAM 功能 |
| DrugCLIP 路由 | 分子与 M4/GPCR 口袋是否匹配 | 分子和口袋表征 | 靶点内结合排序 | M4 使用 2023 GPCR-LOTO；其他 GPCR 融合 2023/2026 分支 |
| 静态与 ensemble docking | 姿势和口袋接触是否合理 | 候选、受体构象 | pose、能量、IFP、口袋覆盖 | 可作结构门控；不能直接预测 PAM 效力 |
| PACER-XR 六路级联 | 如何联合六种 docking 排序 | Glide/Vina 的 PDB、BEmin、BEavg 六路分数 | 六路秩共识与级联排序 | 公共 broad-AM benchmark 上改善全局与早期富集；不等于功能 PAM |
| PACER-FS | 新药化系列效力如何少样本标定 | 分子结构与少量功能锚点 | 系列内相对效力与校准值 | 三锚点协议下优于同架构绝对 QSAR；需要实验锚点 |
| 四上下文 MD | 候选与 ACh 是否产生非加和动态变化 | A、P、C、CP 四套匹配轨迹 | Delta_PAM、Delta_AGO、Delta_INT | MD 是动态数据来源，不是单独的功能分类器 |
| PACER-FKG | 动态变化是否跨 replica 同向重现 | 四上下文轨迹表征 | 区域级方向一致性和幅度 | 已知 PAM/ago-PAM 为正、inactive 近邻为负的回顾性初步证据 |
| 候选决策 | 哪些分子值得进入昂贵验证 | 上述各模块证据 | Pareto 候选、拒绝原因、证据卡 | 不使用未经验证的单一加权总分 |

## 4. 当前最重要的结果

### 4.1 DrugCLIP 结合检索

在统一的 4-GPCR 开发基准上，M4-safe 路由达到 ROC-AUC 0.735、PR-AUC 0.323、BEDROC20 0.398、EF5% 4.474。相对 2023 GPCR-LOTO，ROC、PR、BEDROC20 和 EF5% 的配对置信区间为正；EF1% 仍不确定。该结果支持结合排序改进，不支持功能性 PAM 或效力预测结论。

### 4.2 六路 docking 复现与 PACER-XR

六路分别为 Glide PDB、Glide BEmin、Glide BEavg、Vina PDB、Vina BEmin 和 Vina BEavg。PACER-XR 将每一路转换为同靶点秩百分位，再进行等权共识；Cascade 先保留 Glide-BEmin 前 1%，再按六路共识排列其余分子。在 M4 broad allosteric-modulator 公共基准上，Cascade 的 AUC 为 0.7775、EF0.5% 为 20.40、EF1% 为 13.73。该任务评价广义变构配体检索，不评价 PAM 功能。

### 4.3 PACER-FS 功能效力建模

PACER-FS 将跨来源/系列的基线偏移与系列内 SAR 分开学习。在共同 11 系列三锚点协议中，macro Spearman 从 0.199 提高到 0.263，增益的系列 bootstrap 95% CI 为 `[+0.017,+0.113]`。它适合有少量功能锚点的新系列，不是零样本 PAM 判别器。

### 4.4 四上下文 PACER-FKG

四个体系为 apo（A）、ACh-only（P）、candidate-only（C）和 candidate+ACh（CP）：

```text
Delta_PAM = CP - P
Delta_AGO = C - A
Delta_INT = CP - P - C + A
```

在冻结表示、区域和统计规则后，主区域的跨 replica 方向一致性为：已知 PAM LY2119620 `+0.634`，compound110 `+0.584`，实验 inactive 高相似近邻 CM00734 `-0.327`；CM00734 在 9 个报告区域均为负。该结果是回顾性、预冻结、零结果驱动调参的初步功能特异性证据，但样本量不足以建立通用 PAM 分类器。

## 5. 合作者建议阅读顺序

1. `START_HERE_FOR_COLLABORATORS_V01.md`：本文件，10 分钟掌握全貌。
2. `PACER_M4_PROJECT_MAINLINE_DEVELOPMENT_AND_INFERENCE_V01.md`：开发与推理的完整主线。
3. `PACER_M4_INTEGRATED_ALGORITHM_AND_BENCHMARK_V02.md`：当前主要数字、基线和结论边界。
4. `PACER_DRUGCLIP_ROUTER_INTERFACE_V01.md`：结合筛选模型与部署路由。
5. `THOMPSON_MIAO_2026_EXACT_REPRODUCTION.md` 与 `PACER_XR_CASCADE_METHOD.md`：六路 docking 的来源、复现和级联创新。
6. `PACER_FACTORIAL_KERNEL_GRAPH_METHOD.md` 与 `PACER_FKG_IMPLEMENTATION_AUDIT_20260927.md`：四上下文算法定义和实现审计。
7. `PACER_DC_CLOSE_LOOP_20NS_FINAL_REPORT_v01.md`：已知 PAM 与 inactive 困难负样本闭环。
8. `PACER_M4_LEAD_PAM_HYPOTHESES.md`：候选证据包和湿实验判据。
9. `DOCUMENT_AND_EVIDENCE_INDEX_V01.md`：需要追溯全部方法、支线和原始结果时再查阅。

## 6. 哪些旧文档不能单独代表当前项目

- `project/PACER_M4_FINAL_REPORT.md` 冻结于 2026-08-30，记录 PACER-FS、早期生成和静态结构阶段，缺少后续 DrugCLIP 路由与 PACER-FKG 闭环。
- `project/PIPELINE.md` 中“唯一入口为 OneProt-MD”的表述已经过时。OneProt-MD 后来未成为最终主表示。
- `PACER_M4_RESEARCH_CHARTER.md` 是早期研究约束与理想架构，不等于所有模块均已实现或通过验证。
- 名称含 `PREREGISTRATION`、`PLAN`、`CONTRACT`、`RUNBOOK`、`HANDOFF` 的文件主要记录预设协议或执行方式，不能当作完成结果。
- OneProt、Geom2Vec、CGDA、State-MoE 等失败或关闭支线具有消融价值，但不应放在当前主线首页。

## 7. 对外结论的边界

当前最稳妥的整体表述是：PACER-M4 已建立从候选产生、结合检索、结构门控到四上下文动态功能复核的分层计算框架，并在公开结合基准和少量已知功能分子上获得可复现的模块级证据。项目尚未完成新候选的湿实验确认，也没有一个统一外部数据集能够为整个端到端流程给出单一“SOTA AUC”。最终候选必须以分层证据和不确定性报告，而不是包装成已验证 PAM。

