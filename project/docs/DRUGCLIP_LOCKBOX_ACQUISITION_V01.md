# DrugCLIP 外部锁箱获取方案 v01

## 首选来源

2025 年发表的 M4 PAM biphenyl 系列：

- Luo Y. 等，*European Journal of Medicinal Chemistry* 298 (2025) 117993；
- DOI：<https://doi.org/10.1016/j.ejmech.2025.117993>；
- 文中明确报告 A1、A9 具有 M4 PAM 活性，A9 EC50 为 513 nM；
- 该系列时间较新、化学骨架明确，适合作为 M4 定向外部压力测试的候选来源。

但目前公开摘要只确认少数阳性分子，不能据此计算可靠 AUC。必须通过学校订阅、文章补充材料或联系作者取得完整 SAR 表，包括低活性/无活性同系列分子。

## 获取后如何保持真正盲测

由一名不参与模型选择的队员完成：

1. 从正文和 SI 提取结构、实验条件、EC50/Emax/fold-shift 与 inactive 标记；
2. 保存两份文件：`lockbox_inputs.csv` 不含标签，`lockbox_labels.csv` 只含 `pair_id` 与标签；
3. 将 `lockbox_inputs.csv` 交给建模侧，标签文件由该队员单独保管；
4. 建模侧固定使用本次冻结的 M4-safe 路由、聚合方式和阈值，输出预测并记录 SHA256；
5. 预测冻结后再揭盲，一次性计算结果；
6. 无论成功或失败均归档，不再针对该集合调参。

## 纳入门槛

- 至少包含两个类别，不能只有阳性 PAM；
- 阳性与阴性使用可比的实验体系和 ACh 条件；
- 结构可唯一解析，盐型和重复物完成标准化；
- 与既有训练/开发数据做 exact-SMILES、InChIKey、Murcko scaffold 三层重叠审计；
- 若大部分分子与训练集重复，只能作为回顾性验证，不能称锁箱。

## 备用数据

- PubChem AID 588743 含 M4 PAM confirmatory series，但年代较早，极可能已进入公开训练语料，只适合补充回顾性功能基准；
- PubChem AID 2216033 目前只有 1 个测试化合物，不能单独支持分类跑分；
- 已使用过的 Acadia、Monash、DUD-E 和现有 4-GPCR 集均不得再次命名为盲测锁箱。

## 锁箱的两个终点

1. **DrugCLIP 层**：结合/活性检索，主要看 BEDROC20、EF1%、EF5%。
2. **PACER-FKG 层**：在确认能结合的同系列分子中，区分 functional PAM 与 inactive binder；该终点需要功能标签，且不能由 DrugCLIP 分数替代。
