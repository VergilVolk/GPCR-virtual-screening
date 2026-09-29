# CGDA 外部验证协议

更新日期：2026-09-26

## 目的

Full LIT-PCBA 的15个靶点不足以充分检验 CGDA 的跨靶点稳定性。本阶段使用与源训练数据完全分开的两个外部证据层：DUD-E 的外部 GPCR 子集，以及 MF-PCBA 的真实多阶段 HTS 数据。

## 源模型冻结

CGDA 源模型只使用官方 LIT-PCBA 15 靶点训练。训练配置在外部数据接触前冻结：两个 rank-2 专家、目标上下文门控、BCE、跨靶点检索、困难负样本 ranking、轻量表示保持；使用三个独立随机种子。

外部候选标签不得用于：

- 模型训练；
- epoch/超参数选择；
- 专家数量或rank选择；
- 融合权重选择。

## 外部 DUD-E GPCR 面板

初始面板为 ADRB1、AA2AR、DRD3、CXCR4 和 M3R（DUD-E 名称 `mcr`）。ADRB2 与 PPARG 因在 LIT-PCBA 源训练靶点中有相同目标而排除。

每个靶点提供：

1. DUD-E 官方 `mols.lmdb` 内的活性物/decoy；
2. 官方 `pocket.lmdb`；
3. 官方 `crystal_ligand.mol2` 转换的一个参考配体；
4. 原始 DrugCLIP、reference retrieval 和三个CGDA源模型均值的逐候选分数。

标签与embedding必须按LMDB数值键顺序对齐，并逐条核对SMILES。主要指标为 ROC-AUC、BEDROC80.5、EF0.5%、EF1%、EF2%、EF5%；统计使用目标层级配对bootstrap。

## DUD-E解释边界

DUD-E含102个蛋白、22,886个聚类活性物，且为每个活性物生成property-matched ZINC decoy。它适合大规模压力测试，但生成decoy并不等价于实验inactive。因此，DUD-E上的正结果只能支持跨靶点筛选一致性，不能单独支持真实HTS命中率或PAM功能。

## MF-PCBA验证

MF-PCBA包含60个由PubChem多阶段HTS组成的数据集，primary和confirmatory读出均可获得。Boltzina公开的8靶点测试集被作为独立实验HTS来源。当前下载的测试包只有SMILES、实验标签和Boltzina/GNINA/Vina等既有分数，未携带可复用的蛋白结构、口袋定义或共晶参考配体。因此它尚不能作为CGDA的外部结构评测输入；必须逐一建立 assay→靶标构象→口袋→参考配体的可追溯映射后才可启动。其意义是为后续“少量真实标签/多保真”功能模块提供独立评估，而不是拿现成Boltz分数冒充CGDA结果。

## 当前主张门槛

外部结论必须同时满足：

1. 在源模型冻结后得到；
2. 不使用外部标签选模；
3. 报告逐靶点而非只报全体合并；
4. 与原始DrugCLIP和reference retrieval比较；
5. 将DUD-E与实验HTS结论分开陈述。
