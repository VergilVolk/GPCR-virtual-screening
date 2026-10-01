# DrugCLIP 后续方法调研与决策

更新：2026-09-26。目的：为 PACER 的 Binding 模块确定可复现、可比较且不夸大的下一步，而非继续随机堆叠 score fusion。

## 已核实的相关方法

### 1. DrugCLIP：双塔检索基线

DrugCLIP 将虚拟筛选表述为 pocket-to-molecule dense retrieval，用 Uni-Mol 口袋/分子双编码器及对比学习共享空间；候选分子可离线编码。原论文强调 LIT-PCBA 的实验 active/inactive 比 DUD-E 更严格，并以 BEDROC/EF 作为筛选核心指标。它不是针对目标标签的 post-hoc score stacker。

**对 PACER 的含义**：任何新模块应比较并保留 pocket DrugCLIP，而不能只超过弱的 reference-ligand retrieval。

### 2. AANet（NeurIPS 2025）：结构不确定性下的对齐与聚合

AANet 是目前最直接相关的已发表方法，公开代码可用。其两阶段设计为：（i）tri-modal contrastive alignment；（ii）以 cross-attention 对多个候选口袋聚合，并以口袋置信度监督稳定该聚合。它在 DUD-E/LIT-PCBA 的 holo、apo 与 blind-predicted pocket 条件下比较 DrugCLIP，明确把多口袋/结构不确定性作为问题，而非给两个最终分数加权。

**对 PACER 的含义**：AANet 是我们的多构象/ensemble 结构模块必须复现的强基线。仅仅将 MD snapshot max/mean pooling、reference retrieval 或线性 blend 称为创新都不成立。

### 3. LigUnity（Patterns 2025）：粗粒度筛选 + 细粒度靶点内排序

LigUnity 在共享 pocket-ligand 空间上联合训练 active/inactive 区分与 pocket-specific ligand ranking，并加入 scaffold discrimination、pharmacophore ranking 和异构 pocket-ligand 图传播。它的启示不是“套一个 LoRA”，而是：VS 的类别区分与 hit-to-lead 的同靶点相对排序需要不同监督；若只有 active/inactive 标签，不能伪造 potency ranking。

**对 PACER 的含义**：可借鉴为“筛选损失 + 有真实同靶点剂量/效力数据时的独立排序损失”。M4 PAM 的效力/协同性数据不足时，后者必须不训练或只作为独立 pilot。

### 4. 近期 GenDrugCLIP / CausalBind 等

这些工作提示 pose-level、因果/概念分解或生成式监督可能改善筛选，但部分仍为预印本/在审，且需要大规模训练数据、严格去重和 GPU。它们不适合作为当前 CPU 条件下可直接宣称的“快速微调方案”。

## 已用本项目数据否决的方案

- reference retrieval 与 CGDA 的嵌套线性混合：早期富集低于纯 CGDA；
- pocket DrugCLIP 与 CGDA 的嵌套线性混合：BEDROC 低于纯 CGDA；
- 跨靶点监督的 score router：Full LIT-PCBA target-LOO BEDROC `0.0211`，显著低于 pocket 与 CGDA。

因此不能再使用“不同分数的后处理路由/融合”作为主算法路线。

## 建议的、按证据排序的后续路线

1. **先复现 AANet 官方 LIT-PCBA/DrugCLIP 对照**：固定官方脚本、checkpoint、12-target 子集和指标，确认环境、输入口袋和结果量级。它是对老师所说“docking 太古早、要处理多构象”的直接现代答案。
2. **再做 PACER 的最小差异扩展**：用 M4 的 MD/GaMD ensemble 替换 AANet 的候选口袋来源，保持其聚合器与训练目标不变；先比较 single crystal、mean/max ensemble、AANet aggregation。这个消融能回答“动态构象是否真的增益”。
3. **功能层完全分开**：Binding 模块只能输出 binding-ranking。PAM functionality 仍需 ACh 条件下的独立药理真值或四上下文动力学证据；不得将 AANet/DrugCLIP 分数说成 PAM efficacy。
4. **外部评估**：主结果以 LIT-PCBA target-LOO 和 matched DrugCLIP/AANet 为准；DUD-E 仅压力测试；MF-PCBA 必须完成 assay-to-structure 映射后才可作真实 HTS 外部证据。

## 当前决策

停止继续开发无文献依据的 DrugCLIP score router。下一项实际开发是 AANet 官方复现审计（代码、权重、数据、硬件门槛与官方 12-target LIT-PCBA 结果），再决定能否在 GPU 上执行 PACER ensemble 扩展。
