# M4 PAM 变构耦合假设与可证伪检验

## 为什么不能把 docking affinity 当 PAM 效力

M4 PAM 的功能输出由至少四部分共同决定：PAM 结合、正构探针结合、受体构象动力学、受体—G 蛋白复合物稳定性。7TRQ 只提供 `VU0467154 + iperoxo + M4 + Gi1` 的一个活化态快照；它能约束合理结合模式，但不能单独证明功能协同性。

主要结构依据：Vuckovic 等同时量化 LY2033298/VU0467154 与 ACh/iperoxo 的 affinity、efficacy、cooperativity 和 probe dependence，并结合四套 cryo-EM 结构与 GaMD，指出最终药理由正构—别构—受体—转导蛋白复合物的动态稳定性决定（eLife 2023, DOI: 10.7554/eLife.83477；PDB 7TRQ）。

## 三个预注册结构模块

1. **ECV anchor 模块**：Y89、Q184、F186、Y439。论文 GaMD 直接跟踪 PAM 与这些残基的距离；假设它们主要反映别构配体在 extracellular vestibule 的锚定稳定性。
2. **W435 gate 模块**：W435。论文比较不同配体状态下 W435 chi2，并将其作为别构口袋构象动力学的重要读出；假设它比单纯 TYR439 接触更接近功能耦合门控。
3. **species/probe context 模块**：D432、T433。论文的人/鼠差异和突变实验表明这些位点参与 VU0467154 的物种差异及复合物动力学；假设其作用依赖具体 chemotype/probe，不能形成普适单调分数。

## 当前可计算检验

对 430 个具备 pEC50 的分子，在统一 7TRQ 受体中生成无标签 pose，并提取：

- Vina affinity（仅作基线）；
- 共晶接触指纹 Jaccard、质心偏移、口袋覆盖；
- 13 个共晶口袋残基的最短距离和原子对接触数；
- 上述三个文献模块的消融特征。

在 12 个完整药化系列 LOSO 中比较 `2D RF`、`Vina`、`IFP`、`TYR439-only`、`literature coupling`、`full structure`、`2D + structure`。主要指标只采用系列内 Spearman 的 macro 和 worst-group。

## 可证伪判据

- 若结构融合不能同时改善 macro 和 worst-series Spearman，则不宣称“耦合图提高 PAM 效力预测”；
- 若 W435/D432/T433 的效应方向在系列间反转，则结论应为 chemotype-dependent coupling，而不是全局药效团；
- 若 Vina 优于结构网络，只能说明静态结合排序有信息，不能等同于 PAM 功能机制；
- 任何结构模块必须在冻结测试系列之外定义，禁止用测试标签反向选择 hub。

## 生物学结论的最高可达等级

无湿实验时，最终只能提出“被跨系列回顾性数据支持、可实验检验的 M4 PAM 耦合假设”。真正的 PAM 身份仍需至少在 ACh 条件下测定浓度—响应、alpha/beta 或 alpha-beta，并排除显著内在激动活性；若研究 probe dependence，还需更换正构激动剂重复测试。
