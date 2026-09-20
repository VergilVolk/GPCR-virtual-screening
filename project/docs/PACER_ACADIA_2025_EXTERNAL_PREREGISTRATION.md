# WO2025122811 外部功能效力验证预注册

> 登记日期：2026-08-30  
> 状态：已查看专利元数据与 assay 描述，尚未人工读取 Table 2 数值标签。

## 1. 外部来源

- 来源：WO2025122811A1，`M4 positive allosteric modulators`；
- 申请人：Acadia Pharmaceuticals Inc；
- 优先权日：2023-12-07；公开日：2025-06-12；
- 预期实例：Table 1 的 Examples 1–49，另有少量立体异构体/参考化合物；
- 数据源：WIPO 专利 PDF 与 Google Patents 镜像。

该来源未被用于 PACER-M4 当前 430 分子训练、方法选择或超参数调节。完成 PDF 哈希与结构恢复后，必须做 canonical-SMILES、InChIKey、ECFP4 近邻和专利族重叠审计。

## 2. 实验端点

人 M4 PAM assay 使用 CHO-K1/M4/Gα15 稳定细胞、ACh EC20（M4 final 15 nM）、FDSS μCELL 钙流。化合物最高 10 µM，四参数曲线拟合；报告 relative efficacy（RE），仅当 RE≥50% 时报告 EC50，否则为 ND。

主要任务预注册为：

1. `PAM_RE >= 50%` 的二分类/删失排序；
2. 对有 EC50 的分子预测 `pEC50`；
3. 同时报告 M4 agonist RE/EC50，作为内在激动倾向的负向约束；
4. 若 M2 数据完整，报告 M4/M2 功能选择性。

## 3. 与 PACER-ACM v2 的边界

本专利没有公开 operational-model 的 `log(αβ)` 与 `log(τB)`，所以不能作为 PACER-ACM v2 机制分解的外部验证。它是 PACER-M4 **功能效力与效能门控** 的第四套外部 campaign。禁止把 PAM EC50 当作 `log(αβ)`，也禁止从单一 EC20 曲线反演 affinity、cooperativity 或 intrinsic efficacy。

## 4. 冻结切分与 baseline

- Table 1/2 中所有可唯一恢复分子构成一次性外测；
- 三个锚点按 ECFP4 medoid + max-min、完全不看标签选取；
- 其余分子为查询，锚点只用于 assay-level 校准；
- 结构/标签无法唯一对齐的实例进入隔离表，不参与主指标；
- censored/ND 不按任意固定 pEC50 冒充精确标签。

Baseline：ECFP4-1NN、ECFP4-Ridge、RF/ExtraTrees、冻结 PACER ligand-only、单构象 docking、ensemble docking、三锚点常数偏移。主模型不得在读取 Table 2 后改变表示、阈值或锚点。

## 5. 评估

- 精确 pEC50：Spearman、MAE、pairwise accuracy；
- RE/ND：PR-AUC、ROC-AUC、balanced accuracy；
- 早期筛选：EF@10%、top-k hit rate；
- 每项报告 95% CI、覆盖率、结构恢复率和拒绝率；
- scaffold/series family 为 bootstrap 单位；
- 任何 AUC=1 必须同时报告正负样本数、置信区间和逐分子预测，禁止作为小样本成功结论。

## 6. 成功与失败口径

只有相对所有预注册强 baseline 的 paired CI 下界大于 0、worst-family 不恶化、且没有数据泄漏，才称为该冻结数据上的领先结果。无湿实验时，候选仍只是可检验的 PAM 假设；本外测不能确认新分子有效。

