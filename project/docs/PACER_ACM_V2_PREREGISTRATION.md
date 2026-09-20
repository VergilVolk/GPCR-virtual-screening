# PACER-ACM v2 预注册协议

> 冻结日期：2026-08-30  
> 冻结时点：第四套外部 M4 PAM campaign 的数值标签提取之前。

## 1. 研究假设

跨论文的 M4 PAM 功能数据同时包含两类信息：同一药化系列内部的相对 SAR，以及由正构探针、细胞系统、信号通路和 operational-model 拟合造成的 assay 坐标偏移。PACER-ACM v2 检验：在不查看查询标签的条件下，三个同系列、同实验上下文的功能锚点能否校准冻结基础模型，使其更可靠地预测其余分子的 `log(αβ)` 和 `log(τB)`。

该假设不等价于“用三个点训练一个新 QSAR”。锚点只估计每个 endpoint 的常数残差，主要排序必须来自冻结基础模型。

## 2. 冻结算法

1. 以 Morgan/ECFP4（radius=2，2048 bits）表示分子。
2. 基础预测器必须在读取目标 campaign 数值标签前冻结。
3. 从目标 campaign 的可测分子中，以 ECFP4 medoid 起始，再用 max-min 选出三个结构多样锚点。选择过程不得读取活性值。
4. 对 endpoint `e`，计算三个锚点残差：`r_i,e = y_i,e - f_base(x_i,e)`。
5. 主预测为：`f_ACM(x,e) = f_base(x,e) + mean_i(r_i,e)`。
6. 主分析禁用可改变排序的相似度加权残差。kernel residual 只作为 baseline，不能在看到结果后替代主方法。
7. 不输出 `pKB`。已有开发数据表明，三锚点校准没有为功能亲和力提供稳定外推证据。

## 3. 拒绝规则

任一条件不满足即返回 `REFUSE`，而不是强行输出高置信度预测：

- 恰有三个通过 QC 的锚点；
- 锚点与查询属于同一参考 chemotype/series；
- species、orthosteric probe、probe regime、assay family、readout、cell system 一致；
- 三个锚点的 operational-model 拟合合格，且 `log(αβ)`/`log(τB)` 标签完整；
- 查询分子到至少一个锚点的 ECFP4 Tanimoto 不低于 0.30；
- 三个锚点残差的样本标准差不高于 0.75 log unit；
- 基础预测为有限数值。

拒绝率与适用域覆盖率必须随最终结果一起报告。

## 4. 不确定性

区间由三个锚点的 leave-one-anchor-out 残差构建，并设置 0.30 log unit 的最小半宽，防止小样本产生虚假窄区间。由于只有三个锚点，这个区间是工程性不确定性界限，不解释为严格总体置信区间。

## 5. 外部评估合同

第四套 campaign 的来源与纳排标准先登记，再抽取结构与标签。所有分子按 paper/patent series 和 matched-pair family 组织。主要指标为 Spearman、MAE、worst-series Spearman、适用域 coverage；同时报告区间覆盖率和区间宽度。

强 baseline：冻结基础模型、ECFP4-1NN、ECFP4-Ridge、三锚点常数偏移、三锚点相似度加权残差。所有增益用 series/family 为单位做 5000 次 paired bootstrap。

只有同时满足以下条件，才把 PACER-ACM v2 称为当前冻结 benchmark 上的领先方法：

- 对每个预注册 baseline 的主要指标增益，其 95% CI 下界均大于 0；
- 最差系列不恶化；
- 拒绝覆盖率完整披露；
- 没有在外测后改锚点、阈值或主指标。

若不能满足，只报告负结果或限定范围内的优势；不得声称 SOTA。

## 6. 生物学解释边界

`log(αβ)` 表示功能协同性，`log(τB)` 表示 PAM 自身内在激动/效能成分；二者不能被单一 docking 分数代替。PACER-ACM v2 输出的是同一 assay 坐标系中的功能排序假设，不是已验证 PAM，也不能替代 BRET/cAMP/β-arrestin 或 operational-model 湿实验。

