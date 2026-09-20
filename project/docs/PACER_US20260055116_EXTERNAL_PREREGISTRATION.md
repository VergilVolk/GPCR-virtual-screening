# US20260055116A1 多通路外部压力测试冻结协议

冻结日期：2026-08-30。

## 透明性说明

PACER-FunctionSpace v3 的算法和拒绝规则已先冻结。随后在来源资格筛查时，Google Patents 页面暴露了 Table 23 的 A/B/C/D 活性分档，因此本数据不是完全盲的最终验证；从本文件冻结起，不允许根据这些标签修改模型、特征、锚点策略或主指标。该数据只称为“冻结模型后的外部多通路压力测试”，不能单独建立 SOTA。

## 来源与 endpoint

- 来源：US20260055116A1，substituted tetrahydropyrrolo-pyridinone compounds；
- primary：human M4 PAM pERK EC50 四档；
- secondary：rat M4 PAM pERK EC50、human M4 PAM GTPγS EC50；
- A：`<100 nM`；B：`100-500 nM`；C：`>500-2000 nM`；D：`>2000 nM`；
- N/A 不插补，不当作 inactive。

## 冻结评估

1. 从专利结构图/名称恢复结构，RDKit 规范化，并按 active moiety 去重；
2. 与历史 430、Suven、Acadia、VU/Monash 数据做 exact-SMILES 与 ECFP4 重叠审计；
3. primary metric：不同真实活性档之间的 pairwise ordinal concordance，真实同档对忽略；
4. secondary：Spearman（A=3、B=2、C=1、D=0）、四档 balanced accuracy、A-vs-rest PR-AUC；
5. zero-shot baselines：constant、ECFP4-1NN、similarity-kNN、Ridge、ExtraTrees、frozen LightGBM/PACER centered；
6. seven-anchor：仅使用预先按 ECFP4 选择的 7 个功能锚点；比较 constant offset、Tanimoto GP、DeltaSAR 与 DeltaSARHybrid；
7. human pERK 模型不得借用 rat pERK 或 GTPγS 标签；跨 readout 只在主分析完成后作为机制审计。

对需要连续支持值的 seven-anchor 适配器，冻结区间代表值为 A=7.30、B=6.65、C=6.00、D=5.30 pEC50；主评价仍是只依赖顺序、不依赖这些代表值的 ordinal concordance。

## 晋级门槛

候选方法必须在唯一 active-moiety 集上相对每个强 baseline 的 paired bootstrap concordance 差值 95% CI 下界大于 0，并且 human pERK 的最差核心家族不恶化。若结构无法可靠恢复、精确历史重叠过高或只剩单一活性档，则拒绝主分析。

## 声明边界

分档 EC50 不能验证精确效力、cooperativity、intrinsic agonism 或 probe dependence；该专利也不能确认本项目生成分子为 PAM。
