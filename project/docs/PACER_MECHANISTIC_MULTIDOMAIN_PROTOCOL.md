# PACER-M4 多通路机制模型开发协议

冻结日期：2026-08-30。

## 目标

开发并验证一个不把 M4 PAM 压缩为单一结合能或单一 EC50 的模型。算法把以下量分开：

- 变构配体亲和力 `pKB`；
- 直接激动性 `log τB`；
- ACh 条件下功能协同性 `log αβ`；
- 通路/检测域：GoB、cAMP、β-arrestin、pERK。

## 开发数据

- Jörg 2023：9 个 PAM，binding、cAMP、β-arrestin 多终点；
- Liu 2024：18 个 PAM，GoB、cAMP、β-arrestin 多终点；
- 两篇均来自 ChEMBL 标准化原始记录，按 assay description 显式映射终点；
- 跨文献完全重叠结构从查询指标中排除。

## 预先指定验证

1. `Jörg2023 -> Liu2024` 文献外推；
2. `Liu2024 -> Jörg2023` 文献外推；
3. 对共同通路 cAMP 与 β-arrestin，逐一评价 `logτB` 和 `logαβ`；
4. 对 `pKB` 做单终点文献外推；
5. 2015 LY2033298/pERK 只作已经被前轮结果启发后的 exploratory endpoint transfer，不再称 untouched external。

## 方法与基线

- Molecular-1NN：最相似训练分子的同终点值；
- Pooled-Ridge：忽略通路、合并同一参数；
- AssayConditional-Ridge：共享分子支路 + 通路特异残差支路；
- PACER-Triad：低秩联合预测 `pKB / consensus logτB / consensus logαβ`，并引入冻结历史 pEC50 作为外部功能先验；
- Historical-Absolute：冻结 430 分子 pEC50 模型，作为错误但重要的单终点迁移基线。

所有超参数只在训练文献内部用分子级交叉验证选择。查询文献的任何标签不得进入标准化、选择或拟合。

## 指标和成功条件

- 主指标：跨文献 pooled Spearman；次指标：MAE、逐方向 Spearman、worst direction；
- 同时报告 paired molecule bootstrap 的增益 CI；
- 只有主方法高于全部基线、增益 CI 下界大于 0、worst direction 不恶化，才支持方法学优越性；
- 小样本正结果只称跨文献机制迁移证据，不称通用 SOTA。

## 生物学边界

模型预测的是特定正构探针 ACh 和特定信号读出的 operational-model 参数。它不能脱离探针、受体状态和信号通路被解释为分子固有常数，也不能替代功能实验确认 PAM。
