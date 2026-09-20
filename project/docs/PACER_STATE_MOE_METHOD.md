# PACER-StateMoE：探针状态条件化的 M4 PAM 效力模型

> 状态：**已完成并被证伪，不进入最终效力模型。** StateResidualMoE v0.1 仅由 0.199 提至 0.203，delta 95% CI 跨零；后续嵌套四专家稳健融合由同轮 ligand-only 0.207 降至 0.151，delta 95% CI `[-0.115,-0.003]`。本文件保留为预注册设计和负结果记录。

## 1. 由数据逼出的科学问题

单一 7TRQ 静态结构对不同药化系列的效应方向相反；7TRP 又呈现不同的 macro/worst 权衡。M4 PAM 效力不能被“一个 docking score”或“一个共享药效团”充分描述。核心假设改写为：

> 分子 chemotype、正构探针状态与 extracellular vestibule 微状态共同决定可转移的功能 SAR；模型应学习何时信任哪一种结构状态，并允许退回 ligand-only null expert。

## 2. 模型结构

### 2.1 四个专家

1. `E0 ligand-only`：D-MPNN/ECFP-LightGBM，防止错误结构强迫进入预测；
2. `EQ 7TRQ`：VU0467154 + iperoxo + Gi1 状态；
3. `EP 7TRP`：LY2033298 + iperoxo + Gi1 状态；
4. `ES 7TRS`：ACh + Gi1 状态。

结构专家读取逐残基距离/接触、IFP、口袋覆盖及 state-to-state delta；不把 Vina affinity 当功能标签。

### 2.2 条件门控

门控器输入分子表示、assay family、orthosteric probe、cell system 和数据缺失 mask，输出四专家概率。主要约束：

- entropy floor：防止小数据下所有样本塌缩到单专家；
- null-expert prior：结构证据不可靠时允许回退 E0；
- source-balanced loss：每个训练来源等权，避免大系列支配；
- leave-one-source inner validation：状态权重和正则不得查看外部测试系列。

### 2.3 三类联合损失

`L = L_potency + lambda_pair L_pair + lambda_domain L_source + lambda_cal L_uncertainty`

- `L_potency`：masked robust pEC50 regression；
- `L_pair`：同实验上下文近邻 pair 的 delta/ranking，强调 activity cliffs；
- `L_source`：Group-balanced worst-risk 或 variance penalty，但仅在 inner CV 证明有效时保留；
- `L_uncertainty`：deep ensemble/conformal calibration 所需的不确定性头。

## 3. 与普通“多构象 docking + 加权分数”的区别

- 状态不是无条件平均，而由分子与实验上下文条件化；
- 包含 ligand-only null expert，可证伪结构是否有用；
- 训练目标是功能 pEC50 与局部 SAR，不是结合 affinity；
- 结构状态来自正构探针/PAM 共晶机制，而非随机 receptor snapshots；
- 外层整药化系列完全封存，杜绝事后选择最有利状态。

## 4. 必须击败的冻结 baseline

- ECFP RF：macro Spearman 0.170；
- ECFP LightGBM：0.199，worst -0.130；
- ChemBERTa Ridge：0.148；
- 单结构 7TRQ 融合：0.208，worst -0.152（不显著）；
- 单结构 7TRP 融合：0.148，worst -0.051；
- Chemprop D-MPNN：待完成；
- 无条件三状态拼接：待完成。

## 5. 晋级标准

1. 12 系列 macro Spearman 超过最强 baseline；
2. worst-series 不恶化；
3. series bootstrap 95% CI 的 delta 下界大于 0，或至少在更多独立系列上方向一致；
4. 结构/门控/pair loss 消融各自可解释；
5. 不依赖 pooled 指标、单系列或测试集选择状态；
6. 对生成候选输出 calibrated interval 和适用域，不输出“必然 PAM”。

## 6. 生物学可检验预测

如果模型成立，应观察到：

- 不同 chemotype 的门控权重偏向不同 ECV 微状态；
- W435 gate 与 D432/T433 context 的效应随 probe/series 改变，而非全局单调；
- activity cliffs 更常伴随状态权重切换或关键接触网络重排；
- ACh 数据应比 iperoxo 状态在 7TRS 专家上获得更稳定的跨系列贡献。

这些是计算产生的机制假设，最终需以 ACh 条件的浓度—响应、协同性参数和内在激动活性实验验证。
