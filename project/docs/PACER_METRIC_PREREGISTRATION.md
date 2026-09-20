# PACER-MetaMetric v0.1 预注册

冻结日期：2026-08-30

## 假设

固定Tanimoto核无法知道“哪些局部结构差异对M4 PAM系列内SAR重要”。PACER-MetaMetric在外层训练系列上进行episodic训练，学习一个仅用于三个锚点残差传播的低维化学度量；绝对QSAR仍提供全局基线。

## 无泄漏协议

- 外层：11个`source_component` LOSO系列；
- 每个外层训练分子的absolute base来自不含其自身系列的inner-LOSO模型；
- meta训练每轮在各训练系列内采3个support，其余为query；每5轮使用predicted-span support，其余为固定种子的随机support；
- 输入：外层训练集方差最高的256个ECFP4 bit、9个标准化描述符和无标签base prediction；
- 网络：265→128→32归一化embedding，学习温度；
- query residual为三个support residual的softmax距离加权；
- 损失：query MSE + 0.2系列内pairwise ranking loss；固定220轮，不看外层结果早停；
- 外层测试固定3个predicted-span anchors。

## 主比较

主终点为equal-series macro Spearman；比较Absolute、固定Tanimoto AnchorResidual、冻结PACER-FS 0.263。只有相对Absolute的series-bootstrap 95% CI下界大于0且不低于0.263才晋级。失败结果完整保留。

## 冻结扩展：PACER-Hybrid

MetaMetric单模型首次运行得到macro Spearman 0.254，未达到晋级线。结果查看后只允许一个互补性检验：用同一3个support、同一absolute base，将PACER-FS DeltaSARHybrid与MetaMetric的pEC50预测作固定50:50算术平均。不得搜索权重、按系列选专家或改变网络；主比较仍使用11系列配对bootstrap。若固定融合不超过0.263或CI下界不大于0，则作为负结果结束该路线。
