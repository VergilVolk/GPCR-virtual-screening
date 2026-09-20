# PACER-FunctionSpace v3 冻结协议

冻结日期：2026-08-30。该版本在 Acadia 分子级去重外测失败和锚点敏感性诊断之后开发，因此 Acadia 只能作为开发/压力测试集，不能再作为独立验证。

## 科学问题

固定三个“化学最分散”锚点可能覆盖结构空间，却不能保证覆盖功能响应空间。M4 PAM 的 potency、PAM efficacy、intrinsic agonism 与 M2 selectivity 也不是一个标量。算法必须能够拒绝信息不足的预测，并把额外实验预算用在最能扩展功能范围的分子上。

## 冻结算法

1. 在候选系列中预先划分 anchor pool 与完全隐藏的 query；
2. 先按 ECFP4 medoid + max-min 选择 3 个锚点；
3. 同时测量 M4 PAM EC50/RE、候选单独激动 EC50/RE、匹配条件 M2 反筛；
4. 若 M4 PAM pEC50 实测跨度 `<0.75 log`，不输出正式排名；
5. 用当前锚点拟合 ΔSAR，在未测 anchor pool 中选择预测最低与最高的两个分子；
6. 重复一次，最多 7 个锚点；若跨度仍不足则拒绝该系列的绝对功能预测；
7. 分终点给出预测及不确定性，不构造 affinity/docking/QED 人为总分；
8. 最终按高 PAM potency/RE、低 intrinsic agonism、高 M4/M2 selectivity 做 Pareto，并保留探索性分子。

`0.75 log` 相当于约 5.6 倍效力范围，是预先冻结的功能夹逼门槛，不根据后续测试集调参。

## 开发评估

- 历史 430 分子按来源系列外留；每个测试系列在任何标签读取前固定分成 anchor pool 与 query；
- query 永远不参与锚点选择、停止决策或模型选择；
- 比较固定 3、5、7 个多样性锚点、随机三锚点和自适应预算；
- 主指标为系列宏平均 Spearman、worst-series、MAE、平均锚点预算与拒绝率；
- Acadia 只允许做 post-hoc 开发压力测试；最终证据必须来自第五套未触碰 campaign。

## 晋级与声明边界

自适应策略只有在相同固定 query 上同时改善宏平均排序、不恶化 worst-series，并以少于固定 7 点的平均实验预算完成时才晋级。开发集成功不能称外部 SOTA；没有新功能实验不能称候选为真实 PAM。
