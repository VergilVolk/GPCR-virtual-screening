# PACER-FKG 实现审计与纠偏

日期：2026-09-27

## 结论

队友已正确完成 Geom2Vec GPU 提取、输入哈希核验、atom14 映射和 PBC 屏幕；相关 52 个单元测试通过。当前低复现性不是简单的软件报错，而是第一版 FKG 判据存在两个统计定义错误，且现有表征仍没有显式利用时间顺序。

因此：旧报告中的 `compound110_extension` 和 `pam_contact_union` “passing synergy regions”作废；原始结果保留用于追溯，不删除、不覆盖。

## 两个必须纠正的错误

### 1. 不同核空间的数值不能相减

第一版对每个区域分别拟合 median bandwidth；区域残基数也不同。区域核统计与 `stable_core_control` 不在同一特征空间和量尺上，直接计算“区域 U² − control U²”没有严格可比性。

### 2. 平方核距离没有方向

`kernel_u2`只表示分布差异大小。两个 replica 即使发生方向相反的变化，也可能同时得到较大的正数。因此，两个 replica 的标量幅度同时较高不能证明它们复现了同一个动态响应。

修正后的 G2-B 使用共同 RFF 坐标比较带符号的四上下文对比向量，才可以讨论跨 replica 方向一致性。

## 当前真正支持的结果

### 数据和工程链

- 40/40 Geom2Vec 输入完成，共4,000帧；
- checkpoint、序列、形状、有限性和哈希均通过；
- 八条原始DCD的PBC屏幕未发现受体断裂或明显wrap事件；
- 这些结果支持“编码与数据链可运行”，不支持PAM预测。

### 共同核方向审计

- synergy轴在PAM共识口袋方向不一致：`pam_contact_consensus=-0.135`，`pam_contact_union=-0.228`；
- orthosteric activation core为`0.467`，但stable comparator也为`0.238`，不能据此认定PAM特异性；
- conditional candidate effect在PAM接触区较稳定（consensus `0.514`，union `0.414`），说明候选物存在可复现的条件性构象影响，但不等于正向PAM作用。

### 图消融中最有生物学意义的线索

compound-110的intrinsic-agonism轴在胞内microswitch区域经冻结图传播后，跨replica cosine由`0.053`升至`0.388`，变化`+0.335`，描述性block区间`[0.167,0.417]`。该方向与compound-110具有变构激动活性的已知属性一致，是目前最值得保留的机制线索。

但只有一个化合物、两个replica，不能宣称分类性能或因果通路。

## 另外两个架构缺口

1. `stable_core_control`是“跨结构低位移”面板，不是经过药理学证明的中性残基集合；其中包含靠近TM3/TM5胞内激活区域的残基，不能作为唯一功能阴性对照。
2. 当前方法把每个上下文当作构象分布，打乱帧顺序后结果不变。它利用了MD采样，却没有建模转移、时间相关或慢动力学，严格说是ensemble模型，不是trajectory-dynamics模型。

## 立即执行的修正路线

1. 废止跨区域kernel magnitude减control的功能通过门禁；
2. 保留G2-B共同核带符号方向作为基础表示；
3. control改成多个按区域大小、图度数和结构方差匹配的随机面板，报告相对null分位数，不依赖单一stable control；
4. 新增time-lagged transition kernel：比较`(z_t,z_{t+lag})`的四上下文分布，而非只比较`z_t`；
5. 以state-only共同核、transition kernel、无图和置乱图为消融；
6. 在LY2119620纯PAM、compound-110变构激动剂和CM00734 inactive上冻结验证，之后才允许训练triplet/projector。

## 主张边界

目前只能说：Geom2Vec能稳定编码现有轨迹；compound-110在候选物条件效应与胞内激动传播上出现部分可复现线索。尚不能说PACER-FKG预测了PAM、达到SOTA或优于完整baseline。
