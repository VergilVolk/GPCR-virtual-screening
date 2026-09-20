# PACER-DC v1 开发合同

冻结日期：2026-09-05

## 1. 可检验假设

在已经通过别构口袋结合门控的 M4 配体中，匹配的有/无 ACh 轨迹差分，比静态结合分数或单一复合物轨迹更能区分功能性 PAM 与内在别构激动剂。

这是一项待验证的方法学假设，不是现有结论，也不是 SOTA 声明。

## 2. 最小数据单元

每个分子需要：

- `candidate_probe`：M4 + 候选物 + ACh；
- `candidate_no_probe`：M4 + 候选物；
- 匹配的 `probe_only`：M4 + ACh；
- 匹配的 `apo`：M4；
- 每种条件至少三个配对 replica。

独立统计单位是 replica 或独立体系，不是轨迹帧。训练/测试按分子、chemotype 或独立研究体系切分。

## 3. 模型层次

### B0：非动态基线

- 2D ECFP + PACER-FS/LightGBM；
- 单结构 docking；
- ensemble docking 的最佳分和平均分。

### B1：简单动态基线

- ACh RMSD 和原生接触占有率；
- 别构配体口袋接触及漂移；
- 冻结内部距离耦合坐标；
- 四上下文的 replica 级均值差分。

### B2：无监督轨迹表征

- PCA/tICA；
- VAMP 或 SPIB；
- 所有超参数只在训练体系中选择。

### M1：PACER-DC

对轨迹表示计算两项对照差分：

`Delta_PAM = E(candidate + ACh) - E(ACh-only)`

`Delta_AGO = E(candidate without ACh) - E(apo)`

输出协同变化、ACh 稳定作用、内在激动风险、结合兼容性和不确定性，不输出未经校准的 PAM 概率。

### M2：药理 Triplet 微调

只有达到数据闸门后才启用。Anchor/positive 必须来自不同 replica，外部验证时优先来自不同 chemotype；hard negative 必须包含别构激动剂和无功能结合物。Triplet 是 M1 的消融增量，不是默认成功模块。

## 4. 数据闸门

进入 Triplet 开发至少需要：

- 三种或以上 PAM chemotype；
- 一种或以上明确内在别构激动剂；
- 一种或以上可结合但无 PAM 功能的对照；
- 上述类别具有可比较的四上下文、多 replica 轨迹；
- 标签尽量包含协同性、效力和内在激动，而不只是一列 EC50。

未满足时，只开发数据流程和无监督/简单差分基线。

## 5. 主验证与判废

主验证采用 leave-system-out 或 leave-chemotype-out。报告逐体系效应、replica 一致性、置信区间和置换检验；样本量不足时不报告容易误导的 AUC。

下列任一情况发生，核心假设不成立：

1. M1 不能超过 B1；
2. M2 不能在外部体系稳定超过 M1；
3. PAM 和别构激动剂的差异只在 `candidate+ACh` 中出现，而无 ACh 风险轴不能区分；
4. 结果依赖 frame split、单个 chemotype 或单个 replica；
5. 对适用域外或缺失上下文分子仍强制排序。

## 6. 近期交付顺序

1. 自动审计四上下文清单、软件环境和证据覆盖；
2. 补齐可重复的配体参数化环境；
3. 先运行已知 PAM、别构激动剂及共享对照的小规模短轨迹 QC；
4. QC 通过后运行冻结的三个 replica 试验；
5. 完成 B1，再进入 B2；
6. 数据闸门通过后才启动 M2；
7. 所有候选结论最后由功能实验验证。

## 7. 已冻结的困难对照扩展

v2 清单增加同一实验来源的近邻分子对：CM00717（PAM）与 CM00734（实验非活性）。二者 ECFP4 Tanimoto 约为 0.955，适合检验模型是否只记忆 chemotype。CM00734 在通过结合/MD 门控前只称为“实验非活性困难负样本”，不得称为已确认的中性别构结合物。
