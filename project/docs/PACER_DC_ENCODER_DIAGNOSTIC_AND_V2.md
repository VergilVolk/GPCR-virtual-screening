# PACER-DC 四上下文表示诊断与 v2 路线

更新日期：2026-09-25

## 一、当前结论

现有 MD、轨迹恢复、atom14 转换和冻结 checkpoint forward 已通过工程审计；当前失败点是**表示和差分的生物学稳定性**。compound110 的五个匹配窗口中，R2/R3 `dPAM` cosine 为 `0.503, 0.712, -0.354, 0.262, 0.085`，方向未跨 replica 稳定。继续增加训练样本或直接训练 triplet 会学习 replica/window 噪声，训练门应保持关闭。

## 二、为什么当前表示可能不合适

### 1. OneProt-MD 最终头过度全局化

本仓库 `TrajectoryEncoder.forward` 对 MDGen 输出执行：

1. 沿时间维求均值；
2. 沿全部残基求均值；
3. 将仅 21 维的 pooled vector 经 MLP 扩展到 1024 维；
4. L2 归一化。

因此，最终 1024 维向量的有效信息在投影前已经被压缩为 21 维全局平均。M4 PAM 可能只改变别构口袋、正构口袋、激活微开关和 G 蛋白界面的少数残基；这些局部变化会被约 270 个残基和 100 帧的双重平均稀释。MLP 扩维不能恢复已经丢失的信息。

OneProt 的多模态对比目标还倾向于让同一蛋白不同模态表示接近，适合学习蛋白身份和总体动态，但未必保留同一受体不同配体状态之间的微小差异。当前观察到原始 context embedding 跨 replica cosine 较高、差分方向却不稳定，与这一风险一致。

### 2. 归一化向量相减是病态操作

最终 embedding 均被归一化到单位球面。两个接近向量相减时，差分方向由很小的残差决定；replica 噪声会被相对放大。原始向量相似并不保证差分向量稳定。

### 3. 当前 dPAM 不是纯协同项

四上下文是一个 2×2 因子设计：候选物有/无 × ACh 有/无。严格的候选物-ACh交互项应为：

\[
d_{INT}=(z_{C+A}-z_A)-(z_C-z_0)
=z_{C+A}-z_A-z_C+z_0.
\]

当前 `dPAM=z_CA-z_A` 同时包含候选物一般作用和候选物-ACh交互；`dAGO=z_C-z_0` 虽单独计算，却没有从 `dPAM` 中扣除。v2 必须显式报告 `dINT`，并将 `dAGO` 作为独立风险轴。

### 4. 当前估计量方差过大

四个 context 是独立 MD 实现。差分的方差是多个轨迹方差之和；再做差分之差时更明显。单个1 ns窗口不能作为独立生物学样本，五个连续窗口也不能代替五个独立分子或replica。

## 三、立即执行的 encoder qualification（不新增 MD）

先使用现有 compound110 的40个窗口完成以下实验，任何模型训练之前必须交付结果。

### G1：保存三个表示层

对每条轨迹同时保存：

1. `H[t,residue,21]`：MDGen最后隐藏状态；
2. `h_global_preproj[21]`：时间/残基全局平均、投影前；
3. `z_oneprot[1024]`：当前归一化最终向量。

目的是判断不稳定来自MDGen本体、全局池化、MLP投影还是单位球归一化。

### G2：无标签敏感性测试

在同一现有窗口上生成确定性变换，不重新跑MD：

- identical repeat：必须逐位一致；
- rigid translation/rotation：应基本不变；
- frame reversal、固定种子frame permutation：判断编码器是否保留时间顺序；
- 前50帧重复、后50帧重复：判断表示是否只反映平均构象；
- 仅在已知口袋/微开关残基做小幅扰动：判断局部变化能否被最终向量观察到。

若frame permutation几乎不改变最终向量，或局部扰动远小于replica噪声，则最终OneProt向量不适合作为PACER-DC主表征。

### G3：方差分解

分别在 `H`区域池化、21维pre-projection和1024维最终向量上计算：

- context方差；
- replica方差；
- window方差；
- context×replica交互；
- 同context跨replica距离与不同context距离的分布。

主判据是 context 信号必须大于 replica/window 噪声，而不是只看PCA图是否分开。

### G4：共享基底的非神经基线

将所有训练轨迹共同拟合一个共享的PCA/tICA/VAMP基底，再冻结后投影每个context。禁止每个context独立拟合基底。比较其 `dINT`/`dAGO` 跨replica稳定性与OneProt。

如果简单的共享tICA/VAMP比OneProt稳定，说明瓶颈在预训练表示或池化，而不是MD完全没有信号。

## 四、推荐的 PACER-DC v2 表示

暂不抛弃MDGen；弃用其当前“全时空平均→21维→1024维”的最终OneProt头作为唯一主表征。

### 1. 区域感知池化

从 `H[t,residue,21]` 分别池化：

- 正构ACh口袋；
- M4别构PAM口袋；
- PIF/DRY/NPxxY等激活微开关；
- G蛋白耦联界面；
- 其余受体作为全局背景。

每个区域保留mean、standard deviation、early-late change和time-lagged statistic，避免时间平均抹除动态信息。

### 2. 四上下文交互表示

对每个区域 (R) 计算：

\[
d_{INT}^{R}=h_{C+A}^{R}-h_A^{R}-h_C^{R}+h_0^{R},
\qquad
d_{AGO}^{R}=h_C^{R}-h_0^{R}.
\]

最终输入为多个区域的 `dINT`、`dAGO` 及其不确定度，而非一个全局1024维差分。

### 3. 小型可训练头

冻结MDGen主干，只训练：

- 区域attention/weighted pooling；
- 低维projection head；
- 两个输出轴：ACh-dependent interaction与intrinsic agonism。

该设计保留预训练动态知识，同时把学习容量集中在M4药理相关区域，所需标签远少于从零训练轨迹模型。

## 五、通过/否决门槛

在进入任何功能分类训练前，至少满足：

1. 对预定义阳性/激动/阴性控制，区域级 `dINT` 或 `dAGO` 的跨replica方向中位cosine > 0.5；
2. 分子成簇或replica bootstrap的95%区间下界 > 0；
3. context方差/replica方差比 > 1；
4. 至少一个区域级或共享tICA/VAMP表示稳定优于最终OneProt全局向量；
5. window仅作为增强，统计单位仍是molecule×replica；
6. 训练、阈值和表示选择不能读取blind candidates。

若MDGen隐藏层、区域表示和共享动力学基线全部失败，再替换encoder；在完成这些诊断前直接更换另一个foundation model无法定位瓶颈。

## 六、队友下一步任务

1. 暂停新增候选和无目的延长轨迹；保留现有GPU用于LY2119620匹配replica。
2. 修改extractor，保存 `H[t,residue,21]` 和21维pre-projection向量，不改变现有1024维输出。
3. 对现有40个窗口完成G2敏感性面板。
4. 在三个表示层完成G3方差分解。
5. 用同一批轨迹建立共享PCA/tICA/VAMP基线。
6. 显式增加 `dINT=z_CA-z_A-z_C+z_0`，保留当前 `dPAM` 仅作历史对照。
7. 交付一张表：每种表示的跨replica cosine、context/replica方差比、局部扰动效应与frame-order敏感性；通过后才讨论微调。

