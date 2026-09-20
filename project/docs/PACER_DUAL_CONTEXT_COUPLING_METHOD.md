# PACER-DC：Dual-Context Coupling Metric 方法

## 1. 科学问题

M4 PAM不是“在变构口袋结合得更强”的分子。真正需要区分的是：

1. PAM：在正构探针存在时增强正构-变构-受体耦合；
2. allosteric agonist / ago-PAM：没有正构探针也能推动激活；
3. neutral allosteric ligand：结合变构口袋但不产生正协同性；
4. inactive/decoy：既不稳定正构探针，也不形成目标耦合状态。

静态结构检验已经证明单一“PAM-like geometry”会把compound-110 allosteric agonist
排在真实PAM之前，因此算法必须显式建模正构探针依赖。

## 2. PACER-DC表示

### 2.1 Probe-matched coupling coordinate

使用相同正构探针背景下的成对差分定义结构轴：

- `7TRQ(VU0467154+iperoxo) - 7TRK(iperoxo)`；
- `7TRP(LY2033298+iperoxo) - 7TRK(iperoxo)`。

只保留两种PAM同向的内部距离变化。该轴在独立
`7V68(LY2119620+iperoxo) - 7V69(iperoxo)`中复现约+0.54增量，但仍不能单独
排除allosteric agonist。

### 2.2 Functional triplet

Triplet作用于snapshot/trajectory embedding，而不是直接把所有PAM帧标为阳性：

- Anchor：某PAM在`+probe`上下文的一条trajectory bag；
- Positive：不同replica或不同chemotype、相近cooperativity/orthosteric-stabilization
  表型的bag；
- Hard negative：几何相近但属于probe-only、allosteric-agonist-without-probe、弱PAM
  或neutral allosteric ligand的bag。

现有anchor-guided snapshot triplet可跨replica恢复连续结构轴，但hard cluster种子与
时间稳定性不足；因此不再把KMeans cluster id作为模型输出。

### 2.3 Sticky soft prototypes

连续轴被投影到有序soft prototypes；驻留概率由训练轨迹自相关时间决定，通过
forward-backward得到posterior occupancy。prototype只用于表示亚稳态占比、转移与
attention，不预先命名为“PAM状态”。

当前单MK-97体系中posterior连续性0.928、ACh稳定性6/6折单调，但effective
prototype count=2.979，未过3.0冻结门槛，故仍是开发模块。

## 3. 双上下文候选评价

每个候选至少需要两个匹配模拟/实验上下文：

- `candidate + orthosteric probe`；
- `candidate without orthosteric probe`；
- 并各自配套probe-only与apo receptor controls。

输出保持为Pareto向量，不压成未经校准的总分：

1. `CoupledShift`：+probe时coupling coordinate/prototype occupancy相对probe-only的变化；
2. `OrthostericStabilization`：正构配体pose RMSD、contact persistence或residence proxy改善；
3. `IntrinsicActivationRisk`：无probe时相对apo的激活/耦合状态位移；
4. `BindingCompatibility`：姿势、口袋接触、应变与基础affinity门控；
5. `Uncertainty/AD`：chemotype、动态采样与assay域外程度。

PAM假设要求`CoupledShift`与`OrthostericStabilization`高、`IntrinsicActivationRisk`低；
allosteric agonist则在无probe上下文也产生高位移。未获得多上下文数据时禁止输出
“PAM概率”。

## 4. 现有证据

- DeepRLI/静态binding embedding不能预测M4 PAM potency；
- 普通ensemble docking在27 PAM/27 inactive上AUC约0.50；
- 第三PAM静态结构支持跨PAM共同结构增量，但compound-110证明其不具PAM特异性；
- Wang 2022三条1 us公开MD中，LY2119620使iperoxo RMSD平均降低0.651 A，
  replica级精确单侧p=0.05，50 ns block-bootstrap 95% CI 0.527--0.762 A；
- MK-97六副本中，蛋白coupling coordinate与ACh稳定性6/6同向，mean rho=0.196，
  循环移位p<0.0001；
- hard cluster和当前soft prototype尚不能作为效力模型。

## 5. 达到可声称算法性能所需验证

1. 至少3种PAM chemotype、probe-only、apo、allosteric agonist、weak/neutral ligand
   的匹配多replica轨迹；
2. trajectory-bag级leave-chemotype/system-out，不允许frame split；
3. 终点分别使用`log alpha`、`log alpha beta`、`log tauB`、intrinsic agonism，不用
   单一EC50替代；
4. baseline包括PCA/tICA/HMM、无tripletStateMIL、静态QSAR、单结构/ensemble docking；
5. 最终仍需前瞻功能实验确认。没有湿实验前，候选只能称计算假设。
