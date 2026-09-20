# M4R 大规模早期富集路线审计

> 冻结日期：2026-08-29  
> 决策：**Miao/DUD-E active–decoy 数据可作为 docking 失败对照，不能作为 PACER-M4 生物学 SOTA 的主要监督 benchmark。**

## 1. 原始假设

利用约 2300 个 M4R active 和 11.5 万个 decoy，开发直接优化 EF(0.5%)/BEDROC 的多状态排序器，再以 430 个功能效力分子精排。该路线的吸引力是样本量大、Vina 早期富集弱、表面提升空间充足。

## 2. 防泄漏 benchmark

- 117,568 条 PDB docking 记录；
- canonical 去重并排除 8 个活性/decoy 标签冲突后，117,552 个分子：2306 active、115,246 decoy；
- 57,854 个 Murcko scaffold；833 个 active-bearing scaffold；
- 活性和 decoy 的整个 Murcko scaffold 均只进入一个折；
- 五折各约 23,510 个分子、461 个 active，骨架零跨折；
- 主指标为 EF(0.5%)、EF(1%)、BEDROC20/80，ROC-AUC 只作辅助。

## 3. 原始 decoy benchmark

| 方法 | ROC-AUC | BEDROC80 | EF(0.5%) | worst-fold EF(0.5%) |
|---|---:|---:|---:|---:|
| PDB Vina | 0.725 | 0.047 | 0.864 | 0.432 |
| Ensemble BEmin | 0.734 | 0.032 | 0.346 | 0.000 |
| Ensemble BEavg | 0.696 | 0.044 | 0.950 | 0.000 |
| 9 descriptors LightGBM | 0.980 | 0.787 | 48.472 | 46.575 |
| ECFP classifier | 0.9997 | 0.998 | 50.890 | 50.566 |
| ECFP LambdaRank | 0.906 | 0.693 | 50.977 | 50.887 |

该数据比例下 EF(0.5%) 理论上限约 51。仅 9 个理化描述符已接近上限，说明模型主要识别 active/decoy 生成来源，而不是 M4 生物学。LambdaRank 可以把每折 Top-100 全部排成 active，但这不能被解释为真实虚拟筛选 SOTA。

三随机种子的训练标签置乱对照得到 macro ROC-AUC `0.468 ± 0.076`、EF(0.5%) `2.42 ± 1.57`。AUC 总体接近随机但方差很大，极早期 EF 甚至可因折间分布漂移和极小 top-k 偶然高于 1；因此任何单种子、小 top-k 的 EF 都必须和置乱分布比较。真实标签模型 EF≈51 远高于置乱，但该差异仍主要反映 decoy-origin shortcut，而非外部功能泛化。

## 4. hard-decoy 与 adversarial-decoy

我们依次构造：

1. 每个 active 在同一 scaffold fold 内匹配 10 个唯一理化近邻 decoy；
2. 用其他四折训练的描述符模型挑选“最像 active”的 decoy，再为每个 active 分配同折 ECFP4 最近的 10 个唯一 decoy。

即使在第二种对抗集合中，active–decoy 平均最近 ECFP4 Tanimoto 仍只有 0.211，说明原始 decoy 库缺乏真正的近邻实验阴性。ECFP 分类器 EF(0.5%) 仍达到 11.0（该 1:10 数据比例下的理论上限），描述符分类器也达到 11.0。通过重采样无法消除 decoy-origin shortcut。

## 5. 外部实验阴性检验

为测试这种“完美模型”能否转移到真实功能标签，我们构建 positive–unlabeled 训练：

- 1813 个 M4R active；
- 32,179 个类药未标注分子；
- 训练前删除全部外部测试分子及其 Murcko scaffold；
- 完全封存 463 PAM / 66 实验非活性；
- 只在三个同时含 PAM 与非活性的独立药化系列上计算 macro/worst AUC。

| 方法 | pooled AUC | macro series AUC | worst series AUC |
|---|---:|---:|---:|
| Positive-neighborhood similarity | 0.496 | 0.506 | 0.177 |
| Descriptor PU | 0.649 | 0.528 | 0.400 |
| ECFP PU bagging | 0.727 | 0.530 | 0.264 |
| DUD-E supervised | 0.730 | 0.527 | 0.329 |

pooled AUC 再次达到约 0.73，但系列内 macro 只有约 0.53；不同系列的标签比例和预测尺度制造了 Simpson 型假提升。大规模 decoy 监督没有转化为真实实验阴性识别能力。

## 6. 方法学结论

1. 数据量大不等于 ground truth 强。十万人工 decoy 不能替代几十个实验非活性；
2. scaffold split 只能阻止骨架记忆，不能消除 active/decoy 数据来源偏差；
3. 直接优化 EF 可以把 benchmark 做到近乎完美，同时完全不改善真实系列内 PAM 判别；
4. Miao 数据仍可严谨支持“Vina 全局 AUC 与早期富集脱钩”，但不能训练或证明功能 PAM SOTA；
5. PACER-M4 不把该路线晋级为主算法，也不把 EF≈51 写进摘要。

## 7. 后续唯一可信的提升路径

真正的方法学跃迁需要新增**实验条件明确的 hard negatives 和近邻 SAR**：至少跨多个 chemotype 的 PAM/inactive 对、ACh 浓度–响应、候选单独激动活性及最好具有第二正构探针。没有这些数据时，计算模型应输出候选与不确定性，不能承诺“真正有效的 PAM”。

在现有数据上，可继续研究的合理任务是：系列内 lead optimization、少样本新系列适配、activity-cliff 预测和前瞻候选选择；不能再以 DUD-E 式随机 decoy AUC/EF 冒充生物学 SOTA。
