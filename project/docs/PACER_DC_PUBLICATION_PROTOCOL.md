# PACER-DC 可发表性验证协议 v1.0

冻结日期：2026-09-15

## 1. 研究问题与主张边界

主假设：对于已经通过 M4 别构口袋结合门控的分子，有、无 ACh 的匹配轨迹差分，比二维分子模型、静态对接、ensemble docking 和单上下文 MD 更有助于区分功能性 PAM、内在别构激动剂及实验非活性近邻。

本研究评价的是计算判别和候选优先级，不把未经功能实验确认的候选称为 PAM。PACER-DC 在完成预注册比较前不声明 SOTA。

## 2. 冻结输入单位

每个候选需要四个上下文：

1. `candidate_probe`：M4 + candidate + ACh；
2. `candidate_no_probe`：M4 + candidate；
3. `probe_only`：M4 + ACh；
4. `apo`：M4。

每个上下文至少三个独立 replica。独立统计单位是分子、体系或 replica，不是轨迹帧。共享对照必须与候选使用相同受体模板、力场、膜组成、平衡流程、生产长度和 paired seed group。

## 3. 药理角色与最小数据闸门

- 至少三个独立 chemotype 的功能性 PAM 阳性对照；
- 至少一个明确的内在别构激动剂或 ago-PAM 对照；
- 至少一个实验非活性近邻；
- 若声称区分 binder 与 functional PAM，还必须有至少一个实验确认结合但无 PAM 功能的分子；
- 计算候选只能在对照阈值冻结后分析。

未达到上述闸门时，只报告轨迹工程、机制端点和个体效应，不报告功能分类 SOTA。

## 4. 比较方法

所有方法使用同一外层切分和同一候选集合。

| 编号 | 方法 | 输入 | 作用 |
|---|---|---|---|
| B0-1 | ECFP + Ridge/LightGBM | 二维结构 | 化学结构基线 |
| B0-2 | 单结构 docking | 一个实验结构 | 静态结合基线 |
| B0-3 | ensemble BEmin/BEavg | 多受体构象 | 受体柔性结合基线 |
| B1-1 | 单上下文 MD | `candidate_probe` | 检验仅增加动态是否有效 |
| B1-2 | 四上下文手工差分 | 四上下文、可解释端点 | PACER-DC 最小动态基线 |
| B2-1 | PCA/tICA | 四上下文轨迹特征 | 线性/慢模态表示基线 |
| B2-2 | VAMP 或 SPIB | 四上下文轨迹特征 | 非线性动力学表示基线 |
| M1 | PACER-DC | 双差分、多端点、不确定性 | 主方法 |
| M2 | PACER-DC + pharmacology triplet | M1 表示与药理 hard negative | 数据闸门通过后的增量模块 |

M2 只有在完全留出的 chemotype 上稳定优于 M1 才能作为算法贡献。训练帧不得随机分到测试集。

## 5. 冻结动态端点

第一版只使用预先定义且可解释的端点：

- ACh 重原子 RMSD；
- ACh 原生接触占有率；
- allosteric ligand pose retention 与口袋接触覆盖；
- 冻结的受体内部距离 coupling coordinate；
- TM6 胞内端、DRY、PIF、NPxxY 等激活相关几何；
- replica 间方向一致性、有效样本量和自相关校正的不确定性。

不得查看候选功能标签后选择残基、时间窗、聚类数或端点符号。

## 6. PACER-DC 输出

- `CoupledShift`：`candidate_probe - probe_only`；
- `OrthostericStabilization`：候选对 ACh 构象和接触稳定性的改善；
- `IntrinsicActivationRisk`：`candidate_no_probe - apo` 的激活样变化；
- `BindingCompatibility`：两个配体上下文中的口袋滞留和接触；
- `Uncertainty/AD`：replica 不确定性、上下文缺失和化学适用域。

主输出为多维向量和 Pareto 关系，不使用事后调权得到的单一总分。缺少任一必要上下文时拒绝完整 PACER-DC 判定。

## 7. 验证设计

- 主切分：leave-chemotype-out；
- 辅助切分：leave-system-out 或 leave-study-out；
- 超参数只在外层训练体系中选择；
- 指标：逐分子效应、replica 一致性、macro balanced accuracy/AUC（样本量足够时）、PR-AUC、worst-group、拒绝率；
- 不确定性：分层 bootstrap，以 chemotype/分子和 replica 为层级；
- 显著性：相对 baseline 的配对 bootstrap 差值和标签置换；
- 多重比较：对多个主要端点进行预注册校正；
- 所有点估计同时报告置信区间，不以单个小样本 AUC 作为成功依据。

## 8. 主方法晋级标准

PACER-DC 只有同时满足以下条件才可称为有意义的算法增量：

1. M1 在完全留出的 chemotype 上优于 B1-2，且配对置信区间下界大于零；
2. B1-2 优于单上下文 MD，证明无 ACh 对照具有增量；
3. PAM 与内在别构激动剂的分离不能仅由二维相似度解释；
4. 结果方向在多数独立 replica 中一致；
5. 去除任何一个上下文、换随机标签或进行 frame-level 泄漏控制后，优势应消失或显著下降；
6. 至少一个外部 chemotype、公开独立体系或前瞻候选不参与任何特征和阈值选择。

若仅满足部分条件，论文主张降级为机制性 proof-of-concept 或 benchmark，不声明通用功能预测。

## 9. 候选分子报告规则

候选物仅分为：

- 通过结合门控；
- 获得动态支持；
- 存在内在激动风险；
- 证据不足/拒绝判断。

在 ACh 浓度响应、PAM Emax/EC50、候选单独给药和必要反筛完成前，统一称为“计算优先候选”，不得称为“发现的新 PAM”。

## 10. 论文最小交付物

1. 冻结数据字典、来源和药理角色表；
2. 可重建的四上下文系统及参数审计；
3. 全部 baseline 的统一运行入口；
4. 外层留出预测、逐 replica 端点及置信区间；
5. 消融、置换、适用域和失败结果；
6. 计算候选证据表及明确的功能确认边界；
7. 环境锁定、随机种子、清单和完整复现命令。
