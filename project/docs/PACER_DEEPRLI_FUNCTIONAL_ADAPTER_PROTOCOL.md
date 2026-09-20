# PACER-DeepRLI Functional Adapter 协议

日期：2026-08-30

## 1. 方法纠正

先前 `PACER-FM-HardTriplet v0.1` 是从 ECFP4+描述符训练 MLP，并非对蛋白–配体预训练模型微调。它只否定“从零训练的小样本分子triplet”，不能否定预训练复合物表示迁移。

本路线使用 DeepRLI 官方多目标权重。DeepRLI 在蛋白–配体复合物图上联合学习 affinity scoring、pose discrimination 与 virtual screening；十层 graph transformer 后得到64维复合物向量 `z`，再进入三个任务readout。官方权重SHA256为 `94E24A2DA92FFBEEFDEFBFE20CDA8542132EDFBFF8E1CF2EB19A790F957CBA3B`，599402参数已在本机CPU环境完整加载。

## 2. 核心模型

对每个分子在7TRQ、7TRP、7TRS及GaMD代表状态上的复合物图提取预训练向量：

`z_pretrained(receptor_state, ligand_pose) -> adapter -> z_functional`

adapter初始为残差瓶颈 `64 -> 32 -> 64`，主encoder冻结。功能头分别预测：

- 系列内相对 PAM potency；
- `log alpha beta` cooperativity（有标签子集）；
- `log tauB` intrinsic agonism（有标签子集）；
- PAM/inactive margin（严格语义标签子集）。

Triplet作用在 `z_functional`，不是原始ECFP：正样本要求同assay context且功能谱接近；困难负样本优先选择DeepRLI预训练空间或ECFP中接近、但potency/cooperativity/agonism/PAM身份出现cliff的复合物。标签缺失的终点不参与相应loss。

## 3. 防止遗忘结合知识

总损失为：

`L = L_relative + lambda_t L_triplet + lambda_d L_distill + lambda_m L_multitask`

- `L_distill`保持微调前后的DeepRLI affinity/docking/screening输出与复合物距离关系；
- 第一阶段只训练adapter和功能头；
- 只有adapter阶段在内层验证通过，才允许以更低学习率解冻最后两个graph-transformer层；
- 不在同一外层系列结果上选择margin、解冻层数或loss权重。

## 4. 冻结baseline与消融

1. Vina单结构/ensemble；
2. DeepRLI官方零样本三个输出；
3. DeepRLI冻结 `z` + 线性/LightGBM readout；
4. DeepRLI adapter，仅回归；
5. DeepRLI adapter + 随机triplet；
6. DeepRLI adapter + functional hard-triplet；
7. 去掉蒸馏约束；
8. 单7TRQ vs 三静态状态 vs GaMD状态聚合；
9. PACER-FS ligand-only主基线；
10. RTMScore冻结外部评分（环境可用时）。

## 5. 评价合同

- 外层完整 `source_component/chemotype` 留出，划分后再构造triplet；
- 主要指标：macro within-series Spearman及对PACER-FS的paired series-bootstrap CI；
- 次要指标：worst-series、potency-cliff方向、多终点MAE、选择性拒绝覆盖率；
- 三个受体状态或多个MD snapshot按bag聚合，任何snapshot不得被当作独立功能标签样本；
- DeepRLI affinity改善不等于PAM效力改善；最终功能结论必须在cooperativity、intrinsic agonism及PAM/inactive终点分别验证。

内部晋级要求：hard-triplet adapter相对冻结DeepRLI readout及PACER-FS均有正增益，主要CI下界大于0，worst-series不恶化，且至少一个功能终点改善。外部SOTA仍需未参与开发的独立campaign或前瞻湿实验。

## 6. 当前执行状态

- DeepRLI源码固定于commit `9d6e73993edded01de57b765da42fe58463c670a`；
- 官方v1.0.2权重已下载并通过严格state-dict加载；
- 独立Python 3.11 CPU环境已建立，PyTorch 2.0.1 CPU、DGL 1.1.2可用；
- 430分子的7TRQ/7TRP/7TRS对接姿势已存在；
- 下一步为姿势格式转换、复合物图编译、零样本DeepRLI审计与64维 `z` 导出。

该路线当前是预注册迁移实验，不是已完成SOTA，也没有确认任何新PAM。

## 7. v0.1 冻结实现参数（结果产生前登记）

- 每个分子固定使用7TRQ/7TRP/7TRS三个seed-42 pose组成一个bag；
- 官方DeepRLI encoder完全冻结，每状态输入64维 `z`；
- 共享残差adapter：`64 -> 32 -> 64`，最后一层零初始化；
- 三状态attention pooling后接32维归一化metric head与中心化potency回归头；
- AdamW 1e-3、weight decay 1e-4、160 epochs、种子17/29/43；
- 回归损失为训练系列内中心化MSE；adapter残差L2权重0.1；triplet cosine margin 0.2、权重0.5；
- positive：同训练系列且 `abs(delta pEC50)<=0.30`；
- hard negative：同训练系列且 `abs(delta pEC50)>=1.0`，并从冻结DeepRLI mean-z cosine距离最近者中选择；
- 每个可用anchor每epoch最多采样一个positive和一个negative，各系列先等权采anchor，避免大系列triplet数量支配；
- random-triplet对照从同系列非positive中随机取negative；
- 外层完整source_component留出；三个标签盲predicted-span anchor只用于MAE offset并从排序查询中移除；
- 主比较为HardTriplet Adapter对PACER-FS、Adapter-Reg和RandomTriplet的macro Spearman配对series-bootstrap；
- 不在该外层结果上调整margin、adapter宽度、loss权重或状态数。

## 8. v0.1 冻结结果（2026-09-01）

1290个M4复合物图（430分子×7TRQ/7TRP/7TRS）全部成功，官方DeepRLI权重导出 `1290×64` 复合物表示且分数全部有限。严格11个整系列留出结果如下：

| 方法 | macro Spearman | worst-series | macro MAE |
|---|---:|---:|---:|
| PACER-FS | 0.263 | -0.005 | 0.443 |
| DeepRLI-Scoring-Mean | -0.016 | -0.455 | 0.686 |
| DeepRLI-Screening-Mean | -0.093 | -0.453 | 0.951 |
| FrozenZ-Ridge | 0.020 | -0.615 | 0.541 |
| Adapter-Reg | 0.003 | -0.511 | 0.587 |
| RandomTriplet-Adapter | 0.012 | -0.522 | 0.598 |
| HardTriplet-Adapter | 0.005 | -0.533 | 0.582 |

HardTriplet相对PACER-FS差值为 -0.258，series-bootstrap 95% CI `[-0.424,-0.095]`；相对Adapter-Reg仅+0.001，95% CI `[-0.036,+0.037]`。内部Go失败。

因此已正确检验并否定“通用结合预训练空间经小样本functional triplet即可恢复M4 PAM potency”的v0.1假设。结合/pose预训练仍可作结构质控，但不进入功能效力排序。按照预注册，不在该外层结果上解冻encoder或调margin追分；若未来重新开放，只能由新增独立功能数据或预先冻结的新迁移合同触发。
