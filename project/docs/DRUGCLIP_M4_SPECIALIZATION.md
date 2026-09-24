# DrugCLIP 的 M4 靶点特化实验

## 为什么这样迁移

DrugCLIP（NeurIPS 2023）的核心不是亲和力回归，而是把虚拟筛选改写为口袋—分子密集检索：分别编码蛋白口袋和分子，再用对比学习将真实配对拉近。论文还通过同源口袋扩增提高表征泛化。PACER 对这一优势的迁移应当发生在口袋—分子联合空间，而不是只把128维分子向量当作普通描述符。

## M4-specific re-alignment

- M4多正口袋视图：7TRQ、7TRP、7TRS；
- 正样本：有功能标签的M4 PAM；
- 困难负样本：DrugCLIP结合分数高、但实验上inactive的M4分子；
- 模型：冻结DrugCLIP，只训练rank-16 molecule/pocket residual adapters；
- 损失：不平衡BCE + PAM高于hard inactive的pairwise ranking + 预训练空间保持；
- 可选去偏：source-adversarial head；
- 验证：固定source holdout与完整series holdout，三个随机种子。

这种设计避免原始in-batch InfoNCE把不同M4 PAM错误地互当负样本，并允许不同PAM偏好不同M4构象。

## 初步结果

| 方法 | Source holdout AUC | Series holdout AUC |
|---|---:|---:|
| DrugCLIP molecule-only triplet | 0.660 | 0.515 |
| M4 pocket BCE re-alignment | 0.553 | 0.629 |
| M4 pocket hard-negative ranking | 0.564 | **0.637** |
| 既有2D RandomForest | 0.630 | 0.622 |
| 50:50 rank fusion（探索性） | — | **0.658** |

M4 pocket-ranking首次在series holdout上略高于2D baseline（+0.015），但配对分子bootstrap 95% CI为−0.057–0.087，不能确认单模型优势。固定50:50秩融合的表面增益为+0.036，分子bootstrap CI为−0.003–0.075，来源成簇CI为−0.023–0.089，接近但仍未跨过确认阈值。

## 当前判断

1. 靶点特化是可行方向；它把series AUC从molecule-only的约0.52提高到约0.64。
2. DrugCLIP提供的价值与2D模型互补，而不是完全取代2D模型。
3. Source与series结果方向不一致，提示来源混杂仍未解决。
4. 50:50融合属于看到结果后的探索，下一轮必须冻结后在新的外部数据上验证。
5. 本结果仍只预测历史功能标签，不能证明任何候选物是真实PAM。

## 后续消融：哪些路线没有成立

### 连续效力微调

使用430个M4 PAM的pEC50训练口袋回归与跨来源排序。DrugCLIP在source/scaffold holdout的Spearman分别约为0.321/0.529，低于2D RandomForest的0.472/0.637；固定50:50融合也未改善。说明预训练口袋匹配空间不适合直接承担精细效力标定，该支线关闭。

### 毒蕈碱受体家族预适配

从ChEMBL M1–M5共15,419条原始活性中整理出5,970个分子。严格要求同一分子在一个非M4亚型pChEMBL≥6、另一亚型≤5，并排除PACER全部训练/外部结构后，仅得到33个分子、54个显式选择性对。通过同源序列映射构建M1/M2/M3/M5对应的胞外别构口袋后，家族adapter的scaffold-held-out pair accuracy由DrugCLIP零样本59.3%提高到77.2%。但将其用于M4初始化后，series AUC由0.637降至0.628，source AUC也下降。说明正构/混合ChEMBL亚型选择性知识不能直接迁移为M4 PAM功能，该支线关闭。

### 分系列稳健性

50:50融合在三个series folds上的AUC分别为0.599、0.553、0.691；对应2D baseline为0.622、0.539、0.602。融合改善两个较小系列，但损害最大系列。其macro AUC约0.614、worst-series约0.553，属于互补性信号而非稳定优势，下一步需要冻结后新增独立系列，不能继续在现有三折上选权重。
