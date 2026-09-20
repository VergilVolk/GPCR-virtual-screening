# PACER-M4 首轮 PAM 候选证据包

状态：**前瞻计算假设，未经过新湿实验确认**。本文件冻结首轮优先级，后续不得根据同一批实验结果回头修改排序后再声称前瞻成功。

## Lead 1：PACER0076

SMILES：`Cc1cc2nncn2c2sc(C(=O)NC3CN(c4cncc(F)c4)C3)c(N)c12`

- 最近已知 PAM：CM00240，功能 pEC50 `7.9208`；ECFP4 Tanimoto `0.754`；
- local kNN pEC50参考 `7.61`，但不作为校准效力；
- strict inactive risk `0.026`；
- 7TRS ACh状态口袋覆盖 `0.923`；
- 十构象 dynamic Pareto：通过；BEavg `-7.087`；
- 角色：同系列 exploitation、首要功能验证分子。

## Lead 2：PACER0057

SMILES：`Cc1c(Cl)c2nncn2c2sc(C(=O)NC3CN(C4(F)CCCC4)C3)c(N)c12`

- 最近已知 PAM：CM00202，功能 pEC50 `7.7959`；Tanimoto `0.726`；
- local kNN参考 `7.71`；strict inactive risk `0.035`；
- 口袋覆盖 `0.923`；dynamic Pareto通过；BEavg `-7.031`；
- 与Lead 1提供末端疏水基变化，用于检验同系列SAR转移。

## Lead 3：PACER0026

SMILES：`Cc1c(Cl)c2nnc(C)n2c2sc(C(=O)NC3CN(c4cncc(F)c4)C3)c(N)c12`

- 最近已知 PAM：CM00594，功能 pEC50 `7.6576`；Tanimoto `0.662`；
- local kNN参考 `7.83`；strict inactive risk `0.034`；
- 口袋覆盖 `0.923`；dynamic Pareto通过；BEavg `-7.246`；
- 在三者中动态平均能量最有利，但该值不等于功能效力。

## 首轮实验与成功判据

三枚校准锚点为 PACER0024、PACER0153、PACER0039；上述三枚Lead作为同系列封存查询，另有两枚跨系列探索和两对历史功能对照。批次见 `results/pacer_assay_closed_loop_v01/round1_assay_batch.csv`。

### 生物学命中

1. ACh EC20条件下至少两个相邻浓度呈可重复增强；
2. 至少一个浓度增强不低于1.5倍；
3. 候选单独加药响应不高于同板ACh Emax的20%；
4. 全矩阵拟合的 `log(alpha*beta)` 95% CI下界高于0；
5. 无细胞毒性、聚集或读出干扰。

只有满足第4项的确认实验才称功能PAM；前三项只称初筛命中。

### 算法前瞻终点

- 三锚点仅校准SCOMP0003系列偏移；
- 其余同系列分子保持封存，测定后报告全部查询的Spearman和MAE；
- 不删除inactive，不按结果重新定义系列；
- 当前只有11个候选，低于原benchmark默认12个系列规模，因此结果须报告精确置信区间并视为小样本前瞻证据；
- 跨系列探索分子不使用SCOMP0003偏移，也不混入同系列相关性。

## 结论边界

PACER0076是当前证据最完整的首要候选，不是已发现或已验证的PAM。计算闭环已经完成；“真正找到PAM”的最后一跳是含正构探针的功能实验，任何对接、MD、XR共识或QSAR都不能替代。
