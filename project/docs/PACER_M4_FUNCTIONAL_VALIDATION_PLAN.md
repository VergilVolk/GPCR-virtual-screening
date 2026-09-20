# PACER-M4 候选功能验证计划

## 最小可执行版本

目标不是验证“是否结合”，而是区分 PAM、PAM-agonist、直接 agonist 和非活性分子。

### 第一轮：8 个候选 + 4 个历史功能对照

首轮不是从 24 个分子中任意挑高分，而是冻结为：

- 3 个 `SCOMP0003` 参考系列多样性校准锚点：PACER0024、PACER0153、PACER0039；
- 3 个同系列 exploitation 查询：PACER0076、PACER0057、PACER0026；
- 2 个跨系列探索假设：PACER0166、PACER0172；
- 两组性质匹配的历史 functional PAM/inactive 对照。历史标签不能代替本板复测。

完整批次见 `results/pacer_assay_closed_loop_v01/round1_assay_batch.csv`。三锚点偏移只允许用于 `SCOMP0003` 参考系列，绝不跨 chemotype 传播。

- 人源 CHRM4 稳定表达细胞，使用与 benchmark 尽可能一致的 calcium/FLIPR readout；
- 每个候选建议 0.03、0.1、0.3、1、3、10 µM；
- 条件 A：候选 + 预先实测的 ACh EC20；
- 条件 B：候选单独加药，检测内在激动活性；
- 阳性：VU0467154、LY2033298；阴性：vehicle 和至少两个实验非活性近邻；
- 每条件至少 3 个独立生物重复，板内技术重复；随机化孔位并对分析者隐藏候选角色。

预注册初筛命中：在无明显细胞毒性/检测干扰时，至少两个相邻浓度对 ACh EC20 呈方向一致且可重复增强，其中一个浓度的均值增强至少 1.5 倍；候选单独最大响应不超过同板 ACh Emax 的 20%。这是晋级标准，不是最终 PAM 定义。

### 第二轮：首轮命中 + 同系列剩余候选全矩阵

- ACh 8 点浓度—响应 × PAM 6–8 点浓度；
- 同时测定候选单独响应；
- 联合拟合 operational model of allosterism，报告 affinity、功能 cooperativity、efficacy 和参数置信区间；
- 不只报告“EC50 shift”，同时报告 Emax、曲线斜率、重复间方差和模型拟合诊断。

确认 PAM 的主判据为 `log(alpha*beta)` 的 95% CI 下界高于 0，同时候选单独加药不满足直接激动剂标准。若只有 Emax 改变或模型不可辨识，报告“功能增强待确认”，不强行标为 PAM。

### 第三轮：选择性与 probe dependence

- M1/M2/M3/M5 同源受体反筛；
- 至少在 ACh 与 iperoxo 两种正构探针下重复全矩阵；
- 若某候选只在特定探针下增强，应报告 probe-dependent PAM，而不是算法失败；
- 补充细胞毒性、聚集、荧光/发光干扰和基础安全性反筛。

## 预注册主终点

1. `log(alpha*beta)` 或等价功能协同性参数及 95% CI；
2. 候选单独加药的最大响应；
3. ACh Emax 与 EC50 的联合变化；
4. M4 相对 M1/M2/M3/M5 的选择性；
5. ACh 与 iperoxo 条件下参数差异。

不以 docking affinity、预测 pEC50 或单一浓度的增强百分比作为最终 PAM 证据。

## 计算—实验闭环

实验结果写入 `results/pacer_assay_closed_loop_v01/round1_results_template.csv` 后运行：

```powershell
python project/scripts/apply_pacer_fs_candidates.py --results project/results/pacer_assay_closed_loop_v01/round1_results_template.csv
```

只有 3 个锚点均通过 QC 且有有效 PAM EC50 时才生成系列 pEC50 绝对偏移；否则维持 relative-only。若要启用机制校准，三个锚点还必须在同一 probe、同一 readout 下完成全浓度矩阵，并得到通过拟合 QC 的 `functional pKB/log(alpha*beta)/log(tauB)` 及置信区间。当前证据只允许校准 `log(alpha*beta)` 与 `log(tauB)`；`pKB` 跨文献迁移失败，保持不预测。

```powershell
python project/scripts/apply_pacer_mechanism_calibration.py --results project/results/pacer_assay_closed_loop_v01/round1_results_template.csv
```

实验结果回流时必须保留完整条件和阴性数据。首轮不能在重训后对同一批报告 AUC；同系列剩余 8 个分子是封存查询集，新的 chemotype、probe 或 readout 需要自己的 3 个锚点。
