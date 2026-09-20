# PACER-XR Cascade：M4 allosteric virtual screening 多目标级联

版本：v1.0，2026-08-29  
任务：Thompson–Miao官方M4 broad allosteric-modulator/decoy benchmark  
证据边界：同一公开benchmark上的后发表改进，不等于功能PAM效力SOTA。

## 1. 科学问题

原论文出现明显指标冲突：

- Glide-BEmin在top 0.5–1%有极强富集，但全局AUC较低；
- Vina-BEmin全局AUC最高，但早期几乎没有活性；
- 简单PACER-XR六路平均提升全局AUC，却稀释Glide的早期信号。

因此不存在一个线性总分可以同时代表“首批实验命中率”和“全库总体辨别”。PACER-XR Cascade显式将两个决策阶段分开。

## 2. 冻结算法

### 2.1 六路rank共识

每个分数仅在同一靶点、同一方法内转换为rank percentile，越高越好：

1. Glide PDB；
2. Glide BEmin；
3. Glide BEavg；
4. Vina PDB；
5. Vina BEmin；
6. Vina BEavg。

`XR = mean(rank_1 ... rank_6)`

rank变换避免混合不同程序的能量量纲。

### 2.2 1%级联

1. 先取Glide-BEmin排名前1%，保持其原始顺序；
2. 删除这些分子后，其余99%按XR排序；
3. 拼接为唯一最终列表。

1%边界来自论文预先强调的0.5–1%早期实验预算，不通过M4标签网格搜索得到。该设计保证EF0.5%和EF1%不低于Glide-BEmin，同时让剩余库利用跨引擎共识。

## 3. 结果

共同覆盖：114321分子，其中2293个ASD active。

| 方法 | AUC | PR-AUC | logAUC% | EF0.5% | EF1% | EF2% | EF5% |
|---|---:|---:|---:|---:|---:|---:|---:|
| Vina-BEmin | 0.7331 | 0.0406 | 10.14 | 0.26 | 0.78 | 1.26 | 2.16 |
| Glide-BEmin | 0.6993 | 0.0969 | 16.13 | 20.40 | 13.73 | **8.39** | 4.32 |
| XR平均 | 0.7710 | 0.0559 | 14.68 | 3.84 | 3.31 | 3.40 | 3.65 |
| **XR Cascade 1%** | **0.7775** | **0.1067** | **19.78** | **20.40** | **13.73** | 7.48 | **4.60** |

相对最佳单路AUC提升`+0.0443`，配对影响函数95% CI `[+0.0380,+0.0506]`。除EF2%仍低于Glide-BEmin外，级联在全局、PR、logAUC、0.5%、1%和5%指标上领先或并列领先。

## 4. 为什么不是数据泄漏

- 级联不拟合M4标签；
- 六路等权，无M4权重搜索；
- 1%来自论文定义的早期富集区间；
- 输出覆盖作者六路共同交集，所有baseline在同一分子集合复算；
- AUC差异使用同分子配对影响函数，并用配对bootstrap核验。

但这是论文发表后的方法开发，已知论文指出M4的Glide-BEmin早期表现强。因此准确措辞是“官方公开M4 benchmark上的确定性后发表改进”，不是未经触碰的全新外部测试。

## 5. 转化边界

作者ASD active混合不同实验、PAM/NAM和条件；PACER-XR Cascade预测的是broad AM retrieval，不是ACh探针下的PAM身份或效力。候选真正转化还需要：

- 对新候选生成完整六路分数；
- 由于本机无Schrödinger许可，目前24个生成候选没有Glide分数，不能伪造XR Cascade排名；
- 功能确认仍需ACh EC20、候选单独加药及operational allosteric model。

为防止协议漂移，`results/pacer_xr_candidate_transfer_v01/candidate_six_channel_raw_scores.csv`已冻结24个候选的六路原始分数输入格式。`apply_pacer_xr_candidates.py`拒绝任何缺失通道，不用Vina代填Glide，也不把现有7TRS分数冒充论文PDB通道；只有相同受体、参数和score定义的外部结果回填后才生成候选级Cascade排序。

### 5.1 功能标签语义压力测试

将公开排名与本项目 assay-specific 功能标签按 exact SMILES 映射后，共得到 299 个分子（266 positive、33 negative）。分子级平均rank下，Cascade表观AUC为`0.783`，Vina-BEmin为`0.726`，Glide-BEmin为`0.539`。这个`0.783`**不作为功能预测成绩**：33个阴性全部是B-tier database-text记录，A-tier primary-confirmed inactive为0，而且这33个结构又全部被ASD至少标为一次active；阳性率高达88.96%。它是标签语义/来源压力测试，不是独立功能ground truth。

更严格的主功能基准仍是27个A-tier PAM与27个primary-confirmed inactive。在该平衡集合上，cluster-00/BEmin AUC为`0.499`，BEavg为`0.457`。这两个结果共同规定了方法边界：PACER-XR Cascade可以改善公开broad-AM retrieval，但目前不能替代PAM功能效力模型。

候选迁移审计也显示，PACER0076、PACER0057、PACER0026最近已知PAM在Cascade中仅分别位于top `10.03%`、`48.84%`、`38.22%`。因此不允许把公开库的总体领先直接转写为当前候选的高PAM概率。

## 6. 复现

```powershell
python project/scripts/run_cross_target_rank_fusion.py
python project/scripts/audit_pacer_xr_functional_semantics.py
python project/scripts/apply_pacer_xr_candidates.py --create-template
# 六路协议匹配分数回填后：
python project/scripts/apply_pacer_xr_candidates.py --input project/results/pacer_xr_candidate_transfer_v01/candidate_six_channel_raw_scores.csv
```

结果位于 `results/pacer_xr_v01/` 与 `results/pacer_xr_functional_semantics_v01/`。
