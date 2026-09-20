# Thompson–Miao 2026 GPCR-AM ensemble docking：复现与创新审计

冻结日期：2026-08-29  
论文：*Benchmarking Docking Protocols for GPCR Allosteric Modulators*，bioRxiv DOI `10.64898/2026.08.12.744492`  
作者仓库 commit：`44798c841ee77230b1f89fe41074e41b458c7070`

## 1. 已完成的原文级复现

### 1.1 官方分数与官方分析代码

直接执行作者 `auc_calculator.py`、`enrichment_PDB.py`，输入作者发布的 M4R Glide/Vina PDB、BEmin、BEavg 六张分数表。没有重排标签或更换 decoy。

| 方法 | n | actives | EF0.5% | EF1% | AUC% | logAUC% |
|---|---:|---:|---:|---:|---:|---:|
| Glide PDB | 116392 | 2313 | 11.499 | 8.646 | 67.400 | 11.612 |
| Glide BEmin | 117642 | 2314 | 20.370 | 13.692 | 70.042 | 16.097 |
| Glide BEavg | 117642 | 2314 | 6.128 | 5.097 | 71.384 | 11.625 |
| Vina PDB | 117568 | 2314 | 0.778 | 1.512 | 72.527 | 10.618 |
| Vina BEmin | 115771 | 2294 | 0.349 | 0.784 | 73.377 | 10.180 |
| Vina BEavg | 115771 | 2294 | 0.959 | 1.482 | 69.577 | 8.972 |

当前仓库与论文表4有三处超过舍入误差的差异：Glide-BEavg AUC `71.384 vs 71.83`；Vina-BEmin EF0.5% `0.349 vs 0.261`；对应 EF′ `0.235 vs 0.199`。其余 M4 主指标复现到舍入精度。最可能原因是投稿表格与当前 GitHub score table 存在版本漂移或表格录入差异；不能自行选择更好看的值。

逐文件 SHA256、论文值、代码值和 delta 见 `results/m4_official_exact_reproduction/m4_table4_exact_reproduction.csv`。

### 1.2 GaMD 公共轨迹独立聚类

Figshare `33283491` 的17个文件全部下载并按字节数校验，总计 `561,993,763` bytes。公开坐标包含六条500 ns轨迹的1 ns stride版本，共3000个蛋白帧。

按作者 M4 CPPTRAJ 输入中的23个 residue、239个重原子做质量加权参考对齐，再用 average-linkage、`k=10` 聚类：

- 公共stride主簇占比 `82.77%`；作者全流统计为 `79.5%`；
- rank-matched cluster fraction L1差异 `0.081`；
- 十个公共stride代表构象与十个作者全流代表构象 Hungarian 一一匹配；
- 口袋RMSD均 `<2 Å`，均值 `1.12 Å`，最大 `1.73 Å`。

这独立支持“单一主导低能构象 + 少量次级状态”的M4结构结论。它不是全分辨率CPPTRAJ逐帧完全重跑：作者全流summary约299万帧，公开下载是3000个1 ns stride帧。

### 1.3 本地 Vina 协议

已按论文执行 Vina 1.2.7、molscrub 0.1.1 pH7状态枚举、Meeko、30 Å盒、exhaustiveness 8、9 poses、energy range 3，并在作者十个M4 GaMD cluster上完成：

- 27 primary-confirmed PAM + 27 primary-confirmed inactive：810个状态/cluster任务，零失败；
- cluster-00与BEmin AUC均约0.499，BEavg 0.457；
- 24个候选：300个状态/cluster任务，零失败。

该功能真值结果与ASD broad-AM benchmark不是同一任务，不能混合报告。

## 2. 尚未假装完成的部分

- 未从头运行3/6×500 ns全分辨率GaMD生产模拟；当前完成的是公共stride独立聚类。
- 未从头运行11.7万×10的Glide HTVS；没有Schrödinger商业许可。已精确复算作者发布的Glide分数。
- 论文Methods概述写“三条500 ns”，M4 Figshare与作者M4 CPPTRAJ输入均为六条轨迹；版本差异已冻结。

## 3. PACER-XR与PACER-XR Cascade：跨引擎创新

每个靶点内，将六路分数转换成“越高越好”的rank percentile，消除Glide/Vina量纲差异。最稳健方法不拟合参数，直接等权平均六路排名：

`XR = mean(Glide-PDB, Glide-BEmin, Glide-BEavg, Vina-PDB, Vina-BEmin, Vina-BEavg)`

在六路都有结果的114321个M4分子交集上：

| 方法 | ROC-AUC | PR-AUC | logAUC% | EF0.5% | EF1% |
|---|---:|---:|---:|---:|---:|
| 最佳单路 Vina-BEmin | 0.7331 | 0.0406 | 10.142 | 0.261 | 0.784 |
| Glide-BEmin | 0.6993 | 0.0969 | 16.127 | **20.396** | **13.728** |
| **PACER-XR 等权共识** | **0.7710** | 0.0559 | 14.677 | 3.835 | 3.312 |

PACER-XR相对最佳单路的全局AUC提升`+0.0379`，配对影响函数95% CI `[+0.0320,+0.0438]`；500次分层配对bootstrap核验为约`[+0.0316,+0.0439]`。M4标签没有参与规则设计或拟合；等权共识在看M4结果前已定义。一个只在M2R/CCR2/β2AR训练的logistic版本得到0.7695，但跨开发靶点LOTO很差，故不作为主方法。

结论边界：PACER-XR是当前作者官方M4 broad-AM交集库上更强的**全局辨别**方法，不是早期富集全面SOTA。Glide-BEmin仍显著主导top 0.5–1%。更不能把ASD AM标签解释为ACh条件下的功能PAM效力。

进一步采用固定1%级联：前1%保持Glide-BEmin，剩余99%按XR排序。该方法取得AUC `0.7775`、PR-AUC `0.1067`、logAUC `19.78%`、EF0.5% `20.40`、EF1% `13.73`、EF5% `4.60`；相对最佳单路AUC提升`+0.0443`，95% CI `[+0.0380,+0.0506]`。除EF2% `7.48`仍低于Glide-BEmin `8.39`外，其余主指标领先或并列领先。详见 `docs/PACER_XR_CASCADE_METHOD.md`。

### 3.1 功能语义审计

公开XR排名与assay-specific功能集的exact-SMILES交集为299个分子。Cascade表观AUC为`0.783`，但该集合有266 positive和33 negative；33个negative全部是B-tier database-text，A-tier confirmed negative为0，而且33个又全部是ASD active。这个结果在程序与报告中强制标记为`promotable_functional_metric=False`，不能用来声称功能PAM预测SOTA。

真正的功能压力测试仍是平衡的27+27 A-tier集合，结果约随机。两条证据不矛盾：XR解决的是公开broad-AM排序，A-tier实验检验的是指定探针下的功能增强。候选最近邻在公开Cascade中的top比例为PACER0076 10.03%、PACER0057 48.84%、PACER0026 38.22%，说明总体benchmark优势尚未覆盖所有候选chemotype。

## 4. 当前分子发现结果

首要可检验假设为 `PACER0076`：

- 与已知 `CM00240` Tanimoto `0.754`，后者功能pEC50 `7.9208`；
- local kNN参考pEC50 `7.61`，strict inactive risk `0.026`；
- 7TRS口袋覆盖 `0.923`，十构象dynamic Pareto通过；
- SMILES：`Cc1cc2nncn2c2sc(C(=O)NC3CN(c4cncc(F)c4)C3)c(N)c12`。

备选为 `PACER0057`、`PACER0026`。这些是经过化学域、功能SAR、静态/动态口袋和风险审计的候选，不是已确认PAM。确认必须依赖ACh EC20增强、候选单独加药、全浓度矩阵与operational allosteric model。

六路候选迁移合同位于`results/pacer_xr_candidate_transfer_v01/`。当前三路Glide为空，因此程序不会输出伪造的候选XR名次；外部许可节点按论文协议回填全部六路后方可应用冻结级联。

## 5. 一键复现

```powershell
python project/scripts/run_thompson_miao_reproduction.py
```

加 `--recluster` 可重新执行3000帧聚类和代表构象匹配；默认复用已校验输出。
