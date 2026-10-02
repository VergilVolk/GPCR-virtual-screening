# PACER-M4 外部功能 Benchmark v02

## 1. 为什么需要两层评测

“外部 functional PAM benchmark”不能只写成一个 AUC。高质量确认性药理数据规模较小，公开初筛数据规模很大但噪声较高。因此，本项目冻结两层互补证据：

1. **定量功能层**：评价已知 M4 PAM 的效力、相对效能、通路读出及内在激动风险；
2. **大规模筛选层**：评价模型能否从真实高通量实验库中早期富集 M4 PAM 初筛阳性物。

两层数据不混合，结论也不互相替代。

## 2. 定量功能层

### 2.1 数据规模

- 5 个公开来源；
- 147 个唯一分子，32 个 Murcko scaffold；
- 398 条 endpoint 记录，去除精确训练重叠后为 384 个“分子—终点”单元；
- 13 个独立终点，包括 PAM/非 PAM、PAM pEC50、相对效能、pERK/GTPγS 分档、内在激动和 M4/M2 选择性；
- 原始标签不跨实验体系合并。

主宏评估只纳入样本数不少于 15、且类别分布可用的主要功能终点：Acadia PAM pEC50、Acadia PAM 相对效能、Suven CRE-luc pEC50、VU6025733 calcium pEC50、专利 human pERK 和 rat pERK 分档。小样本二分类仍报告，但不进入主宏分数。

### 2.2 当前冻结结果

宏分数使用 endpoint 内的方向性能力：二分类为 `2×AUC−1`，连续/序数终点为 `2×pairwise concordance−1`，随后在终点间等权平均。

| 方法 | 6 终点宏方向分数 |
|---|---:|
| official DrugCLIP | +0.008 |
| GPCR-triplet DrugCLIP | −0.054 |
| M4-matched DrugCLIP | −0.068 |
| ECFP-RF 功能模型 | −0.136 |
| 最大训练集相似度 | +0.056 |

这些结果不是“模型成功”，而是重要的边界证据：现有静态结合表示在不同化学系列和读出体系之间不能稳定预测 PAM 功能效力。个别来源内部存在信号，例如 official DrugCLIP 在 Acadia PAM pEC50 和相对效能上的 pairwise concordance 分别为 0.663 和 0.642；但在 VU6025733 系列上方向反转，因此不能用单一来源得出普遍结论。

## 3. 大规模筛选层

### 3.1 数据来源与冻结规则

PubChem AID 624126 是人 CHRM4 PAM 的细胞钙流功能初筛。公开记录包含 1,450 个 active 和 362,354 个 inactive。当前冻结子集保留全部 active，并按 `SHA256(AID624126:CID)` 排序确定性抽取 20,000 个 inactive；选择过程不读取结构或模型分数。

结构获取后共有 21,450 个唯一结构。与 529 分子历史训练表精确去重后，正式评测为：

- 21,445 个分子；
- 1,445 active、20,000 inactive；
- 15,421 个 Murcko scaffold。

该数据是单次 primary screen，只能作为筛选级标签，不能称为确认性 PAM 真值。

### 3.2 首轮二维基线

| 方法 | ROC-AUC | 校正后 AP | EF@0.5% | EF@1% |
|---|---:|---:|---:|---:|
| ECFP-Ridge | 0.535 | 0.0050 | 3.44 | 2.69 |
| ECFP-RF | 0.508 | 0.0089 | 3.44 | 2.62 |
| 最大训练集相似度 | 0.528 | 0.0102 | **4.12** | **2.83** |
| 理化性质 Logistic | 0.522 | 0.0043 | 0.41 | 0.97 |

ROC-AUC 接近随机，说明该集合不能被普通二维模型轻易解决；同时早期富集高于随机，说明历史 M4 化学知识仍能在最前端提供有限筛选价值。后续核心比较应优先看 PR-AUC、EF@0.5%、EF@1% 和 scaffold-stratified bootstrap，而不是只看 ROC-AUC。

## 4. 下一步冻结评测

对 21,445 个分子按同一顺序输出以下分数：

1. official DrugCLIP；
2. GPCR-triplet / famaug 路由 DrugCLIP；
3. 单结构 docking；
4. 六路跨构象 docking 共识；
5. DrugCLIP + docking 冻结融合；
6. 可计算子集上的四上下文动态复核。

所有模型禁止读取 AID 624126 标签调权。DrugCLIP 与 docking 评价大规模筛选和早期富集；四上下文只在预先冻结的小规模复核子集上评价，不能把未模拟分子补成动态结果。

## 5. 可成立的主张边界

- 可以称：建立了覆盖确认性药理与大规模真实初筛的多层外部 M4 PAM benchmark；
- 可以称：普通二维模型在大规模初筛上总体区分能力有限，但存在早期富集；
- 可以称：静态结合表示不能稳定预测跨体系 PAM 效力；
- 不能称：AID 624126 的 active 均为确认性 PAM；
- 不能称：现有模型已达到 functional PAM SOTA；
- 不能用该 benchmark 替代最终湿实验验证。

## 6. 复现入口

```bash
python project/scripts/build_pacer_m4_external_functional_benchmark_v02.py \
  --benchmark-root project/data/benchmarks/m4_pam_v1 \
  --training project/data/benchmarks/m4_pam_v1/pam_vs_inactive.csv \
  --output-dir project/data/benchmarks/pacer_m4_external_functional_v02

python project/scripts/acquire_pubchem_m4_pam_screen.py \
  --output-dir project/data/benchmarks/pubchem_aid624126_m4_pam_v01 \
  --inactive 20000

python project/scripts/evaluate_pubchem_m4_pam_screen_baselines.py \
  --screen project/data/benchmarks/pubchem_aid624126_m4_pam_v01/aid624126_screen_subset.csv \
  --training project/data/benchmarks/m4_pam_v1/pam_vs_inactive.csv \
  --output-dir project/results/pubchem_aid624126_m4_pam_v01
```

数据来源：PubChem [AID 624126](https://pubchem.ncbi.nlm.nih.gov/bioassay/624126)、[AID 624135](https://pubchem.ncbi.nlm.nih.gov/bioassay/624135)；确认性数据的专利、论文和哈希见各来源目录中的审计文件。
