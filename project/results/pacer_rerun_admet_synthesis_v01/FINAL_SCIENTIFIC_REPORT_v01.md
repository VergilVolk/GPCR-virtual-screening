# PACER-200 药化性质、ADMET与合成可行性综合评价 v01

分析日期：2026-10-12（Asia/Shanghai）。这是独立事后评价，不是2026年10月原始筛选门槛。冻结Top3、PACER-XR、DrugCLIP、Stage4和PACER-FKG均不因本报告改变；没有执行分子生成、Docking或MD。全部数值来自本次保存的实际计算或推理，无本次实验数据。

## 1. 输入与环境（PHASE A）

唯一分子输入：`project/results/pacer_xr_rerun_v01/pacer_xr_rerun_full200_ranking.csv`。SHA256：`723830461f517315172a47686b941b468da12d23cdb60b8e4e041cdc86eef944`。起始Git提交：`9bb03d6dfa7ef29f1da877ccdfd6426198b2ac0a`；起始结果分支：`results/pacer-stage4-pbcfix01-fkg-v02-v01`；分析分支：`analysis/pacer-rerun-admet-synthesis-v01`。初始工作区干净。

200行、200个唯一candidate_id、200个唯一原始canonical SMILES及200个RDKit规范异构SMILES；全部成功解析和价态消毒检查，candidate_id与SMILES一一对应；原始SMILES列与canonical列一致，规范化序列化差异0。Top3均从此文件读取；未使用或映射历史PACER编号。没有身份冲突。

隔离环境`.venv-pacer-admet-v01`采用Python 3.11.9，RDKit 2026.3.6；完整依赖见两个lock文件与ENVIRONMENT。原默认Python 3.14.7缺少RDKit/pandas/matplotlib；未改变项目既有环境。仓库检索未发现可直接复用的专用ADMET或逆合成权重/库存；原药化脚本涉及其他分析背景，因此本批独立实现，不借用历史分子。

## 2. 全200分子基础性质与类药性（PHASE B）

| 指标               |     min |   median |     max |
|:-----------------|--------:|---------:|--------:|
| molecular_weight | 324.771 |  466.951 | 547.422 |
| MolLogP          |   2.464 |    4.513 |   5.486 |
| TPSA             |  41.050 |   98.200 | 123.980 |
| rotatable_bonds  |   3.000 |    3.000 |   4.000 |
| FractionCSP3     |   0.042 |    0.241 |   0.619 |

| 规则       |   零违反PASS |   FAIL |
|:---------|----------:|-------:|
| Lipinski |       123 |     77 |
| Veber    |       200 |      0 |
| Ghose    |       114 |     86 |
| Egan     |       200 |      0 |
| Muegge   |       133 |     67 |

五条规则全部零违反：91/200。Lipinski严格PASS指0违反，另保存“至多1违反”字段。规则每个具体阈值、违反条件与数量均保留；QED是原始连续分数，未赋予淘汰阈值或综合权重。

QED中位数0.4115，均值0.4383，范围0.2876–0.8527。PAINS命中0/200；Brenk命中47/200，单独保存规则名称、catalog标识、输入分子0起始原子索引和片段SMILES。PAINS零命中不证明安全。

实现约定：RDKit MolLogP一致用于规则；Muegge采用MolLogP替代原始/其他实现中的脂溶性算法，Egan采用常见矩形阈值近似，不能称为原版判别器完全复现；Ghose原子数包括隐式H。QED为RDKit默认加权几何平均。规则是经验标准，通过不证明ADMET优良，违反也不等于不能成药。[规则实现背景](https://pmc.ncbi.nlm.nih.gov/articles/PMC5335600/)

## 3. ADMET模型、证据等级与不确定性（PHASE C）

实际运行ADMET-AI **2.0.1** / Chemprop 2.3.1，分类及回归各5个本地权重，CPU推理覆盖200×41=8200个原始端点输出，逐成员预测均保存。权重、元数据、SA文件SHA256及端点训练来源、单位、报告误差全部记录于MODEL_PROVENANCE。包为MIT许可证；训练来源为TDC，底层数据集许可证各异。2024论文描述v1，v2重新训练后输出不同；本报告不将v1验证指标自动视为v2或PACERGEN验证。[官方版本说明](https://github.com/swansonk14/admet_ai)

证据等级：A为可复现描述符/结构计算（不是实验验证）；B为经验规则/结构警示；C为可追溯模型的探索性预测，适用域未验证；D为缺失、失效或不可靠。所有模型的逐分子适用域为**UNKNOWN**：权重包未提供端点训练结构或经验证的适用域估计器，不能声称IN_DOMAIN或任意设置OUT_OF_DOMAIN阈值。五成员样本标准差仅表示集成分歧，不能量化全部模型偏差，不能作为置信区间或保证低误差。分类保留P(class1)，没有应用安全阈值。

在52个记录端点（41个模型+11个明确缺项）中，每个分子都有记录；总体状态计数为 `{'PREDICTED': 7766, 'NOT_ASSESSED': 2600, 'MODEL_FAILED': 34}`。PREDICTED仅表示数值推理成功且通过物理边界检查，不表示模型在本批化学空间已被验证。无实验参考数据被当成本次实验，也未计算综合成药概率。

### 水溶性

200/200 PREDICTED；中位数 -5.365，范围 -6.785–-3.975，单位log10(mol/L)。使用AqSolDB训练的回归模型；已另外保存代数转换的mol/L及mg/L，明确不是新实验。pH、盐型、晶型和温度未指定，不能视为某固定条件下的本征溶解度。[训练端点](https://tdcommons.ai/single_pred_tasks/adme/#solubility-aqsoldb)

### BBB、CNS与通透性

BBB：200/200 PREDICTED；中位数 0.531，范围 0.161–0.980。这是穿透分类概率，不是logBB、脑组织浓度或Kp,uu。实际脑暴露200/200 NOT_ASSESSED。CNS MPO仅保存MW、MolLogP代理、TPSA、HBD四个单项分量及明确标记的0–4部分和；缺乏pKa与pH7.4 logD，完整0–6 CNS MPO未计算。[CNS MPO原始方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC3368654/)

Caco2：200/200 PREDICTED；中位数 -5.075，范围 -5.356–-4.757。模型包声明`log(10^-6 cm/s)`，TDC榜单标注`cm/s`，两者不一致；本次保留原始模型尺度并明确标记单位冲突，**没有换算绝对通透率或作高/低通透判定**。[TDC端点单位](https://tdcommons.ai/benchmark/admet_group/01caco2/)

PAMPA为分类概率，200/200完成；它是人工膜体系，不与Caco2或BBB等价。P-gp抑制：200/200 PREDICTED；中位数 0.886，范围 0.324–0.984。P-gp底物200/200 NOT_ASSESSED；不能从抑制概率推出外排底物性，CNS外排风险仍是关键证据缺口。

### CYP、微粒体与蛋白结合

五个CYP抑制端点及CYP2C9/2D6/3A4底物均完成200/200预测；底物与抑制端点独立保存。代谢位点200/200 NOT_ASSESSED，无原子位点或确定代谢途径结论。

人微粒体CLint：200/200 PREDICTED；中位数 52.501，范围 18.549–95.469，µL/min/mg；随包报告R²≈0.277、MAE≈26.74，精度有限，适用域未知。肝细胞CLint也完成200/200，保留其原始µL/min/10^6 cells单位及混合物种来源说明。微粒体半衰期与稳定性分类均未评价，没有进行体内外缩放。

PPBR：166/200 PREDICTED；中位数 97.061，范围 87.760–99.980。另34个输出超过100%，保留raw_prediction并将prediction置空，状态MODEL_FAILED；不截断为100%。Fraction unbound没有独立可靠模型，未伪装为另一个预测结果。

### 清除率、半衰期与口服生物利用度

系统清除率0/200可靠预测；体外CLint没有被称作全身清除率。Half_Life_Obach虽实际推理200次，但随包元数据R²≈−2.386、MAE≈31.85 h，缺乏可靠支持，全部NOT_ASSESSED，原始模型值仅归档。VDss模型同样因负R²未采纳。Bioavailability_Ma完成200个分类概率，不能解释为口服F百分比；连续F为NOT_ASSESSED。未拼接独立模型计算任何PK参数。

### 毒性及反应性代谢物

hERG：200/200 PREDICTED；中位数 0.701，范围 0.259–0.983；Ames：200/200 PREDICTED；中位数 0.575，范围 0.112–0.903；DILI：200/200 PREDICTED；中位数 0.992，范围 0.567–0.999。另保留ClinTox、Tox21等模型输出，不能将受体/应激端点活性等同于已确认系统毒性。DILI是特定数据集端点，未据此生成独立泛肝毒性结论。任何低概率都不是“安全”结论。

非穷尽性结构匹配发现190/200命中反应性或生物活化相关motif（具体规则、SMARTS、原子及片段单独保存）。所有此类结果均标为STRUCTURAL_ALERT；包含噻吩、呋喃、硝基芳环、酚及少数反应性官能团，不能推断这些分子必然形成反应性代谢物。环融合、取代位置、竞争代谢与解毒都会改变意义。专门反应性代谢物模型未运行。[结构警示研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC5871347/)

## 4. 逆合成、采购、收率与成本（PHASE D）

在另一隔离环境安装AiZynthFinder 4.4.1，修复PyArrow 26与NumPy 1.26不兼容问题（固定PyArrow 18.1.0），完成实际导入及对象初始化。正式搜索前检查策略模型、模板、库存及许可证。公开模型Zenodo记录为CC-BY-4.0；代码MIT。Zenodo首次预检超时，随后沙箱外复核模型URL返回HTTP200；官方ZINC库存标准GET仍返回403，本地未找到所需完整资产，因此没有继续下载策略权重并声称可搜索。没有绕过限制，也没有把安装包当作完成逆合成。[官方工具与资产](https://github.com/MolecularAI/aizynthfinder)

因此200个分子均为NOT_ASSESSED，真正运行路线搜索0个、获得路线0个；**NO_ROUTE_FOUND为0个**。这里的0路线是未评价结果，不能说明分子不可合成。没有路线中间体、起始原料、步骤或路线评分可提供。SA Score全部计算，只用于相对结构复杂度描述，不能替代路线可行性、实验收率或成本。

起始原料现实采购核实0个，AVAILABILITY_NOT_VERIFIED；无供应商、CAS、产品链接、库存或报价被虚构。反应步数0个可靠证据，收率和总成本均200/200 NOT_ESTIMATED；币种、日期、纯化/人工/耗材及收率考虑字段留空，因为尚无路线和真实报价基础。保护/脱保护、关键转化和特殊条件同样不作无证据推断。

## 5. 立体化学与盐型

RDKit识别的已指定四面体中心分子0个、未指定四面体中心分子0个、未指定双键立体分子0个。所有分子为中性单组分。按保留已指定构型、仅枚举未指定构型的方法，理论异构体计数与枚举上限均已保存。零RDKit中心不排除阻转异构、构象或实验稳定立体形式。

结构存在可能可电离的胺/芳香氮等motif；匹配不意味着这些位点都在生理pH质子化。未运行pKa或pH依赖质子化预测，未确定最佳盐型，原始SMILES从未被改变。

## 6. 冻结Top3完整比较

下列列名中的预测值均为模型原尺度；所有概率适用域未知。

| candidate_id   |   molecular_weight |   MolLogP |    TPSA |   FractionCSP3 |   QED |   Solubility_AqSolDB__prediction |   BBB_Martins__prediction |   hERG__prediction |   AMES__prediction |   DILI__prediction |   SA_score |
|:---------------|-------------------:|----------:|--------:|---------------:|------:|---------------------------------:|--------------------------:|-------------------:|-------------------:|-------------------:|-----------:|
| PACERGEN02432  |            516.034 |     5.262 | 103.880 |          0.200 | 0.332 |                           -5.888 |                     0.447 |              0.846 |              0.479 |              0.997 |      3.575 |
| PACERGEN00123  |            480.041 |     3.118 |  74.720 |          0.417 | 0.468 |                           -5.194 |                     0.862 |              0.975 |              0.794 |              0.851 |      3.472 |
| PACERGEN00462  |            494.068 |     3.426 |  74.720 |          0.440 | 0.453 |                           -5.054 |                     0.905 |              0.973 |              0.836 |              0.800 |      3.461 |

### PACERGEN02432

A. MW 516.034; MolLogP 5.262; TPSA 103.88; 6 aromatic rings; lower saturation.

B. Lipinski 2 violations (MW/logP), Ghose 2, Muegge 1; Veber/Egan pass.

C. Exploratory 5-member ADMET-AI v2 outputs; per-endpoint status/SD retained; AD unknown. Solubility and permeability require assay confirmation.

D. BBB classification only; no brain exposure, Kp,uu, Pgp substrate, full CNS MPO, logD or pKa evidence.

E. hERG, DILI and CYP3A4 inhibition model probabilities warrant experimental confirmation; no safe label.

F. SA heuristic only. Route NOT_ASSESSED; no reaction step count, verified starting materials, yield or cost.

G. No RDKit potential tetrahedral/double-bond stereocenters found; neutral single component. Ionizable motifs present, no pKa/speciation or best-salt evidence.

H. Unknown AD and uncalibrated probabilities; no experimental ADMET; Caco2 units unresolved; microsomal/systemic PK gap; retrosynthesis asset/stock unavailable.

| endpoint               |   prediction |   raw_prediction |   ensemble_sample_std | unit            | status       |
|:-----------------------|-------------:|-----------------:|----------------------:|:----------------|:-------------|
| BBB_Martins            |        0.447 |            0.447 |                 0.176 | probability     | PREDICTED    |
| Bioavailability_Ma     |        0.753 |            0.753 |                 0.107 | probability     | PREDICTED    |
| CYP1A2_Veith           |        0.439 |            0.439 |                 0.219 | probability     | PREDICTED    |
| CYP2C19_Veith          |        0.856 |            0.856 |                 0.085 | probability     | PREDICTED    |
| CYP2C9_Veith           |        0.911 |            0.911 |                 0.034 | probability     | PREDICTED    |
| CYP2D6_Veith           |        0.084 |            0.084 |                 0.040 | probability     | PREDICTED    |
| CYP3A4_Veith           |        0.924 |            0.924 |                 0.062 | probability     | PREDICTED    |
| Pgp_Broccatelli        |        0.975 |            0.975 |                 0.009 | probability     | PREDICTED    |
| Caco2_Wang             |       -5.308 |           -5.308 |                 0.088 | log(10^-6 cm/s) | PREDICTED    |
| Clearance_Microsome_AZ |       73.637 |           73.637 |                 4.860 | uL/min/mg       | PREDICTED    |
| PPBR_AZ                |      nan     |          100.618 |                 2.764 | %               | MODEL_FAILED |
| Solubility_AqSolDB     |       -5.888 |           -5.888 |                 0.181 | log(mol/L)      | PREDICTED    |
| AMES                   |        0.479 |            0.479 |                 0.072 | probability     | PREDICTED    |
| DILI                   |        0.997 |            0.997 |                 0.001 | probability     | PREDICTED    |
| hERG                   |        0.846 |            0.846 |                 0.064 | probability     | PREDICTED    |

预测溶解浓度：1.29e-06 mol/L，0.667 mg/L（仅模型logS的代数换算）。

### PACERGEN00123

A. MW 480.041; MolLogP 3.118; TPSA 74.72; FractionCSP3 0.417.

B. Lipinski/Veber/Egan/Muegge pass, Ghose 2 violations (MW just above 480 and MR).

C. Exploratory 5-member ADMET-AI v2 outputs; per-endpoint status/SD retained; AD unknown. Solubility and permeability require assay confirmation.

D. BBB classification only; no brain exposure, Kp,uu, Pgp substrate, full CNS MPO, logD or pKa evidence.

E. hERG, DILI and CYP3A4 inhibition model probabilities warrant experimental confirmation; no safe label.

F. SA heuristic only. Route NOT_ASSESSED; no reaction step count, verified starting materials, yield or cost.

G. No RDKit potential tetrahedral/double-bond stereocenters found; neutral single component. Ionizable motifs present, no pKa/speciation or best-salt evidence.

H. Unknown AD and uncalibrated probabilities; no experimental ADMET; Caco2 units unresolved; microsomal/systemic PK gap; retrosynthesis asset/stock unavailable.

| endpoint               |   prediction |   raw_prediction |   ensemble_sample_std | unit            | status       |
|:-----------------------|-------------:|-----------------:|----------------------:|:----------------|:-------------|
| BBB_Martins            |        0.862 |            0.862 |                 0.136 | probability     | PREDICTED    |
| Bioavailability_Ma     |        0.676 |            0.676 |                 0.069 | probability     | PREDICTED    |
| CYP1A2_Veith           |        0.512 |            0.512 |                 0.199 | probability     | PREDICTED    |
| CYP2C19_Veith          |        0.626 |            0.626 |                 0.218 | probability     | PREDICTED    |
| CYP2C9_Veith           |        0.541 |            0.541 |                 0.127 | probability     | PREDICTED    |
| CYP2D6_Veith           |        0.648 |            0.648 |                 0.133 | probability     | PREDICTED    |
| CYP3A4_Veith           |        0.856 |            0.856 |                 0.063 | probability     | PREDICTED    |
| Pgp_Broccatelli        |        0.882 |            0.882 |                 0.076 | probability     | PREDICTED    |
| Caco2_Wang             |       -5.158 |           -5.158 |                 0.131 | log(10^-6 cm/s) | PREDICTED    |
| Clearance_Microsome_AZ |       48.518 |           48.518 |                 7.840 | uL/min/mg       | PREDICTED    |
| PPBR_AZ                |      nan     |          100.147 |                 1.401 | %               | MODEL_FAILED |
| Solubility_AqSolDB     |       -5.194 |           -5.194 |                 0.400 | log(mol/L)      | PREDICTED    |
| AMES                   |        0.794 |            0.794 |                 0.181 | probability     | PREDICTED    |
| DILI                   |        0.851 |            0.851 |                 0.117 | probability     | PREDICTED    |
| hERG                   |        0.975 |            0.975 |                 0.020 | probability     | PREDICTED    |

预测溶解浓度：6.4e-06 mol/L，3.07 mg/L（仅模型logS的代数换算）。

### PACERGEN00462

A. MW 494.068; MolLogP 3.426; TPSA 74.72; FractionCSP3 0.440.

B. Lipinski/Veber/Egan/Muegge pass, Ghose 2 violations (MW/MR).

C. Exploratory 5-member ADMET-AI v2 outputs; per-endpoint status/SD retained; AD unknown. Solubility and permeability require assay confirmation.

D. BBB classification only; no brain exposure, Kp,uu, Pgp substrate, full CNS MPO, logD or pKa evidence.

E. hERG, DILI and CYP3A4 inhibition model probabilities warrant experimental confirmation; no safe label.

F. SA heuristic only. Route NOT_ASSESSED; no reaction step count, verified starting materials, yield or cost.

G. No RDKit potential tetrahedral/double-bond stereocenters found; neutral single component. Ionizable motifs present, no pKa/speciation or best-salt evidence.

H. Unknown AD and uncalibrated probabilities; no experimental ADMET; Caco2 units unresolved; microsomal/systemic PK gap; retrosynthesis asset/stock unavailable.

| endpoint               |   prediction |   raw_prediction |   ensemble_sample_std | unit            | status    |
|:-----------------------|-------------:|-----------------:|----------------------:|:----------------|:----------|
| BBB_Martins            |        0.905 |            0.905 |                 0.096 | probability     | PREDICTED |
| Bioavailability_Ma     |        0.721 |            0.721 |                 0.072 | probability     | PREDICTED |
| CYP1A2_Veith           |        0.442 |            0.442 |                 0.212 | probability     | PREDICTED |
| CYP2C19_Veith          |        0.584 |            0.584 |                 0.226 | probability     | PREDICTED |
| CYP2C9_Veith           |        0.399 |            0.399 |                 0.112 | probability     | PREDICTED |
| CYP2D6_Veith           |        0.653 |            0.653 |                 0.156 | probability     | PREDICTED |
| CYP3A4_Veith           |        0.820 |            0.820 |                 0.074 | probability     | PREDICTED |
| Pgp_Broccatelli        |        0.870 |            0.870 |                 0.076 | probability     | PREDICTED |
| Caco2_Wang             |       -5.151 |           -5.151 |                 0.100 | log(10^-6 cm/s) | PREDICTED |
| Clearance_Microsome_AZ |       45.470 |           45.470 |                 4.469 | uL/min/mg       | PREDICTED |
| PPBR_AZ                |       98.885 |           98.885 |                 1.563 | %               | PREDICTED |
| Solubility_AqSolDB     |       -5.054 |           -5.054 |                 0.287 | log(mol/L)      | PREDICTED |
| AMES                   |        0.836 |            0.836 |                 0.132 | probability     | PREDICTED |
| DILI                   |        0.800 |            0.800 |                 0.133 | probability     | PREDICTED |
| hERG                   |        0.973 |            0.973 |                 0.027 | probability     | PREDICTED |

预测溶解浓度：8.83e-06 mol/L，4.37 mg/L（仅模型logS的代数换算）。

PACERGEN02432的高MW/脂溶性、较高TPSA及较低预测logS构成更明显的药化暴露挑战；其hERG概率较另外两者低，但仍不能归为低风险。PACERGEN00123与PACERGEN00462具有更低TPSA、更低MolLogP及更高FractionCSP3，BBB模型倾向更高；两者hERG及Ames概率提示优先验证需求。00462的BBB/logS模型均值略优于00123，但五模型分歧存在，不能把微小差异当作实验优越性。三者的P-gp抑制、CYP3A4抑制和DILI输出均值得确认；P-gp底物性、脑游离暴露及真实合成路线仍未解决。这里不计算候选总分、不排序、不重新定义PAM活性。

## 7. 证据缺口、后续实验与比赛展示

首要证据缺口：模型逐分子适用域与概率校准、Caco2单位冲突、pKa/logD、P-gp底物、真实Kp,uu/脑暴露、代谢位点、可靠半衰期/全身清除率、连续口服F、专门反应性代谢物模型、完整计算路线及真实采购/收率/成本。PPBR超界尤其说明不能盲信输出。可考虑针对冻结Top3补充实际溶解度（明确pH/盐/晶型）、MDR1外排、hERG、Ames、微粒体稳定性/代谢物与蛋白结合试验，并开展药化人员路线评审；这些仅为建议，没有在本次执行实验。

适合比赛展示：输入身份审计；RDKit性质分布、QED及透明规则计数；Top3结构及实算描述符；带“模型预测、AD未知、非实验”标记的logS/BBB/CYP/hERG等；未评价项与失败项组成的证据矩阵。不得展示为“已证实脑暴露良好”“安全无毒”“已完成合成路线”“低成本可采购”或“10月原始筛选已通过的门槛”。Figure 1–6均从落盘CSV生成，PNG 300 dpi及矢量PDF；图中的灰色为未评价，不是PASS。

## 8. 复现与文件

从仓库根目录运行（环境依赖精确版本见lock；需可访问公开PyPI）：

```sh
/usr/local/bin/python3 -m venv .venv-pacer-admet-v01
.venv-pacer-admet-v01/bin/pip install -r project/results/pacer_rerun_admet_synthesis_v01/requirements_admet_lock_v01.txt
/usr/local/bin/python3 -m venv .venv-pacer-retro-v01
.venv-pacer-retro-v01/bin/pip install -r project/results/pacer_rerun_admet_synthesis_v01/requirements_retro_lock_v01.txt
MPLCONFIGDIR=/tmp/pacer_mpl .venv-pacer-admet-v01/bin/python project/scripts/pacer_rerun_medchem_v01.py --output /tmp/pacer_admet_reproduction
```

后续阶段以相同`--output`运行pacer_rerun_admet_v01.py、pacer_rerun_synthesis_v01.py（retro环境）、pacer_rerun_report_v01.py（admet环境）。D阶段需要把本次EXTERNAL_ASSET_PREFLIGHT_v01.json审计快照复制至新目录，复现的是本次资产不可用条件；若未来资产可用，应另建版本、实现并验证正式搜索，不可把旧失败状态当作新网络状态。报告与图表阶段只读取保存的计算结果，不再次推理。脚本以独占创建方式防止覆盖已有结果。

输出包括规定的13个CSV/JSON/报告文件、补充端点级原始模型均值/五成员输出、结构警示、反应性motif、CNS单项、溶解浓度换算、模型元数据、阶段完成记录、两个依赖lock、Top3结构图、6组PNG/PDF。ALL200_ADMET与TOXICITY为长表，每个candidate_id×endpoint一行，prediction空值代表缺失/失败，raw_prediction保留已运行但不可采用的模型原值；0从不替代失败。

Git仅创建指定分析分支；未commit、push或创建PR。所有输出原本被.gitignore排除，版本管理采用仅针对本分析目录的精确例外；冻结目录与既有排序未修改。最终Git状态、文件SHA256和验证见ANALYSIS_VALIDATION与DELIVERY_MANIFEST。
