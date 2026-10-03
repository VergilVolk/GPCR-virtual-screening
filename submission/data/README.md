# 数据说明与获取方式（附件5 第三节）

提交包自带数据仅 `data/demo/`（50 分子 + 10 口袋的 512 维表征，用于演示；
提取自主仓库全量表征，种子 20261002）。**不含受限/大体量原始数据**，全量
数据按下列来源获取：

| 数据 | 来源与版本 | 获取方式 | 许可 | 用途 |
|---|---|---|---|---|
| ChEMBL 毒蕈碱族标签 | ChEMBL 33 | https://www.ebi.ac.uk/chembl/ | CC-BY-4.0 | 微调 |
| LIT-PCBA | v1.0（15 靶） | https://www.esat.kuleuven.be/it/lit-pcba/ | 学术使用 | 外部评测 |
| DUD-E GPCR 5 靶 | DUD-E | https://dude.docking.org/ | 学术使用 | 外部评测 |
| M4 GaMD 构象 | figshare 33283491（6×500ns） | figshare 公开下载 | CC-BY | 口袋集合 |
| DrugCLIP 2023 ckpt | DrugCLIP 官方 | GitHub@partridger/DrugCLIP 或主仓库快照 | Apache-2.0 | 基座 |
| Science-2026 ckpt | 合作方提供（LIT-PCBA identity-90 重训） | 主仓库 data/external/drugclip_science2026/ | 内部 | 对照基座 |
| 28,519 筛选库 | 内部构建（ChEMBL+市售可及性过滤），分子/口袋表征 npz | 主仓库 results/gpcr_drugclip_screening_v01/ensemble_embeddings.npz | 内部 | 筛选输入 |

## 数据清洗与防泄漏记录（摘要）

- 全部分子 RDKit 标准化 + 最大有机片段盐清洗（2,503→2,328，记录于主仓库 prefilter 日志）；
- 标签冲突行剔除（audit.json 有逐靶计数）；
- 训练/测试隔离：靶留出（LOSO）+ **Murcko scaffold 全局 purge**（测试靶 scaffold 一律不得入训练）；
- 家族增强行仅作训练、永不入测试；留出同靶时同亚型增强行强制排除；
- 所有评测结果带 frozen SHA256 / tag（主仓库 PROJECT_STATUS.md）。

## 划分方式

- 主基准：13 靶（B2AR/CCR2/M2R/M4R + 9 个 LIT-PCBA 靶）留一靶out ×3 种子；
- M4 决胜对照：GaMD M4R 筛选对 25,249 行（模型未见），见 MODEL_CARD。
