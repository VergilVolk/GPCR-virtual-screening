# 冻结复现包：CHRM4 PAM 动态—药理协同智能筛选

本目录是答辩与评审的最小可核验入口，对应四个模块：

1. 候选生成与药物化学过滤；
2. GPCR 迁移后的 DrugCLIP 口袋匹配；
3. Glide/Vina 六路构象共识排序；
4. apo、ACh、候选、候选+ACh 四条件受体动力学分析。

前三个模块用于优先排序可能的 M4 结合候选；第四模块寻找可重复的条件依赖构象变化。任何单一计算步骤都不能证明 PAM 活性。

## 快速核验

已在 Windows 11、Python 3.11、CPU 环境测试：

```bash
python -m pip install -r requirements.txt
python predict.py
python verify_submission.py
```

`verify_submission.py` 会对 2,605 个分子完成 Module 2 的真实前向计算，核对三个模型权重哈希、Top-200 顺序、分子身份、六路融合、四条件记录和最终结果表。成功时输出：

```text
PASS: frozen four-module replay, model hashes and final table
```

## 冻结结果层级

| 文件 | 含义 | 数量 |
| --- | --- | ---: |
| `results/01_generated_candidates_2605.csv` | 生成且通过初筛的分子 | 2,605 |
| `results/01_candidate_audit_2592.csv` | 去除已知参考分子后的审计集合 | 2,592 |
| `results/02_drugclip_top200.csv` | DrugCLIP 排序前 200 | 200 |
| `results/03_six_channel_8.csv` | 5 个新候选 + 3 个历史锚点 | 8 |
| `results/04_four_context_3.csv` | 进入四条件动力学的历史锚点 | 3 |
| `results/results.csv` | 标准化的八候选证据表 | 8 |

三名历史锚点 PACER0010、PACER0027、PACER0073 的 `functional_prediction` 均为 `unresolved`，`experimental_pam_validation` 均为 `not_performed`。

## 可追溯的 Module 4 计算样例

`data/examples/module4_pacer0073_r1/` 提供 PACER0073 第一重复每个条件 20 帧的真实生产轨迹样例及对应 PDB 拓扑。`run_module4_example.py` 会执行冻结编码、数值状态应用、四条件对比和图/区域分析；它不是完整三重复 Stage 4 结果的重跑。

该样例使用独立的 PyTorch 2.6 环境，不能与 PyTorch 2.8 的结果重放环境混用。安装与核验命令见 `docs/MODULE4_REPRODUCIBILITY.md`。

## 能复现什么

- Module 2 冻结表示上的模型前向计算与 Top-200；
- 六路结构评分的融合与排序；
- 四条件动力学冻结摘要的导出；
- 最终表及模型、输入文件的完整性校验。

## 不能从本包推断什么

- 不能把 DrugCLIP 或 docking 分数解释为 PAM 活性；
- 不能仅凭动力学方向余弦宣称机制成立；
- 不能在缺少湿实验时宣称候选已验证；
- 不含完整 Glide、膜体系和原始轨迹，因此不声称一键重跑昂贵生产计算。

进一步参见 [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md)、[`docs/REPRODUCIBILITY_MATRIX.md`](docs/REPRODUCIBILITY_MATRIX.md)、[`docs/SUBMISSION_CLAIMS.md`](docs/SUBMISSION_CLAIMS.md) 和 [`docs/DATA_PROVENANCE.md`](docs/DATA_PROVENANCE.md)。
