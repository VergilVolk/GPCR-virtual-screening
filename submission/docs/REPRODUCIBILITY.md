# Reproducibility

## 环境

- Python 3.11
- 普通 CPU；冻结重放不需要 GPU、Glide 或分子动力学环境
- 依赖版本见 `requirements.txt`

## 从零核验

在本目录运行：

```bash
python -m pip install -r requirements.txt
python predict.py
python verify_submission.py
```

`predict.py` 依次校验冻结输入、重新计算六路融合、导出四条件记录并写入 `results/results.csv`。`verify_submission.py` 额外完成：

1. 对 2,605 个候选使用三份已打包投影权重做真实 Module 2 前向计算；
2. 核对模型 SHA-256、分数误差和 Top-200 次序；
3. 以规范分子图核对模块连接，不盲信展示编号；
4. 对最终八候选表逐列比较。

预期末行：

```text
PASS: frozen four-module replay, model hashes and final table
```

## 单模块入口

```bash
python module1_generate_filter.py --input data/frozen/generated_candidates_2605.csv
python screen.py --molecules MOLECULES.npz --pockets M4_POCKETS.npz
python module3_ensemble_dock.py --input data/frozen/six_channel_top5_plus_anchors.csv
python module4_four_context.py --summary data/frozen/four_context_results.json
```

Module 2 单独运行需要相应表示文件；完整训练还需要未随包分发的上游大模型和训练表示。冻结核验无需这些外部大文件。

## 完整生产重跑边界

完整重跑还需要已准备的 M4 结构、AutoDock Vina、合法授权的 Schrödinger Glide、膜体系和原始 MD 轨迹。仓库保留生产脚本与允许再分发的派生结果，但不以替代引擎伪造缺失步骤。

## 历史记录限制

已保存模型和输入哈希、随机种子、验证报告与冻结输出。部分早期实验没有保存原始控制台日志、硬件和墙钟时间；这些缺口作为限制披露，不事后编造。
