# DrugCLIP 非对称适配器 CPU 复现说明

## 环境

Windows/Linux 均可，GPU 非必需。建议 Python 3.10–3.11。

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install rdkit pandas numpy scipy scikit-learn
```

DrugCLIP embedding 已预先保存在 `project/results/*/embeddings.npz`。重跑本页实验不需要重新执行 Uni-Mol，也不需要下载大型检索库。

## 1. 冻结选择并评价 M4

```bash
python project/scripts/select_molecular_adapter_on_non_m4.py \
  --train-representations project/results/muscarinic_triplet_expansion_v01/embeddings.npz \
  --train-triplets project/data/benchmarks/muscarinic_triplet_expansion_v01/train_triplets.csv \
  --strict-test project/data/benchmarks/muscarinic_triplet_expansion_v01/strict_external_test.csv \
  --zero-representations project/results/muscarinic_m4_zero_shot_v01/embeddings.npz \
  --zero-test project/data/benchmarks/muscarinic_m4_zero_shot_v01/m4_zero_shot_pairs.csv \
  --output project/results/muscarinic_triplet_expansion_v01/molecular_adapter_selection.json
```

预期选择 `ecfp3_bit`。三种子权重自动保存为：

```text
molecular_adapter_selection.ecfp3_bit.seed20260924.pt
molecular_adapter_selection.ecfp3_bit.seed20260925.pt
molecular_adapter_selection.ecfp3_bit.seed20260926.pt
```

部署规则为三种子分数差的无权平均。

## 2. 五亚型严格 target-LOSO

```bash
python project/scripts/evaluate_muscarinic_ecfp_target_loso.py \
  --representations project/results/muscarinic_target_loso_v01/embeddings.npz \
  --benchmark-dir project/data/benchmarks/muscarinic_target_loso_v01 \
  --output project/results/muscarinic_target_loso_v01/ecfp3_pocket_loso_result.json \
  --epochs 100 --fingerprint-radius 3 --skip-random-control
```

预期宏准确率约 `0.730`。该结果是开发 benchmark；冻结 M4 面板仍是更重要的外推检查。

## 3. 汇总配对不确定性

```bash
python project/scripts/summarize_drugclip_adapter_evidence.py \
  --drugclip-loso project/results/muscarinic_target_loso_v01/target_loso_result.predictions.csv \
  --ecfp2-loso project/results/muscarinic_target_loso_v01/ecfp_pocket_loso_result.predictions.csv \
  --ecfp3-loso project/results/muscarinic_target_loso_v01/ecfp3_pocket_loso_result.predictions.csv \
  --concat-loso project/results/muscarinic_target_loso_v01/concat_adapter_loso_result.predictions.csv \
  --drugclip-m4 project/results/muscarinic_m4_zero_shot_v01/zero_shot_result.predictions.csv \
  --ecfp2-m4 project/results/muscarinic_triplet_expansion_v01/ecfp_to_frozen_pocket_baseline.m4_predictions.csv \
  --selected-m4 project/results/muscarinic_triplet_expansion_v01/molecular_adapter_selection.m4_predictions.csv \
  --output project/results/muscarinic_target_loso_v01/molecular_adapter_evidence_summary.json
```

## 结果解释

- 输出分数只用于毒蕈碱受体亚型活性排序。
- 不得将其解释为 M4 PAM 概率、PAM EC50、协同性或内在激动预测。
- 禁止继续使用冻结 M4 面板调超参数；后续修改必须使用新验证集。
- 跨家族 GPCR 适配实验为负结果，详见 `DRUGCLIP_ASYMMETRIC_ADAPTER_RESULT.md`。

