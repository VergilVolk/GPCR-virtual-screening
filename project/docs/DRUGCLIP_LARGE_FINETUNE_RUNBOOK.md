# DrugCLIP 大样本微调复现入口

## 方法边界

该流程冻结两个 Uni-Mol encoder，直接更新 DrugCLIP 内部 `mol_project` 与 `pocket_project` 的 rank-8 LoRA。它是 DrugCLIP 参数高效微调，不是全 encoder 微调，也不是外接 ECFP 分类器。

## 环境

```bash
micromamba create -f project/environment_drugclip_cpu.yml
micromamba activate drugclip-cpu
```

官方权重放在：

```text
project/tools/DrugCLIP/artifacts/checkpoint_best.pt
```

若使用 PyTorch 2.6 及以上，只有在核验该 checkpoint 来源后，给 embedding 命令增加 `--trusted-checkpoint`。

## 主训练与五折 OOF

```bash
python project/scripts/finetune_drugclip_large_multitask.py \
  --screen-representations project/results/gpcr_drugclip_screening_v01/ensemble_embeddings.npz \
  --projection project/results/gpcr_drugclip_screening_v01/ensemble_embeddings.projection.pt \
  --pairs project/data/benchmarks/gpcr_drugclip_screening_v01/screening_pairs.csv \
  --triplet-representations project/results/muscarinic_triplet_expansion_v01/embeddings.npz \
  --triplets project/data/benchmarks/muscarinic_triplet_expansion_v01/train_triplets.csv \
  --output project/results/gpcr_drugclip_screening_v01/large_multitask_controlled_result.json \
  --epochs 60
```

该命令同时运行：BCE-only、BCE+retrieval、完整多任务、随机 triplet 和全随机标签对照，并保存 `bce_retrieval` 与 `full_multitask` 两套可部署投影权重。

## 严格 target-LOSO

```bash
python project/scripts/evaluate_drugclip_large_target_loso.py \
  --representations project/results/gpcr_drugclip_screening_v01/ensemble_embeddings.npz \
  --projection project/results/gpcr_drugclip_screening_v01/ensemble_embeddings.projection.pt \
  --pairs project/data/benchmarks/gpcr_drugclip_screening_v01/screening_pairs.csv \
  --output project/results/gpcr_drugclip_screening_v01/large_target_loso_3seed_result.json \
  --epochs 60 --bootstrap 500
```

验证要求：每个留出靶点的 `held_target_in_loss=false`、`scaffold_overlap=0`；主张必须使用宏平均、逐靶点结果和 scaffold-bootstrap CI。

## 当前推荐权重

```text
project/results/gpcr_drugclip_screening_v01/
  large_multitask_controlled_result.bce_retrieval.projection.pt
```

原因：triplet 分支尚未通过随机 triplet 消融。该权重用于 GPCR active/decoy 检索，不用于 PAM 效力预测。

## 13 靶点扩展验证

13-target LOSO 使用：

```bash
python project/scripts/evaluate_drugclip_13target_loso.py \
  --gpcr-representations project/results/gpcr_drugclip_screening_v01/ensemble_embeddings.npz \
  --gpcr-pairs project/data/benchmarks/gpcr_drugclip_screening_v01/screening_pairs.csv \
  --external-representations project/results/litpcba_drugclip_external_v01/embeddings.npz \
  --external-pairs project/data/benchmarks/litpcba_drugclip_external_v01/pairs.csv \
  --projection project/results/gpcr_drugclip_screening_v01/ensemble_embeddings.projection.pt \
  --output project/results/litpcba_drugclip_external_v01/combined_13target_loso_seed20260926.json \
  --epochs 40 --seed 20260926
```

正式结果需运行三个 seed，并用 `summarize_drugclip_13target_multiseed.py` 合并；二维基线由 `evaluate_drugclip_13target_ecfp_loso.py` 在完全相同切分下计算。

## 适配器与损失消融

```bash
python project/scripts/evaluate_drugclip_13target_ablation.py \
  --gpcr-representations project/results/gpcr_drugclip_screening_v01/ensemble_embeddings.npz \
  --gpcr-pairs project/data/benchmarks/gpcr_drugclip_screening_v01/screening_pairs.csv \
  --external-representations project/results/litpcba_drugclip_external_v01/embeddings.npz \
  --external-pairs project/data/benchmarks/litpcba_drugclip_external_v01/pairs.csv \
  --projection project/results/gpcr_drugclip_screening_v01/ensemble_embeddings.projection.pt \
  --output project/results/litpcba_drugclip_external_v01/combined_13target_adapter_ablation_seed20260926.json \
  --epochs 40 --seed 20260926
```

当前证据支持 `dual_no_preserve`。训练全数据部署权重时增加：

```text
--modes dual_no_preserve --full-only \
--save-full-checkpoint <output.projection.pt>
```

部署模型不能用于报告训练性能；性能数字必须来自 LOSO 预测文件。
