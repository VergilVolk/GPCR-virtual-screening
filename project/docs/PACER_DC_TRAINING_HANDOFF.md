# PACER-DC 最终训练部署说明（CPU 版）

> 当前协作任务以仓库根目录 `START_HERE_ONEPROT_PACER_DC.md` 为准。
> OneProt-MD 是本轮必须复现和比较的候选预训练 encoder；只有在留出验证中超过
> tICA/VAMP 后，才会晋级为最终 backbone。

## 1. 这一步训练什么

最终阶段不重新训练 MD 或蛋白基础模型，而是在冻结的轨迹表征上训练两个轻量功能头：

- `d_PAM = E(candidate + ACh) - E(ACh only)`：学习 ACh 存在时由候选物带来的额外动态变化；
- `d_AGO = E(candidate only) - E(apo)`：学习候选物在无 ACh 时引起的内在激动风险。

两个头分别输出 `PAM prioritization probability` 与 `intrinsic agonism risk`。这些值是回顾性模型输出，不等同于实验确认的 PAM 活性。

## 2. 计算资源

- **训练双头模型不需要 GPU**，普通 8 核 CPU、16 GB 内存即可；
- MD 轨迹生成仍然昂贵。本仓库建议复用已完成或外部计算得到的轨迹；
- 双头训练可在 CPU 上进行；OneProt-MD embedding 提取通常需要官方 CUDA/Apptainer 环境；
- tICA/VAMP 是必须保留的强 baseline，不能用 OneProt-MD 的模型名称替代实际比较。

## 3. 建立环境

```bash
conda env create -f project/environment_pacer_dc_train_cpu.yml
conda activate pacer-dc-train-cpu
```

## 4. 输入文件

输入为一个 CSV。每行对应一个候选分子、独立 replica、trajectory window 和 context。必须包含：

| 字段 | 含义 |
|---|---|
| `candidate_id` | 分子唯一编号 |
| `chemotype` | 化学系列或 Murcko scaffold 分组 |
| `replicate_id` | 独立 MD replica，不能用 window 代替 |
| `window_id` | 轨迹窗口，仅用于数据增强 |
| `context` | `candidate_probe`、`candidate_no_probe`、`probe_only`、`apo` |
| `pam_label` | 功能性 PAM 标签，0/1 |
| `agonism_label` | 内在激动/ago-PAM 标签，0/1 |
| `split` | `train`、`val` 或 `test` |
| `f_000...` | 冻结编码器输出的数值特征 |

示例见 `project/pacer_dc_training/example_four_context_features.csv`。

硬性要求：

1. 同一候选的四个 context 必须使用相同的 replica 和 window 对齐；
2. 同一分子及同一 chemotype 不得跨 train/val/test；
3. window 不计作独立样本，指标必须按分子聚合；
4. 训练集每个功能类别至少包含两个独立分子。这个下限只保证程序可运行，不足以支持发表级结论；
5. blind candidates 只能放入 test，不能用于调参。

## 5. 启动训练

```bash
python project/pacer_dc_training/train_dual_context_heads.py \
  --features /path/to/four_context_features.csv \
  --config project/config/pacer_dc_training_v01.json \
  --output runs/pacer_dc_dual_head_v01
```

输出：

- `model.pt`：模型权重和特征列定义；
- `candidate_predictions.csv`：按分子聚合的两个概率；
- `training_history.csv`：训练曲线；
- `audit.json`：数据量、切分、指标和主张边界。

## 6. 正式比较协议

必须同时报告以下 baseline，并使用完全相同的 molecule/chemotype-held-out split：

1. 2D ECFP4 模型；
2. 静态单结构 docking；
3. ensemble docking 的 min/mean score；
4. 不含 triplet loss 的双差分头；
5. 只使用 `d_PAM` 或只使用 `d_AGO` 的消融；
6. 完整 PACER-DC 双头模型。

主指标为外部或 chemotype-held-out ROC-AUC、PR-AUC、balanced accuracy，并按独立分子或化学系列 bootstrap 置信区间。不得按 trajectory frame 计算置信区间。

## 7. 当前证据边界

当前仓库已经完成四上下文体系构建、短平衡和公开已知 PAM 的动态稳定化验证，但尚没有足够的多 chemotype 四上下文生产轨迹。因此：

- 本训练包是可部署的最终算法接口；
- 目前不能宣称已完成 PACER-DC 正式训练；
- 目前不能把候选物称为“已发现的 PAM”；
- 只有在独立外部系列优于 2D 与 docking baseline 后，才能主张算法具有预测增益。
