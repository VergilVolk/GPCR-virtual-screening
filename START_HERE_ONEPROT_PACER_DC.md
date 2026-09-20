# 当前队员任务：OneProt-MD → PACER-DC Adapter

> 本文件是当前协作任务的唯一入口。它不是整个仓库的项目综述。

## 一句话目标

复现 **OneProt-MD** 预训练轨迹编码器，用它提取 M4 四上下文轨迹 embedding，并在冻结 encoder 的前提下训练 **PACER-DC 双差分 Adapter + chemotype-aware triplet head**。

这项工作不是 PACER-FS，不是单纯 docking，也不是“把六个 OpenMM 体系构建出来”即完成。

## 方法定义

对每个候选分子构造四种匹配上下文：

1. `candidate_probe`：M4 + candidate + ACh；
2. `candidate_no_probe`：M4 + candidate；
3. `probe_only`：M4 + ACh；
4. `apo`：apo M4。

由同一个冻结 OneProt-MD encoder 得到轨迹窗口表征：

```text
z_CP = E(T_candidate_probe)       z_C = E(T_candidate_no_probe)
z_P  = E(T_probe_only)            z_0 = E(T_apo)

d_PAM = z_CP - z_P
d_AGO = z_C  - z_0
```

- `d_PAM` 学习 ACh 存在时的额外协同变化；
- `d_AGO` 学习无 ACh 时的内在激动风险；
- adapter/triplet 只在训练集内部构造，外层按 molecule/chemotype 留出；
- trajectory window 是数据增强，不是独立分子样本。

## 当前交付物

### G0：OneProt-MD 技术复现

- 使用官方仓库 `oneprot-models/oneprot-embeddings`；
- 固定版本：`53fa9c0f9150aa44a9258db6948277a4602036ab`；
- 初始化 `external/mdgen` 子模块；
- 获取官方公开 checkpoint，记录来源、文件名、SHA256；
- 删除/覆盖源码中的作者集群硬编码路径，不得伪造 checkpoint；
- 在一条公开 M4 trajectory 上输出有限、非全零且可重复的 embedding。

G0 输出：`oneprot_environment_audit.json`、`checkpoint_sha256.txt`、一份 embedding smoke sample 和运行日志。

若官方 OneProt checkpoint 下载后仍缺少源码默认的 `forward_sim.ckpt`，先运行：

```bash
python project/pacer_dc_training/inspect_oneprot_md_checkpoint.py \
  /path/to/epoch_012_01100-v1.ckpt \
  --output runs/oneprot_checkpoint_audit.json
```

- 若不存在 `network.md.transformer.*` 权重：立即停止，报告硬阻塞；
- 若存在：只能说明外层 checkpoint **可能**已经内嵌 MDGen 权重；
- 随后用 `pretrained=false`、`model_path=null` 实例化，再加载外层 checkpoint；
- 必须检查 `load_state_dict` 返回值，要求 `network.md.transformer.*` 下零 missing key、零 shape mismatch；
- 不得仅使用 `strict=false` 后看到程序不报错，就认定权重加载成功。

### G1：四上下文轨迹输入

- 六个参考体系全部构建并完成数值 QC；
- 缺 compound-110 位姿时的惰性加载修复可以保留，但 `4/6` 不是任务完成；
- 每个用于功能比较的候选必须有四上下文和匹配 replica；
- 短平衡只标记为 `equilibration_pilot`，不得当作正式功能训练数据。

### G2：OneProt-MD embedding 导出

输出长表 CSV，每行一个 candidate/replica/window/context，至少包含：

```text
candidate_id, chemotype, replicate_id, window_id, context,
pam_label, agonism_label, split, f_000 ... f_N
```

四上下文必须按照相同 replica 与 window 对齐。格式示例：
`project/pacer_dc_training/example_four_context_features.csv`。

### G3：PACER-DC Adapter 训练

```bash
conda env create -f project/environment_pacer_dc_train_cpu.yml
conda activate pacer-dc-train-cpu

python project/pacer_dc_training/audit_oneprot_pacer_dc_delivery.py \
  --oneprot-root /path/to/oneprot-embeddings \
  --checkpoint /path/to/official_checkpoint.ckpt \
  --features /path/to/four_context_features.csv \
  --output runs/oneprot_delivery_audit.json

python project/pacer_dc_training/train_dual_context_heads.py \
  --features /path/to/four_context_features.csv \
  --config project/config/pacer_dc_training_v01.json \
  --output runs/pacer_dc_oneprot_adapter_v01
```

第一阶段冻结 OneProt-MD encoder，只训练 projection/functional heads。完整解冻或 LoRA 只能作为后续消融，不能在小样本上直接宣称更优。

### G4：必须完成的 baseline

使用完全相同的 molecule/chemotype-held-out split 比较：

1. 手工动态端点；
2. PCA/tICA；
3. VAMP；
4. 冻结 OneProt-MD embedding + 线性头；
5. OneProt-MD + Adapter，但无 triplet；
6. OneProt-MD + random triplet；
7. OneProt-MD + functional hard triplet。

只有 OneProt-MD Adapter 在留出 chemotype 上稳定超过简单动态 baseline，且置信区间支持增益，才算方法学结果。

## 以下情况不算完成

- 只完成蛋白准备、重对接、膜体系或 OpenMM QC；
- 只报告原子数、碰撞距离、Vina score 或某个 seed；
- 只在同一个 PAM 的不同 trajectory window 上训练和测试；
- 把 window 数量当作独立样本量；
- 没有 OneProt-MD embedding；
- 没有冻结 baseline 和 chemotype-held-out 结果；
- 把计算优先候选称为实验确认的 PAM。

## 当前真实状态

- 上游六体系构建和短平衡已有基础；
- OneProt-MD 官方代码已审计，但尚未完成可复现推理；
- PACER-DC 双头/triplet 训练入口已经提交；
- 正式四上下文生产轨迹与 OneProt-MD embedding 尚未产生；
- 因此当前任务仍处于 G0/G1，不是训练完成状态。
