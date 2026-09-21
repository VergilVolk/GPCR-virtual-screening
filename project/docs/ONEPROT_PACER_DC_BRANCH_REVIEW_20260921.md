# OneProt-MD / PACER-DC 远端分支审查

审查对象：`origin/audit/oneprot-pacer-dc-g0`，commit `2e232714`。

## 结论

该分支完成的是 **G0 权重装载审计**，不是四上下文微调，也没有产生功能预测结果。
外层 checkpoint 看起来足以替代缺失的 `forward_sim.ckpt`，但必须完成真实轨迹 forward
后，才允许进入四上下文特征提取。

## 已通过

- OneProt 仓库和 MDGen 子模块固定到指定 commit；
- checkpoint 大小、SHA256 和 state-dict 键已记录；
- 125 个 MD transformer 张量全部找到；
- `pretrained=false`、`model_path=null` 初始化后，MD transformer 报告零 missing key、
  零 shape mismatch；
- 4-byte `log_logit_scale` 差异已合理归因于未启用的对比温度缓冲区。

## 未通过或尚未完成

1. 生成 `runs/oneprot_md_branch_runtime_audit.json` 的脚本未提交，结果不可独立复现。
2. 没有真实 M4 轨迹的 forward smoke；尚未证明预处理、单位、atom14 映射、序列对齐和
   mask 均正确。
3. 没有四上下文 embedding CSV，因此没有 `d_PAM` 或 `d_AGO`。
4. 没有 PACER-DC Adapter/Triplet 的训练日志、模型权重、候选级预测或 held-out 指标。
5. 当前 18 个 MD job 只有两个独立有标签分子。trajectory window 和 replica 不能冒充独立
   分子，因此正式训练门禁不可能通过。
6. 当前训练脚本只训练冻结 embedding 上的 MLP projection/classification heads；它不是对
   OneProt-MD backbone 的参数微调，也没有实现 cluster-attention head。
7. 文档要求 bootstrap CI 和完整 baseline/ablation，但当前训练脚本尚未实现这些评估。

## 下一步唯一允许的执行顺序

### G0-forward

- 提交 runtime audit 生成脚本；
- 将一条真实 M4 蛋白轨迹转换成 OneProt-MD 训练时一致的 39-frame atom14/sequence/mask 输入；
- 明确坐标单位为 Å，核对残基顺序和缺失原子 mask；
- `model.eval()` 下对同一输入运行两次；要求输出 finite、非零、维度正确且重复误差接近零；
- 保存输入清单、checkpoint SHA256、embedding、运行日志和审计 JSON。

### G1-four-context

对同一 candidate、replica 和 window 生成：

- `z_CA = E(candidate_probe)`
- `z_C = E(candidate_no_probe)`
- `z_A = E(probe_only)`
- `z_0 = E(apo)`
- `d_PAM = z_CA - z_A`
- `d_AGO = z_C - z_0`

四个 context 必须使用相同的蛋白原子/残基选择、窗口长度和归一化规则。小分子原子本身不应
直接进入蛋白 trajectory encoder；候选物通过其诱导的蛋白动态变化产生信号。

### G2-pilot

LY2119620 和 compound-110 只能用于检查：

- replica 内重复性；
- context 差分是否非零；
- d_PAM 与 d_AGO 是否表现出预期的定性分离趋势；
- OneProt-MD 是否优于 PCA/tICA/VAMP 的同协议表示。

该阶段禁止报告 AUC、SOTA 或“发现 PAM”。

### G3-formal training

只有补齐多个独立 PAM、ago-PAM/allosteric agonist、inactive/binder controls，并按分子和
chemotype 隔离 train/validation/test 后，才运行 `train_dual_context_heads.py`。正式比较必须
包含 2D、静态 docking、ensemble docking、PCA/tICA/VAMP、无 triplet、单差分和完整双差分。

## 分支合并建议

- 可保留惰性 ligand loader 和 Pareto copy 修复；
- G0 两份 JSON 在补交生成脚本及 forward smoke 前，只能作为暂定审计证据；
- `cuda-version=13.2` 是特定主机兼容修复，不宜作为所有机器共享环境的无条件固定项；
- OneProt G0、OpenMM 环境修复和膜体系构建应拆成独立提交，避免算法审计与基础设施混合。
