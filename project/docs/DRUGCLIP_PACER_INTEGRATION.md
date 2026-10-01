# DrugCLIP × PACER-DC 接入方案

## 方法定位

DrugCLIP 将口袋和分子编码到共享的 128 维空间，用余弦相似度完成快速结合检索。它比单一
docking score 更适合作为大规模静态结合筛选器，但仍不能直接判断一个别构配体是否增强 ACh
介导的 M4 功能响应。

因此本项目不把 DrugCLIP 当作 PAM 分类器，而把它接入 PACER-DC 的 Binding 层：

```text
DrugCLIP pocket encoder:  receptor snapshot -> p_context
DrugCLIP molecule encoder: candidate conformer -> m
OneProt-MD encoder: matched trajectory window -> z_context

F_context = [z_context, p_context, p_context * m, cosine(p_context, m)]
d_PAM = F_candidate+probe - F_probe-only
d_AGO = F_candidate-only - F_apo
```

这里的逐元素交互项使分子身份不会在上下文相减时完全抵消；DrugCLIP cosine 保留可解释的
静态结合轴；OneProt-MD 保留轨迹动态轴。

## Triplet 微调

不建议用当前小样本更新完整 DrugCLIP/Uni-Mol backbone。第一阶段冻结两个 DrugCLIP encoder，
只训练 PACER 的 projection heads：

- anchor/positive：不同 chemotype、但均为 functional PAM 的候选级动态差分；
- hard negative：DrugCLIP 结合分数高、但实验上为 non-PAM/inactive 的分子；
- AGO 分支负担独立激动风险，不能与 PAM 分支混为一个标签。

这一困难负样本定义直接对应核心科学问题：区分“可能结合”与“具有 PAM 功能”。如果没有实验
确认的 binder-but-non-PAM，不能把普通 decoy 冒充功能困难负样本。

## 最小验证矩阵

所有方法必须使用相同的 molecule/chemotype-held-out split：

1. Vina/ensemble docking；
2. DrugCLIP 原始 cosine；
3. DrugCLIP + 普通分类头；
4. OneProt-MD 四上下文差分；
5. DrugCLIP 四上下文交互差分；
6. OneProt-MD + DrugCLIP，不含 triplet；
7. 完整融合模型 + chemotype-aware hard triplet。

只有第 7 项在候选级外部测试中稳定优于第 2、4、6 项，才能把融合和 triplet 作为算法贡献。

## 当前工程状态与边界

- 官方代码固定在 `bowen-gao/DrugClip@7a3a3fa3`；代码为 Apache-2.0；模型权重和输出为
  CC BY-NC 4.0。
- 官方 checkpoint 约 1.18 GB，示例 molecule LMDB 约 4.68 GB；不应提交到 Git。
- 官方实现依赖旧版 Uni-Mol/Uni-Core、RDKit 2022.9.5，并有写死的 CUDA 调用，需要单独兼容层。
- 当前新增的 `build_drugclip_context_features.py` 只负责冻结 embedding 的有机融合，不产生
  PAM 性能结论。
- DrugCLIP输入是局部口袋静态构象；四上下文中必须采用同一口袋原子定义、坐标处理和窗口抽样。

## 已完成与下一门禁

- 已冻结官方 checkpoint 大小与 SHA256，并完成真实 7TRS M4 口袋和 3 个真实分子的 CPU
  前向测试；输出为 128 维、分数有限，权重缺失键与异常键均为空。
- 尚未完成 529 分子 × 3 个 M4 状态的正式 benchmark，因此当前不能报告 DrugCLIP 的
  M4 PAM 区分性能。
- 下一步先完整运行冻结 baseline，再检查重复预处理和分子构象敏感性；只有这些门禁通过，
  才批量提取四上下文 snapshot embedding。
- `build_drugclip_context_features.py` 会在导出前强制阻断 candidate 或 chemotype 跨拆分泄漏。
