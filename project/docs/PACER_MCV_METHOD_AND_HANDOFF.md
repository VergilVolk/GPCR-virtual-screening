# PACER-MCV：无编码器四上下文物理构象基线

冻结版本：v0.1（2026-09-27）

## 为什么做

OneProt-MD 和 Geom2Vec 能生成数值表示，但目前没有在 compound-110 的独立 replica 中给出稳定的 PAM 协同方向。继续更换黑箱编码器，无法回答失败来自表示、采样还是候选物本身。

PACER-MCV 不训练 encoder，直接测量 M4 上有明确几何含义的量，用作四上下文方法的物理基线：

- PAM 接触区、正构核心和胞内微开关的 Cα 距离；
- 相同位置的侧链质心距离；
- 各功能区域的紧致程度；
- PAM 口袋—正构核心—胞内微开关之间的区域距离。

特征清单在 `project/config/pacer_m4_mechanism_cv_v01.json` 中预先冻结，不能看过 R3 结果后删改。

## 四上下文与三个输出轴

对每个候选物保留四个体系：apo（0）、ACh（A）、候选物（C）、候选物+ACh（C+A）。先分别汇总每条轨迹，再计算：

```text
条件性候选物效应  ΔC|A = f(C+A) - f(A)
内在激动响应      ΔAGO = f(C) - f(0)
交互/非加性响应   ΔINT = f(C+A) - f(A) - f(C) + f(0)
```

`ΔINT` 只表示候选物与 ACh 的非加性构象响应，不自动等于正协同性，更不自动等于 PAM 药效。

## 冻结的第一判据

主终点为逐轨迹均值特征上的 `ΔINT`。只有同时满足以下条件，才记为“可重复的物理构象信号”：

1. R2/R3 整体方向余弦不低于 0.50；
2. 单特征方向一致率不低于 0.65；
3. 以独立窗口为单位重采样，方向余弦 95% 描述区间下限大于 0。

这是机制信号门禁，不是显著性检验或 PAM 分类性能。要声称 PAM 判别，还必须加入已知 PAM、ago-PAM 与实验无功能对照，并按分子或 chemotype 留出验证。

## 队友一条命令运行

队友设备上需要保留 R2/R3 的 40 个原始 atom14 `.npy` 和对应序列 `.csv`，目录形式与 Geom2Vec 批任务一致，但输入必须是 `.npy`，不能拿 `.geom2vec.npz` 代替：

```bash
python project/pacer_dc_training/run_pacer_mechanism_cv.py \
  --input-root project/results/pacer_dc_four_context_v01/compound110 \
  --mapping project/results/pacer_dc_four_context_v01/compound110/G2_RESIDUE_MAPPING_v01.csv \
  --config project/config/pacer_m4_mechanism_cv_v01.json \
  --output project/results/pacer_mcv_compound110_R2R3_v01
```

输出：

- `PACER_MCV_AUDIT_v01.json`：主门禁、全部轴和来源哈希；
- `PACER_MCV_FEATURE_REPLICATION_v01.csv`：每个物理特征在 R2/R3 的方向与效应；
- `PACER_MCV_WINDOW_FEATURES_v01.csv`：窗口级审计数据。

## 当前已完成的验证

- 单元测试证明刚体平移不改变特征；
- 纯加性的四上下文构造严格给出 `ΔINT=0`；
- 注入非加性变化时准确恢复该变化；
- 已在仓库内 R1/W0 的四条真实 atom14 轨迹上完成 53 个特征的提取烟雾测试。

## 当前缺口

仓库只含 R1/W0 原始 atom14；R2/R3 仅提交了 Geom2Vec 输出及审计，没有原始坐标。因此本地不能诚实计算 R2/R3 的 PACER-MCV 结果。该任务不需要 GPU；有原始 atom14 的队友在 CPU 上即可运行。
