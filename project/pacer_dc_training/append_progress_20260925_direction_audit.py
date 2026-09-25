from pathlib import Path

doc = Path("project/docs/PACER_DC_PROGRESS_20260923.md")
if not doc.is_file():
    raise SystemExit(f"MISSING_DOC: {doc}")

text = doc.read_text(encoding="utf-8")
header = "## 2026-09-25 统一QC、PCA与 dPAM 方向敏感性更新"

if header in text:
    raise SystemExit("SECTION_ALREADY_EXISTS: refusing to append duplicate progress section")

section = """

## 2026-09-25 统一QC、PCA与 dPAM 方向敏感性更新

### 一、compound110 三 replica × 四 context 统一 QC（共同 1 ns 窗口）

已对 compound110 的 12 组轨迹（4 context × 3 replica）完成统一 QC 与共同降维预处理。

共同分析范围与结论如下：

- 统一时间范围：每条轨迹取前 100 帧（约 1 ns，10 ps/frame）；
- 总帧数：1200；
- 统一受体选择：`chainID E and protein and name CA`；
- 共同受体映射检查通过：270 个残基、270 个 Cα、原子顺序一致；
- 全部 1200 帧坐标有限性检查通过；
- 12 条轨迹逐帧连续性检查通过；
- R3 apo 在共同 1 ns 分析中使用修复轨迹 `trajectory_corrected.dcd`。

当前结论：compound110 的三 replica 均可在统一 1 ns 范围内进行共同 PCA 与轨迹可视化分析。

### 二、共同 PCA 结果（compound110，12 组轨迹）

对 1200 帧共同对齐后的受体 Cα 坐标进行了确定性 PCA（full SVD）。

解释方差：

- PC1 = 0.286183
- PC2 = 0.256741
- PC3 = 0.184046
- PC1 + PC2 = 0.542924
- PC1 + PC2 + PC3 = 0.726970

阶段性观察：

- 二维 PCA（PC1/PC2）中，四个 context 与三个 replica 均存在分布差异；
- 12 面板时间演化图显示，部分轨迹存在明显时间漂移；
- R3 apo 在 PC3 方向上漂移最大，但前 1 ns 内未观察到异常逐帧跳变；
- R3 的存在会显著影响二维 PCA 子空间，但前三维子空间整体仍较稳定。

留一 replica 分析（相对于完整 1200 帧 PCA）：

- 排除 R1：PC1–PC2 overlap = 0.980830；PC1–PC3 overlap = 0.986651；
- 排除 R2：PC1–PC2 overlap = 0.973557；PC1–PC3 overlap = 0.989586；
- 排除 R3：PC1–PC2 overlap = 0.842806；PC1–PC3 overlap = 0.950084。

进一步的 R3 留一轨迹分析显示：

- 排除 R3 candidate + ACh：top2 overlap = 0.996180；top3 overlap = 0.998083；
- 排除 R3 candidate only：top2 overlap = 0.935915；top3 overlap = 0.994576；
- 排除 R3 ACh-only：top2 overlap = 0.567767；top3 overlap = 0.986661；
- 排除 R3 apo：top2 overlap = 0.962165；top3 overlap = 0.988111。

组件旋转诊断表明：排除 R3 ACh-only 后，二维 PCA 的变化主要来源于前三个主成分内部的旋转，尤其表现为原 PC2/PC3 的显著重排，而不是前三维子空间整体失稳。

当前结论：二维 PCA 图对个别轨迹（尤其 R3 ACh-only）较敏感；前三维子空间整体仍保持较高一致性。二维图适合做可视化提示，但不宜单独作为构象稳定性的定论依据。

### 三、R2/R3 五窗口 dPAM 方向敏感性审计

已基于现有冻结 OneProt-MD embedding，对 compound110 的 R2/R3 五个匹配窗口（W0–W4）完成 dPAM 方向敏感性审计。

定义：

- dPAM = candidate_probe - probe_only
- 跨 replica 变化定义为：Δ = R3 - R2

对每个窗口均验证：

- dPAM 数学一致性残差 < 1e-5；
- 差分分解平方残差 < 1e-5。

五窗口主要结果如下：

| Window | CA cosine | A cosine | dPAM cosine | dPAM angle (deg) | cos(ΔzCA, ΔzA) | cross term |
| --- | --- | --- | --- | --- | --- | --- |
| W0 | 0.895527 | 0.847766 | 0.502744 | 59.818 | -0.310129 | 0.156445 |
| W1 | 0.915691 | 0.918110 | 0.712309 | 44.577 | 0.169672 | -0.056393 |
| W2 | 0.832024 | 0.868390 | -0.353644 | 110.710 | -0.203894 | 0.121264 |
| W3 | 0.886168 | 0.614202 | 0.262047 | 74.808 | 0.476550 | -0.399468 |
| W4 | 0.882035 | 0.785358 | 0.085430 | 85.099 | -0.008958 | 0.005702 |

阶段性结论：

1. 原始 embedding 的跨 replica 相似度较高，并不保证 dPAM 方向一致；
2. 五个窗口中，两个原始变化向量既可能相互增强（W0/W2）、也可能相互抵消（W1/W3），W4 近乎正交；
3. W2 的 dPAM cosine 为负（-0.353644），显示最明显的方向反转；
4. W3 虽然存在较强抵消（cross term = -0.399468），但 dPAM 方向仍不稳定；
5. 当前结果支持“差分几何关系是 dPAM 方向波动的重要因素”，但尚不能区分其根源来自 MD 采样、编码器表示还是差分操作对两者相对方向的敏感性。

### 四、当前阶段状态

- 统一 QC：已完成（compound110，三 replica，四 context，共同 1 ns）；
- 轨迹降维可视化：已完成第一轮 PCA、12 面板时间图、PC1/PC3 诊断与稳定性分析；
- dPAM 方向敏感性：已完成 R2/R3 五窗口初步审计；
- TRAINING_GATE：继续保持 `CLOSED`。

### 五、已新增的本地结果与脚本

新增脚本：

- `project/pacer_dc_training/plot_dpam_direction_summary_R2_R3.py`
- `project/pacer_dc_training/audit_dpam_direction_R2_R3.py`

新增结果：

- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_decomposition_R2_R3_v01.json`
- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_decomposition_R2_R3_v01.csv`
- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_summary_R2_R3_v01.png`
- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_summary_R2_R3_v01.svg`
- `project/results/pacer_dc_four_context_v01/compound110/direction_audit_R2_R3_v01/dpam_direction_summary_R2_R3_v01.json`
"""

doc.write_text(text.rstrip() + "\n" + section + "\n", encoding="utf-8")
print("PROGRESS_DOC_UPDATE_PASS:", doc)
