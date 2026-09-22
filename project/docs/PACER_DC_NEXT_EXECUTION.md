# PACER-DC 下一步执行顺序（2026-09-22冻结）

## 当前证据

当前只完成了 `compound110 / replica_01 / window_000` 的四上下文技术测试。每个上下文
100帧、间隔10 ps，即约1 ns；已得到四个1024维 OneProt-MD embedding，以及数学上正确的
`dPAM`和`dAGO`。这证明提取链路可运行，不证明轨迹收敛、药理功能或模型性能。

对外准确表述：

> compound-110 replica 1 的首个约1 ns四上下文窗口已完成冻结表征提取。

不得表述为“compound-110完整MD已完成”或“模型已识别其激动/PAM功能”。

## 立即执行：先完成可重复性门禁

1. 用修复后的提取器重新生成当前unit；确认四个context分别保留独立拓扑文件。
2. 验收 checkpoint：`missing_keys=[]`；只允许已知的
   `unexpected_keys=[norm.1.log_logit_scale]`。
3. 完成 compound-110 的 replica 2、3，并让 apo、probe-only 使用相同replica编号和时间窗。
4. 每个体系先推进到5 ns技术pilot，形成5个不重叠的1 ns窗口；不要把窗口当作独立分子。
5. 运行稳定性审计，报告同一候选物的跨replica与跨window `dPAM/dAGO`余弦和L2距离。

若差分方向在replica之间不稳定，停止延长轨迹并检查初始构象、周期边界、配体稳定性和编码器
敏感性；不能用挑选“最好窗口”的方式补救。

## 第二阶段：形成最小功能对照矩阵

在同一协议下至少补齐：

- functional PAM：LY2119620、CM00717；
- allosteric agonist / ago-PAM对照：compound-110；
- experimental inactive / non-PAM：CM00734；
- shared controls：ACh-only、apo。

所有候选必须使用3个独立replica、匹配时间窗和相同预处理。只有当5 ns pilot通过数值稳定、
配体未明显逸出口袋、跨replica差分具有可重复方向后，才把通过的体系延长到100 ns；不要在
QC失败前盲目消耗100 ns算力。

## 训练门禁

在下列条件满足前，`TRAINING_GATE`保持关闭：

- 至少两个PAM正例、两个PAM负例和两个独立激动性阳性分子；
- 具有可审计的chemotype字段；
- candidate和chemotype均不跨train/validation/test；
- 每个分子至少有完整四上下文和多个独立replica；
- 统计单位是分子，window只能作为数据增强；
- baseline使用相同拆分，至少包含静态/ensemble docking、DrugCLIP、OneProt-MD无triplet、
  融合无triplet和完整模型。

当前仅一个compound-110 unit，不得训练分类器、计算AUC或报告SOTA。

## 必须随结果提交的文件

每次推送结果时同时提交小型审计文件：

- `progress.json`：system、replica、target_ns、completed_ns、status；
- `unit_summary.json`；
- `four_context_long.csv`与`four_context_differential.csv`；
- 稳定性审计JSON；
- checkpoint、拓扑和轨迹的来源及哈希。

大型DCD、checkpoint和第三方仓库不进入Git。只有审计文件显示完成，才允许在汇报中称为完成。
