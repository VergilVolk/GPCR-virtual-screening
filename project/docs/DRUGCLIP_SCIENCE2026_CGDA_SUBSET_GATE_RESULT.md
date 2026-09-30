# CGDA 双 checkpoint 对照：分层子集临时门（v01）

日期：2026-09-30。执行链：`cgda_subset_gate.ps1`（提交 0a08f224）。

## 协议

- 15 靶 = 8 个已完成全量靶 + 7 个分层子集靶
  （ADRB2/FEN1/GBA/IDH1/KAT2A/OPRK1/VDR：全部活性物 + 每靶 2 万固定种子诱饵，
  行索引由旧 checkpoint 标签确定性导出，新旧两侧完全同行——行对齐已在
  ALDH1 上逐行验证）；
- CGDA 旗舰配方（rank2/experts2/epochs5/hardrank）× 3 种子 × 两侧；
  配对靶点 bootstrap 95% CI（5000–20000 次）；
- 全量 15 靶正式门继续后台运行，本结果标注为 interim。

## 结果（3 种子均值，15 靶宏平均）

| | pocket | reference | CGDA |
|---|---:|---:|---:|
| **旧 2023 ckpt** | ROC 0.5674 / BEDROC 0.0658 | 0.5697 / 0.0793 | **0.5815 / 0.0821** |
| **新 Science-2026** | 0.5932 / 0.0793 | **0.5941 / 0.0889** | 0.5975 / 0.0866 |

CGDA − reference 的 CI：

| 指标 | 旧 ckpt | 新 ckpt |
|---|---|---|
| ROC-AUC | **[+0.0090, +0.0149] 正** | **[+0.0013, +0.0067] 正**（幅度缩为 1/3） |
| BEDROC80.5 | **[+0.0011, +0.0039] 正** | **[−0.0051, −0.0007] 负（翻转）** |
| EF0.5% | n.s. | **[+0.029, +0.415] 正** |
| EF5% | n.s. | [−0.688, −0.059] 负 |

## 三条结论

1. **Science-2026 checkpoint 本身是更强的 LIT-PCBA 基线**：同一批行上
   reference 检索 ROC 0.5697→0.5941（+0.024）、BEDROC 0.0793→0.0889。
   若目标是该 benchmark 家族的筛选性能，checkpoint 迁移在基线层面直接获益。
2. **CGDA 的相对增益随基线变强而收缩**：旧权重上 ROC+BEDROC 双正的格局，
   在新权重上变为 ROC 小幅为正（+0.004）+ BEDROC 翻负 + EF5% 翻负。
   "强基线吃掉适配器余量"的经典形态。
3. **配置建议（interim 证据）**：LIT-PCBA 型任务直接用 Science-2026
   reference/CGDA；CGDA 的方法学主张仍以旧权重全量 15 靶的已冻结结果
   （0.5693→0.5788，CI 为正）为准，新权重上的正式表述等全量门落地。

## 边界

- 子集靶的指标定义在全量标签比例下（EF/BEDROC 对子集密度敏感），
  两侧同行对照有效，绝对值不可与全量表直接比较；
- 7 小时临时门，全量 15 靶（含 ADRB2 补录）继续运行中；
- 不改变结合/功能模块边界。

## 复现

- 运行器：`cgda_subset_gate.ps1`；子集编码 `run_drugclip_science2026_subset.py`
- 输出：`results/drugclip_science2026/cgda_subset_gate_v01/`
  （双侧各 3 个 seed JSON + 双侧 3seed summary + chain.log）
