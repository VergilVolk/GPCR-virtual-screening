# PACER DrugCLIP 路由器：算法说明与代码接口 v01

## 1. 它到底用了哪两个模型

不是“2023 原模型和 2026 原模型直接平均”，而是两个经过我们适配的模型：

### 分支 A：2023 GPCR 专化模型

- backbone：2023 DrugCLIP checkpoint；
- 我们的训练：GPCR 整靶点留出（LOTO）迁移；
- 训练时，被测试靶点及其重叠 Murcko scaffold 不进入损失；
- 当前特点：整体 ROC 不一定最高，但 BEDROC、EF1、EF5 等头部富集较强；
- 接口名称：`drugclip2023_gpcr_loto`。

### 分支 B：2026 family-aware 模型

- backbone：Science 2026 DrugCLIP checkpoint；
- 我们的训练：13-target 微调，并加入 family augmentation；
- 当前特点：跨靶点整体排序较好，但 M4 出现负迁移，头部富集也不稳定；
- 接口名称：`drugclip2026_13t_family_aug`。

## 2. 路由算法

设两个模型对同一靶点候选库输出：

- `s23(i,t)`：2023 GPCR LOTO 分数；
- `s26(i,t)`：2026 13T family-aug 分数。

因为两个 checkpoint 的原始分数尺度不同，不能直接相加。先在每个靶点内部转换为百分位秩：

```text
r23(i,t) = percentile_rank_t(s23(i,t))
r26(i,t) = percentile_rank_t(s26(i,t))
```

然后执行冻结规则：

```text
如果 target ∈ {M4, M4R, CHRM4}：
    PACER_binding_score = r23
否则：
    PACER_binding_score = 0.5 × r23 + 0.5 × r26
```

M4 不采用 2026 分支，是因为此前独立审计已经发现该分支在 M4 上负迁移；这不是本次跑分后临时挑选的单分子规则。

## 3. 在全项目中的位置

```text
分子生成/商业库
      ↓
理化性质、PAINS、反应性过滤
      ↓
多构象 docking：排除明显不合理姿势
      ↓
2023 GPCR LOTO ─┐
                 ├─ PACER DrugCLIP M4-safe router
2026 family-aug ─┘
      ↓
结合候选排序与化学多样性选择
      ↓
PACER-FKG 四上下文动态差分
      ↓
区分 functional PAM、ago-PAM 风险和 inactive binder
      ↓
最终实验候选
```

DrugCLIP 路由器只负责回答“哪些分子更值得作为结合候选继续分析”。它不输出 PAM 标签，也不能代替 PACER-FKG 或湿实验。

## 4. 输入接口

CSV 一行对应一个“靶点—分子”配对。必需字段：

| 字段 | 类型 | 含义 |
|---|---|---|
| `pair_id` | string | 唯一配对编号，不允许重复 |
| `target` | string | 靶点名称，如 `M4R`、`B2AR` |
| `score_2023_gpcr_loto` | float | 2023 GPCR 适配模型分数，越大越好 |
| `score_2026_13t_famaug` | float | 2026 family-aug 模型分数，越大越好 |

可附加 `canonical_smiles`、`compound_id`、`docking_score` 等字段，路由器会原样保留。

示例：

```csv
pair_id,target,compound_id,score_2023_gpcr_loto,score_2026_13t_famaug
M4_0001,M4R,PACER0001,0.183,-0.047
M4_0002,M4R,PACER0002,0.121,0.206
B2_0001,B2AR,PACER0001,0.092,0.318
B2_0002,B2AR,PACER0002,-0.044,0.267
```

每个靶点至少需要两个候选分子，否则无法计算靶点内百分位秩，程序会拒绝运行。

## 5. 命令行接口

```powershell
python project/scripts/pacer_drugclip_router.py `
  --input project/results/my_screen/two_model_scores.csv `
  --output project/results/my_screen/pacer_routed_ranking.csv
```

默认同时生成：

```text
pacer_routed_ranking.audit.json
```

审计文件记录输入/输出 SHA256、行数、靶点数、路由数量、模型标识和主张边界。

如果修改融合权重：

```powershell
python project/scripts/pacer_drugclip_router.py `
  --input two_model_scores.csv `
  --output routed.csv `
  --old-weight 0.7 `
  --new-weight 0.3
```

程序仍会运行，但 `frozen_protocol=false`。这种输出属于新实验，不能与冻结 v01 的跑分混用。

## 6. Python API

```python
import sys
import pandas as pd

# 从仓库根目录运行时
sys.path.insert(0, "project/scripts")
from pacer_drugclip_router import RouterConfig, route_scores

scores = pd.read_csv("two_model_scores.csv")

# 正式冻结协议
ranking = route_scores(scores)

# 只允许用于新消融实验
custom = route_scores(
    scores,
    RouterConfig(old_weight=0.7, new_weight=0.3),
)
```

正式函数签名：

```python
route_scores(
    frame: pandas.DataFrame,
    config: RouterConfig | None = None,
) -> pandas.DataFrame
```

## 7. 输出接口

在保留输入字段的基础上增加：

| 字段 | 含义 |
|---|---|
| `rankpct_2023_gpcr_loto` | 2023 分支的靶点内百分位秩 |
| `rankpct_2026_13t_famaug` | 2026 分支的靶点内百分位秩 |
| `pacer_binding_score` | 路由后的最终结合候选分数 |
| `route` | `m4_2023_gpcr_loto` 或 `gpcr_fixed_rank_fusion` |
| `final_rank` | 每个靶点内部最终排名，1为最高 |
| `router_protocol` | 固定协议版本 |
| `frozen_protocol` | 是否严格使用冻结的 0.5/0.5 协议 |
| `eligible_for_pam_claim` | 始终为 `false` |
| `claim_boundary` | 结合候选检索，不代表 PAM 功能或效力 |

下游应优先读取：

```text
target + pair_id + pacer_binding_score + final_rank + route
```

## 8. benchmark 接口

统一评测入口：

```text
project/scripts/run_drugclip_generation_bakeoff.py
```

它接收两代 embedding、各微调模型预测、配对表、口袋群体比例和 docking baseline，输出：

- macro/per-target 指标；
- BEDROC20、EF1、EF5、PR-AUC、ROC-AUC；
- scaffold bootstrap 95% CI；
- 与 Glide/Vina 的共同样本比较；
- 图表和输入 SHA256。

该脚本用于方法评估；`pacer_drugclip_router.py` 用于真实候选库推理。两者不能混为一个脚本。

## 9. 当前性能与部署结论

在 4-GPCR 回顾性开发集上：

| 方法 | ROC | PR | BEDROC20 | EF1 | EF5 |
|---|---:|---:|---:|---:|---:|
| 2023 GPCR LOTO | 0.655 | 0.218 | 0.313 | **4.476** | 3.568 |
| 2026 family-aug | 0.718 | 0.215 | 0.246 | 1.733 | 2.420 |
| M4-safe 路由融合 | **0.735** | **0.323** | **0.398** | 3.786 | **4.474** |

相对 2023 GPCR LOTO，路由融合的 ROC、PR、BEDROC20、EF5 配对 bootstrap 95% CI 均为正；EF1 仍跨零。

因此当前部署策略是：

- 主排序：M4-safe routed score；
- 最前 1%：同时保留 2023 GPCR 分数和 Glide 结果进行复核；
- 最终 PAM 判断：进入 PACER-FKG 四上下文层，不使用 DrugCLIP 分数直接下结论。

## 10. 主张边界

该接口已经可以稳定用于候选排序和全流程联调，但当前数字来自开发/选择基准，不是最终外部盲测。锁箱揭盲前不得更改模型、权重、M4 路由规则或评价指标。
