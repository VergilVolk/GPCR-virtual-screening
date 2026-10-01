# DrugCLIP 代码、权重与结果账本

审计日期：2026-09-26。此文件冻结本项目当前实际使用的 DrugCLIP 工件身份；未完成的来源验证不得被写成已确认事实。

## 结论先行

1. **DrugCLIP 方法并非“2026 才出现”**：最初版本为 NeurIPS 2023；同一作者团队的扩展工作 *Deep contrastive learning enables genome-wide virtual screening* 于 2026-01-08 发表在 Science。
2. 本项目当前的代码来自官方 GitHub 的**当前 HEAD**，但当前推理 `.pt` 的内部训练元数据指向 **2023-05-06**；因此不能把本项目结果称为“对 2026 Science DrugCLIP 最新权重的微调/超越”。
3. 当前 `.pt` **不是裸 Uni-Mol 预训练权重**。它包含完整的 molecule encoder、pocket encoder、128维投影层和 temperature，并通过 `in_batch_softmax` 任务训练；它是 DrugCLIP 的任务训练 checkpoint。
4. 当前 `.pt` 的文件路径来自官方 Drive 说明的手工下载流程，但仓库中没有官方 checksum、下载 URL 记录或签名。其“官方文件一致性”状态为 **待独立验证**，不能只凭文件名或代码仓库推断。

## 工件清单

| 项目 | 已验证事实 | 状态 |
|---|---|---|
| 官方代码远端 | `https://github.com/bowen-gao/DrugClip.git` | 已验证 |
| 本地代码 commit | `7a3a3fa33673f8668c811790f2e4681c98af44ef`（2026-01-26，`save results`） | 已验证；与远端 HEAD 一致 |
| 本地权重 | `project/tools/DrugCLIP/artifacts/checkpoint_best.pt` | 可读取 |
| 权重 SHA256 | `dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e` | 已计算 |
| 权重大小 | 1,183,713,459 bytes | 已计算 |
| checkpoint 内训练目录 | `/data/protein/local_data/drug_clip_pdb_general/train_no_test_af/` | 已读取 |
| checkpoint 内训练时间 | TensorBoard 路径含 `affinity_2023-05-06_22-08-56` | 已读取 |
| 训练状态 | epoch 43，14,449 updates，`in_batch_softmax`，最佳验证指标 `valid_bedroc` | 已读取 |
| 模型初始权重 | `mol_pre_no_h_220816.pt`、`pocket_pre_220816.pt` | 已读取：说明先有 Uni-Mol 预训练，再有 DrugCLIP 任务训练 |
| 与官方 Drive 二进制一致 | 无官方 hash/下载审计可比对 | **未验证** |
| 是否为 Science 2026 配套最新权重 | 无明确版本映射或权重 hash | **未验证，当前不得假定为是** |

## 当前 checkpoint 的内部结构证据

checkpoint 顶层键为：`args`、`model`、`loss`、`optimizer_history`、`task_state`、`extra_state`、`last_optimizer_state`。`model` 含 419 项权重，包括 `mol_model.*`、`pocket_model.*`、投影层及 `logit_scale`；这排除了“只有通用 embedding、尚未训练 DrugCLIP”的解释。

checkpoint 的 `args` 明确记录：

- `task=binding_affinity`，`arch=binding_affinity`（本项目 CPU loader 将 registry 名映射为可执行的 `drugclip`，不改权重张量）；
- `loss=in_batch_softmax`；
- `finetune_mol_model=mol_pre_no_h_220816.pt`；
- `finetune_pocket_model=pocket_pre_220816.pt`；
- `best_checkpoint_metric=valid_bedroc`；
- 4 GPU/NCCL 训练元数据。

因此，所有目前所谓“官方 DrugCLIP baseline”必须准确表述为：

> 使用本地保存、内部元数据为 2023 DrugCLIP contrastive task checkpoint 的冻结双塔 baseline；其二进制与官方 Drive 文件的一致性尚待以官方 hash 或重新下载比对确认。

## 我们的“提升”到底相对什么

| 实验线 | 基线身份 | 是否可称为超过 2026 Science 最新 DrugCLIP |
|---|---|---|
| 13-target LoRA/CGDA | 上述 2023 task-trained checkpoint 的冻结输出 | 否 |
| Full LIT-PCBA CGDA | 同上，reference retrieval 的上下文适配 | 否 |
| DUD-E 外部 GPCR | 同上 | 否 |
| M4 PAM/function probe | 同上 | 否 |

也就是说，现有结果仅能说“在固定的本地 DrugCLIP task checkpoint 上，某些严格协议下的适配改善/不改善”。它们**不能**作为“我们超过 DrugCLIP 2026 SOTA”的证据。

## 必须补齐的验证闸门

1. 从 DrugCLIP 官方 Drive 重新下载权重到临时路径，记录下载 URL、大小、SHA256；与上述 SHA256 比较；
2. 确认 2026 Science 工作是否发布了不同版本的代码或 checkpoint，并记录其唯一文件名/hash；
3. 若有新权重：在不改变 LIT-PCBA/DUD-E 输入与指标的条件下重跑冻结 baseline，再重跑我们的适配；
4. 若无新公开权重：明确声明我们复现的是官方公开 2023 checkpoint，2026 Science 论文仅为同方法的后续发表和应用扩展，不能拿其标题替代权重版本。

在闸门 1--3 前，停止使用“最新 DrugCLIP”“超过 DrugCLIP SOTA”“2026 DrugCLIP 微调”等表述。

## Science 2026 工件与复现边界（2026-09-26 新核验）

Science 论文配套的 `THU-ATOM/Drug-The-Whole-Genome` 仓库指向独立的 Hugging Face 数据集 `bgao95/DrugCLIP_data`。其公开文件清单中至少包括：

| 工件 | 大小 | 用途 | 与本项目旧 checkpoint 的关系 |
|---|---:|---|---|
| `benchmark_weights/litpcba_identity_90.pt` | 1,183,713,651 bytes | LIT-PCBA 90% identity benchmark 单模型权重 | 独立文件，下载中、尚未比较 |
| `model_weights.zip` | 13,007,383,205 bytes | 6-fold/8-fold screening 权重 | 与旧单 checkpoint 不是同一发布组织方式 |
| `encoded_mol_embs.zip` | 5,911,416,616 bytes | 对应的预编码分子库 | 不能与旧模型 embedding 混用 |

论文/配套仓库描述的升级包括 ProFSA 口袋预训练、BioLip2 微调、配体构象增强以及多 fold 集成。这表明它不应被视为“旧 checkpoint 换一个论文标题”；但核心仍是 pocket--molecule 对比双塔，而非完全换成不可衔接的新范式。

### 当前公开脚本的一个可复现性缺口

`Drug-The-Whole-Genome/test.sh` 传入 `--use-folds True`，但 `unimol/tasks/drugclip.py` 的 `test_pcba()` 函数开头将 `use_folds` 无条件重设为 `False`。因此，当前 checkout 按原脚本运行时**不会**执行其中定义的 6-fold LIT-PCBA 分支。它不能直接被用来复现论文声称的 6-fold 数字。

后续比较分两层进行，并分别报告：

1. 下载并校验 `litpcba_identity_90.pt`，在完全固定的 LIT-PCBA 输入/指标上运行其单模型基线；
2. 下载并校验完整 fold 权重，以明确记录的最小补丁恢复真正的 6-fold 分支，复现集成基线；
3. 仅在每一层固定的 Science baseline 上，重训/重评 PACER 适配器。绝不把“旧单模型上的改进”外推成“超过 Science 2026”。

下载审计：2026-09-26 已从官方 Hugging Face resolve URL 启动 `litpcba_identity_90.pt` 下载到 `project/data/external/drugclip_science2026/`，不覆盖旧权重。官方 LFS 对象 SHA256 应为 `997d343b0ba4b9436a30d76b7b3a0281f99d7e7f78f2e1fbdef5dc3bb0c67a31`；下载完成前不得读取或执行该不完整文件。
