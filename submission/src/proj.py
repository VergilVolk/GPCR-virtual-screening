"""DrugCLIP 双侧投影头（LoRA 微调后物化的线性适配层）。

结构（逐字核对自主仓库 finetune_drugclip_muscarinic_triplet.Proj）：
    linear1 : 512 -> hidden
    linear2 : hidden -> 256
    forward : L2 归一化的嵌入（余弦空间）

权重文件（.projection.pt）为 dict：
    {'mol_project': state_dict, 'pocket_project': state_dict}
打分语义：score(mol, pocket) = <mol_embed, pocket_embed>（余弦）。
"""
from __future__ import annotations
import torch
from torch import nn
from torch.nn import functional as F


class Proj(nn.Module):
    def __init__(self, state: dict):
        super().__init__()
        self.linear1 = nn.Linear(state['linear1.weight'].shape[1], state['linear1.weight'].shape[0])
        self.linear2 = nn.Linear(state['linear2.weight'].shape[1], state['linear2.weight'].shape[0])
        self.load_state_dict({k: torch.as_tensor(v) for k, v in state.items()})

    def forward(self, x):
        return F.normalize(self.linear2(F.relu(self.linear1(x))), dim=-1)
