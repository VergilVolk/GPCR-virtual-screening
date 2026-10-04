from __future__ import annotations

import torch
from torch import nn
from torch.nn import functional as F


class ProjectionHead(nn.Module):
    """Two-layer projection head used by the DrugCLIP adapters."""

    def __init__(self, state: dict):
        super().__init__()
        self.linear1 = nn.Linear(
            state["linear1.weight"].shape[1], state["linear1.weight"].shape[0]
        )
        self.linear2 = nn.Linear(
            state["linear2.weight"].shape[1], state["linear2.weight"].shape[0]
        )
        self.load_state_dict({key: torch.as_tensor(value) for key, value in state.items()})

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        hidden = F.relu(self.linear1(values))
        return F.normalize(self.linear2(hidden), dim=-1)


Proj = ProjectionHead
