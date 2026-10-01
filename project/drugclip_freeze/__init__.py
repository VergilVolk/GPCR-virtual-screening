"""Frozen, apply-only DrugCLIP production path (Model A + Model B).

Hard prohibitions enforced by this package:

  * no training, no fine-tuning, no adapter fitting, no PocketRealign;
  * no parameter update of any kind;
  * no score fusion invented -- the dual-model decision is exactly
    `old_top50  INTERSECT  new_top25`.

Model A  = official frozen DrugCLIP backbone + official checkpoint, embeddings
           and cosine scores only.
Model B  = the frozen, already-trained projection heads applied to Model A's
           frozen 512-D representations.  Projection-only: the backbone is not
           re-run and no gradient is ever computed.
"""

from __future__ import annotations

VERSION = "PACER_DRUGCLIP_FREEZE_v01"

__all__ = ["projection", "model_a", "model_b", "dual_model", "apply"]
