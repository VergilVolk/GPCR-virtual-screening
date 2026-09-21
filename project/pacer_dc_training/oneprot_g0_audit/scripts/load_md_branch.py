#!/usr/bin/env python
"""Step 4 (take 2): load the OneProt MD branch without importing the Lightning stack.

`src/models/__init__.py` pulls in OneProtLitModule, which requires
pytorch_lightning.  PACER-DC never touches the multimodal Lightning module, so
that import is pure overhead here.  Registering a stub ``src.models`` package
whose __path__ points at the real directory lets ``src.models.components.*``
resolve while skipping the __init__ that needs Lightning.
"""
from __future__ import annotations

import sys
import types
from pathlib import Path

CKPT = Path(
    "project/tools/oneprot-embeddings/artifacts/"
    "Pocket_Text_ST_SG_MD/epoch_012_01100-v1.ckpt"
)
REPO = Path("project/tools/oneprot-embeddings").resolve()
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "external" / "mdgen"))  # provides the `mdgen` package

# --- stub packages so __init__.py files are not executed -------------------
for name, rel in [
    ("src", "src"),
    ("src.models", "src/models"),
    ("src.models.components", "src/models/components"),
    ("src.data", "src/data"),
    ("src.data.datasets", "src/data/datasets"),
]:
    mod = types.ModuleType(name)
    mod.__path__ = [str(REPO / rel)]
    sys.modules[name] = mod
# ---------------------------------------------------------------------------

import torch  # noqa: E402


def main() -> None:
    obj = torch.load(CKPT, map_location="cpu", weights_only=False)
    sd = obj["state_dict"]
    md = {k: v for k, v in sd.items() if k.startswith("network.md.")}
    print(f"checkpoint keys : {len(sd):,} total, {len(md)} under network.md.")

    print()
    print("=== construct TrajectoryEncoder(pretrained=False) ===", flush=True)
    print("   forcing the faked cluster path does not matter; the module's own")
    print("   sys.path.append('/p/project1/...') is a harmless no-op here.")
    from src.models.components.md_encoder import TrajectoryEncoder

    encoder = TrajectoryEncoder(
        output_dim=1024,
        hidden_size=21,
        num_layers=4,
        num_heads=8,
        pretrained=False,
        frozen=True,
        proj_type="mlp",
        num_frames=100,
        suffix="_i1",
    )
    print("   constructed OK WITHOUT forward_sim.ckpt")

    own = encoder.state_dict()
    print(f"   encoder parameters: {len(own):,}")

    print()
    print("=== load OneProt md weights, strict=False ===", flush=True)
    remapped = {k[len("network.md."):]: v for k, v in md.items()}
    result = encoder.load_state_dict(remapped, strict=False)

    print(f"   missing_keys    : {len(result.missing_keys)}")
    for k in sorted(result.missing_keys):
        print(f"      MISSING    {k}")
    print(f"   unexpected_keys : {len(result.unexpected_keys)}")
    for k in sorted(result.unexpected_keys):
        print(f"      UNEXPECTED {k}")

    covered = len(own) - len(result.missing_keys)
    print()
    print(f"   coverage: {covered}/{len(own)} of the encoder's tensors were supplied")
    print()
    print("STEP4_DONE")


if __name__ == "__main__":
    main()
