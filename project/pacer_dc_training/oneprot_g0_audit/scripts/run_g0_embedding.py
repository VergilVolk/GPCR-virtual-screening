#!/usr/bin/env python
"""G0 smoke test: run the frozen OneProt MD encoder on a public M4 trajectory.

Replicates MDDataset.collate_fn + NewMDGenWrapper.prep_batch in plain torch so
that pytorch_lightning (needed only by the multimodal wrapper) is not required.
Loads TrajectoryEncoder(pretrained=False) and injects the OneProt
``network.md.*`` weights, then runs the forward pass twice to confirm the
embedding is finite, non-zero, and bit-reproducible.
"""
from __future__ import annotations

import argparse
import sys
import types
from pathlib import Path

import numpy as np
import torch

REPO = Path("project/tools/oneprot-embeddings").resolve()
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "external" / "mdgen"))

# Skip src/*/__init__.py which import pytorch_lightning; PACER-DC only needs the
# MD encoder and none of the Lightning plumbing.
for name, rel in [
    ("src", "src"),
    ("src.models", "src/models"),
    ("src.models.components", "src/models/components"),
]:
    mod = types.ModuleType(name)
    mod.__path__ = [str(REPO / rel)]
    sys.modules[name] = mod

from mdgen import residue_constants as rc  # noqa: E402
from mdgen import geometry as geo  # noqa: E402
from mdgen.rigid_utils import Rigid, Rotation  # noqa: E402
from mdgen.utils import get_offsets  # noqa: E402
from src.models.components.md_encoder import TrajectoryEncoder  # noqa: E402

CKPT = REPO / "artifacts" / "Pocket_Text_ST_SG_MD" / "epoch_012_01100-v1.ckpt"
NPY = REPO / "artifacts" / "m4_atom14" / "m4xan_i1.npy"
CSV = REPO / "artifacts" / "m4_atom14" / "m4xan_test.csv"


def build_batch(arr: np.ndarray, seqres: str):
    """Mirror MDDataset.collate_fn (atlas=False, frame_start=0) + prep_batch."""
    T, L, _, _ = arr.shape
    aa = np.array([rc.restype_order[c] for c in seqres])  # [L]
    aatype = torch.from_numpy(aa)[None].expand(T, -1)      # [T, L]

    arr_t = torch.from_numpy(arr.astype(np.float32))       # [T, L, 14, 3]
    frames = geo.atom14_to_frames(arr_t)                   # Rigid [T, L]

    # atom14 -> atom37 (pure torch, equivalent to geometry.atom14_to_atom37)
    atom37_map = torch.from_numpy(rc.RESTYPE_ATOM37_TO_ATOM14).long()   # [21, 37]
    atom37_mask = torch.from_numpy(rc.RESTYPE_ATOM37_MASK).float()      # [21, 37]
    idx = atom37_map[aatype]                                # [T, L, 37]
    atom37 = arr_t.gather(2, idx.unsqueeze(-1).expand(T, L, 37, 3))
    atom37 = atom37 * atom37_mask[aatype].unsqueeze(-1)

    torsions, torsion_mask = geo.atom37_to_torsions(atom37, aatype)  # [T,L,7,2],[T,L,7]
    torsion_mask = torsion_mask[0]                                    # [L, 7]

    traj = {
        "torsions": torsions.unsqueeze(0),                     # [1,T,L,7,2]
        "torsion_mask": torsion_mask.unsqueeze(0),             # [1,L,7]
        "trans": frames._trans.unsqueeze(0),                   # [1,T,L,3]
        "rots": frames._rots._rot_mats.unsqueeze(0),           # [1,T,L,3,3]
        "seqres": torch.from_numpy(aa).unsqueeze(0),           # [1,L]
        "mask": torch.ones(1, L, dtype=torch.float32),         # [1,L]
    }

    # --- replicate prep_batch (standard path: no no_frames, no hyena, no tps) ---
    B = 1
    rigids = Rigid(trans=traj["trans"], rots=Rotation(rot_mats=traj["rots"]))
    offsets = get_offsets(rigids[:, 0:1], rigids)
    offsets[..., :4] *= torch.where(offsets[:, :, :, 0:1] < 0, -1, 1)
    frame_loss_mask = traj["mask"].unsqueeze(-1).expand(-1, -1, 7)       # [1,L,7]
    torsion_loss_mask = (
        traj["torsion_mask"].unsqueeze(-1).expand(-1, -1, -1, 2).reshape(B, L, 14)
    )
    latents = torch.cat([offsets, traj["torsions"].view(B, T, L, 14)], -1)  # [1,T,L,21]
    loss_mask = torch.cat([frame_loss_mask, torsion_loss_mask], -1)

    cond_mask = torch.zeros(B, T, L, dtype=int)
    aatype_mask = torch.ones_like(traj["seqres"])
    model_kwargs = {
        "start_frames": rigids[:, 0],
        "end_frames": rigids[:, -1],
        "mask": traj["mask"].unsqueeze(1).expand(-1, T, -1),              # [1,T,L]
        "aatype": torch.where(aatype_mask.bool(), traj["seqres"], 20),
        "x_cond": torch.where(cond_mask.unsqueeze(-1).bool(), latents, 0.0),
        "x_cond_mask": cond_mask,
    }
    return latents, model_kwargs


def parse_args(argv=None):
    """Defaults reproduce the original 100-frame invocation exactly."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npy", type=Path, default=NPY, help="atom14 input array")
    parser.add_argument("--csv", type=Path, default=CSV, help="name,seqres CSV for the same window")
    parser.add_argument("--checkpoint", type=Path, default=CKPT)
    parser.add_argument("--output", type=Path,
                        default=REPO / "artifacts" / "g0_outputs" / "embedding_smoke.npy")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    arr = np.load(args.npy)
    seqres = args.csv.read_text().strip().splitlines()[1].split(",")[1]
    T, L = arr.shape[0], arr.shape[1]
    print(f"atom14     : {arr.shape}  seqres={len(seqres)}  frames={T}")

    latents, model_kwargs = build_batch(arr, seqres)
    print(f"latents    : {tuple(latents.shape)}  (expect [1, {T}, {L}, 21])")

    obj = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    md = {k: v for k, v in obj["state_dict"].items() if k.startswith("network.md.")}
    encoder = TrajectoryEncoder(
        output_dim=1024, hidden_size=21, num_layers=4, num_heads=8,
        pretrained=False, frozen=True, proj_type="mlp", num_frames=T, suffix="_i1",
    )
    encoder.load_state_dict({k[len("network.md."):]: v for k, v in md.items()}, strict=False)
    encoder.eval()

    with torch.no_grad():
        out1 = encoder(latents, 0, **model_kwargs)
        out2 = encoder(latents, 0, **model_kwargs)

    print()
    print("=== G0 embedding results ===")
    print(f"  output shape    : {tuple(out1.shape)}")
    print(f"  finite          : {bool(torch.isfinite(out1).all())}")
    print(f"  non-zero        : {bool(out1.abs().max() > 0)}")
    print(f"  bit-reproducible: {bool(torch.equal(out1, out2))}")
    print(f"  L2 norm         : {float(out1.norm()):.6f}")
    print(f"  mean/abs-max    : {float(out1.mean()):.6f} / {float(out1.abs().max()):.6f}")
    print(f"  head (first 8)  : {[round(float(x), 4) for x in out1[0, :8]]}")

    out_np = out1[0].cpu().numpy()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.output, out_np)
    print(f"  saved smoke     : {args.output}  shape={tuple(out_np.shape)}")

    assert torch.isfinite(out1).all(), "non-finite embedding"
    assert out1.abs().max() > 0, "all-zero embedding"
    assert torch.equal(out1, out2), "not reproducible"
    print()
    print("G0_EMBEDDING_OK")


if __name__ == "__main__":
    main()
