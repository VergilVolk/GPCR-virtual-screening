from pathlib import Path
import math

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path("project/results/pacer_dc_four_context_v01/compound110")
CKPT = Path(
    "project/tools/oneprot-embeddings/artifacts/"
    "Pocket_Text_ST_SG_MD/epoch_012_01100-v1.ckpt"
)
PREFIX = "network.md.transformer."

folder = ROOT / "five_layer_diagnostic_R2_W0_v01"
context = "candidate_probe"

l0 = np.load(folder / f"L0_{context}_R2_W0.npy")
l1 = np.load(folder / f"L1_{context}_R2_W0.npy")

assert l0.shape == (100, 270, 384)
assert l1.shape == (100, 270, 21)
assert np.isfinite(l0).all()
assert np.isfinite(l1).all()

state = torch.load(
    CKPT, map_location="cpu", weights_only=False
)["state_dict"]

def param(name, shape):
    x = state[PREFIX + name]
    assert tuple(x.shape) == shape
    assert torch.isfinite(x).all()
    return x.detach()

tw0 = param("t_embedder.mlp.0.weight", (384, 256))
tb0 = param("t_embedder.mlp.0.bias", (384,))
tw2 = param("t_embedder.mlp.2.weight", (384, 384))
tb2 = param("t_embedder.mlp.2.bias", (384,))

aw = param("emb_to_latent.adaLN_modulation.1.weight", (768, 384))
ab = param("emb_to_latent.adaLN_modulation.1.bias", (768,))
lw = param("emb_to_latent.linear.weight", (21, 384))
lb = param("emb_to_latent.linear.bias", (21,))

with torch.no_grad():
    # Reproduce the original timestep embedding at t=0.
    t = torch.zeros(1, dtype=torch.float32)
    half = 128
    freqs = torch.exp(
        -math.log(10000)
        * torch.arange(half, dtype=torch.float32)
        / half
    )
    phase = t[:, None] * freqs[None]
    t_freq = torch.cat(
        [torch.cos(phase), torch.sin(phase)], dim=-1
    )

    c = F.linear(t_freq, tw0, tb0)
    c = F.silu(c)
    c = F.linear(c, tw2, tb2)

    shift, scale = F.linear(
        F.silu(c), aw, ab
    ).chunk(2, dim=-1)

    x = torch.from_numpy(
        np.array(l0, dtype=np.float32, copy=True)
    )

    # Reproduce the original FinalLayer.
    x_norm = F.layer_norm(x, (384,), eps=1e-6)
    x_mod = (
        x_norm * (1 + scale.reshape(1, 1, 384))
        + shift.reshape(1, 1, 384)
    )
    reconstructed = F.linear(x_mod, lw, lb)

    reference = torch.from_numpy(
        np.array(l1, dtype=np.float32, copy=True)
    )

    error = torch.abs(reconstructed - reference)
    max_error = float(error.max())
    mean_error = float(error.mean())

print("MAX_ABS_ERROR:", max_error)
print("MEAN_ABS_ERROR:", mean_error)

assert max_error < 1e-4
assert torch.allclose(
    reconstructed, reference,
    rtol=1e-5, atol=1e-4
)

print("L1_RECONSTRUCTION_QC: PASS")
