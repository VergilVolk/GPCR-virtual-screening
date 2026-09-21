#!/usr/bin/env python
"""Independent re-derivation of the G0 embedding from this package's own inputs.

Purpose: prove that `outputs/embedding_smoke.npy` is reproducible from
`input/m4xan_i1.npy` + `input/m4xan_test.csv` and the published checkpoint, by
re-running the forward in a fresh process and comparing bit-for-bit.

This script deliberately reuses the ORIGINAL preprocessing code
(`scripts/run_g0_embedding.py::build_batch`) instead of re-implementing it, so a
disagreement would point at the model/checkpoint rather than at a re-written
data path. It writes nothing; its report is captured to
`outputs/forward_reverification.txt` by the caller.

Target: OneProt-MD encoder (TrajectoryEncoder -> LatentMDGenModel, the MD branch),
instantiated with pretrained=False and then loaded from the outer checkpoint,
exactly as START_HERE_ONEPROT_PACER_DC.md prescribes.

Prerequisites (not shipped in this repository, see README):
  * project/tools/oneprot-embeddings       - upstream OneProt checkout @53fa9c0
    with its external/mdgen submodule @1b82cc1
  * the 3.6 GB checkpoint at
    project/tools/oneprot-embeddings/artifacts/Pocket_Text_ST_SG_MD/epoch_012_01100-v1.ckpt
  * a torch CPU environment (validated with torch 2.4.1, python 3.11, numpy 2.4.6)

Run from the repository root:
  python project/pacer_dc_training/oneprot_g0_audit/scripts/reverify_forward.py
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import torch

PKG = Path(__file__).resolve().parents[1]
REPO_ROOT = PKG.parents[2]
NESTED = REPO_ROOT / "project" / "tools" / "oneprot-embeddings"
CKPT = NESTED / "artifacts" / "Pocket_Text_ST_SG_MD" / "epoch_012_01100-v1.ckpt"


def parse_args(argv=None):
    """Defaults reproduce the original 100-frame check exactly."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npy", type=Path, default=PKG / "input" / "m4xan_i1.npy")
    parser.add_argument("--csv", type=Path, default=PKG / "input" / "m4xan_test.csv")
    parser.add_argument("--saved", type=Path,
                        default=PKG / "outputs" / "embedding_smoke.npy")
    parser.add_argument("--label", default="100-frame",
                        help="free-text label for the report header")
    return parser.parse_args(argv)


ARGS = parse_args()

print(f"python      : {sys.version.split()[0]}")
print(f"torch       : {torch.__version__}")
print(f"numpy       : {np.__version__}")
print(f"variant     : {ARGS.label}")
print(f"package     : {PKG}")
print(f"input array : {ARGS.npy}")
print(f"saved embed : {ARGS.saved}")
print(f"nested repo : {NESTED}  exists={NESTED.is_dir()}")
print(f"checkpoint  : {CKPT}  exists={CKPT.is_file()}")
if not CKPT.is_file():
    raise SystemExit("checkpoint absent - this verification cannot run without it")

# --- reuse the original preprocessing verbatim -------------------------------
spec = importlib.util.spec_from_file_location("g0", PKG / "scripts" / "run_g0_embedding.py")
g0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g0)          # registers src/mdgen stubs, imports TrajectoryEncoder

arr = np.load(ARGS.npy)
seqres = ARGS.csv.read_text().strip().splitlines()[1].split(",")[1]
T, L = arr.shape[0], arr.shape[1]
print(f"\natom14 input: shape={arr.shape} dtype={arr.dtype} L={L} seqres_len={len(seqres)} T={T}")
print(f"finite      : {bool(np.isfinite(arr).all())}")

latents, model_kwargs = g0.build_batch(arr, seqres)
print(f"latents     : {tuple(latents.shape)}  (expect [1, {T}, {L}, 21])")
print(f"finite      : {bool(torch.isfinite(latents).all())}")

t0 = time.time()
obj = torch.load(CKPT, map_location="cpu", weights_only=False)
md = {k: v for k, v in obj["state_dict"].items() if k.startswith("network.md.")}
print(f"checkpoint  : {len(obj['state_dict'])} top-level keys, {len(md)} under network.md. "
      f"(loaded in {time.time() - t0:.1f}s)")

encoder = g0.TrajectoryEncoder(
    output_dim=1024, hidden_size=21, num_layers=4, num_heads=8,
    pretrained=False, frozen=True, proj_type="mlp", num_frames=T, suffix="_i1",
)
load_result = encoder.load_state_dict(
    {k[len("network.md."):]: v for k, v in md.items()}, strict=False
)
encoder.eval()
print(f"load_state_dict: missing={len(load_result.missing_keys)} "
      f"unexpected={len(load_result.unexpected_keys)} "
      f"{load_result.unexpected_keys}")

with torch.no_grad():
    first = encoder(latents, 0, **model_kwargs)
    second = encoder(latents, 0, **model_kwargs)

saved = np.load(ARGS.saved)
out = first[0].cpu().numpy()

print("\n=== re-derivation result ===")
print(f"  output shape          : {tuple(first.shape)}   saved shape: {tuple(saved.shape)}")
print(f"  finite                : {bool(torch.isfinite(first).all())}")
print(f"  non-zero              : {bool(first.abs().max() > 0)}")
print(f"  same-process repeat   : {bool(torch.equal(first, second))}")
print(f"  identical_to_package  : {bool(np.array_equal(out, saved))}")
print(f"  max_abs_diff_vs_saved : {float(np.abs(out - saved).max())}")
print(f"  L2 norm               : {float(first.norm()):.6f}")
print(f"  nonzero values        : {int((out != 0).sum())}/{out.size}")
print(f"  sha256(re-derived)    : {hashlib.sha256(out.tobytes()).hexdigest()}")
print(f"  sha256(package copy)  : {hashlib.sha256(saved.tobytes()).hexdigest()}")
print(f"  saved file sha256     : {hashlib.sha256(ARGS.saved.read_bytes()).hexdigest()}")

ok = bool(np.array_equal(out, saved) and torch.equal(first, second)
          and torch.isfinite(first).all() and first.abs().max() > 0)
print("\nREVERIFY_PASSED" if ok else "\nREVERIFY_FAILED")
raise SystemExit(0 if ok else 2)
