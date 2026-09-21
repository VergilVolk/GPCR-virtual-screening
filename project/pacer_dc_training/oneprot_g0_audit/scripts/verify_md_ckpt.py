#!/usr/bin/env python
"""Verify the extracted OneProt MD checkpoint and decide whether the separate
MDGen ``forward_sim.ckpt`` is actually required.

Run this after fetch_md_ckpt.py reports MD_CHECKPOINT_READY.

Checks, in order:
  1. size + CRC32 of the inflated member, plus locally computed sha256/md5.
     NOTE: Zenodo publishes only the md5 of the *whole* 14.73 GB checkpoints.zip
     (cfda36445f6c21d74058c3a13154d430).  That value CANNOT be reproduced from a
     single extracted member, so the authoritative per-file guarantee here is the
     ZIP member CRC32 recorded in the central directory (2f07819f).  We report the
     local md5/sha256 for the record, not as a cross-check against Zenodo.
  2. state_dict key inventory: does ``model.components.md.transformer.*`` exist?
     If yes, forward_sim.ckpt is not needed and TrajectoryEncoder(pretrained=False)
     is a valid drop-in.
  3. what the OneProt checkpoint would leave missing/unexpected when loaded into
     a freshly constructed (pretrained=False) OneProtLitModule-like namespace.
"""
from __future__ import annotations

import hashlib
import sys
import zipfile
import zlib
from pathlib import Path

CKPT = Path(
    "project/tools/oneprot-embeddings/artifacts/"
    "Pocket_Text_ST_SG_MD/epoch_012_01100-v1.ckpt"
)
EXPECTED_SIZE = 3_614_862_456
EXPECTED_CRC32 = 0x2F07819F


def hash_file(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    print("=" * 72)
    print("STEP 1  file integrity")
    print("=" * 72)
    if not CKPT.is_file():
        raise SystemExit(f"missing: {CKPT}")
    size = CKPT.stat().st_size
    print(f"  path        : {CKPT}")
    print(f"  size        : {size:,}")
    print(f"  expected    : {EXPECTED_SIZE:,}   match={size == EXPECTED_SIZE}")

    # Recompute the member CRC32 straight from the inflated bytes.
    crc = 0
    with CKPT.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            crc = zlib.crc32(block, crc)
    crc &= 0xFFFFFFFF
    print(f"  crc32       : {crc:08x}   expected {EXPECTED_CRC32:08x}"
          f"   match={crc == EXPECTED_CRC32}")
    print(f"  md5 (local) : {hash_file(CKPT, 'md5')}")
    print(f"  sha256      : {hash_file(CKPT, 'sha256')}")
    print("  zip md5 from Zenodo cfda36445f6c21d74058c3a13154d430 "
          "(whole 14.73 GB archive; NOT reproducible from one member)")

    if size != EXPECTED_SIZE or crc != EXPECTED_CRC32:
        raise SystemExit("integrity check FAILED - stopping")

    print()
    print("=" * 72)
    print("STEP 2  state_dict key inventory")
    print("=" * 72)
    try:
        import torch
    except ImportError:
        print("  torch not importable in this interpreter; skipping")
        print("  -> rerun with the pacer-dc-train-cpu environment")
        return

    obj = torch.load(CKPT, map_location="cpu", weights_only=False)
    print(f"  top-level type : {type(obj).__name__}")
    if isinstance(obj, dict):
        print(f"  top-level keys : {sorted(obj.keys())[:20]}")

    sd = obj.get("state_dict", obj) if isinstance(obj, dict) else obj
    if not isinstance(sd, dict):
        raise SystemExit(f"cannot locate a state_dict in {type(obj)}")
    keys = list(sd.keys())
    print(f"  parameter count: {len(keys):,}")

    groups: dict[str, int] = {}
    for k in keys:
        parts = k.split(".")
        g = ".".join(parts[:4]) if len(parts) >= 4 else k
        groups[g] = groups.get(g, 0) + 1
    print("  --- top key groups ---")
    for g, n in sorted(groups.items(), key=lambda x: -x[1])[:25]:
        print(f"    {n:6d}  {g}")

    md_keys = [k for k in keys if ".md." in k or k.startswith("md.")]
    trans_keys = [k for k in keys if "md.transformer" in k]
    print()
    print(f"  keys containing '.md.'            : {len(md_keys)}")
    print(f"  keys containing 'md.transformer'  : {len(trans_keys)}")
    print()
    if trans_keys:
        print("  *** md.transformer.* IS PRESENT ***")
        print("  -> forward_sim.ckpt is NOT required;")
        print("     TrajectoryEncoder(pretrained=False) + OneProt ckpt is sufficient.")
        print("  sample keys:")
        for k in trans_keys[:12]:
            print(f"     {k}   {tuple(sd[k].shape) if hasattr(sd[k], 'shape') else type(sd[k])}")
    else:
        print("  *** md.transformer.* IS ABSENT ***")
        print("  -> the separate MDGen forward_sim.ckpt IS required, and its public")
        print("     download URL currently 404s.  STOPPING as instructed.")
    print()
    print("STEP2_DONE")


if __name__ == "__main__":
    main()
