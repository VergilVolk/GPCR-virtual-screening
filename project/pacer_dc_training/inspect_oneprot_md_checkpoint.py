#!/usr/bin/env python
"""Inspect whether a OneProt checkpoint embeds the full MD transformer weights.

This is a necessary pre-instantiation audit.  A final runtime audit must still
verify that loading into the instantiated model leaves no missing MD keys.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import torch


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tensor_bytes(value: object) -> int:
    return int(value.numel() * value.element_size()) if isinstance(value, torch.Tensor) else 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.checkpoint.is_file():
        raise FileNotFoundError(args.checkpoint)

    try:
        payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False, mmap=True)
    except TypeError:
        payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    state = payload.get("state_dict", payload) if isinstance(payload, dict) else payload
    if not isinstance(state, dict):
        raise TypeError(f"Unsupported checkpoint payload: {type(state).__name__}")

    keys = [str(k) for k in state]
    md_keys = [k for k in keys if k.startswith("network.md.") or ".network.md." in k]
    transformer_keys = [k for k in md_keys if ".transformer." in k]
    md_prefixes = Counter(".".join(k.split(".")[:3]) for k in md_keys)
    transformer_bytes = sum(tensor_bytes(state[k]) for k in transformer_keys)
    report = {
        "checkpoint": str(args.checkpoint),
        "file_size_bytes": args.checkpoint.stat().st_size,
        "sha256": sha256(args.checkpoint),
        "state_dict_key_count": len(keys),
        "network_md_key_count": len(md_keys),
        "md_transformer_key_count": len(transformer_keys),
        "md_transformer_tensor_bytes": transformer_bytes,
        "md_prefix_counts": dict(md_prefixes),
        "sample_md_transformer_keys": transformer_keys[:30],
        "embedded_md_transformer_candidate": bool(transformer_keys and transformer_bytes > 0),
        "decision": (
            "candidate_embedded_weights_present_runtime_coverage_audit_required"
            if transformer_keys and transformer_bytes > 0
            else "hard_block_no_embedded_md_transformer_weights"
        ),
        "runtime_requirement": (
            "Instantiate with pretrained=false/model_path=null, load the outer checkpoint, "
            "and require zero missing or shape-mismatched keys under network.md.transformer. "
            "Do not accept strict=false without checking its returned incompatibility report."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["embedded_md_transformer_candidate"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
