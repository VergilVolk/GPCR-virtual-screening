#!/usr/bin/env python3
"""Apply a small, audited CPU compatibility patch to an official DrugCLIP checkout.

The official repository is kept outside Git. This script refuses unknown source
text and records every replacement, so CPU inference cannot silently drift from
the pinned implementation.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


EXPECTED_COMMIT = "7a3a3fa33673f8668c811790f2e4681c98af44ef"
REPLACEMENTS = {
    'torch.ones([1], device="cuda")': 'torch.ones([1])',
    'torch.tensor(pocket_duplicate_matrix, dtype=ba_predict.dtype).cuda()': (
        'torch.tensor(pocket_duplicate_matrix, dtype=ba_predict.dtype, device=ba_predict.device)'
    ),
    'torch.tensor(mol_duplicate_matrix, dtype=ba_predict.dtype).cuda()': (
        'torch.tensor(mol_duplicate_matrix, dtype=ba_predict.dtype, device=ba_predict.device)'
    ),
    'torch.eye(bsz).cuda()': 'torch.eye(bsz, device=ba_predict.device)',
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()

    head = subprocess.check_output(
        ["git", "-C", str(args.checkout), "rev-parse", "HEAD"], text=True
    ).strip()
    if head != EXPECTED_COMMIT:
        raise RuntimeError(f"DrugCLIP commit mismatch: {head} != {EXPECTED_COMMIT}")

    target = args.checkout / "unimol" / "models" / "drugclip.py"
    text = target.read_text(encoding="utf-8")
    changed, already = [], []
    for old, new in REPLACEMENTS.items():
        if old in text:
            if text.count(old) != 1:
                raise RuntimeError(f"Ambiguous source occurrence: {old}")
            text = text.replace(old, new)
            changed.append({"from": old, "to": new})
        elif new in text:
            already.append(new)
        else:
            raise RuntimeError(f"Expected source text not found: {old}")
    target.write_text(text, encoding="utf-8")

    payload = {
        "official_commit": head,
        "target": str(target),
        "changed": changed,
        "already_patched": already,
        "scope": "device placement only; architecture and checkpoint tensors unchanged",
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
