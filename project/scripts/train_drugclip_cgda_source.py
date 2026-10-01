#!/usr/bin/env python3
"""Train a frozen-protocol CGDA checkpoint on all source benchmark targets."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

from evaluate_drugclip_cgda_loso import fit_fold, load_data
from evaluate_drugclip_full_litpcba_reference_fusion import pocket_map


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pocket-root", type=Path, required=True)
    parser.add_argument("--pocket-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--rank", type=int, default=2)
    parser.add_argument("--experts", type=int, default=2)
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--max-negatives-per-target", type=int, default=1000)
    parser.add_argument("--hard-negative-fraction", type=float, default=0.5)
    parser.add_argument("--ranking-weight", type=float, default=0.2)
    parser.add_argument("--ranking-margin", type=float, default=0.1)
    parser.add_argument("--lr", type=float, default=2e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.12)
    parser.add_argument("--retrieval-temperature", type=float, default=0.07)
    parser.add_argument("--retrieval-weight", type=float, default=0.2)
    parser.add_argument("--preserve-weight", type=float, default=0.5)
    parser.add_argument("--validation-fraction", type=float, default=0.2)
    parser.add_argument("--chunk-size", type=int, default=65536)
    args = parser.parse_args()

    pockets = pocket_map(args.pocket_root, args.pocket_archive)
    data = load_data(args.root, pockets)
    model, audit = fit_fold(data, "__external_target__", args, args.seed)
    payload = {
        "model_state": model.state_dict(),
        "architecture": {"dimension": 128, "rank": args.rank, "experts": args.experts},
        "training_targets": sorted(data),
        "training_audit": audit,
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "claim_boundary": "Source-trained checkpoint; external target labels must never be used for selection or fitting.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    torch.save(payload, args.output)
    audit_path = args.output.with_suffix(".audit.json")
    audit_path.write_text(json.dumps({key: value for key, value in payload.items() if key != "model_state"}, indent=2),
                          encoding="utf-8")
    print(json.dumps({"output": str(args.output), "audit": str(audit_path),
                      "targets": len(data), "training_audit": audit}, indent=2))


if __name__ == "__main__":
    main()
