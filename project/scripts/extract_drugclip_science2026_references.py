#!/usr/bin/env python3
"""Encode crystallographic reference ligands with the Science-2026 DrugCLIP checkpoint.

This is intentionally separate from the multi-million-molecule extraction so
existing target archives do not need to be recomputed merely to add the
reference-ligand context required by CGDA/PACER.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from run_drugclip_science2026_full_litpcba import encode, load_model, move, sha256


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--drugclip", type=Path, required=True)
    parser.add_argument("--unicore", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--trusted-checkpoint", action="store_true")
    args = parser.parse_args()

    device = torch.device(args.device)
    task, model = load_model(args, device)
    args.output.mkdir(parents=True, exist_ok=True)
    audit = {}
    for folder in sorted(p for p in args.data_root.iterdir() if p.is_dir()):
        source = folder / "ligand.lmdb"
        if not source.exists():
            continue
        dataset = task.load_retrieval_mols_dataset(str(source), "atoms", "coordinates")
        loader = torch.utils.data.DataLoader(
            dataset, batch_size=args.batch_size, shuffle=False, collate_fn=dataset.collater
        )
        representations, embeddings, identifiers = [], [], []
        with torch.inference_mode():
            for raw in loader:
                identifiers.extend(map(str, raw["smi_name"]))
                representation, embedding = encode(model, move(raw, device), "molecule")
                representations.append(representation)
                embeddings.append(embedding)
        output = args.output / f"{folder.name}.references.npz"
        np.savez_compressed(
            output,
            reference_ids=np.asarray(identifiers),
            reference_representations=np.concatenate(representations).astype(np.float32),
            reference_embeddings=np.concatenate(embeddings).astype(np.float32),
        )
        audit[folder.name] = {"n_references": len(identifiers), "output": str(output)}
        print(f"target={folder.name} references={len(identifiers)}", flush=True)

    report = {
        "method": "Science-2026 DrugCLIP crystallographic reference-ligand encoding",
        "checkpoint": str(args.checkpoint),
        "checkpoint_sha256": sha256(args.checkpoint),
        "targets": audit,
        "claim_boundary": "Unlabelled reference context only; no activity labels or screening result.",
    }
    (args.output / "references.audit.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
