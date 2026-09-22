#!/usr/bin/env python3
"""Extract frozen DrugCLIP molecule and pocket embeddings on CPU.

This bypasses the official retrieval helper's unconditional CUDA transfer but
uses the official dataset transforms, model architecture, and checkpoint.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch


def encode(model, sample: dict, kind: str) -> torch.Tensor:
    net = sample["net_input"]
    prefix = "mol" if kind == "molecule" else "pocket"
    component = model.mol_model if kind == "molecule" else model.pocket_model
    project = model.mol_project if kind == "molecule" else model.pocket_project
    tokens = net[f"{prefix}_src_tokens"]
    distance = net[f"{prefix}_src_distance"]
    edge_type = net[f"{prefix}_src_edge_type"]
    padding_mask = tokens.eq(component.padding_idx)
    x = component.embed_tokens(tokens)
    n_node = distance.size(-1)
    bias = component.gbf_proj(component.gbf(distance, edge_type))
    bias = bias.permute(0, 3, 1, 2).contiguous().view(-1, n_node, n_node)
    representation = component.encoder(x, padding_mask=padding_mask, attn_mask=bias)[0][:, 0, :]
    embedding = project(representation)
    return torch.nn.functional.normalize(embedding, dim=1)


def batches(dataset, batch_size: int):
    return torch.utils.data.DataLoader(
        dataset, batch_size=batch_size, shuffle=False, collate_fn=dataset.collater
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--drugclip", type=Path, required=True)
    parser.add_argument("--unicore", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--molecules", type=Path, required=True)
    parser.add_argument("--pockets", type=Path, required=True, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=8)
    args = parser.parse_args()

    sys.path[:0] = [str(args.drugclip), str(args.unicore)]
    import unimol  # noqa: F401  registers DrugCLIP model/task
    from unicore import checkpoint_utils, tasks

    state = checkpoint_utils.load_checkpoint_to_cpu(str(args.checkpoint))
    model_args = state["args"]
    original_registry = {"task": model_args.task, "arch": model_args.arch}
    # Public checkpoint kept the pre-release registry name while the public
    # repository renamed both registrations to ``drugclip``.
    if model_args.task == "binding_affinity":
        model_args.task = "drugclip"
    if model_args.arch == "binding_affinity":
        model_args.arch = "drugclip"
    model_args.data = str(args.drugclip / "data")
    model_args.cpu = True
    model_args.fp16 = False
    # These private paths initialized the two encoders during original training;
    # every resulting tensor is already present in the full DrugCLIP checkpoint.
    model_args.finetune_mol_model = None
    model_args.finetune_pocket_model = None
    task = tasks.setup_task(model_args)
    model = task.build_model(model_args)
    missing, unexpected = model.load_state_dict(state["model"], strict=False)
    model.cpu().float().eval()

    mol_dataset = task.load_retrieval_mols_dataset(str(args.molecules), "atoms", "coordinates")
    molecule_embeddings, molecule_ids = [], []
    pocket_embeddings, pocket_ids = [], []
    with torch.inference_mode():
        for sample in batches(mol_dataset, args.batch_size):
            molecule_embeddings.append(encode(model, sample, "molecule").cpu().numpy())
            molecule_ids.extend(map(str, sample["smi_name"]))
        for pocket_path in args.pockets:
            pocket_dataset = task.load_pockets_dataset(str(pocket_path))
            for sample in batches(pocket_dataset, args.batch_size):
                pocket_embeddings.append(encode(model, sample, "pocket").cpu().numpy())
                pocket_ids.extend(map(str, sample["pocket_name"]))
    mol = np.concatenate(molecule_embeddings).astype(np.float32)
    poc = np.concatenate(pocket_embeddings).astype(np.float32)
    scores = poc @ mol.T
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, molecule_ids=np.asarray(molecule_ids), molecule_embeddings=mol,
                        pocket_ids=np.asarray(pocket_ids), pocket_embeddings=poc, scores=scores)
    audit = {
        "checkpoint": str(args.checkpoint), "pocket_inputs": list(map(str, args.pockets)),
        "device": "cpu", "molecule_shape": list(mol.shape),
        "checkpoint_registry_remap": {"from": original_registry, "to": {
            "task": model_args.task, "arch": model_args.arch}},
        "pocket_shape": list(poc.shape), "score_shape": list(scores.shape),
        "finite": bool(np.isfinite(scores).all()), "score_range": [float(scores.min()), float(scores.max())],
        "missing_checkpoint_keys": list(missing), "unexpected_checkpoint_keys": list(unexpected),
        "claim_boundary": "Frozen DrugCLIP binding-compatibility embeddings; not PAM efficacy.",
    }
    args.output.with_suffix(".audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
