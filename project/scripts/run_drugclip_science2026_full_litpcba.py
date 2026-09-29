#!/usr/bin/env python3
"""Resumable Science-2026 DrugCLIP evaluation on the official full LIT-PCBA archive."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import torch
from rdkit.ML.Scoring.Scoring import CalcBEDROC
from sklearn.metrics import average_precision_score, roc_auc_score


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def move(value, device: torch.device):
    if torch.is_tensor(value):
        return value.to(device)
    if isinstance(value, dict):
        return {key: move(item, device) for key, item in value.items()}
    if isinstance(value, list):
        return [move(item, device) for item in value]
    if isinstance(value, tuple):
        return tuple(move(item, device) for item in value)
    return value


def encode(model, sample: dict, kind: str) -> tuple[np.ndarray, np.ndarray]:
    prefix = "mol" if kind == "molecule" else "pocket"
    component = model.mol_model if kind == "molecule" else model.pocket_model
    projection = model.mol_project if kind == "molecule" else model.pocket_project
    net = sample["net_input"]
    tokens = net[f"{prefix}_src_tokens"]
    distance = net[f"{prefix}_src_distance"]
    edge_type = net[f"{prefix}_src_edge_type"]
    padding = tokens.eq(component.padding_idx)
    x = component.embed_tokens(tokens)
    nodes = distance.size(-1)
    bias = component.gbf_proj(component.gbf(distance, edge_type))
    bias = bias.permute(0, 3, 1, 2).contiguous().view(-1, nodes, nodes)
    representation = component.encoder(x, padding_mask=padding, attn_mask=bias)[0][:, 0, :]
    embedding = torch.nn.functional.normalize(projection(representation), dim=-1)
    return representation.detach().float().cpu().numpy(), embedding.detach().float().cpu().numpy()


def metrics(labels: np.ndarray, scores: np.ndarray) -> dict:
    ranked = [[float(score), int(label)] for score, label in sorted(
        zip(scores, labels), key=lambda item: item[0], reverse=True
    )]
    result = {
        "n": int(len(labels)),
        "actives": int(labels.sum()),
        "active_fraction": float(labels.mean()),
        "roc_auc": float(roc_auc_score(labels, scores)),
        "pr_auc": float(average_precision_score(labels, scores)),
        "bedroc_alpha80_5": float(CalcBEDROC(ranked, 1, 80.5)),
    }
    order = np.argsort(-scores, kind="mergesort")
    for fraction in (0.005, 0.01, 0.02, 0.05):
        count = max(1, int(np.ceil(fraction * len(labels))))
        result[f"ef{fraction:g}"] = float(labels[order[:count]].mean() / labels.mean())
    return result


def find_targets(root: Path) -> list[Path]:
    targets = []
    for mols in root.rglob("mols.lmdb"):
        folder = mols.parent
        if (folder / "pockets.lmdb").exists():
            targets.append(folder)
    # Smallest LMDBs first so a long CPU run produces auditable checkpoints
    # early instead of blocking for hours on the largest alphabetical target.
    return sorted(targets, key=lambda value: ((value / "mols.lmdb").stat().st_size, value.name))


def load_model(args: argparse.Namespace, device: torch.device):
    if args.trusted_checkpoint:
        os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
    sys.path[:0] = [str(args.drugclip), str(args.unicore)]
    import unimol  # noqa: F401
    from unicore import checkpoint_utils, tasks

    state = checkpoint_utils.load_checkpoint_to_cpu(str(args.checkpoint))
    model_args = state["args"]
    if model_args.task == "binding_affinity":
        model_args.task = "drugclip"
    if model_args.arch == "binding_affinity":
        model_args.arch = "drugclip"
    model_args.data = str(args.drugclip / "data")
    model_args.cpu = device.type == "cpu"
    model_args.fp16 = False
    model_args.finetune_mol_model = None
    model_args.finetune_pocket_model = None
    task = tasks.setup_task(model_args)
    model = task.build_model(model_args)
    missing, unexpected = model.load_state_dict(state["model"], strict=False)
    if missing or unexpected:
        raise RuntimeError(f"Checkpoint mismatch: missing={missing}, unexpected={unexpected}")
    return task, model.to(device).float().eval()


def run_target(task, model, folder: Path, output: Path, device: torch.device,
               batch_size: int, save_representations: bool) -> dict:
    mol_loader_fn = getattr(task, "load_mols_dataset", None) or task.load_retrieval_mols_dataset
    mol_dataset = mol_loader_fn(str(folder / "mols.lmdb"), "atoms", "coordinates")
    pocket_dataset = task.load_pockets_dataset(str(folder / "pockets.lmdb"))
    molecule_embeddings, molecule_representations, molecule_ids, labels = [], [], [], []
    pocket_embeddings, pocket_representations, pocket_ids = [], [], []
    mol_loader = torch.utils.data.DataLoader(
        mol_dataset, batch_size=batch_size, shuffle=False, collate_fn=mol_dataset.collater
    )
    pocket_loader = torch.utils.data.DataLoader(
        pocket_dataset, batch_size=batch_size, shuffle=False, collate_fn=pocket_dataset.collater
    )
    with torch.inference_mode():
        for batch_index, raw in enumerate(mol_loader, 1):
            labels.extend(torch.as_tensor(raw["target"]).cpu().numpy().astype(int).tolist())
            molecule_ids.extend(map(str, raw["smi_name"]))
            representation, embedding = encode(model, move(raw, device), "molecule")
            molecule_embeddings.append(embedding)
            if save_representations:
                molecule_representations.append(representation)
            if batch_index % 100 == 0:
                print(f"target={folder.name} molecule_batches={batch_index} rows={len(labels)}", flush=True)
        for raw in pocket_loader:
            pocket_ids.extend(map(str, raw["pocket_name"]))
            representation, embedding = encode(model, move(raw, device), "pocket")
            pocket_embeddings.append(embedding)
            if save_representations:
                pocket_representations.append(representation)

    mol = np.concatenate(molecule_embeddings).astype(np.float32)
    pocket = np.concatenate(pocket_embeddings).astype(np.float32)
    labels_array = np.asarray(labels, dtype=np.int8)
    score_matrix = pocket @ mol.T
    scores = score_matrix.max(axis=0)
    result = metrics(labels_array, scores)
    result.update({"target": folder.name, "pockets": int(len(pocket)), "device": str(device)})
    payload = {
        "molecule_ids": np.asarray(molecule_ids),
        "molecule_embeddings": mol,
        "pocket_ids": np.asarray(pocket_ids),
        "pocket_embeddings": pocket,
        "labels": labels_array,
        "scores": scores.astype(np.float32),
    }
    if save_representations:
        payload["molecule_representations"] = np.concatenate(molecule_representations).astype(np.float32)
        payload["pocket_representations"] = np.concatenate(pocket_representations).astype(np.float32)
    np.savez_compressed(output / f"{folder.name}.npz", **payload)
    (output / f"{folder.name}.metrics.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result), flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--drugclip", type=Path, required=True)
    parser.add_argument("--unicore", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--source-zip", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--targets", default="", help="Optional comma-separated target subset")
    parser.add_argument("--save-representations", action="store_true")
    parser.add_argument("--trusted-checkpoint", action="store_true")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but torch.cuda.is_available() is false")
    args.output.mkdir(parents=True, exist_ok=True)
    targets = find_targets(args.data_root)
    if args.targets:
        requested = {value.strip() for value in args.targets.split(",") if value.strip()}
        targets = [folder for folder in targets if folder.name in requested]
        missing = requested - {folder.name for folder in targets}
        if missing:
            raise ValueError(f"Requested targets not found: {sorted(missing)}")
    if not targets:
        raise FileNotFoundError(f"No target folders with mols.lmdb and pockets.lmdb under {args.data_root}")
    task, model = load_model(args, device)
    per_target = {}
    for folder in targets:
        metric_path = args.output / f"{folder.name}.metrics.json"
        archive_path = args.output / f"{folder.name}.npz"
        if metric_path.exists() and archive_path.exists():
            per_target[folder.name] = json.loads(metric_path.read_text(encoding="utf-8"))
            print(f"resume_skip={folder.name}", flush=True)
            continue
        per_target[folder.name] = run_target(
            task, model, folder, args.output, device, args.batch_size, args.save_representations
        )
    metric_names = ["roc_auc", "pr_auc", "bedroc_alpha80_5", "ef0.005", "ef0.01", "ef0.02", "ef0.05"]
    summary = {
        "method": "Verified Science-2026 DrugCLIP single checkpoint",
        "checkpoint": str(args.checkpoint),
        "checkpoint_sha256": sha256(args.checkpoint),
        "source_zip": str(args.source_zip),
        "source_zip_sha256": sha256(args.source_zip),
        "targets": list(per_target),
        "macro": {name: float(np.mean([value[name] for value in per_target.values()])) for name in metric_names},
        "per_target": per_target,
        "claim_boundary": "Official full LIT-PCBA single-checkpoint reproduction; not six-fold ensemble or PAM efficacy.",
    }
    summary_name = "summary.json" if not args.targets else f"summary.{args.targets.replace(',', '_')}.json"
    (args.output / summary_name).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
