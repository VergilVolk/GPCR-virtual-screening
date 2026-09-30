#!/usr/bin/env python3
"""Stratified-subset encoding of remaining LIT-PCBA targets with the
Science-2026 checkpoint (7-hour interim gate).

Row indices are derived deterministically from the OLD-checkpoint pickle
labels (all actives + N seeded random decoys per target). The LMDB iteration
order matches the old pickle order (verified on ALDH1: identical labels).
Writes the same npz/metrics layout as run_drugclip_science2026_full_litpcba.py
plus a row_indices array. Full-scale run continues independently.
"""
from __future__ import annotations
import argparse, json, os, pickle, sys
from pathlib import Path
import numpy as np
import torch

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import run_drugclip_science2026_full_litpcba as full_runner


def subset_indices(labels: np.ndarray, decoys: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    active = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    rng.shuffle(negative)
    keep = np.sort(np.concatenate([active, negative[:decoys]]))
    return keep


def run_subset_target(task, model, folder, output, device, batch_size, labels_old, decoys, seed):
    mol_loader_fn = getattr(task, "load_mols_dataset", None) or task.load_retrieval_mols_dataset
    mol_dataset = mol_loader_fn(str(folder / "mols.lmdb"), "atoms", "coordinates")
    keep = subset_indices(labels_old, decoys, seed)
    subset = torch.utils.data.Subset(mol_dataset, keep.tolist())
    loader = torch.utils.data.DataLoader(subset, batch_size=batch_size, shuffle=False,
                                         collate_fn=mol_dataset.collater)
    pocket_dataset = task.load_pockets_dataset(str(folder / "pockets.lmdb"))
    pocket_loader = torch.utils.data.DataLoader(pocket_dataset, batch_size=batch_size,
                                                shuffle=False, collate_fn=pocket_dataset.collater)
    molecule_embeddings, molecule_ids, labels = [], [], []
    molecule_representations, pocket_representations = [], []
    pocket_embeddings, pocket_ids = [], []
    with torch.inference_mode():
        for batch_index, raw in enumerate(loader, 1):
            labels.extend(torch.as_tensor(raw["target"]).cpu().numpy().astype(int).tolist())
            molecule_ids.extend(map(str, raw["smi_name"]))
            representation, embedding = full_runner.encode(model, full_runner.move(raw, device), "molecule")
            molecule_embeddings.append(embedding)
            molecule_representations.append(representation)
            if batch_index % 100 == 0:
                print(f"target={folder.name} subset_batches={batch_index} rows={len(labels)}", flush=True)
        for raw in pocket_loader:
            pocket_ids.extend(map(str, raw["pocket_name"]))
            representation, embedding = full_runner.encode(model, full_runner.move(raw, device), "pocket")
            pocket_embeddings.append(embedding)
            pocket_representations.append(representation)
    mol = np.concatenate(molecule_embeddings).astype(np.float32)
    pocket = np.concatenate(pocket_embeddings).astype(np.float32)
    labels_array = np.asarray(labels, dtype=np.int8)
    scores = (pocket @ mol.T).max(axis=0)
    result = full_runner.metrics(labels_array, scores)
    result.update({"target": folder.name, "pockets": int(len(pocket)),
                   "subset_rows": int(len(labels_array)),
                   "subset_actives": int(labels_array.sum()), "device": str(device)})
    np.savez_compressed(output / f"{folder.name}.npz",
                        molecule_ids=np.asarray(molecule_ids),
                        molecule_embeddings=mol,
                        molecule_representations=np.concatenate(molecule_representations).astype(np.float32),
                        pocket_ids=np.asarray(pocket_ids),
                        pocket_embeddings=pocket,
                        pocket_representations=np.concatenate(pocket_representations).astype(np.float32),
                        labels=labels_array,
                        scores=scores.astype(np.float32),
                        row_indices=keep.astype(np.int64))
    (output / f"{folder.name}.metrics.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result), flush=True)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--drugclip", type=Path, required=True)
    p.add_argument("--unicore", type=Path, required=True)
    p.add_argument("--checkpoint", type=Path, required=True)
    p.add_argument("--data-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--old-root", type=Path, required=True)
    p.add_argument("--targets", required=True)
    p.add_argument("--decoys-per-target", type=int, default=20000)
    p.add_argument("--subset-seed", type=int, default=20261001)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--trusted-checkpoint", action="store_true")
    a = p.parse_args()
    a.device = "cpu"
    if a.trusted_checkpoint:
        os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
    a.output.mkdir(parents=True, exist_ok=True)
    folders = {f.name: f for f in full_runner.find_targets(a.data_root)}
    requested = [v.strip() for v in a.targets.split(",") if v.strip()]
    task, model = full_runner.load_model(a, torch.device("cpu"))
    for name in requested:
        if (a.output / f"{name}.metrics.json").exists() and (a.output / f"{name}.npz").exists():
            print(f"resume_skip={name}", flush=True)
            continue
        m, _, lab = pickle.load((a.old_root / name / "drugclip_emb" / "mols.lmdb.pkl").open("rb"))
        run_subset_target(task, model, folders[name], a.output, torch.device("cpu"),
                          a.batch_size, np.asarray(lab), a.decoys_per_target, a.subset_seed)


if __name__ == "__main__":
    main()
