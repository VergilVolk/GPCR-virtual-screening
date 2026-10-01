#!/usr/bin/env python3
"""Encode a molecules-only LMDB with the Science-2026 checkpoint:
outputs npz with molecule_ids / embeddings / REPRESENTATIONS (for LoRA
fine-tuning) keyed by canonical SMILES (RDKit-canonicalized)."""
from __future__ import annotations
import argparse, os, sys
from pathlib import Path
import numpy as np
import torch
from rdkit import Chem

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import run_drugclip_science2026_full_litpcba as full_runner


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--drugclip", type=Path, required=True)
    p.add_argument("--unicore", type=Path, required=True)
    p.add_argument("--checkpoint", type=Path, required=True)
    p.add_argument("--mols-lmdb", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--trusted-checkpoint", action="store_true")
    a = p.parse_args()
    a.device = "cpu"
    if a.trusted_checkpoint:
        os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
    a.output.parent.mkdir(parents=True, exist_ok=True)
    task, model = full_runner.load_model(a, torch.device("cpu"))
    mol_loader_fn = getattr(task, "load_mols_dataset", None) or task.load_retrieval_mols_dataset
    ds = mol_loader_fn(str(a.mols_lmdb), "atoms", "coordinates")
    loader = torch.utils.data.DataLoader(ds, batch_size=a.batch_size, shuffle=False, collate_fn=ds.collater)
    embs, reps, ids = [], [], []
    with torch.inference_mode():
        for i, raw in enumerate(loader, 1):
            ids.extend(map(str, raw["smi_name"]))
            rep, emb = full_runner.encode(model, full_runner.move(raw, torch.device("cpu")), "molecule")
            reps.append(rep); embs.append(emb)
            if i % 100 == 0:
                print(f"batches={i} rows={len(ids)}", flush=True)
    canon = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) or s for s in ids]
    np.savez_compressed(a.output,
                        molecule_ids=np.asarray(canon),
                        molecule_embeddings=np.concatenate(embs).astype(np.float32),
                        molecule_representations=np.concatenate(reps).astype(np.float32))
    print(f"saved {len(canon)} molecules -> {a.output}")


if __name__ == "__main__":
    main()
