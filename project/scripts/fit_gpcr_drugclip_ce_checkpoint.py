#!/usr/bin/env python3
"""Fit deployable balanced-CE DrugCLIP projection adapters after OOF selection."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from finetune_drugclip_gpcr_retrieval import TARGETS, fit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seeds", default="20260925,20260926,20260927")
    parser.add_argument("--rank", type=int, default=4)
    parser.add_argument("--epochs", type=int, default=250)
    parser.add_argument("--lr", type=float, default=5e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-3)
    parser.add_argument("--temperature", type=float, default=0.07)
    parser.add_argument("--margin", type=float, default=0.1)
    parser.add_argument("--hard-weight", type=float, default=0.5)
    parser.add_argument("--preserve-weight", type=float, default=0.2)
    args = parser.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    archive = np.load(args.representations, allow_pickle=False)
    initial = torch.load(args.projection, map_location="cpu")
    table_raw = pd.read_csv(args.benchmark)
    molecule_index = {value: i for i, value in enumerate(map(str, archive["molecule_ids"]))}
    table = table_raw[table_raw.canonical_smiles.astype(str).isin(molecule_index)].copy().reset_index(drop=True)
    order = np.asarray([molecule_index[value] for value in table.canonical_smiles.astype(str)])
    molecule_rep = torch.as_tensor(archive["molecule_representations"].astype(np.float32)[order])
    pocket_rep = torch.as_tensor(archive["pocket_representations"].astype(np.float32))
    labels = torch.as_tensor([TARGETS.index(value) for value in table.target], dtype=torch.long)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    checkpoints = []
    for seed in map(int, args.seeds.split(",")):
        model = fit(molecule_rep, pocket_rep, labels, initial, seed, "ce", args)
        path = args.output_dir / f"drugclip_gpcr_balanced_ce_seed{seed}.projection.pt"
        torch.save({
            "mol_project": model.mol.materialized_state(),
            "pocket_project": model.pocket.materialized_state(),
            "logit_scale": copy.deepcopy(initial.get("logit_scale")),
            "seed": seed,
            "training_task": "four-GPCR target-balanced CE",
            "benchmark_rows": int(len(table)),
            "targets": TARGETS,
            "selection_basis": "five-fold Murcko-scaffold OOF; CE selected over CE+triplet",
        }, path)
        checkpoints.append(str(path))
    manifest = {
        "selected_method": "target-balanced CE dual-projection rank-4 LoRA",
        "rejected_add_on": "hardest-wrong-pocket triplet did not improve OOF macro Recall@1",
        "checkpoints": checkpoints,
        "claim_boundary": "Deployable target-retrieval adapter; not a PAM efficacy predictor.",
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
