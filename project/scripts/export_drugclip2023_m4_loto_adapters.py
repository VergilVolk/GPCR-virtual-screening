#!/usr/bin/env python3
"""Rebuild, validate and export the frozen 2023 M4-held-out adapter ensemble."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import torch

from evaluate_drugclip_large_target_loso import fit
from finetune_drugclip_gpcr_screening import TARGETS


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--representations", type=Path, required=True)
    p.add_argument("--projection", type=Path, required=True)
    p.add_argument("--pairs", type=Path, required=True)
    p.add_argument("--reference-predictions", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--seeds", default="20260925,20260926,20260927")
    p.add_argument("--rank", type=int, default=8)
    p.add_argument("--epochs", type=int, default=60)
    p.add_argument("--lr", type=float, default=5e-3)
    p.add_argument("--weight-decay", type=float, default=1e-3)
    p.add_argument("--temperature", type=float, default=.07)
    p.add_argument("--retrieval-weight", type=float, default=.25)
    p.add_argument("--preserve-weight", type=float, default=.2)
    p.add_argument("--max-reproduction-error", type=float, default=2e-6)
    a = p.parse_args()
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    arc = np.load(a.representations, allow_pickle=False)
    initial = torch.load(a.projection, map_location="cpu", weights_only=False)
    pairs0 = pd.read_csv(a.pairs)
    mol_index = {str(v): i for i, v in enumerate(arc["molecule_ids"])}
    pairs = pairs0[pairs0.canonical_smiles.astype(str).isin(mol_index)].copy().reset_index(drop=True)
    pair_mol = np.asarray([mol_index[str(v)] for v in pairs.canonical_smiles], dtype=int)
    pair_target = np.asarray([TARGETS.index(v) for v in pairs.target], dtype=int)
    labels = pairs.label.to_numpy(int)
    scaffolds = pairs.murcko_scaffold.astype(str).to_numpy()
    ids = list(map(str, arc["pocket_ids"]))
    pocket_order = [ids.index(f"{target}_cluster0") for target in TARGETS]
    molecule_rep = torch.as_tensor(arc["molecule_representations"].astype(np.float32))
    pocket_rep = torch.as_tensor(arc["pocket_representations"][pocket_order].astype(np.float32))

    held = TARGETS.index("M4R")
    test = pair_target == held
    held_scaffolds = set(scaffolds[test])
    train = (pair_target != held) & ~np.asarray([s in held_scaffolds for s in scaffolds])
    seen = [i for i in range(len(TARGETS)) if i != held]
    seeds = [int(v) for v in a.seeds.split(",")]
    fit_args = SimpleNamespace(**{k: getattr(a, k) for k in [
        "rank", "epochs", "lr", "weight_decay", "temperature",
        "retrieval_weight", "preserve_weight",
    ]})
    a.output_dir.mkdir(parents=True, exist_ok=True)
    test_predictions, files = [], []
    for base_seed in seeds:
        effective_seed = base_seed + held * 1009
        model = fit(
            molecule_rep, pocket_rep, pair_mol[train], pair_target[train], labels[train],
            seen, initial, effective_seed, fit_args, random_labels=False,
        )
        with torch.inference_mode():
            zm = model.mol.from_hidden(model.mol.hidden(molecule_rep))
            zp = model.pocket.from_hidden(model.pocket.hidden(pocket_rep))
            matrix = (zm @ zp.T).cpu().numpy()
        prediction = matrix[pair_mol[test], held]
        test_predictions.append(prediction)
        payload = {
            "protocol_id": "drugclip2023_m4_heldout_loto_3seed_v01",
            "base_seed": base_seed,
            "effective_seed": effective_seed,
            "held_target": "M4R",
            "seen_targets": [TARGETS[i] for i in seen],
            "mol_project": model.mol.materialized_state(),
            "pocket_project": model.pocket.materialized_state(),
            "hyperparameters": vars(fit_args),
            "claim_boundary": "Binding-candidate retrieval only; not PAM function or potency.",
        }
        path = a.output_dir / f"drugclip2023_m4_loto_seed{base_seed}.projection.pt"
        torch.save(payload, path)
        files.append(path)

    mean_prediction = np.mean(test_predictions, axis=0)
    reference = pd.read_csv(a.reference_predictions).set_index("pair_id")
    expected = reference.loc[pairs.loc[test, "pair_id"], "tuned"].to_numpy(float)
    error = np.abs(mean_prediction - expected)
    max_error = float(error.max())
    if max_error > a.max_reproduction_error:
        raise RuntimeError(f"Exported adapter ensemble failed reproduction: max error {max_error}")
    manifest = {
        "protocol_id": "drugclip2023_m4_heldout_loto_3seed_v01",
        "held_target": "M4R",
        "seen_targets": [TARGETS[i] for i in seen],
        "base_seeds": seeds,
        "training_pairs": int(train.sum()),
        "m4_validation_pairs": int(test.sum()),
        "held_target_in_training": False,
        "held_scaffold_overlap": 0,
        "reference_reproduction": {
            "max_abs_error": max_error,
            "mean_abs_error": float(error.mean()),
            "threshold": a.max_reproduction_error,
            "passed": True,
        },
        "artifacts": {path.name: {"sha256": sha256(path), "bytes": path.stat().st_size} for path in files},
        "source_sha256": {
            "representations": sha256(a.representations), "projection": sha256(a.projection),
            "pairs": sha256(a.pairs), "reference_predictions": sha256(a.reference_predictions),
        },
        "claim_boundary": "M4-held-out GPCR transfer model for binding-candidate retrieval; not a PAM classifier.",
    }
    (a.output_dir / "adapter_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
