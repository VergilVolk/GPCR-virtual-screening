from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import torch

from src.proj import ProjectionHead


def load_npz(path: Path) -> dict[str, np.ndarray]:
    archive = np.load(path, allow_pickle=False)
    return {key: archive[key] for key in archive.files}


def screen(
    molecule_npz: Path,
    pocket_npz: Path,
    model_paths: list[Path],
    pooling: str = "mean",
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Score molecules against pocket conformations with a seed ensemble."""
    molecules = load_npz(molecule_npz)
    pockets = load_npz(pocket_npz)
    molecule_ids = [str(value) for value in molecules["molecule_ids"]]
    molecule_repr = torch.as_tensor(
        molecules["molecule_representations"].astype(np.float32)
    )
    pocket_repr = torch.as_tensor(
        pockets["pocket_representations"].astype(np.float32)
    )

    seed_scores = []
    with torch.inference_mode():
        for model_path in model_paths:
            checkpoint = torch.load(model_path, map_location="cpu", weights_only=False)
            molecule_head = ProjectionHead(checkpoint["mol_project"]).eval()
            pocket_head = ProjectionHead(checkpoint["pocket_project"]).eval()
            scores = (molecule_head(molecule_repr) @ pocket_head(pocket_repr).T).numpy()
            if pooling == "max":
                seed_scores.append(scores.max(axis=1))
            elif pooling == "mean":
                seed_scores.append(scores.mean(axis=1))
            else:
                raise ValueError("pooling must be 'mean' or 'max'")

    per_seed = np.stack(seed_scores)
    return molecule_ids, per_seed.mean(axis=0), per_seed


def write_score_table(
    output_path: Path,
    molecule_ids: list[str],
    scores: np.ndarray,
    model_tag: str,
) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    order = np.argsort(-scores)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["molecule_id", "binding_score", "binding_rank", "model_version"])
        for rank, index in enumerate(order, 1):
            writer.writerow(
                [molecule_ids[index], f"{scores[index]:.6f}", rank, model_tag]
            )
    return output_path
