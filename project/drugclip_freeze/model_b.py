"""Model B: frozen projection-only apply path.

Model B never re-runs the DrugCLIP backbone.  It consumes Model A's frozen
512-D `molecule_representations` / `pocket_representations` and applies the
already-trained projection heads, then takes the official max-pooled pocket
score.  Three deployed seeds are averaged exactly as the frozen pipeline did.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from project.drugclip_freeze import model_a, projection

DEPLOYED_SEEDS = (20260925, 20260926, 20260927)
FAMILY_AUG_DIR = "project/results/drugclip_science2026/family_aug_v01"
FAMILY_AUG_GLOB = "science2026_13target_ep80_famaug.seed{seed}.projection.pt"


class ModelBError(ValueError):
    """Raised when the frozen Model B apply contract is not satisfied."""


def default_projection_paths(repo_root: Path) -> list[Path]:
    directory = repo_root / FAMILY_AUG_DIR
    paths = [directory / FAMILY_AUG_GLOB.format(seed=seed) for seed in DEPLOYED_SEEDS]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise ModelBError("deployed Model B projection weights are missing: " + ", ".join(missing))
    return paths


def score(archive_path: Path, projection_paths: Sequence[Path], *,
          pocket_ids: Sequence[str] = model_a.M4_POCKET_IDS) -> dict[str, Any]:
    """Apply the frozen projections and return per-seed and ensemble scores."""
    if len(projection_paths) != len(DEPLOYED_SEEDS):
        raise ModelBError(f"expected {len(DEPLOYED_SEEDS)} deployed projection files, got {len(projection_paths)}")
    if not Path(archive_path).is_file():
        raise ModelBError(f"frozen Model A archive not found: {archive_path}")
    archive = np.load(Path(archive_path), allow_pickle=False)
    molecule_ids = [str(x) for x in archive["molecule_ids"]]
    pocket_index = {str(x): i for i, x in enumerate(archive["pocket_ids"])}
    missing = [p for p in pocket_ids if p not in pocket_index]
    if missing:
        raise ModelBError(f"archive lacks declared M4 pockets: {missing}")
    rows = [pocket_index[p] for p in pocket_ids]
    mol_rep = archive["molecule_representations"].astype(np.float32)
    poc_rep = archive["pocket_representations"][rows].astype(np.float32)

    per_seed: dict[str, np.ndarray] = {}
    metadata = None
    for seed, path in zip(DEPLOYED_SEEDS, projection_paths):
        payload = projection.load_projection_file(Path(path))
        if metadata is None:
            metadata = {k: v for k, v in payload["metadata"].items() if isinstance(v, (str, int, float))}
        zm = projection.forward(mol_rep, payload["mol_project"])
        zp = projection.forward(poc_rep, payload["pocket_project"])
        matrix = zm @ zp.T
        per_seed[str(seed)] = matrix.max(axis=1).astype(np.float64)

    ensemble = np.mean(np.stack([per_seed[str(s)] for s in DEPLOYED_SEEDS]), axis=0)
    order = np.argsort(-ensemble, kind="mergesort")
    ranks = np.empty(len(ensemble), dtype=np.int64)
    ranks[order] = np.arange(1, len(ensemble) + 1)

    return {
        "schema": "pacer.drugclip_freeze.model_b_scores.v1",
        "archive": str(archive_path),
        "projection_files": [str(p) for p in projection_paths],
        "projection_metadata": metadata,
        "seeds": list(DEPLOYED_SEEDS),
        "pocket_ids": list(pocket_ids),
        "molecule_ids": molecule_ids,
        "per_seed_scores": {k: [float(x) for x in v] for k, v in per_seed.items()},
        "ensemble_score": [float(x) for x in ensemble],
        "new_rank": [int(x) for x in ranks],
        "pooling": "max over declared M4 pockets (official DrugCLIP semantics)",
        "claim_boundary": ("Frozen projection-only re-scoring of frozen DrugCLIP representations; "
                           "second-opinion evidence, not a PAM classifier and not a fused score."),
    }


def write_csv(report: dict[str, Any], path: Path) -> Path:
    import csv

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    seeds = [str(s) for s in report["seeds"]]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["candidate_smiles", "famaug_m4_score", "new_rank"] + [f"seed{s}" for s in seeds])
        for index, smiles in enumerate(report["molecule_ids"]):
            writer.writerow([smiles, report["ensemble_score"][index], report["new_rank"][index]]
                            + [report["per_seed_scores"][s][index] for s in seeds])
    return path
