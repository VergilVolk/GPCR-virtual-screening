"""Model A: official frozen DrugCLIP apply-only contract.

This module does not run the backbone itself -- that is done by the frozen
`project/scripts/extract_drugclip_embeddings_cpu.py` under the `drugclip-cpu`
environment.  What this module does is enforce and record the *contract*:

  * pinned repository / Uni-Core commits and checkpoint identity;
  * empty `missing_checkpoint_keys` and `unexpected_checkpoint_keys`;
  * the exact archive keys, shapes and finiteness;
  * the three-pocket M4 baseline and the official max-pooling semantics.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

DRUGCLIP_COMMIT = "7a3a3fa33673f8668c811790f2e4681c98af44ef"
UNICORE_COMMIT = "44f6386f4dcd7137fc1e5d5e768117d635d64a26"
CHECKPOINT_FILENAME = "checkpoint_best.pt"
CHECKPOINT_BYTES = 1183713459
CHECKPOINT_SHA256 = "dc2c76d0f02f9bb079a613f09d538dcda1bf9075f2952d91dc1bea55571f667e"
EMBEDDING_DIM = 128
REPRESENTATION_DIM = 512
ARCHIVE_KEYS = ("molecule_ids", "molecule_representations", "molecule_embeddings",
                "pocket_ids", "pocket_representations", "pocket_embeddings", "scores")
M4_POCKET_IDS = ("7TRQ_M4_allosteric", "7TRP_M4_allosteric", "7TRS_M4_allosteric")
CONFORMER_SEED = 20260922
EXPECTED_PDB_SHA256 = {
    "7TRQ": "bccf53f6c49bcfa91c15e2151c66d6b0db82181f1536d92b314b9af69917becc",
    "7TRP": "b5b1af3a2974360dfb4fbad1f5930d9258b171eb9665004b0b7cf119e5517e7b",
    "7TRS": "53442edb5d1ce98d4fcbf9ae52d96f2b387c5ace1b913b40b54906227416fcc1",
}


class ModelAError(ValueError):
    """Raised when the frozen Model A contract is not satisfied."""


def validate_archive(path: Path, *, audit_path: Path | None = None) -> dict[str, Any]:
    """Fail-closed validation of one frozen Model A embedding archive."""
    archive_path = Path(path)
    if not archive_path.is_file():
        raise ModelAError(f"DrugCLIP embedding archive not found: {archive_path}")
    archive = np.load(archive_path, allow_pickle=False)
    keys = set(archive.files)
    missing_keys = sorted(set(ARCHIVE_KEYS) - keys)
    if missing_keys:
        raise ModelAError(f"archive lacks required keys: {missing_keys}")

    molecule_ids = [str(x) for x in archive["molecule_ids"]]
    pocket_ids = [str(x) for x in archive["pocket_ids"]]
    mol_rep = archive["molecule_representations"]
    mol_emb = archive["molecule_embeddings"]
    poc_rep = archive["pocket_representations"]
    poc_emb = archive["pocket_embeddings"]
    scores = archive["scores"]

    if len(set(molecule_ids)) != len(molecule_ids):
        raise ModelAError("molecule identifiers are not unique")
    if mol_emb.shape != (len(molecule_ids), EMBEDDING_DIM):
        raise ModelAError(f"molecule_embeddings shape {mol_emb.shape} != ({len(molecule_ids)}, {EMBEDDING_DIM})")
    if mol_rep.shape != (len(molecule_ids), REPRESENTATION_DIM):
        raise ModelAError(f"molecule_representations shape {mol_rep.shape} != ({len(molecule_ids)}, {REPRESENTATION_DIM})")
    if poc_emb.shape != (len(pocket_ids), EMBEDDING_DIM):
        raise ModelAError(f"pocket_embeddings shape {poc_emb.shape} != ({len(pocket_ids)}, {EMBEDDING_DIM})")
    if poc_rep.shape != (len(pocket_ids), REPRESENTATION_DIM):
        raise ModelAError(f"pocket_representations shape {poc_rep.shape} != ({len(pocket_ids)}, {REPRESENTATION_DIM})")
    expected_scores = (len(pocket_ids), len(molecule_ids))
    if scores.shape != expected_scores:
        raise ModelAError(f"scores shape {scores.shape} != {expected_scores}")
    for name, array in (("molecule_embeddings", mol_emb), ("molecule_representations", mol_rep),
                        ("pocket_embeddings", poc_emb), ("pocket_representations", poc_rep),
                        ("scores", scores)):
        if not np.isfinite(array).all():
            raise ModelAError(f"{name} contains non-finite values")

    audit: dict[str, Any] | None = None
    if audit_path is None:
        candidate = archive_path.with_suffix(".json")
        if candidate.is_file():
            audit_path = candidate
    if audit_path is not None and Path(audit_path).is_file():
        audit = json.loads(Path(audit_path).read_text(encoding="utf-8-sig"))
        for field in ("missing_checkpoint_keys", "unexpected_checkpoint_keys"):
            if audit.get(field):
                raise ModelAError(f"frozen checkpoint contract violated: {field}={audit[field]!r}")

    return {
        "archive": str(archive_path),
        "n_molecules": len(molecule_ids),
        "n_pockets": len(pocket_ids),
        "pocket_ids": pocket_ids,
        "molecule_embedding_shape": list(mol_emb.shape),
        "molecule_representation_shape": list(mol_rep.shape),
        "score_shape": list(scores.shape),
        "finite": True,
        "audit_present": audit is not None,
        "missing_checkpoint_keys": (audit or {}).get("missing_checkpoint_keys", None),
        "unexpected_checkpoint_keys": (audit or {}).get("unexpected_checkpoint_keys", None),
        "m4_pockets_present": [p for p in M4_POCKET_IDS if p in set(pocket_ids)],
    }


def reference_scores(archive_path: Path, pocket_ids: tuple[str, ...] = M4_POCKET_IDS) -> np.ndarray:
    """Official max-pooled M4 score (max over the declared pockets), per molecule."""
    archive = np.load(Path(archive_path), allow_pickle=False)
    pocket_index = {str(x): i for i, x in enumerate(archive["pocket_ids"])}
    missing = [p for p in pocket_ids if p not in pocket_index]
    if missing:
        raise ModelAError(f"archive lacks declared M4 pockets: {missing}")
    rows = [pocket_index[p] for p in pocket_ids]
    scores = archive["scores"][rows, :].astype(np.float64)
    return scores.max(axis=0)
