#!/usr/bin/env python3
"""Offline scoring for the packaged PACER-200 2023 M4-LOTO bundle.

No DrugCLIP checkout or base checkpoint is required.  Only PyTorch, NumPy and
pandas are needed because the bundle contains frozen 512-D representations and
three materialized projection adapters.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def project(representation: torch.Tensor, state: dict[str, torch.Tensor]) -> torch.Tensor:
    hidden = F.relu(F.linear(representation, state["linear1.weight"], state["linear1.bias"]))
    return F.normalize(F.linear(hidden, state["linear2.weight"], state["linear2.bias"]), dim=-1)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    manifest_path = args.bundle_dir / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    representation_path = args.bundle_dir / manifest["representations_file"]
    if sha256(representation_path) != manifest["representations_sha256"]:
        raise ValueError(f"Representation checksum mismatch: {representation_path.name}")
    archive = np.load(representation_path, allow_pickle=False)
    required = {"candidate_ids", "canonical_smiles", "molecule_representations", "m4_pocket_representation"}
    if required - set(archive.files):
        raise ValueError(f"Incomplete inference archive: {sorted(required - set(archive.files))}")
    molecule_rep = torch.as_tensor(archive["molecule_representations"].astype(np.float32))
    pocket_rep = torch.as_tensor(archive["m4_pocket_representation"].astype(np.float32)).reshape(1, -1)
    if molecule_rep.ndim != 2 or molecule_rep.shape[1] != 512 or pocket_rep.shape != (1, 512):
        raise ValueError("Expected molecule N x 512 and pocket 1 x 512 representations")
    scores, seed_columns = [], []
    artifact_hashes = {representation_path.name: sha256(representation_path)}
    for record in manifest["adapters"]:
        path = args.bundle_dir / record["file"]
        if sha256(path) != record["sha256"]:
            raise ValueError(f"Adapter checksum mismatch: {path.name}")
        artifact_hashes[path.name] = record["sha256"]
        payload = torch.load(path, map_location="cpu", weights_only=False)
        with torch.inference_mode():
            zm = project(molecule_rep, payload["mol_project"])
            zp = project(pocket_rep, payload["pocket_project"])
            value = (zm @ zp.T).squeeze(1).numpy()
        column = f"seed{payload['base_seed']}"
        seed_columns.append(column); scores.append(value)
    matrix = np.stack(scores, axis=1)
    output = pd.DataFrame({
        "pair_id": archive["candidate_ids"].astype(str),
        "target": "M4R",
        "canonical_smiles": archive["canonical_smiles"].astype(str),
    })
    for i, column in enumerate(seed_columns): output[column] = matrix[:, i]
    output["score_2023_gpcr_loto"] = matrix.mean(axis=1)
    output["score_2023_seed_sd"] = matrix.std(axis=1, ddof=0)
    output = output.sort_values("score_2023_gpcr_loto", ascending=False, kind="mergesort").reset_index(drop=True)
    output["rank_2023_gpcr_loto"] = np.arange(1, len(output) + 1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    audit = {
        "protocol_id": manifest["protocol_id"], "rows": int(len(output)),
        "target": "M4R", "output": str(args.output), "output_sha256": sha256(args.output),
        "artifact_sha256": artifact_hashes,
        "finite": bool(np.isfinite(matrix).all()),
        "score_range": [float(matrix.mean(1).min()), float(matrix.mean(1).max())],
        "claim_boundary": "Binding-candidate retrieval only; not PAM function or potency.",
    }
    args.output.with_suffix(".audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
