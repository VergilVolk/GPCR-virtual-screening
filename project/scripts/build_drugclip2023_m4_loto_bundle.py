#!/usr/bin/env python3
"""Build a base-checkpoint-free PACER-200 M4-LOTO inference bundle."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--candidate-representations", type=Path, required=True)
    p.add_argument("--candidate-table", type=Path, required=True)
    p.add_argument("--reference-representations", type=Path, required=True)
    p.add_argument("--adapters-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    candidates = np.load(a.candidate_representations, allow_pickle=False)
    reference = np.load(a.reference_representations, allow_pickle=False)
    table = pd.read_csv(a.candidate_table)
    if len(table) != 200 or table.candidate_id.nunique() != 200:
        raise ValueError("Frozen PACER candidate table must contain 200 unique candidate IDs")
    rep_index = {str(v): i for i, v in enumerate(candidates["molecule_ids"])}
    missing = [s for s in table.canonical_smiles.astype(str) if s not in rep_index]
    if missing:
        raise ValueError(f"{len(missing)} candidates lack 2023 representations")
    order = [rep_index[s] for s in table.canonical_smiles.astype(str)]
    pocket_ids = list(map(str, reference["pocket_ids"]))
    pocket_index = pocket_ids.index("M4R_cluster0")
    a.output_dir.mkdir(parents=True, exist_ok=True)
    representation_name = "pacer200_m4_2023_frozen_representations.npz"
    representation_path = a.output_dir / representation_name
    np.savez_compressed(
        representation_path,
        candidate_ids=np.asarray(table.candidate_id.astype(str).tolist(), dtype=str),
        canonical_smiles=np.asarray(table.canonical_smiles.astype(str).tolist(), dtype=str),
        molecule_representations=candidates["molecule_representations"][order].astype(np.float32),
        m4_pocket_representation=reference["pocket_representations"][pocket_index].astype(np.float32),
    )
    adapter_manifest = json.loads((a.adapters_dir / "adapter_manifest.json").read_text(encoding="utf-8"))
    adapters = []
    for name, record in adapter_manifest["artifacts"].items():
        source = a.adapters_dir / name
        target = a.output_dir / name
        shutil.copy2(source, target)
        actual = sha256(target)
        if actual != record["sha256"]:
            raise ValueError(f"Copied adapter checksum mismatch: {name}")
        adapters.append({"file": name, "sha256": actual, "bytes": target.stat().st_size})
    manifest = {
        "protocol_id": "drugclip2023_m4_heldout_loto_3seed_v01",
        "scope": "frozen PACER-200 M4 candidates only",
        "base_checkpoint_required": False,
        "drugclip_checkout_required": False,
        "representations_file": representation_name,
        "representations_sha256": sha256(representation_path),
        "adapters": adapters,
        "candidate_rows": 200,
        "target": "M4R",
        "molecule_representation_shape": [200, 512],
        "pocket_representation_shape": [512],
        "adapter_validation": adapter_manifest["reference_reproduction"],
        "source_sha256": {
            "candidate_representations": sha256(a.candidate_representations),
            "candidate_table": sha256(a.candidate_table),
            "reference_representations": sha256(a.reference_representations),
        },
        "claim_boundary": "Binding-candidate retrieval only; not PAM function or potency.",
    }
    (a.output_dir / "bundle_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
