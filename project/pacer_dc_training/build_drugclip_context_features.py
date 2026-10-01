#!/usr/bin/env python
"""Build PACER four-context features from frozen DrugCLIP embeddings.

Expected columns are metadata plus dc_p_000... pocket embeddings and
dc_m_000... molecule embeddings. Optional op_000... columns add OneProt-MD
trajectory embeddings. The output uses f_000... columns accepted by
train_dual_context_heads.py.

DrugCLIP is not treated as a PAM classifier. It contributes a pocket state,
a molecule-conditioned interaction feature, and a static binding score. The
downstream trainer forms d_PAM and d_AGO from matched contexts.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


CONTEXTS = {"candidate_probe", "candidate_no_probe", "probe_only", "apo"}
META = [
    "candidate_id", "chemotype", "replicate_id", "window_id", "context",
    "pam_label", "agonism_label", "split",
]


def numbered_columns(frame: pd.DataFrame, prefix: str) -> list[str]:
    cols = [c for c in frame.columns if c.startswith(prefix)]
    try:
        return sorted(cols, key=lambda c: int(c[len(prefix):]))
    except ValueError as exc:
        raise ValueError(f"Non-numeric feature suffix under {prefix}") from exc


def row_normalize(values: np.ndarray, name: str) -> np.ndarray:
    if not np.isfinite(values).all():
        raise ValueError(f"{name} contains non-finite values")
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    if np.any(norms <= 1e-12):
        raise ValueError(f"{name} contains an all-zero row")
    return values / norms


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--molecule-tolerance", type=float, default=1e-6)
    args = parser.parse_args()

    frame = pd.read_csv(args.input)
    missing = set(META) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing metadata columns: {sorted(missing)}")
    if not set(frame.context).issubset(CONTEXTS):
        raise ValueError(f"Unknown contexts: {sorted(set(frame.context) - CONTEXTS)}")

    pocket_cols = numbered_columns(frame, "dc_p_")
    molecule_cols = numbered_columns(frame, "dc_m_")
    oneprot_cols = numbered_columns(frame, "op_")
    if not pocket_cols or len(pocket_cols) != len(molecule_cols):
        raise ValueError("DrugCLIP pocket and molecule dimensions must be equal and non-zero")

    keys = ["candidate_id", "chemotype", "replicate_id", "window_id", "split"]
    incomplete: list[dict] = []
    molecule_drift: list[dict] = []
    for key, part in frame.groupby(keys, sort=False):
        seen = set(part.context)
        if seen != CONTEXTS or part.context.duplicated().any():
            incomplete.append({"key": list(map(str, key)), "contexts": sorted(seen)})
        mol = part[molecule_cols].to_numpy(float)
        drift = float(np.max(np.abs(mol - mol[:1])))
        if drift > args.molecule_tolerance:
            molecule_drift.append({"key": list(map(str, key)), "max_abs_drift": drift})
    if incomplete:
        raise ValueError(f"Incomplete four-context units: {incomplete[:3]}")
    if molecule_drift:
        raise ValueError(f"Molecule embedding changes across contexts: {molecule_drift[:3]}")

    pocket = row_normalize(frame[pocket_cols].to_numpy(np.float32), "DrugCLIP pocket")
    molecule = row_normalize(frame[molecule_cols].to_numpy(np.float32), "DrugCLIP molecule")
    interaction = pocket * molecule
    similarity = np.sum(pocket * molecule, axis=1, keepdims=True)

    blocks = []
    block_dims: dict[str, int] = {}
    if oneprot_cols:
        oneprot = row_normalize(frame[oneprot_cols].to_numpy(np.float32), "OneProt-MD")
        blocks.append(oneprot)
        block_dims["oneprot_md"] = oneprot.shape[1]
    blocks.extend([pocket, interaction, similarity])
    block_dims.update({
        "drugclip_pocket": pocket.shape[1],
        "drugclip_pocket_times_molecule": interaction.shape[1],
        "drugclip_cosine": 1,
    })
    features = np.concatenate(blocks, axis=1).astype(np.float32)

    candidate_leak = frame.groupby("candidate_id").split.nunique()
    chemotype_leak = frame.groupby("chemotype").split.nunique()
    leaking_candidates = sorted(map(str, candidate_leak[candidate_leak > 1].index))
    leaking_chemotypes = sorted(map(str, chemotype_leak[chemotype_leak > 1].index))
    if leaking_candidates or leaking_chemotypes:
        raise ValueError(
            "Split leakage detected before feature export: "
            f"candidates={leaking_candidates[:10]}, chemotypes={leaking_chemotypes[:10]}"
        )

    output = frame[META].copy()
    output["drugclip_binding_score"] = similarity[:, 0]
    feature_frame = pd.DataFrame(
        features,
        columns=[f"f_{i:04d}" for i in range(features.shape[1])],
        index=output.index,
    )
    output = pd.concat([output, feature_frame], axis=1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)

    audit = {
        "input": str(args.input),
        "output": str(args.output),
        "rows": int(len(output)),
        "four_context_units": int(frame.groupby(keys).ngroups),
        "drugclip_dimension": len(pocket_cols),
        "oneprot_dimension": len(oneprot_cols),
        "output_dimension": int(features.shape[1]),
        "feature_blocks": block_dims,
        "binding_score_range": [float(similarity.min()), float(similarity.max())],
        "candidate_split_leakage": [],
        "chemotype_split_leakage": [],
        "four_context_complete": True,
        "molecule_embedding_context_invariant": True,
        "claim_boundary": (
            "Prepared frozen representations only; DrugCLIP similarity is binding compatibility, "
            "not PAM efficacy, and no performance claim is produced by this script."
        ),
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
