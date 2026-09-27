#!/usr/bin/env python
"""Extract frozen Geom2Vec ViSNet embeddings from PACER-DC atom14 windows.

The input is the existing ``(frames, residues, 14, 3)`` heavy-atom tensor
produced by the OneProt/MDGen preprocessing stage.  No MD is rerun.  The
output retains per-frame, per-residue invariant features so that downstream
four-context contrasts can be pooled over predeclared biological regions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import torch

from geom2vec import create_model_from_checkpoint


ATOM14 = {
    "A": "N CA C O CB".split(),
    "R": "N CA C O CB CG CD NE CZ NH1 NH2".split(),
    "N": "N CA C O CB CG OD1 ND2".split(),
    "D": "N CA C O CB CG OD1 OD2".split(),
    "C": "N CA C O CB SG".split(),
    "Q": "N CA C O CB CG CD OE1 NE2".split(),
    "E": "N CA C O CB CG CD OE1 OE2".split(),
    "G": "N CA C O".split(),
    "H": "N CA C O CB CG ND1 CD2 CE1 NE2".split(),
    "I": "N CA C O CB CG1 CG2 CD1".split(),
    "L": "N CA C O CB CG CD1 CD2".split(),
    "K": "N CA C O CB CG CD CE NZ".split(),
    "M": "N CA C O CB CG SD CE".split(),
    "F": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ".split(),
    "P": "N CA C O CB CG CD".split(),
    "S": "N CA C O CB OG".split(),
    "T": "N CA C O CB OG1 CG2".split(),
    "W": "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2".split(),
    "Y": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH".split(),
    "V": "N CA C O CB CG1 CG2".split(),
}
ATOMIC_NUMBER = {"C": 6, "N": 7, "O": 8, "S": 16}


def install_cpu_radius_graph_fallback() -> bool:
    """Install a deterministic CPU radius-graph backend when PyG lacks wheels.

    Windows CPU environments commonly lack ``pyg-lib``/``torch-cluster`` for
    the newest PyTorch release.  The fallback is mathematically equivalent for
    inference, but slower; Linux/CUDA runs keep the official compiled backend.
    """
    from torch_geometric.nn import radius_graph as pyg_radius_graph
    try:
        pyg_radius_graph(
            torch.zeros((2, 3)), r=1.0, batch=torch.zeros(2, dtype=torch.long),
            max_num_neighbors=1,
        )
        return False
    except ImportError:
        pass

    from geom2vec.models.representation import visnet

    def radius_graph(x, r, batch=None, loop=False, max_num_neighbors=32, flow="source_to_target"):
        del flow
        if batch is None:
            batch = torch.zeros(x.shape[0], dtype=torch.long, device=x.device)
        sources, targets = [], []
        for graph_id in torch.unique(batch, sorted=True):
            idx = torch.nonzero(batch == graph_id, as_tuple=False).flatten()
            local = x[idx]
            distances = torch.cdist(local, local)
            if not loop:
                distances.fill_diagonal_(float("inf"))
            k = min(int(max_num_neighbors), local.shape[0] - (0 if loop else 1))
            values, neighbours = torch.topk(distances, k=k, dim=1, largest=False, sorted=True)
            target = torch.arange(local.shape[0], device=x.device)[:, None].expand_as(neighbours)
            valid = values <= float(r)
            sources.append(idx[neighbours[valid]])
            targets.append(idx[target[valid]])
        return torch.stack((torch.cat(sources), torch.cat(targets)), dim=0)

    visnet.radius_graph = radius_graph
    return True


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_sequence(path: Path) -> str:
    with path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    return row["seqres"].strip()


def flatten_atom14(coords: np.ndarray, sequence: str):
    if coords.ndim != 4 or coords.shape[2:] != (14, 3):
        raise ValueError(f"expected (frames,residues,14,3), got {coords.shape}")
    if coords.shape[1] != len(sequence):
        raise ValueError(f"sequence length {len(sequence)} != residues {coords.shape[1]}")

    atom_z, residue_index, slots = [], [], []
    for ridx, aa in enumerate(sequence):
        names = ATOM14.get(aa)
        if names is None:
            raise ValueError(f"unsupported residue {aa!r} at index {ridx}")
        for slot, name in enumerate(names):
            atom_z.append(ATOMIC_NUMBER[name[0]])
            residue_index.append(ridx)
            slots.append((ridx, slot))

    xyz = np.stack([coords[:, r, s, :] for r, s in slots], axis=1).astype(np.float32)
    if not np.isfinite(xyz).all():
        raise ValueError("non-finite coordinates")
    return xyz, np.asarray(atom_z, dtype=np.int64), np.asarray(residue_index, dtype=np.int64)


def invariant_residue_features(x, v, residue_index, n_residues):
    """Combine scalar channels and vector norms, then mean-pool by residue."""
    inv = torch.cat((x, torch.linalg.vector_norm(v, dim=1)), dim=-1)
    out = torch.zeros((n_residues, inv.shape[-1]), dtype=inv.dtype, device=inv.device)
    count = torch.zeros((n_residues, 1), dtype=inv.dtype, device=inv.device)
    out.index_add_(0, residue_index, inv)
    count.index_add_(0, residue_index, torch.ones((len(residue_index), 1), device=inv.device))
    return out / count.clamp_min(1.0)


def embed(model, xyz, atom_z, residue_index, n_residues, batch_size, device):
    outputs = []
    z0 = torch.as_tensor(atom_z, dtype=torch.long, device=device)
    ridx0 = torch.as_tensor(residue_index, dtype=torch.long, device=device)
    model.eval()
    with torch.inference_mode():
        for start in range(0, xyz.shape[0], batch_size):
            pos = torch.as_tensor(xyz[start : start + batch_size], dtype=torch.float32, device=device)
            bsz, natoms, _ = pos.shape
            z = z0.repeat(bsz)
            batch = torch.arange(bsz, device=device).repeat_interleave(natoms)
            x, v, _ = model(z=z, pos=pos.reshape(-1, 3), batch=batch)
            x = x.reshape(bsz, natoms, -1)
            v = v.reshape(bsz, natoms, 3, -1)
            for i in range(bsz):
                outputs.append(invariant_residue_features(x[i], v[i], ridx0, n_residues).cpu())
    return torch.stack(outputs).numpy().astype(np.float32)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atom14", type=Path, required=True)
    parser.add_argument("--sequence-csv", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()

    sequence = read_sequence(args.sequence_csv)
    original = np.load(args.atom14, mmap_mode="r")
    frame_ids = np.arange(0, original.shape[0], args.stride, dtype=np.int64)
    xyz, atom_z, residue_index = flatten_atom14(np.asarray(original[frame_ids]), sequence)
    used_radius_fallback = install_cpu_radius_graph_fallback()
    model = create_model_from_checkpoint(str(args.checkpoint), device=args.device)
    started = time.time()
    features = embed(model, xyz, atom_z, residue_index, len(sequence), args.batch_size, args.device)
    elapsed = time.time() - started

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        residue_features=features,
        global_features=features.mean(axis=1),
        frame_ids=frame_ids,
        residue_index=np.arange(len(sequence), dtype=np.int64),
        sequence=np.asarray(sequence),
    )
    audit = {
        "method": "frozen Geom2Vec ViSNet atom14 heavy-atom inference",
        "evidence_level": "encoder_feasibility_pilot",
        "atom14": str(args.atom14),
        "atom14_sha256": sha256(args.atom14),
        "checkpoint": str(args.checkpoint),
        "checkpoint_sha256": sha256(args.checkpoint),
        "frames_total": int(original.shape[0]),
        "frames_encoded": int(len(frame_ids)),
        "stride": args.stride,
        "residues": len(sequence),
        "heavy_atoms": int(len(atom_z)),
        "feature_shape": list(features.shape),
        "finite": bool(np.isfinite(features).all()),
        "elapsed_seconds": round(elapsed, 3),
        "seconds_per_frame": round(elapsed / len(frame_ids), 3),
        "cpu_radius_graph_fallback": used_radius_fallback,
        "claim_boundary": "Feasibility only; no PAM efficacy or cross-replica stability claim.",
    }
    args.output.with_suffix(".audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
