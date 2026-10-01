#!/usr/bin/env python3
"""Build the OLD-checkpoint CGDA root for the stratified-subset interim gate:
6 subset targets use the same deterministic row indices as the Science-2026
subset run; 8 completed targets stay full-size. Identical rows both sides."""
from __future__ import annotations
import argparse, pickle
from pathlib import Path
import numpy as np


def subset_indices(labels: np.ndarray, decoys: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    active = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    rng.shuffle(negative)
    return np.sort(np.concatenate([active, negative[:decoys]]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--old-root", type=Path, required=True)
    p.add_argument("--subset-npz-root", type=Path, required=True,
                   help="Science-2026 subset output dir (for row_indices cross-check)")
    p.add_argument("--output-root", type=Path, required=True)
    p.add_argument("--decoys-per-target", type=int, default=20000)
    p.add_argument("--subset-seed", type=int, default=20261001)
    a = p.parse_args()
    subset_targets = sorted(x.name[:-4] for x in a.subset_npz_root.glob("*.npz")
                            if not x.name.endswith(".references.npz"))
    for t in subset_targets:
        row_file = a.subset_npz_root / f"{t}.npz"
        keep_ref = np.load(row_file, allow_pickle=False)["row_indices"]
        folder = a.output_root / t / "drugclip_emb"
        folder.mkdir(parents=True, exist_ok=True)
        m, _, lab = pickle.load((a.old_root / t / "drugclip_emb" / "mols.lmdb.pkl").open("rb"))
        m, lab = np.asarray(m, dtype=np.float32), np.asarray(lab, dtype=np.int64)
        keep = subset_indices(lab, a.decoys_per_target, a.subset_seed)
        assert np.array_equal(keep, keep_ref), f"row index mismatch for {t}"
        with (folder / "mols.lmdb.pkl").open("wb") as f:
            pickle.dump((m[keep], None, lab[keep]), f, protocol=4)
        src = a.old_root / t / "drugclip_emb" / "ligand.lmdb.pkl"
        with src.open("rb") as f, (folder / "ligand.lmdb.pkl").open("wb") as g:
            g.write(f.read())
        print(f"{t}: subset rows={len(keep)} actives={int(lab[keep].sum())}")
    full_targets = ["ALDH1", "ESR1_ago", "ESR1_ant", "MAPK1", "MTORC1", "PKM2", "PPARG", "TP53"]
    for t in full_targets:
        src = a.old_root / t / "drugclip_emb"
        if not src.exists():
            continue
        folder = a.output_root / t / "drugclip_emb"
        folder.mkdir(parents=True, exist_ok=True)
        if not (folder / "mols.lmdb.pkl").exists():
            for name in ("mols.lmdb.pkl", "ligand.lmdb.pkl"):
                with (src / name).open("rb") as f, (folder / name).open("wb") as g:
                    g.write(f.read())
            print(f"{t}: full copy")


if __name__ == "__main__":
    main()
