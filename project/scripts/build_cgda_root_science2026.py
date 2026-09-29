#!/usr/bin/env python3
"""Build a CGDA-compatible pickle root from Science-2026 full_litpcba npz outputs.

For each <target>.npz + <target>.references.npz in the full_litpcba directory,
write <root>/<target>/drugclip_emb/{mols.lmdb.pkl, ligand.lmdb.pkl} matching the
loader in evaluate_drugclip_cgda_loso.py (3-tuples: (array, None, labels) and
(array, None, None)). Also emit a combined pockets pickle mapping each target to
its pocket_embeddings block, consumed by run_cgda_science2026.py via a patched
pocket_map (bypassing the historical LMDB pocket-root layout).
No fitting happens here; this is a pure format converter.
"""
from __future__ import annotations
import argparse, pickle
from pathlib import Path
import numpy as np

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--full-litpcba", type=Path, required=True)
    p.add_argument("--output-root", type=Path, required=True)
    p.add_argument("--pockets-pkl", type=Path, required=True)
    a = p.parse_args()
    targets = sorted(x.name[: -len(".npz")] for x in a.full_litpcba.glob("*.npz")
                     if not x.name.endswith(".references.npz") and x.name[:-4] in
                     {y.name[:-len(".metrics.json")] for y in a.full_litpcba.glob("*.metrics.json")})
    if len(targets) < 15:
        raise SystemExit(f"only {len(targets)} completed targets; need 15")
    pockets = {}
    for t in targets:
        d = np.load(a.full_litpcba / f"{t}.npz", allow_pickle=False)
        r = np.load(a.full_litpcba / f"{t}.references.npz", allow_pickle=False)
        folder = a.output_root / t / "drugclip_emb"
        folder.mkdir(parents=True, exist_ok=True)
        mols = np.asarray(d["molecule_embeddings"], dtype=np.float32)
        labels = np.asarray(d["labels"], dtype=np.int64)
        refs = np.asarray(r["reference_embeddings"], dtype=np.float32)
        with (folder / "mols.lmdb.pkl").open("wb") as f:
            pickle.dump((mols, None, labels), f, protocol=4)
        with (folder / "ligand.lmdb.pkl").open("wb") as f:
            pickle.dump((refs, None, None), f, protocol=4)
        pockets[t] = np.asarray(d["pocket_embeddings"], dtype=np.float32)
    with a.pockets_pkl.open("wb") as f:
        pickle.dump(pockets, f, protocol=4)
    print(json.dumps({"targets": targets, "n_targets": len(targets),
                      "pockets_pkl": str(a.pockets_pkl)}, indent=2))

import json
if __name__ == "__main__":
    main()
