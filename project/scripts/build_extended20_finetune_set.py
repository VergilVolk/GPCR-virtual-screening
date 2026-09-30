#!/usr/bin/env python3
"""Build the extended external npz + pairs for fine-tuning scale-up:
pilot 9-target npz (has representations) + 7 subset-reps npz (ADRB2, ESR1_ago,
ESR1_ant, IDH1, OPRK1, PPARG, TP53). Subset molecule_ids are 'SMILES ZINCID';
canonicalize the SMILES for the join key. Pockets: mean over structures per
target. Pairs: all actives + seeded decoy fill to per-target cap (matching the
pilot's per-target scale)."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def canon(smi: str) -> str:
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol) if mol else smi

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pilot-npz", type=Path, required=True)
    p.add_argument("--subset-root", type=Path, required=True)
    p.add_argument("--targets", required=True)
    p.add_argument("--per-target-cap", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20261001)
    p.add_argument("--output-npz", type=Path, required=True)
    p.add_argument("--output-pairs", type=Path, required=True)
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)
    pilot = np.load(a.pilot_npz, allow_pickle=False)
    mols_rep = [pilot["molecule_representations"].astype(np.float32)]
    mol_ids = [np.asarray([str(x) for x in pilot["molecule_ids"]])]
    pocket_ids = [str(x) for x in pilot["pocket_ids"]]
    pocket_rep = [pilot["pocket_representations"].astype(np.float32)]
    pairs_rows = []
    pilot_labels = None
    targets = [t.strip() for t in a.targets.split(",") if t.strip()]
    for t in targets:
        d = np.load(a.subset_root / f"{t}.npz", allow_pickle=False)
        ids = np.asarray([str(x).split(" ")[0] for x in d["molecule_ids"]])
        canon_ids = np.asarray([canon(s) for s in ids])
        labels = np.asarray(d["labels"], dtype=np.int64)
        rep = d["molecule_representations"].astype(np.float32)
        # per-target balanced sampling: all actives + decoy fill
        act = np.flatnonzero(labels == 1)
        dec = np.flatnonzero(labels == 0)
        rng.shuffle(dec)
        keep = np.sort(np.concatenate([act, dec[: max(0, a.per_target_cap - len(act))]]))
        mols_rep.append(rep[keep])
        mol_ids.append(canon_ids[keep])
        pocket_ids.append(t)
        pocket_rep.append(d["pocket_representations"].astype(np.float32).mean(axis=0, keepdims=True))
        for i in keep:
            pairs_rows.append((f"LITPCBA:{t}:{i}", t, str(ids[i]), str(canon_ids[i]), int(labels[i])))
    out = {
        "molecule_ids": np.concatenate(mol_ids),
        "molecule_representations": np.concatenate(mols_rep),
        "pocket_ids": np.asarray(pocket_ids),
        "pocket_representations": np.concatenate(pocket_rep),
    }
    np.savez_compressed(a.output_npz, **out)
    with a.output_pairs.open("w", encoding="utf-8") as f:
        f.write("pair_id,target,molecule_name,canonical_smiles,label\n")
        for r in pairs_rows:
            f.write("%s,%s,%s,%s,%d\n" % r)
    print(json.dumps({"n_external_targets": len(pocket_ids), "n_pairs_external": len(pairs_rows),
                      "n_molecules_total": int(len(out["molecule_ids"]))}, indent=2))

if __name__ == "__main__":
    main()
