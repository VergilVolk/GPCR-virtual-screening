#!/usr/bin/env python3
"""Assemble encoded DUD-E target files with auditable label alignment."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import pickle
from pathlib import Path

import lmdb
import numpy as np


def records(path: Path):
    env = lmdb.open(str(path), subdir=False, readonly=True, lock=False, readahead=False)
    with env.begin() as txn:
        result = [pickle.loads(txn.get(str(index).encode("ascii"))) for index in range(txn.stat()["entries"])]
    env.close()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--embedding-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    targets, pocket_embeddings, reference_embeddings = [], [], []
    pocket_offsets, reference_offsets = [0], [0]
    audit = {}
    for folder in sorted(path for path in args.raw_root.iterdir() if path.is_dir()):
        target = folder.name
        candidate = np.load(args.embedding_root / f"{target}_candidates.npz", allow_pickle=False)
        reference = np.load(args.embedding_root / f"{target}_reference.npz", allow_pickle=False)
        source_records = records(folder / "mols.lmdb")
        source_smiles = [str(record["smi"]) for record in source_records]
        extracted_smiles = candidate["molecule_ids"].astype(str).tolist()
        # Official DrugCLIP retrieval datasets may use a deterministic internal
        # ordering that differs from raw LMDB's numeric-key order.  Never
        # presume positional alignment: require exact multiset agreement, then
        # reconstruct labels by the explicit molecule identifier.
        source_counter, extracted_counter = Counter(source_smiles), Counter(extracted_smiles)
        if source_counter != extracted_counter:
            raise RuntimeError(f"Candidate SMILES multiset mismatch for {target}")
        label_lookup = defaultdict(set)
        for record in source_records:
            label_lookup[str(record["smi"])].add(int(record["label"]))
        conflicts = [smiles for smiles, values in label_lookup.items() if len(values) != 1]
        if conflicts:
            raise RuntimeError(f"Ambiguous duplicated-SMILES labels for {target}: {conflicts[:3]}")
        labels = np.asarray([next(iter(label_lookup[smiles])) for smiles in extracted_smiles], dtype=np.int64)
        out = args.output_root / target / "drugclip_emb"
        out.mkdir(parents=True, exist_ok=True)
        with (out / "mols.lmdb.pkl").open("wb") as handle:
            pickle.dump([candidate["molecule_embeddings"].astype(np.float32), extracted_smiles, labels], handle)
        targets.append(target)
        pocket_embeddings.append(candidate["pocket_embeddings"].astype(np.float32))
        reference_embeddings.append(reference["molecule_embeddings"].astype(np.float32))
        pocket_offsets.append(pocket_offsets[-1] + len(pocket_embeddings[-1]))
        reference_offsets.append(reference_offsets[-1] + len(reference_embeddings[-1]))
        audit[target] = {"candidates": int(len(labels)), "actives": int(labels.sum()),
                         "pockets": int(len(pocket_embeddings[-1])),
                         "references": int(len(reference_embeddings[-1])),
                         "label_alignment": "exact SMILES multiset match; labels explicitly reconstructed by SMILES",
                         "duplicate_smiles": int(sum(count - 1 for count in source_counter.values() if count > 1))}
    args.output_root.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output_root / "dude_gpcr_contexts.npz",
                        targets=np.asarray(targets),
                        pocket_embeddings=np.concatenate(pocket_embeddings),
                        pocket_offsets=np.asarray(pocket_offsets),
                        reference_embeddings=np.concatenate(reference_embeddings),
                        reference_offsets=np.asarray(reference_offsets))
    report = {"targets": targets, "audit": audit,
              "claim_boundary": "External DUD-E GPCR panel; candidate labels read only after frozen embedding extraction."}
    (args.output_root / "dude_gpcr_assembly.audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
