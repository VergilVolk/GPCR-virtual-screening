#!/usr/bin/env python3
"""Convert published ten-cluster GaMD GPCR ensembles into DrugCLIP LMDBs."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import lmdb
import numpy as np
import pandas as pd
import pickle


TARGETS = ["B2AR", "CCR2", "M2R", "M4R"]


def parse_resids(path: Path) -> set[int]:
    text = path.read_text(encoding="utf-8")
    matches = re.findall(r"rms reference mass out rmsd-pocket\.dat :([0-9,]+)", text)
    if not matches:
        raise ValueError(f"Could not find pocket residue definition in {path}")
    return {int(value) for value in matches[-1].split(",")}


def parse_models(path: Path, resids: set[int]):
    records, current, model_number = [], None, None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            model_number = int(line.split()[-1]); current = {"atoms": [], "coordinates": []}
        elif line.startswith("ENDMDL"):
            if current is not None:
                records.append((model_number, current))
            current = None
        elif current is not None and line.startswith("ATOM"):
            try:
                resid = int(line[22:26])
                xyz = np.asarray([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=np.float32)
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            element = (line[76:78].strip() or atom_name[:1]).upper()
            if resid in resids and element != "H":
                current["atoms"].append(atom_name)
                current["coordinates"].append(xyz)
    return records


def populations(path: Path):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        rows.append({"cluster": int(parts[0]), "frames": int(parts[1]), "population": float(parts[2]),
                     "avg_dist": float(parts[3]), "centroid_frame": int(parts[5])})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ensembles", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--max-atoms", type=int, default=256)
    args = parser.parse_args()
    records, metadata, audit = [], [], {}
    for target in TARGETS:
        root = args.ensembles / target
        resids = parse_resids(root / f"{target}_cpptraj_cluster.in")
        models = parse_models(root / f"{target}_ensemble_clusters.pdb", resids)
        pop = populations(root / f"{target}_ensemble_clustering_summary.dat")
        if len(models) != 10 or len(pop) != 10:
            raise ValueError(f"{target}: expected 10 models/populations; got {len(models)}/{len(pop)}")
        for (model_number, pocket), meta in zip(models, pop):
            if model_number != meta["cluster"] + 1:
                raise ValueError(f"{target}: model/cluster mismatch {model_number}/{meta['cluster']}")
            if not pocket["atoms"] or len(pocket["atoms"]) > args.max_atoms:
                raise ValueError(f"{target} cluster {meta['cluster']}: atom count {len(pocket['atoms'])}")
            name = f"{target}_cluster{meta['cluster']}"
            records.append({"pocket": name, "pocket_atoms": pocket["atoms"],
                            "pocket_coordinates": pocket["coordinates"]})
            metadata.append({"pocket_id": name, "target": target, "cluster": meta["cluster"],
                             "population": meta["population"], "frames": meta["frames"],
                             "centroid_frame": meta["centroid_frame"], "atom_count": len(pocket["atoms"]),
                             "residue_count": len(resids)})
        audit[target] = {"clusters": len(models), "residues": sorted(resids),
                         "population_sum": float(sum(row["population"] for row in pop)),
                         "atom_counts": [len(pocket["atoms"]) for _, pocket in models]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(args.output), subdir=False, map_size=2**30)
    with env.begin(write=True) as txn:
        for index, record in enumerate(records):
            txn.put(str(index).encode("ascii"), pickle.dumps(record))
    env.sync(); env.close()
    pd.DataFrame(metadata).to_csv(args.metadata, index=False)
    args.output.with_suffix(args.output.suffix + ".audit.json").write_text(
        json.dumps({"records": len(records), "targets": audit,
                    "claim_boundary": "Published GaMD cluster representatives; no learned selection."}, indent=2),
        encoding="utf-8")
    print(json.dumps({"records": len(records), "targets": audit}, indent=2))


if __name__ == "__main__":
    main()
