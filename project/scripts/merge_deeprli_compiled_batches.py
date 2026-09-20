"""Merge the main and isotope-recovery DeepRLI graph ledgers with hard checks."""

import json
import pickle
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
MAIN = PROJECT / "data" / "deeprli_m4_v01"
EXTRA = PROJECT / "data" / "deeprli_m4_missing_v01"


def load(root):
    with (root / "compiled" / "complexes.processed.pkl").open("rb") as fh:
        graphs = pickle.load(fh)
    index = pd.read_csv(root / "index" / "complexes.processed.csv")
    return graphs, index


def main():
    main_graphs, main_index = load(MAIN); extra_graphs, extra_index = load(EXTRA)
    overlap = set(main_graphs) & set(extra_graphs)
    if overlap:
        raise RuntimeError(f"overlapping graph keys: {sorted(overlap)}")
    merged = {**main_graphs, **extra_graphs}
    index = pd.concat([main_index, extra_index], ignore_index=True).drop_duplicates("complex_path")
    index[["state", "canonical_molecule_id"]] = index.complex_path.str.split("/", expand=True)
    counts = index.groupby("canonical_molecule_id").state.nunique()
    if len(merged) != 1290 or len(index) != 1290 or len(counts) != 430 or not (counts == 3).all():
        raise RuntimeError(f"bad merged dimensions graphs={len(merged)} rows={len(index)} molecules={len(counts)}")
    index = index.sort_values(["state", "canonical_molecule_id"]).drop(columns=["state", "canonical_molecule_id"])
    with (MAIN / "compiled" / "complexes.processed.pkl").open("wb") as fh:
        pickle.dump(merged, fh)
    index.to_csv(MAIN / "index" / "complexes.processed.csv", index=False)
    audit = {"main_graphs": len(main_graphs), "recovered_isotope_graphs": len(extra_graphs),
             "merged_graphs": len(merged), "unique_molecules": len(counts),
             "states_per_molecule": 3, "overlap": 0}
    (MAIN / "merge_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
