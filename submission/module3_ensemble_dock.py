from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent


def summarize(input_path: Path, output_path: Path) -> pd.DataFrame:
    table = pd.read_csv(input_path)
    required = {"candidate_id", "canonical_smiles", "cluster_coverage_fraction",
                "best_vina", "mean_vina", "structural_gate_pass"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Missing docking fields: {sorted(missing)}")
    table = table.sort_values(
        ["structural_gate_pass", "mean_vina"], ascending=[False, True]
    ).reset_index(drop=True)
    table.insert(0, "structural_rank", range(1, len(table) + 1))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description="Rank frozen multi-conformation docking evidence")
    parser.add_argument("--input", type=Path, default=ROOT / "data/demo/ensemble_docking_summary.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "results/module3_structure.csv")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    result = summarize(args.input, args.output)
    print(f"Structural gate: {int(result.structural_gate_pass.sum())}/{len(result)} passed")
    print(args.output)


if __name__ == "__main__":
    main()
