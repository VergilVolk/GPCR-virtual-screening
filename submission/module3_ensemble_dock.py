from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
CHANNELS = [
    "Glide_PDB", "Glide_BEmin", "Glide_BEavg",
    "Vina_PDB", "Vina_BEmin", "Vina_BEavg",
]


def fuse_six_channels(input_path: Path, output_path: Path) -> pd.DataFrame:
    table = pd.read_csv(input_path)
    rank_columns = [f"{name}_rankpct" for name in CHANNELS]
    raw_columns = [f"{name}_raw" for name in CHANNELS]
    required = {"candidate_id", *rank_columns, *raw_columns}
    missing_columns = required.difference(table.columns)
    if missing_columns:
        raise ValueError(f"Missing six-channel fields: {sorted(missing_columns)}")

    table[rank_columns + raw_columns] = table[rank_columns + raw_columns].apply(
        pd.to_numeric, errors="coerce"
    )
    missing_scores = table[rank_columns + raw_columns].isna().sum()
    if missing_scores.any():
        raise ValueError(
            "All six protocol-matched channels are required; no imputation is allowed: "
            f"{missing_scores[missing_scores > 0].to_dict()}"
        )

    table["six_channel_complete"] = True
    table["pacer_xr_score"] = table[rank_columns].mean(axis=1)
    table["glide_bemin_top1pct"] = table["Glide_BEmin_rankpct"] >= 0.99
    table["cascade_stage"] = np.where(table["glide_bemin_top1pct"], 1, 2)
    table["cascade_stage_score"] = np.where(
        table["glide_bemin_top1pct"], table["Glide_BEmin_rankpct"], table["pacer_xr_score"]
    )
    table = table.sort_values("pacer_xr_score", ascending=False, kind="stable").reset_index(drop=True)
    table["pacer_xr_rank"] = np.arange(1, len(table) + 1)

    cascade = table.sort_values(
        ["cascade_stage", "cascade_stage_score"],
        ascending=[True, False], kind="stable",
    ).candidate_id.tolist()
    cascade_rank = {candidate_id: rank for rank, candidate_id in enumerate(cascade, 1)}
    table["candidate_cascade_rank"] = table.candidate_id.map(cascade_rank)
    table["claim_boundary"] = "binding retrieval evidence; PAM function requires module 4 and experiment"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description="Fuse Glide and Vina static/ensemble docking ranks")
    parser.add_argument("--input", type=Path, default=ROOT / "data/demo/six_channel_docking_summary.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "results/module3_structure.csv")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    result = fuse_six_channels(args.input, args.output)
    print(result[["pacer_xr_rank", "candidate_id", "pacer_xr_score"]].to_string(index=False))
    print(args.output)


if __name__ == "__main__":
    main()
