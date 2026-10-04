from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
CHANNELS = ["Glide_PDB", "Glide_BEmin", "Glide_BEavg", "Vina_PDB", "Vina_BEmin", "Vina_BEavg"]

def fuse_six_channels(input_path: Path, output_path: Path) -> pd.DataFrame:
    table = pd.read_csv(input_path)
    ranks = [f"{name}_rankpct" for name in CHANNELS]
    raw = [f"{name}_raw" for name in CHANNELS]
    missing = {"candidate_id", *ranks, *raw}.difference(table.columns)
    if missing:
        raise ValueError(f"Missing six-channel fields: {sorted(missing)}")
    table[ranks + raw] = table[ranks + raw].apply(pd.to_numeric, errors="coerce")
    if table[ranks + raw].isna().any().any():
        raise ValueError("All six docking channels are required; no imputation is used")
    table["six_channel_complete"] = True
    table["pacer_xr_score"] = table[ranks].mean(axis=1)
    if "PACER_XR" in table:
        delta = float(np.max(np.abs(table["PACER_XR"] - table["pacer_xr_score"])))
        if delta > 1e-10:
            raise ValueError(f"Frozen score mismatch: {delta}")
    table = table.sort_values("pacer_xr_score", ascending=False, kind="stable").reset_index(drop=True)
    table["pacer_xr_rank"] = np.arange(1, len(table) + 1)
    table["claim_boundary"] = "six-channel binding evidence; not a PAM-function label"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    return table

def main() -> None:
    parser = argparse.ArgumentParser(description="Fuse Glide/Vina PDB, BEmin and BEavg ranks")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "results/module3_structure.csv")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    result = fuse_six_channels(args.input, args.output)
    print(result[["pacer_xr_rank", "candidate_id", "pacer_xr_score", "group"]].to_string(index=False))

if __name__ == "__main__":
    main()
