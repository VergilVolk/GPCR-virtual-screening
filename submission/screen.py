from __future__ import annotations

import argparse
from pathlib import Path

from src.screen_core import screen, write_score_table


ROOT = Path(__file__).resolve().parent
MODEL_TAG = "DrugCLIP-2023 GPCR-adapted projection ensemble, 3 seeds"


def model_paths() -> list[Path]:
    return sorted((ROOT / "models/used_module2").glob("drugclip2023_m4_loto_seed*.projection.pt"))


def run_screen(molecules: Path, pockets: Path, output: Path) -> Path:
    models = model_paths()
    if len(models) != 3:
        raise FileNotFoundError("Expected three packaged M4 adapter checkpoints")
    molecule_ids, scores, _ = screen(molecules, pockets, models, pooling="mean")
    return write_score_table(output, molecule_ids, scores, MODEL_TAG)


def main() -> None:
    parser = argparse.ArgumentParser(description="Score molecules against the M4 pocket")
    parser.add_argument("--molecules", type=Path, required=True)
    parser.add_argument("--pockets", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "results/module2_binding.csv")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    output = run_screen(args.molecules, args.pockets, args.output)
    print(output)


if __name__ == "__main__":
    main()
