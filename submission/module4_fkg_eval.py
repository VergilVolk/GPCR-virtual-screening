from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent


def export_evidence(validation_path: Path, prospective_path: Path, output_path: Path) -> pd.DataFrame:
    validation = pd.read_csv(validation_path)
    prospective = pd.read_csv(prospective_path)
    validation["evidence_set"] = "retrospective_control"
    prospective["evidence_set"] = "prospective_candidate"
    result = pd.concat([validation, prospective], ignore_index=True, sort=False)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Export frozen four-context dynamics evidence")
    parser.add_argument("--validation", type=Path, default=ROOT / "data/demo/four_context_validation.csv")
    parser.add_argument("--prospective", type=Path, default=ROOT / "data/demo/prospective_four_context_summary.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "results/module4_function.csv")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    result = export_evidence(args.validation, args.prospective, args.output)
    print(f"Exported {len(result)} frozen evidence rows")
    print(args.output)


if __name__ == "__main__":
    main()
