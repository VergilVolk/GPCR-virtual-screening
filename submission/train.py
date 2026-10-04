from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description="Rebuild the M4-held-out DrugCLIP adapter ensemble")
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--reference-predictions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "models/retrained")
    parser.add_argument("--seeds", default="20260925,20260926,20260927")
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    script = ROOT / "src/production/export_drugclip2023_m4_loto_adapters.py"
    command = [sys.executable, str(script),
               "--representations", str(args.representations.resolve()),
               "--projection", str(args.projection.resolve()),
               "--pairs", str(args.pairs.resolve()),
               "--reference-predictions", str(args.reference_predictions.resolve()),
               "--output-dir", str(output), "--seeds", args.seeds]
    subprocess.run(command, cwd=script.parent, check=True)


if __name__ == "__main__":
    main()
