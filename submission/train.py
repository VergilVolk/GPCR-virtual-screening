from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description="Train the GPCR-adapted DrugCLIP projection ensemble")
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "models/retrained/model.json")
    parser.add_argument("--seeds", default="20260925,20260926,20260927")
    args = parser.parse_args()
    script = ROOT / "src/production/finetune_drugclip_gpcr_retrieval.py"
    command = [sys.executable, str(script), "--representations", str(args.representations),
               "--projection", str(args.projection), "--benchmark", str(args.benchmark),
               "--output", str(args.output), "--seeds", args.seeds]
    subprocess.run(command, cwd=script.parent, check=True)


if __name__ == "__main__":
    main()
