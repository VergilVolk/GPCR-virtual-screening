from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description="Launch the documented DrugCLIP adapter training")
    parser.add_argument("--repo-root", type=Path, default=ROOT.parent)
    parser.add_argument("--representations", type=Path, default=Path("data/training/representations.pt"))
    parser.add_argument("--projection", type=Path, default=Path("data/training/base_projection.pt"))
    parser.add_argument("--contrasts", type=Path, default=Path("data/training/contrasts.csv"))
    parser.add_argument("--output", type=Path, default=ROOT / "models/retrained")
    parser.add_argument("--seeds", default="20260925,20260926,20260927")
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    script = args.repo_root / "project/scripts/finetune_drugclip_muscarinic_triplet.py"
    if not script.exists():
        raise FileNotFoundError(
            "The full training script is not included in the compact archive. "
            "Run this entry from the complete repository."
        )
    command = [
        sys.executable, str(script),
        "--representations", str(args.representations),
        "--projection", str(args.projection),
        "--contrasts", str(args.contrasts),
        "--output", str(args.output),
        "--seeds", args.seeds,
        "--epochs", str(args.epochs),
    ]
    manifest = {
        "entrypoint": str(script), "seeds": args.seeds, "epochs": args.epochs,
        "note": "Target and scaffold separation must follow the frozen protocol in logs/training_manifest.json",
    }
    print(json.dumps(manifest, indent=2))
    if not args.dry_run:
        subprocess.run(command, cwd=args.repo_root, check=True)


if __name__ == "__main__":
    main()
