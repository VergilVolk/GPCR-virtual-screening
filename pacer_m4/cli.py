"""Command-line entry point for the PACER-M4 source release."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .runner import repository_root, run_stage, stage_capabilities
from .stages import STAGES


def _print_stage_table(root: Path) -> None:
    print(f"PACER-M4 repository: {root}")
    print("\nRegistered stages:")
    for stage_id in sorted(STAGES):
        stage = STAGES[stage_id]
        marker = "ready" if (root / stage.script).is_file() else "missing-code"
        print(f"  {stage_id:27} {marker:12} {stage.description}")


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="pacer-m4",
        description="Discover and run audited PACER-M4 workflow stages.",
    )
    parser.add_argument("--root", type=Path, help="Repository root; defaults to PACER_M4_ROOT or checkout root")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("stages", help="List registered workflow stages")
    capabilities = subparsers.add_parser("capabilities", help="Emit machine-readable stage availability")
    capabilities.add_argument("--output", type=Path)

    run = subparsers.add_parser("run", help="Run one registered stage")
    run.add_argument("stage_id", choices=sorted(STAGES))
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--receipt", type=Path)
    args, stage_args = parser.parse_known_args()
    if args.command != "run" and stage_args:
        parser.error(f"unrecognized arguments: {' '.join(stage_args)}")
    root = repository_root(args.root)

    if args.command == "stages":
        _print_stage_table(root)
        return
    if args.command == "capabilities":
        report = stage_capabilities(root)
        text = json.dumps(report, indent=2, ensure_ascii=False)
        if args.output:
            destination = args.output if args.output.is_absolute() else root / args.output
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(text + "\n", encoding="utf-8")
        print(text)
        return
    receipt = run_stage(
        args.stage_id,
        stage_args,
        root=root,
        dry_run=args.dry_run,
        receipt_path=args.receipt,
    )
    print(json.dumps(receipt, indent=2, ensure_ascii=False))
    if receipt.get("returncode") not in (None, 0):
        raise SystemExit(int(receipt["returncode"]))
