from __future__ import annotations

import argparse
from pathlib import Path

from src.frozen_replay import replay_frozen

ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce the frozen eight-candidate M4 evidence table")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data/frozen")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "run")
    parser.add_argument("--output", type=Path, default=ROOT / "results/results.csv")
    args = parser.parse_args()
    result = replay_frozen(args.data_dir, args.work_dir, args.output)
    print(result[["candidate_rank", "candidate_id", "six_channel_score", "four_context_status"]].to_string(index=False))


if __name__ == "__main__":
    main()
