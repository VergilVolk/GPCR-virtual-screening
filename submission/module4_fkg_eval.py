from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent

def export_evidence(summary_path: Path, output_path: Path) -> pd.DataFrame:
    data = json.loads(summary_path.read_text(encoding="utf-8-sig"))
    if data["source_status"] != "STAGE4_FULL_FROZEN_EVALUATION_COMPLETE":
        raise ValueError("Stage 4 frozen evaluation is not complete")
    rows = []
    for candidate_id, candidate in data["candidates"].items():
        rows.append({"candidate_id": candidate_id, "cluster": candidate["cluster"],
                     "trajectories": 12,
                     "production_ns_per_trajectory": data["dataset"]["production_ns_per_trajectory"],
                     "four_context_status": "complete_unresolved",
                     "robust_cooperative_signal": False,
                     "functional_interpretation": candidate["interpretation"][-1],
                     "claim_boundary": data["claim_boundary"]})
    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)
    return result

def main() -> None:
    parser = argparse.ArgumentParser(description="Export frozen four-context evidence")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "results/module4_function.csv")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    result = export_evidence(args.summary, args.output)
    print(result[["candidate_id", "cluster", "four_context_status"]].to_string(index=False))

if __name__ == "__main__":
    main()
