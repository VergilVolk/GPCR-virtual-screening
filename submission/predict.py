from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from module1_generate_filter import evaluate_library
from module3_ensemble_dock import fuse_six_channels
from module4_fkg_eval import export_evidence
from screen import run_screen


ROOT = Path(__file__).resolve().parent
TRACK = "Track 3: small-molecule virtual screening for GPCR targets"
PIPELINE_VERSION = "PACER-M4 submission 2026-10-04"


def run_pipeline(
    output_path: Path,
    candidate_library: Path = ROOT / "data/demo/candidate_library.csv",
    molecule_representations: Path = ROOT / "data/demo/m4_candidate_representations.npz",
    pocket_representations: Path = ROOT / "data/demo/m4_pocket_representations.npz",
    docking_summary: Path = ROOT / "data/demo/six_channel_docking_summary.csv",
    functional_summary: Path = ROOT / "data/demo/prospective_four_context_summary.csv",
) -> pd.DataFrame:
    results = ROOT / "results"
    chemistry = evaluate_library(
        candidate_library, results / "module1_chemistry.csv"
    )
    run_screen(
        molecule_representations,
        pocket_representations,
        results / "module2_binding.csv",
    )
    binding = pd.read_csv(results / "module2_binding.csv")
    structure = fuse_six_channels(
        docking_summary,
        results / "module3_structure.csv",
    )
    function = export_evidence(
        ROOT / "data/demo/four_context_validation.csv",
        functional_summary,
        results / "module4_function.csv",
    )

    table = chemistry.merge(binding, left_on="smiles", right_on="molecule_id", how="left")
    table = table.merge(
        structure, left_on="candidate_id", right_on="candidate_id", how="inner",
        suffixes=("", "_structure"),
    )
    prospective = function[function.evidence_set == "prospective_candidate"]
    prospective = prospective[["candidate_id", "interpretation"]].rename(
        columns={"interpretation": "functional_evidence"}
    )
    table = table.merge(prospective, on="candidate_id", how="left")
    table["functional_evidence"] = table["functional_evidence"].fillna(
        "not evaluated by prospective four-context MD"
    )
    table["track"] = TRACK
    table["model_version"] = PIPELINE_VERSION
    table["claim_status"] = "high-priority computational candidate; experimental validation required"
    table = table.sort_values(
        ["chemistry_pass", "six_channel_complete", "pacer_xr_score"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    table.insert(0, "final_rank", range(1, len(table) + 1))
    columns = [
        "final_rank", "candidate_id", "track", "smiles", "chemistry_pass", "qed",
        "binding_score", "binding_rank", "pacer_xr_score", "pacer_xr_rank",
        "candidate_cascade_rank", "six_channel_complete", "functional_evidence",
        "model_version", "claim_status",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table[columns].to_csv(output_path, index=False)
    return table[columns]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the packaged M4 screening example")
    parser.add_argument("--output", type=Path, default=ROOT / "results/results.csv")
    parser.add_argument("--candidate-library", type=Path, default=ROOT / "data/demo/candidate_library.csv")
    parser.add_argument("--molecule-representations", type=Path, default=ROOT / "data/demo/m4_candidate_representations.npz")
    parser.add_argument("--pocket-representations", type=Path, default=ROOT / "data/demo/m4_pocket_representations.npz")
    parser.add_argument("--docking-summary", type=Path, default=ROOT / "data/demo/six_channel_docking_summary.csv")
    parser.add_argument("--functional-summary", type=Path, default=ROOT / "data/demo/prospective_four_context_summary.csv")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    result = run_pipeline(
        args.output, args.candidate_library, args.molecule_representations,
        args.pocket_representations, args.docking_summary, args.functional_summary,
    )
    print(f"Wrote {len(result)} candidates to {args.output}")


if __name__ == "__main__":
    main()
